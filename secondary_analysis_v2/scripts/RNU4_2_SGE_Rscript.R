# Load all the libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(factoextra) 
  library(umap)       
  library(tibble)     
  library(scales)      
})



# 1) USER SETTINGS
# Main working directory 
WORKDIR <- "~/Documents/PDF/RNU4_2/secondary_analysis_v2"

# Variant scoring (Figure 2) 
VARIANT_DIR   <-  "metadata/secondary_analysis_v2"
VARIANT_FILE  <- "metadata/RNU4_2_function_scores_for_R_script.csv"

# Clinical PCA (Figure 3 + Supplementary Figure 4) 
CAT_FILE_1    <- "metadata/PCA_variant_labels.csv"
CLIN_FILE_1   <- "metadata/clinical_data_for_PCA.csv"

# New PCA/UMAP dataset (Supplementary Figure 6)
CAT_FILE_2    <- "metadata/category_data_for_PCA_UMAP.csv"
CLIN_FILE_2   <- "metadata/clinical_data_for_PCA_UMAP.csv"


# Variant to optionally exclude from PCA/UMAP 
DROP_VARIANT_FOR_PLOTTING <- "n.64_65insT"

# UMAP settings
RUN_UMAP <- TRUE
UMAP_SEED <- 250

# Output folder (relative to chosen WORKDIR)
OUTDIR <- file.path(WORKDIR, "outputs")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# 2) Global plot styling (consistent look across figures)
# Category colours used throughout
CATEGORY_COLORS <- c(
  "ReNU"         = "#e63946",
  "UKBB/AllofUs" = "#456990",
  "unobserved"   = "#49beaa",
  "strong"       = "#e63946",
  "moderate"     = "#ffbe0b",
  "NA"           = "grey"
)

# Variant type shapes used throughout
VARIANT_SHAPES <- c(
  "SNV"       = 16,  # circle
  "insertion" = 17,  # triangle
  "deletion"  = 15,  # square
  "del/dup"   = 15   # reuse square for CNV-like class
)

# Base theme
theme_pub <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text = element_text(color = "black")
    )
}

# PART A — Variant scoring plots (Figure 2 / Sup Fig 3)
setwd(WORKDIR)

# Read merged variant scoring table
variant_df <- read.csv(VARIANT_FILE, stringsAsFactors = FALSE)

# Replace empty strings with NA
variant_df[variant_df == ""] <- NA

# Minimal cleaning / harmonisation
variant_df <- variant_df %>%
  drop_na(HGVS) %>%
  mutate(
    # Harmonise category naming across plots
    category = case_when(category == "ReNU_syndrome" ~ "ReNU",
                         TRUE ~ category),
    # Enforce consistent type ordering
    variant_type = factor(Type_expanded, levels = c("SNV", "insertion", "deletion"))
  )

# If allele counts are missing, set to 0 (used in Fig 2D)
variant_df <- variant_df %>%
  mutate(
    UKBiobank_AC = ifelse(is.na(UKBiobank_AC), 0, UKBiobank_AC),
    AoU_AC       = ifelse(is.na(AoU_AC), 0, AoU_AC)
  )


# Figure 2A — function score vs transcript position (full length)
p_fig2a <- ggplot(
  variant_df,
  aes(x = plot_position, y = function_score, colour = category, shape = variant_type)
) +
  # Grey highlight blocks (critical regions)
  annotate("rect", fill = "grey", alpha = 0.4, xmin = 61.25, xmax = 70.25, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.4, xmin = 74.75, xmax = 78.25, ymin = -Inf, ymax = Inf) +
  # Vertical boundaries of the region
  geom_vline(xintercept = 61.25, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 79.75, linetype = "dashed", colour = "red") +
  # Threshold lines
  geom_hline(yintercept = -0.3, linetype = "dashed") +
  geom_hline(yintercept = -0.9, linetype = "dashed", colour = "grey") +
  # Points
  geom_point(size = 1.5) +
  # Label one key variant (adjust xlim/ylim for label placement)
  geom_text_repel(
    data = dplyr::filter(variant_df, HGVS == "n.64_65insT"),
    aes(label = HGVS),
    show.legend = FALSE,
    xlim = c(80, 120), ylim = c(-1.5, -1),
    colour = CATEGORY_COLORS["ReNU"]
  ) +
  scale_color_manual(values = CATEGORY_COLORS) +
  scale_shape_manual(values = VARIANT_SHAPES) +
  scale_x_continuous(limits = c(1, 145), breaks = seq(0, 145, 10)) +
  labs(x = "Position", y = "Function score") +
  theme_pub(base_size = 16)

ggsave(
  filename = file.path(OUTDIR, "Fig2A_score_by_position_full.pdf"),
  plot = p_fig2a, width = 10, height = 3, dpi = 400
)


# Sup Fig 3 — function score vs position (critical region only)
p_supfig3 <- ggplot(
  variant_df,
  aes(x = plot_position, y = function_score, colour = category, shape = variant_type)
) +
  annotate("rect", fill = "grey", alpha = 0.4, xmin = 61.25, xmax = 70.25, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.4, xmin = 74.75, xmax = 78.25, ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 61.25, colour = "red", linetype = "dashed") +
  geom_vline(xintercept = 79.75, colour = "red", linetype = "dashed") +
  geom_hline(yintercept = -0.3, linetype = "dashed") +
  geom_hline(yintercept = -0.9, linetype = "dashed", colour = "grey") +
  geom_point(size = 2) +
  # Label all ReNU variants (as in your original intent)
  geom_text_repel(
    data = dplyr::filter(variant_df, category == "ReNU"),
    aes(label = HGVS),
    show.legend = FALSE,
    colour = CATEGORY_COLORS["ReNU"],
    size = 2.5,
    max.overlaps = Inf,
    box.padding = 1.0,
    point.padding = 0.6,
    force = 3
  ) +
  scale_color_manual(values = CATEGORY_COLORS) +
  scale_shape_manual(values = VARIANT_SHAPES) +
  scale_x_continuous(limits = c(60, 81), breaks = seq(60, 80, 2)) +
  labs(x = "Position", y = "Function score") +
  theme_pub(base_size = 14) +
  theme(legend.position = "right")

ggsave(
  filename = file.path(OUTDIR, "Sup_Fig3_score_by_position_critical.pdf"),
  plot = p_supfig3, width = 7, height = 3, dpi = 400
)


# Figure 2B — density/histogram of function scores by category
p_fig2b <- ggplot(variant_df, aes(x = function_score, fill = category, colour = category)) +
  geom_histogram(aes(y = 0.5 * ..density..), alpha = 0.5, binwidth = 0.04, position = "stack") +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = -0.3, linetype = "dashed") +
  geom_vline(xintercept = -0.9, linetype = "dashed", colour = "grey") +
  scale_color_manual(values = CATEGORY_COLORS) +
  scale_fill_manual(values = CATEGORY_COLORS) +
  labs(x = "Function score", y = "Density") +
  theme_pub(base_size = 16) +
  theme(legend.position = c(0.25, 0.85))

ggsave(
  filename = file.path(OUTDIR, "Fig2B_density_plot.pdf"),
  plot = p_fig2b, width = 6.5, height = 4, dpi = 400
)


# A4) Figure 2D — allele count (UKB + AoU) vs function score (SNVs)

#subset do SNVs within transcripts
snv_df <- variant_df %>%
  filter(Type_expanded == "SNV", !is.na(HGVS)) %>%
  mutate(total_AC = UKBiobank_AC + AoU_AC)

p_fig2d <- ggplot(snv_df, aes(x = total_AC, y = function_score)) +
  geom_jitter(
    data = dplyr::filter(snv_df, UKBiobank_AC < 5),
    aes(color = category),
    alpha = 1, height = 0, width = 0.20, size = 0.8
  ) +
  geom_jitter(
    data = dplyr::filter(snv_df, UKBiobank_AC >= 5),
    aes(color = category),
    alpha = 1, height = 0, width = 0.00, size = 0.8
  ) +
  geom_hline(yintercept = -0.3, linetype = 3) +
  geom_vline(xintercept = 91, linetype = 3) +
  geom_vline(xintercept = 0.5, linetype = 3, colour = "#D3D3D3") +
  scale_color_manual(values = CATEGORY_COLORS) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 2),
    breaks = c(0, 1, 10, 100, 1000, 10000)
  ) +
  labs(x = "UKB and AoU allele count", y = "Function score") +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(OUTDIR, "Fig2D_UKBB_AoU_AC_vs_function_score.pdf"),
  plot = p_fig2d, width = 4, height = 3, dpi = 400
)


# Figure 2E — CADD vs function score (SNVs), with medians
# Compute medians for vertical reference lines
pop_con_med_CADD <- median(variant_df$CADD_score[variant_df$variant_type == "SNV" &
                                                   !is.na(variant_df$HGVS) &
                                                   variant_df$category == "UKBB/AllofUs"],
                           na.rm = TRUE)

ReNU_med_CADD <- median(variant_df$CADD_score[variant_df$variant_type == "SNV" &
                                                !is.na(variant_df$HGVS) &
                                                variant_df$category == "ReNU"],
                        na.rm = TRUE)

p_fig2e <- ggplot(snv_df, aes(x = CADD_score, y = function_score)) +
  geom_jitter(aes(color = category), alpha = 1, height = 0, width = 0.00, size = 0.8) +
  geom_hline(yintercept = -0.3, linetype = 3) +
  geom_vline(xintercept = pop_con_med_CADD, linetype = 2, colour = CATEGORY_COLORS["UKBB/AllofUs"]) +
  geom_vline(xintercept = ReNU_med_CADD,    linetype = 2, colour = CATEGORY_COLORS["ReNU"]) +
  scale_color_manual(values = CATEGORY_COLORS) +
  labs(x = "CADD score", y = "Function score") +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(OUTDIR, "Fig2E_CADD_vs_function_score.pdf"),
  plot = p_fig2e, width = 4, height = 3, dpi = 400
)



# PART B — Clinical PCA + phenotype barplot (Figure 3)
# Build the master annotated clinical table
# - Replace NA phenotype entries with 0 for PCA/UMAP
# - Join category/variant_type annotations by nucleotide-change label
# - Optionally hide one label for plotting (n.64_65insT)
categories_test <- read_csv(CAT_FILE_1)

df_PCA <- read_csv(clinical_df)

clinical_4_PCA = df_PCA %>%
  replace(is.na(.), 0) %>%
  mutate(label = case_when(`Nucleotide change` == "n.64_65insT" ~ NA_character_, TRUE ~ `Nucleotide change`)) %>%
  left_join(categories_test, by = c(`Nucleotide change` = "label")) %>%
  filter(category %in% c("strong", "moderate")) %>%
  mutate(category = case_when(is.na(category) ~ "NA", TRUE ~ category), category = factor(category, levels = c("strong", "moderate")), variant_type = case_when(is.na(variant_type) ~ "del/dup", TRUE ~ variant_type), variant_type = factor(variant_type, levels = c("SNV", "insertion")))


clinical_4_PCA_active = clinical_4_PCA %>%
  column_to_rownames("Patient") %>%
  select(-`Nucleotide change`, -Location, -Classification, -label, -category, -variant_type) %>%
  mutate_all(as.numeric)

res.pca = prcomp(clinical_4_PCA_active, scale = TRUE)

# Figure 3B - PCA (all ReNU variants)
p_pca_all<-fviz_pca_ind(res.pca, geom = c("point"), repel = TRUE, pointsize = 0.5) +
  geom_point(aes(colour = factor(clinical_4_PCA$category), shape = factor(clinical_4_PCA$variant_type)), size=2) +
  scale_color_manual(values = c(strong = "#e63946", moderate = "#ffbe0b")) +
  scale_shape_manual(values = c(SNV = 19, insertion = 17)) +
  guides(colour = guide_legend(title = "category"), shape = guide_legend(title = "variant type")) +
  theme_classic() +
  theme(plot.title = element_blank(), axis.title = element_text(size = 15), axis.text = element_text(size = 13, color = "black"), legend.title = element_blank(), legend.text = element_text(size = 13)) +
  geom_text_repel(aes(label = clinical_4_PCA$label), size = 2.5, na.rm = TRUE, min.segment.length = 0)
ggsave(
  filename = file.path(OUTDIR, "Fig3B_PCA_all_variants.pdf"),
  plot = p_pca_all, width = 6.5, height = 5, dpi = 400
)

# Figure 3C - clinical phenotype barplot
pheno<-read.csv("metadata/proportions_category_phenotypes.csv",header=TRUE)
pheno<-pheno %>%
  mutate(description=paste0(feature," - ",feature_categories)) %>%
  rowwise() %>%
  mutate(lower_confint=binom.test(x=patients, n=total_category_category, conf.level = 0.95)$conf.int[1]) %>%
  mutate(upper_confint=binom.test(x=patients, n=total_category_category, conf.level = 0.95)$conf.int[2])


pheno$description = factor(pheno$description, levels=c("ID - severe","ID - moderate","ID - mild","DD - severe","DD - moderate","DD - mild/none","speech - non verbal","speech - few words","speech - simple sentences/normal","epilepsy - yes","epilepsy - 1 episode","epilepsy - no"))
pheno$category = factor(pheno$category, levels=c("strong","moderate"))


bar_plot<- ggplot(pheno,aes(fill=category,x=description,y=proportion)) +
  geom_bar(position="dodge", stat="identity") +
  geom_linerange(aes(ymin=lower_confint, ymax=upper_confint), position=position_dodge(width = 1), alpha=0.9, colour="#888888") +
  geom_vline(xintercept = 3.5,linetype="dashed",colour="grey") +
  geom_vline(xintercept = 6.5,linetype="dashed",colour="grey") +
  geom_vline(xintercept = 9.5,linetype="dashed",colour="grey") +
  theme_classic() +
  scale_fill_manual(values = c(strong = "#e63946", moderate = "#ffbe0b")) +
  scale_x_discrete(labels=c("severe","moderate","mild","severe","moderate","mild/none","non-verbal","few words","simple\nsentences/normal","yes","one episode","no")) +
  xlab("") +
  ylab("proportion") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size = 14),legend.title=element_blank(),axis.text=element_text(size = 8))
ggsave(filename = file.path(OUTDIR, 'Fig3C_barplot.pdf'), plot=bar_plot, width = 5, height = 4, dpi = 400)


# Supplementary Figure 6A drop n.64_65insT
### Import data
categories_test <- read_csv(CAT_FILE_2)

df_PCA <- read_csv(CLIN_FILE_2)


### All variants
clinical_4_PCA = df_PCA %>%
  replace(is.na(.), 0) %>%
  mutate(label = case_when(`Nucleotide change` == "n.64_65insT" ~ NA_character_,
                           TRUE ~ `Nucleotide change`)) %>%
  left_join(categories_test, by = c(`Nucleotide change` = "label")) %>%
  mutate(category = case_when(is.na(category) ~ "NA",
                              TRUE ~ category),
         category = factor(category, levels = c("strong", "moderate", "NA")),
         variant_type = case_when(is.na(variant_type) ~ "del/dup",
                                  TRUE ~ variant_type),
         variant_type = factor(variant_type, levels = c("SNV", "insertion", "del/dup"))) 

clinical_4_PCA_active = clinical_4_PCA %>%
  column_to_rownames("Patient") %>%
  select(-`Nucleotide change`, -Location, -Classification, -label, -category, -variant_type) %>%
  mutate_all(as.numeric)


### Without n.64_65insT
clinical_4_PCA_wo = clinical_4_PCA %>%
  filter(`Nucleotide change` != "n.64_65insT") %>%
  select(-`Fracture/osteopenia`)

clinical_4_PCA_active_wo = clinical_4_PCA_wo %>%
  column_to_rownames("Patient") %>%
  select(-`Nucleotide change`, -Location, -Classification, -label, -category, -variant_type) %>%
  mutate_all(as.numeric)



### Supplementary Figure 6A - PCA Without n.64_65insT
res.pca.wo = prcomp(clinical_4_PCA_active_wo, scale = TRUE)

PCA_test <-fviz_pca_ind(res.pca.wo, 
                               geom = c("point"),
                               repel = TRUE) +
  geom_point(aes(colour = factor(clinical_4_PCA_wo$category),
                 shape = factor(clinical_4_PCA_wo$variant_type)), 
             size = 4) +
  scale_color_manual(values = c(strong = "#e63946", 
                                moderate = "#ffbe0b",
                                `NA` = "grey")) +
  scale_shape_manual(values = c(
    SNV = 16,         # circle
    insertion = 17   
  )) +
  guides(colour = guide_legend(title = "category"),
         shape = guide_legend(title = "variant type")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  geom_text_repel(aes(label = clinical_4_PCA_wo$label), 
                  size = 3, na.rm = TRUE,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  box.padding = 1.0,
                  point.padding  = 0.6,
                  force = 3)

ggsave(filename = file.path(OUTDIR, "Sup_Fig6a_PCA_wo_n64_65insT.pdf"),
       dpi = 600, width = 8, height = 6, 
       plot = PCA_test)
### UMAP

#compute UMAP
set.seed(250)
umap_result <- umap::umap(clinical_4_PCA_active)

# Convert results to data frame
umap_df <- as.data.frame(umap_result$layout) %>%
  rename(UMAP1 = "V1",
         UMAP2 = "V2") %>%
  rownames_to_column("Patient") %>%
  mutate(Patient = as.double(Patient)) %>%
  left_join(clinical_4_PCA_Nicky %>%
              select(Patient, Classification:category))

# Supplementary Figure 6B - UMAP
umap_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = category, shape = variant_type), size = 4) +
  scale_color_manual(values = c(strong = "#e63946", 
                                moderate = "#ffbe0b",
                                `NA` = "grey")) +
  scale_shape_manual(values = c(
    SNV = 16,         # circle
    insertion = 17   
  )) +
  guides(colour = guide_legend(title = "category"),
         shape = guide_legend(title = "variant type")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  geom_text_repel(aes(label = label), 
                  size = 3, na.rm = TRUE,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  box.padding = 1.0,
                  point.padding  = 0.6,
                  force = 3)

ggsave(filename = file.path(OUTDIR, "Sup_Fig6B_umap.pdf"), dpi = 600, width = 8, height = 6, plot = umap_umap)

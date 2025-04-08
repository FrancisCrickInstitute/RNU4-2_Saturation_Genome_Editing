library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)

setwd("~/RNU4_2/")
data <-read.csv("RNU4_2_ST1.csv")

#colour_scheme
category_colors = c('ReNU' = '#e63946','UKBB/AllofUs' = '#456990', 'unobserved' = '#49beaa')

#fill empty entries outside of transcript with NA
data[data == ""] <- NA
# Plot Figure 2A - gene length plot
data<-data %>%
  drop_na(HGVS) %>%
  mutate(category=case_when(category=="ReNU_syndrome" ~ "ReNU", TRUE ~ category))

data$variant_type = factor(data$Type_expanded, levels=c("SNV","insertion","deletion"))

p<-ggplot(data, aes(x=position_oligonucleotide,y=function_score, colour=category, shape=variant_type)) +
  annotate("rect", fill="grey", alpha=0.4, xmin=61.25,xmax=70.25,ymin=-Inf,ymax=Inf) +
  annotate("rect", fill="grey", alpha=0.4, xmin=74.75,xmax=78.25,ymin=-Inf,ymax=Inf) +
  geom_vline(aes(xintercept=61.25),linetype="dashed",colour="red") +
  geom_vline(aes(xintercept=79.75),linetype="dashed",colour="red") +
  geom_hline(aes(yintercept=-0.39),linetype="dashed") +
  geom_hline(aes(yintercept=-1.0),linetype="dashed",colour="grey") +
  geom_point(size=1.5) +
  theme_classic() +
  xlab("Position") +
  ylab("Function score") +
  scale_color_manual(values=c("red","#e63946","#456990","#49beaa")) +
  scale_x_continuous(limits = c(1,145), breaks = seq(0,145,10)) +
  geom_text_repel(data = filter(data, HGVS=='n.64_65insT'), aes(label=HGVS),show.legend=FALSE, xlim=c(80,120), ylim=c(-1.5,-1),colour="#e63946") +
  theme(text = element_text(size = 16),legend.title=element_blank())
ggsave("Fig2A_score_by_position.pdf", plot = p, width = 10, height = 3, dpi = 400)

# Plot Figure 2B - density plots
p<-ggplot(data, aes(function_score, fill=category, colour=category)) +
  geom_histogram(aes(y=0.5*..density..),alpha = 0.5,binwidth=0.05,position="stack") +
  geom_density(alpha = 0.3) +
  geom_vline(aes(xintercept=-0.39),linetype="dashed") +
  geom_vline(aes(xintercept=-1.0),linetype="dashed",colour="grey") +
  theme_classic() +
  xlab("Function score") +
  ylab("Density") +
  scale_color_manual(values=category_colors) +
  scale_fill_manual(values=category_colors) +
  theme(text = element_text(size = 16),legend.title=element_blank(),legend.position=c(0.25,0.85))
ggsave("Fig2B_density_plot_class.pdf", plot = p, width = 6.5, height = 4, dpi = 400)

# Plot Figure 2D - UKBB Allele Count vs function scores
data[is.na(data$UKBiobank_AC),]$UKBiobank_AC <- 0

p<-ggplot(data[which(data$Type_expanded == "SNV" & !is.na(data$HGVS)),], aes(x= UKBiobank_AC , y=function_score))+
  geom_jitter(alpha=1,height=0.00,width=0.20,aes(color=category),size=0.8,data=data[which(data$Type_expanded == "SNV" & !is.na(data$HGVS) & data$UKBiobank_AC < 5),])+
  geom_jitter(alpha=1,height=0.00,width=0.00,aes(color=category),size=0.8,data=data[which(data$Type_expanded == "SNV" & !is.na(data$HGVS) & data$UKBiobank_AC >= 5),])+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.title=element_blank(),legend.position = "None",legend.key = element_blank(),axis.line.x = element_line(colour = "grey50"),axis.line.y = element_line(colour = "grey50"))+
  ylab('Function score')+
  xlab("UKBB allele count")+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(text = element_text(size=12),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))+
  geom_hline(yintercept=-0.39, lty=3)+
  scale_color_manual(values=category_colors)+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 2),breaks = c(0,1,10,100,1000,10000))+geom_vline(xintercept=58.7, lty=3)+
  geom_vline(xintercept=0.5, lty=3, color = "#D3D3D3")
ggsave("Fig2D_UKBB_AC_vs_function_score.pdf", plot = p, width = 4, height = 3, dpi = 400)

# Plot Figure 2E - CADD vs function scores
#calculate mean scores
ReNU_mean_CADD <- mean(data[which(data$variant_type == "SNV" & !is.na(data$HGVS) & data$category == "ReNU" ),]$CADD_score)

pop_con_mean_CADD <- mean(data[which(data$variant_type == "SNV" & !is.na(data$HGVS) & data$category == "UKBB/AllofUs" ),]$CADD_score)

ReNU_med_CADD <- median(data[which(data$variant_type == "SNV" & !is.na(data$HGVS) & data$category == "ReNU" ),]$CADD_score)

pop_con_med_CADD <- median(data[which(data$variant_type == "SNV" & !is.na(data$HGVS) & data$category == "UKBB/AllofUs" ),]$CADD_score)

p <- ggplot(data[which(data$variant_type == "SNV" & !is.na(data$HGVS)),], aes(x=CADD_score, y=function_score))+
  geom_jitter(alpha=1,height=0.00,width=0.00,aes(color=category),size=0.8)+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.title=element_blank(),legend.position = "None",legend.key = element_blank(),axis.line.x = element_line(colour = "grey50"),axis.line.y = element_line(colour = "grey50"))+
  ylab('Function score')+xlab("CADD score")+theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(text = element_text(size=12),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),strip.background = element_blank())+
  geom_hline(yintercept=-0.39, lty=3)+scale_color_manual(values=category_colors)+
  geom_vline(xintercept=pop_con_med_CADD, lty=2, color = "#49698D")+
  geom_vline(xintercept=ReNU_med_CADD, lty=2, color = "#E63B49")
ggsave("Fig2E_CADD_vs_function_score.pdf", plot = p, width = 4, height = 3, dpi = 400)


## Figure 3
file_categories = "heatmap_PCA_test.csv"
file_clinical = "clinical_data_for_heatmap_PCA.csv"

categories_test <- read_csv(file_categories)

df_PCA <- read_csv(file_clinical)

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

# Plot (FIGURE 3A) - PCA
p<-fviz_pca_ind(res.pca, geom = c("point"), repel = TRUE, pointsize = 0.5) +
  geom_point(aes(colour = factor(clinical_4_PCA$category), shape = factor(clinical_4_PCA$variant_type)), size=2) +
  scale_color_manual(values = c(strong = "#e63946", moderate = "#ffbe0b")) +
  scale_shape_manual(values = c(SNV = 19, insertion = 17)) +
  guides(colour = guide_legend(title = "category"), shape = guide_legend(title = "variant type")) +
  theme_classic() +
  theme(plot.title = element_blank(), axis.title = element_text(size = 15), axis.text = element_text(size = 13, color = "black"), legend.title = element_blank(), legend.text = element_text(size = 13)) +
  geom_text_repel(aes(label = clinical_4_PCA$label), size = 2.5, na.rm = TRUE, min.segment.length = 0)
ggsave('Fig3A_PCA.pdf', plot=p, width = 5.5, height = 4, dpi = 400)


# detailed phenotype info by category
pheno<-read.csv("proportions_category_phenotypes_v3.csv",header=TRUE)
pheno<-pheno %>%
  mutate(description=paste0(feature," - ",feature_categories)) %>%
  rowwise() %>%
  mutate(lower_confint=binom.test(x=patients, n=total_category_category, conf.level = 0.95)$conf.int[1]) %>%
  mutate(upper_confint=binom.test(x=patients, n=total_category_category, conf.level = 0.95)$conf.int[2])


pheno$description = factor(pheno$description, levels=c("ID - severe","ID - moderate","ID - mild","DD - severe","DD - moderate","DD - mild/none","speech - non verbal","speech - few words","speech - simple sentences/normal","epilepsy - yes","epilepsy - 1 episode","epilepsy - no"))
pheno$category = factor(pheno$category, levels=c("strong","moderate"))

# FIGURE 3B - barplot
p<- ggplot(pheno,aes(fill=category,x=description,y=proportion)) +
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
ggsave('Fig3B_barplot.pdf', plot=p, width = 5, height = 4, dpi = 400)





library("vegan")
library("tidyverse")
library("ggpubr")
library("RColorBrewer")
library("ComplexHeatmap")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

## Data preparation ----
ra <- "0_Input_abundances.tsv"
sampleData <- "0_Input_Metadata.tsv"

## Load data
adat <- read.delim(ra, row.names = 1, check.names = F, strip.white = T)
sdat <- read.delim(sampleData, row.names = 1, check.names = F)

## Order relative abundance as sample data
adat <- adat[, row.names(sdat)]

## Set the reference covariate
refcov <- "Class"

sdat$Class <- paste0(sdat$Implant, "_", sdat$Sludge)
sdat$Sludge <- factor(sdat$Sludge, levels=c("Primary", "Secondary", "Pre-thickening", "Digestate"))

ch <- colorRampPalette(brewer.pal(9,"RdBu"))(9)
colorh2 <- c(ch[1], ch[2], ch[3], ch[8])

## 2.1 alpha diversity analysis ----

## Species isolation and Richness
genus <- adat[grep("G__", row.names(adat)),]
species <- adat[grep("S__", row.names(adat)),]

## Convert to relative abundance
genus <- prop.table(as.matrix(genus), 2)
species <- prop.table(as.matrix(species), 2)

## Richness estimation
sdat$richness_genus <- colSums(genus>0)
sdat$richness_species <- colSums(species>0)

## Richness plot ----
ggplot(sdat, aes(x=Sludge, y=richness_genus, col=Sludge))+
  geom_violin(trim=FALSE, alpha=1, width=1) +
  geom_boxplot(width = 0.15, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, alpha=0.75, size=1.5, aes(shape=Implant))+
  theme_bw() +
  labs(y="Richness (Genus)", x="Sludge", col="Sludge") +
  scale_color_manual(values=colorh2) +
  facet_wrap(~Implant)+
  geom_pwc(hide.ns = "p", label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=8))

ggsave("2.1_Boxplot_Alfa_diversity_Richness_genus.pdf")

ggplot(sdat, aes(x=Sludge, y=richness_species, col=Sludge))+
  geom_violin(trim=FALSE, alpha=1, width=1) +
  geom_boxplot(width = 0.15, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, alpha=0.75, size=1.5, aes(shape=Implant))+
  theme_bw() +
  labs(y="Richness (Species)", x="Sludge", col="Sludge") +
  scale_color_manual(values=colorh2) +
  facet_wrap(~Implant)+
  geom_pwc(hide.ns = "p", label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=8))

ggsave("2.1_Boxplot_Alfa_diversity_Richness_species.pdf")

## alpha diversity (Shannon and Inverse Simpson diversity index) -----
sdat$diversity_shannon_genus <- diversity(genus, MARGIN = 2, index = "shannon")
sdat$diversity_shannon_species <- diversity(species, MARGIN = 2, index = "shannon")

sdat$diversity_inverse_simpson_genus <- diversity(genus, MARGIN = 2, index = "invsimpson")
sdat$diversity_inverse_simpson_species <- diversity(species, MARGIN = 2, index = "invsimpson")

## Shannon ------
ggplot(sdat, aes(x=Sludge, y=diversity_shannon_genus, col=Sludge))+
  geom_violin(trim=FALSE, alpha=1, width=1) +
  geom_boxplot(width = 0.15, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, alpha=0.75, size=1.5, aes(shape=Implant))+
  theme_bw() +
  labs(y="Shannon index (Genus)", x="Sludge", col="Sludge") +
  scale_color_manual(values=colorh2) +
  facet_wrap(~Implant)+
  geom_pwc(hide.ns = "p", label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=8))

ggsave("2.1_Boxplot_Alfa_diversity_Shannon_genus.pdf")

ggplot(sdat, aes(x=Sludge, y=diversity_shannon_species, col=Sludge))+
  geom_violin(trim=FALSE, alpha=1, width=1) +
  geom_boxplot(width = 0.15, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, alpha=0.75, size=1.5, aes(shape=Implant))+
  theme_bw() +
  labs(y="Shannon index (Species)", x="Sludge", col="Sludge") +
  scale_color_manual(values=colorh2) +
  facet_wrap(~Implant)+
  geom_pwc(hide.ns = "p", label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=8))

ggsave("2.1_Boxplot_Alfa_diversity_Shannon_species.pdf")

## Inverse Simpson -----
ggplot(sdat, aes(x=Sludge, y=diversity_inverse_simpson_genus, col=Sludge))+
  geom_violin(trim=FALSE, alpha=1, width=1) +
  geom_boxplot(width = 0.15, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, alpha=0.75, size=1.5, aes(shape=Implant))+
  theme_bw() +
  labs(y="Inverse Simpson index (Genus)", x="Sludge", col="Sludge") +
  scale_color_manual(values=colorh2) +
  facet_wrap(~Implant)+
  geom_pwc(hide.ns = "p", label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=8))

ggsave("2.1_Boxplot_Alfa_diversity_Simpson_genus.pdf")

ggplot(sdat, aes(x=Sludge, y=diversity_inverse_simpson_species, col=Sludge))+
  geom_violin(trim=FALSE, alpha=1, width=1) +
  geom_boxplot(width = 0.15, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, alpha=0.75, size=1.5, aes(shape=Implant))+
  theme_bw() +
  labs(y="Inverse Simpson index (Species)", x="Sludge", col="Sludge") +
  scale_color_manual(values=colorh2) +
  facet_wrap(~Implant)+
  geom_pwc(hide.ns = "p", label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=8))

ggsave("2.1_Boxplot_Alfa_diversity_Simpson_species.pdf")

## Evenness ------
sdat$evenness_genus<- diversity(genus, MARGIN = 2, index = "shannon") / log(sdat$richness_genus)
sdat$evenness_species <- diversity(species, MARGIN = 2, index = "shannon") / log(sdat$richness_species)

ggplot(sdat, aes(x=Sludge, y=evenness_genus, col=Sludge))+
  geom_violin(trim=FALSE, alpha=1, width=1) +
  geom_boxplot(width = 0.15, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, alpha=0.75, size=1.5, aes(shape=Implant))+
  theme_bw() +
  labs(y="Evenness (Genus)", x="Sludge", col="Sludge") +
  scale_color_manual(values=colorh2) +
  facet_wrap(~Implant)+
  geom_pwc(hide.ns = "p", label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=8))

ggsave("2.1_Boxplot_Alfa_diversity_Richness_genus.pdf")

ggplot(sdat, aes(x=Sludge, y=evenness_species, col=Sludge))+
  geom_violin(trim=FALSE, alpha=1, width=1) +
  geom_boxplot(width = 0.15, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, alpha=0.75, size=1.5, aes(shape=Implant))+
  theme_bw() +
  labs(y="Evenness (Species)", x="Sludge", col="Sludge") +
  scale_color_manual(values=colorh2) +
  facet_wrap(~Implant)+
  geom_pwc(hide.ns = "p", label = "p.signif")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size=8))

ggsave("2.1_Boxplot_Alfa_diversity_Richness_species.pdf")

write.table(sdat, "2.1_Table_Diversity_metrics.txt", sep="\t", quote=F)

## 2.2 Beta diversity ----
library("caret")
beta_dist <- vegdist(t(genus), index = "bray")

mds <- metaMDS(beta_dist)

mds_data <- as.data.frame(mds$points)
mds_data$SampleID <- rownames(mds_data)
mds_data <- cbind(mds_data, sdat)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = factor(Sludge), shape=Implant)) +
  geom_point(cex=2)+
  labs(color="Sludge")+
  scale_color_manual(values=colorh2) +
  theme_bw()

ggsave("2.2_MDS_Beta_diversity_genus.pdf")

## 2.2 PERMANOVA ----
y_permanova <- adonis2(beta_dist ~ Implant + Sludge,
                       data=sdat,
                       permutations=999,
                       method="bray", 
                       parallel=4, by="term")

write.table(y_permanova, "2.2_Table_PERMANOVA_results_genus.tsv", sep="\t", quote=F)

### Phylum
## Barplot distributions
phylum <- adat[grep("P__", row.names(adat)),]

phylum <- t(phylum)

phyla <- data.frame(cbind(ID=row.names(phylum), phylum, Class = as.character(sdat$Class)))

phylum2plot <- gather(phyla, Phylum, Abundance, -c(Class))

phylum2plot$Phylum <- factor(phylum2plot$Phylum, levels=unique(phylum2plot$Phylum))

phylum2plot$Abundance <- as.numeric(phylum2plot$Abundance)

colbar <- colorRampPalette(brewer.pal(11, "Spectral"))(9)

phylum2plot$Class <- factor(phylum2plot$Class, levels=c("Outdoor", "Indoor", "Indoor_VGT", "Indoor_VGT_eff_air"))

ggplot(phylum2plot, aes(x=ID, y=Abundance, fill=Phylum)) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(cols = vars(Class), scales="free", space = "free") +
  theme_bw() +
  scale_fill_manual(values = colbar) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7))+
  labs(y="Relative abundance")

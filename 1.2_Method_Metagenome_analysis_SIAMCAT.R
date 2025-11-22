library("SIAMCAT")
library("tidyverse")
library("RColorBrewer")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

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

### Avg and sd levels - genus
dat <- data.frame(Class=sdat$Class, t(genus), check.names = F)
dat <- gather(dat, Species, Level, -c(Class))
out <- dplyr::group_by(dat, Species, Class) %>% 
  dplyr::summarize(mean=mean(Level, na.rm=T), 
            sd=sd(Level, na.rm=T),
            prev=sum(Level>0)/n())

out <- pivot_wider(out, names_from = Class, values_from = c(mean, sd, prev))

write.table(out, "2.3_Table_Summary_tested_genus.tsv", sep="\t", quote=F, row.names = F)

### Avg and sd levels - species
dat <- data.frame(Class=sdat$Class, t(species), check.names = F)
dat <- gather(dat, Species, Level, -c(Class))
out <- dplyr::group_by(dat, Species, Class) %>% 
  dplyr::summarize(mean=mean(Level, na.rm=T), 
                   sd=sd(Level, na.rm=T),
                   prev=sum(Level>0)/n())

out <- pivot_wider(out, names_from = Class, values_from = c(mean, sd, prev))

write.table(out, "2.3_Table_Summary_tested_species.tsv", sep="\t", quote=F, row.names = F)

## SIAMCAT analysis -----
library("caret")
species_nzv <- nearZeroVar(t(species))
species_filtered <- species[-species_nzv, ]

run_siamcat <- function(case, control){
  
### creation of siamcat object
sc.obj <- siamcat(feat=genus,
                  meta=sdat, 
                  label='Sludge',
                  case=case,
                  control=control
                  )

### filtering
sc.obj <- filter.features(sc.obj, cutoff=0,
                          filter.method = 'abundance')
# Tests between times
sc.obj <- check.associations(sc.obj,
                             log.n0 = 1e-06,
                             alpha=0.001,
                             test = "wilcoxon",
                             mult.corr = "BH"
                             )

# Tests between classes
association.plot(sc.obj, 
                 fn.plot = paste0("2.3_Figure_DA_plot_", 
                 case, "_vs_", control, ".pdf"),
                 panels = c('fc'), color.scheme = c("#B2182B", "#4393C3"))

### Write output
out <- associations(sc.obj)
write.table(out, paste0("2.3_Table_DA_analysis_", 
                        case, "_vs_", control, ".tsv"),
            sep="\t", quote=F)

}

run_siamcat("Post", "Pre")
run_siamcat("Digestate", "Primary")
run_siamcat("Digestate", "Secondary")
run_siamcat("Digestate", "Pre-thickening")

#### Plot DA species -----
toplot <- cbind(sdat, t(species))

# Heatmap log2FC
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")

DA <- read.delim("0_Input_DA_taxa.tsv", row.names=1, check.names=F)
sdat <- read.delim("0_Input_Metadata.tsv", row.names=1, check.names=F)
adat <- read.delim("0_Input_abundances.tsv", row.names=1, check.names=F)

DAgenus <- DA[grep("G__", row.names(DA)),]
DAspecies <- DA[grep("S__", row.names(DA)),]

colorh<-colorRampPalette(brewer.pal(11,"PRGn"))(256)

h1 <- Heatmap(DAgenus,
              row_names_gp = gpar(fontsize = 2),
              column_names_gp = gpar(fontsize = 6),
              cluster_columns = F,
              cluster_rows = T,
              name = "log2FC",
              col=colorRamp2(breaks=c(-4,-2,0,2,4), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
              border = "black",
              column_names_rot = 45,
              heatmap_legend_param = list(border = "black"),
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2")

# Heatmap abundance -----
daspecies <- adat[row.names(DAspecies), row.names(sdat)]
daspecies <- t(scale(t(log10(daspecies+0.1))))

##
colorh<-colorRampPalette(brewer.pal(11,"PRGn"))(256)
implant.c <- colorRampPalette(brewer.pal(9, "YlOrRd"))(9)

ch <- colorRampPalette(brewer.pal(9,"RdBu"))(9)
colorh2 <- c(ch[1], ch[2], ch[3], ch[8])

sludge.cv <- c("Primary"=ch[1], "Secondary"=ch[2], "Pre-thickening"=ch[3], "Digestate"=ch[8])

implant.cv <- c("WWTP-S"=implant.c[1], "WWTP-M"=implant.c[4])

# Definition of the annotations
column_ha = HeatmapAnnotation(Slugde=sdat$Sludge,
                              Implant=sdat$Implant,
                              col = list(Slugde = sludge.cv,
                                         Implant = implant.cv),
                              border = TRUE,
                              simple_anno_size = unit(0.3, "cm"),
                              annotation_name_gp = gpar(fontsize = 5),
                              annotation_legend_param = list(border="black"))

h2 <- Heatmap(daspecies,
        row_names_gp = gpar(fontsize = 2),
        column_names_gp = gpar(fontsize = 4),
        top_annotation = column_ha,
        cluster_columns = T,
        cluster_rows = T,
        name = "Z-score",
        border = "black",
        col=colorRamp2(breaks=c(-4,-2,0,2,4), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
        column_names_rot = 45,
        heatmap_legend_param = list(border = "black"),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        clustering_distance_columns = "spearman",
        clustering_distance_rows = "spearman")
h1+h2

ggplot("2.3_Figure_Heatmap_DA_genus.pdf")
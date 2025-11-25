Data_Heatmap <- read_excel("Data_Heatmap.xlsx")
View(Data_Heatmap)
df <- Data_Heatmap

library(tidyverse)
library(ggplot2)

metadata_cols <- c("Sample_ID", "Sample_name",
                   "Implant_(1=WWTP-S_2=WWTP-M)",
                   "Sludge_type_(1=primary_sludge_2=secondary_sludge_3=pre-thickening_sludge_4=digestate)",
                   "AD_(1=pre_2=post)",
                   "DNA_ng/ul","RNA_ng/ul","pH","T°_(°C)","T°_digestate_(°C)",
                   "SST%","SSV%","%SSV_tot_(da_usare!)",
                   "16S_(Log_gene_copies/kg_of_sludge)",
                   "Reads_Total","Fragments_After_De-hosting","Fragments_Classified",
                   "%_Fragments_Classified","Fragments_Not_Classified",
                   "Shannon_Index","Species_Number","Sample_group")


df_clean <- df[!is.na(df$Sample_group) & df$Sample_group != "NA" & df$Sample_group != "NA_NA", ]
species_df <- df_clean %>% select(-all_of(metadata_cols[metadata_cols %in% colnames(df_clean)]))
species_df <- species_df %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
species_df <- species_df %>% select(where(~ !all(is.na(.))))
species_df[is.na(species_df)] <- 0


species_rel_with_names <- species_df %>%
  mutate(Sample_original = df_clean$Sample_group)

df_long <- species_rel_with_names %>%
  pivot_longer(cols = -Sample_original, names_to = "Species", values_to = "Relative_Abundance") %>%
  rename(Sample = Sample_original)

df_long <- df_long %>% 
  filter(str_detect(Sample, "WWTP-[MS]_(Primary|Secondary|Pre-thickening|Digestate)"))

df_long$Sample_clean <- str_replace(df_long$Sample, "\\.\\d+$", "")


df_long <- df_long %>%
  group_by(Sample_clean, Species) %>%
  summarise(Relative_Abundance = mean(Relative_Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(Sample = Sample_clean)

# DATA ORDERING
valid_samples <- c("WWTP-M_Primary","WWTP-S_Primary",
                   "WWTP-M_Secondary","WWTP-S_Secondary", 
                   "WWTP-M_Pre-thickening","WWTP-S_Pre-thickening",
                   "WWTP-M_Digestate","WWTP-S_Digestate")

df_long$Sample <- factor(df_long$Sample, levels = valid_samples)

# Order species by mean relative abundance
species_means <- df_long %>%
  group_by(Species) %>%
  summarise(mean_abundance = mean(Relative_Abundance, na.rm = TRUE)) %>%
  arrange(mean_abundance)

df_long$Species <- factor(df_long$Species, levels = species_means$Species)

# LOG TRANSFORMATION FUNCTIONS
apply_log_transform <- function(data, constant = 0.001) {
  data %>%
    mutate(
      Log_Abundance = log10(Relative_Abundance + constant),
      # Also create log2 version for comparison
      Log2_Abundance = log2(Relative_Abundance + constant)
    )
}

# HEATMAP CREATION FUNCTIONS
create_heatmap_standard <- function(data, color_high, max_scale = 1, title_suffix = "") {
  ggplot(data, aes(x = Sample, y = Species, fill = Relative_Abundance)) +
    geom_tile(color = "white", size = 0.5) +
    geom_vline(xintercept = c(2.5, 4.5, 6.5), color = "white", size = 3) +
    scale_fill_gradient(low = "#f7f7f7", high = color_high, 
                        limits = c(0, max_scale), oob = scales::squish,
                        name = "Relative\nAbundance (%)",
                        breaks = seq(0, max_scale, by = max_scale/4),
                        labels = function(x) sprintf("%.3f", x)) +
    annotate("text", x = 1.5, y = max(as.numeric(data$Species)) + 3, 
             label = "1", size = 6, fontface = "bold", color = "black") +
    annotate("text", x = 3.5, y = max(as.numeric(data$Species)) + 3, 
             label = "2", size = 6, fontface = "bold", color = "black") +
    annotate("text", x = 5.5, y = max(as.numeric(data$Species)) + 3, 
             label = "3", size = 6, fontface = "bold", color = "black") +
    annotate("text", x = 7.5, y = max(as.numeric(data$Species)) + 3, 
             label = "4", size = 6, fontface = "bold", color = "black") +
    scale_x_discrete(labels = c("WWTP-M", "WWTP-S", "WWTP-M", "WWTP-S", 
                                "WWTP-M", "WWTP-S", "WWTP-M", "WWTP-S")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.key.height = unit(1, "cm"),
      axis.title = element_blank(),
      plot.margin = margin(40, 10, 10, 50),
      axis.ticks.x = element_blank()
    ) +
    coord_cartesian(clip = 'off') +
    labs(x = "", y = "", title = title_suffix)
}


create_heatmap_log10 <- function(data, color_high, title_suffix = "") {
  ggplot(data, aes(x = Sample, y = Species, fill = Log_Abundance)) +
    geom_tile(color = "white", size = 0.5) +
    geom_vline(xintercept = c(2.5, 4.5, 6.5), color = "white", size = 3) +
    scale_fill_gradient(low = "#f7f7f7", high = color_high,
                        name = "Log10\n(Abundance + 0.001)") +
    annotate("text", x = 1.5, y = max(as.numeric(data$Species)) + 3, 
             label = "1", size = 6, fontface = "bold", color = "black") +
    annotate("text", x = 3.5, y = max(as.numeric(data$Species)) + 3, 
             label = "2", size = 6, fontface = "bold", color = "black") +
    annotate("text", x = 5.5, y = max(as.numeric(data$Species)) + 3, 
             label = "3", size = 6, fontface = "bold", color = "black") +
    annotate("text", x = 7.5, y = max(as.numeric(data$Species)) + 3, 
             label = "4", size = 6, fontface = "bold", color = "black") +
    scale_x_discrete(labels = c("WWTP-M", "WWTP-S", "WWTP-M", "WWTP-S", 
                                "WWTP-M", "WWTP-S", "WWTP-M", "WWTP-S")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.key.height = unit(1, "cm"),
      axis.title = element_blank(),
      plot.margin = margin(40, 10, 10, 50),
      axis.ticks.x = element_blank()
    ) +
    coord_cartesian(clip = 'off') +
    labs(x = "", y = "", title = title_suffix)
}

cat("Applying log transformations...\n")

max_abundance <- max(df_long$Relative_Abundance, na.rm = TRUE)
min_abundance <- min(df_long$Relative_Abundance[df_long$Relative_Abundance > 0], na.rm = TRUE)
cat("Original data range: ", sprintf("%.4f", min_abundance), "% to ", sprintf("%.4f", max_abundance), "%\n")

df_long_log <- apply_log_transform(df_long, constant = 0.001)

cat("Log10 range: ", sprintf("%.3f", min(df_long_log$Log_Abundance, na.rm = TRUE)), 
    " to ", sprintf("%.3f", max(df_long_log$Log_Abundance, na.rm = TRUE)), "\n")
cat("Log2 range: ", sprintf("%.3f", min(df_long_log$Log2_Abundance, na.rm = TRUE)), 
    " to ", sprintf("%.3f", max(df_long_log$Log2_Abundance, na.rm = TRUE)), "\n")


if (max_abundance <= 1) {
  scale_max <- 1
} else if (max_abundance <= 5) {
  scale_max <- 5  
} else {
  scale_max <- 20
}

# Create standard heatmaps (for comparison)
cat("Creating standard heatmaps for comparison...\n")
p_red_standard <- create_heatmap_standard(df_long, "darkred", scale_max, "")

# Create log-transformed heatmaps
cat("Creating log10-transformed heatmaps...\n")
p_darkblue_log10 <- create_heatmap_log10(df_long_log, "#003366", "")
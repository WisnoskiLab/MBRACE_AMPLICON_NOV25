# Run from the repository root.
# ---- Core data handling & manipulation ----
library(tidyverse)

# ---- Microbiome data processing ----
library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(MicrobiomeStat)
library(meconetcomp)
library(WGCNA)
library(ape)
library(picante)
library(Biostrings)

# ---- Differential abundance & compositional analysis ----
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)

# ---- Visualization ----
library(ggplot2)
library(ggnested)
library(ggpubr)
library(ggtree)
library(tidygraph)
library(paletteer)
library(colorspace)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggraph)
library(patchwork)
library(linkET)

# ---- Statistical modeling ----
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(multcompView)
library(usethis)
library(nlMS)
library(iCAMP)
library(minpack.lm)
library(Hmisc)
library(ggpmisc)

dir.create("results/figures/alpha_diversity", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/beta_diversity", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/alpha_beta", recursive = TRUE, showWarnings = FALSE)

location_colors = c("Biloxi Bay" = "#0072B2", "Pascagoula Bay" = "#E69F00")

group_order = c("Limnetic (0-0.5)", "Oligohaline (0.5-5)", "Mesohaline (5-18)", "Polyhaline (18-30)")

group_colors = setNames(c("#A6E1D4", "#5FADAE", "#2B7BBA", "#0D326A"), group_order)

biom = import_biom("data/processed/trimmed_200_230/w_nov_2025_200-230_feture.biom")
metadata = import_qiime_sample_data("data/processed/trimmed_200_230/metadata.txt")
tree = read_tree("data/processed/trimmed_200_230/w_nov_2025_200_230_unrooted-tree.nwk")
rep_fasta = readDNAStringSet("data/processed/trimmed_200_230/dna-sequences_230_200.fasta", format = "fasta")

w_biom = merge_phyloseq(biom, rep_fasta, metadata, tree)

colnames(tax_table(w_biom)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

w_nov_2025 = phyloseq2meco(w_biom)
w_nov_2025$tidy_dataset()


sample_depth = colSums(w_nov_2025$otu_table)

summary(sample_depth)
quantile(sample_depth, probs = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 1))

stopifnot(min(sample_depth) >= 40000)

depth_df = tibble(SampleID = names(sample_depth), Reads = as.numeric(sample_depth))

ggplot(depth_df, aes(Reads)) +
  geom_histogram(bins = 20, fill = "#0072B2", color = "white") +
  theme_classic(base_size = 13) +
  labs(x = "Sequencing reads per sample", y = "Number of samples")


tmp = trans_norm$new(dataset = w_nov_2025)

set.seed(123)
w_nov_2025_rarefied = tmp$norm(method = "rarefy", sample.size = 40000)

# Define salinity groups from actual PSU values
w_nov_2025_rarefied$sample_table = w_nov_2025_rarefied$sample_table %>%
  mutate(
    Salinity_PSU = as.numeric(as.character(Salinity_PSU)),
    Group = case_when(
      Salinity_PSU >= 0 & Salinity_PSU < 0.5 ~ group_order[1],
      Salinity_PSU >= 0.5 & Salinity_PSU < 5 ~ group_order[2],
      Salinity_PSU >= 5 & Salinity_PSU < 18 ~ group_order[3],
      Salinity_PSU >= 18 & Salinity_PSU <= 30 ~ group_order[4],
      TRUE ~ NA_character_
    ),
    Group = factor(Group, levels = group_order),
    Location = factor(Location, levels = names(location_colors))
  )

w_nov_2025_rarefied$tidy_dataset()

# Check missing group classifications
w_nov_2025_rarefied$sample_table %>% filter(is.na(Group) | is.na(Salinity_PSU))

# Examine sample distribution
sample_distribution = w_nov_2025_rarefied$sample_table %>%
  as_tibble(rownames = "SampleID") %>%
  count(Location, Group, .drop = FALSE)

print(sample_distribution)
write_csv(sample_distribution, "results/tables/alpha_beta/Sample_distribution.csv")


# Combined
w_combined = w_nov_2025_rarefied$clone(deep = TRUE)

# Biloxi
w_biloxi = w_nov_2025_rarefied$clone(deep = TRUE)
w_biloxi$sample_table = w_biloxi$sample_table %>% filter(Location == "Biloxi Bay")
w_biloxi$tidy_dataset()

# Pascagoula
w_pascagoula = w_nov_2025_rarefied$clone(deep = TRUE)
w_pascagoula$sample_table = w_pascagoula$sample_table %>% filter(Location == "Pascagoula Bay")
w_pascagoula$tidy_dataset()

# Verify sample counts
sapply(list(Combined = w_combined, Biloxi = w_biloxi, Pascagoula = w_pascagoula),
       function(x) nrow(x$sample_table))

datasets = list(Combined = w_combined, Biloxi = w_biloxi, Pascagoula = w_pascagoula)

for (nm in names(datasets)) {
  datasets[[nm]]$cal_alphadiv()
  datasets[[nm]]$cal_betadiv()
}

# Confirm calculations
names(w_combined$alpha_diversity)
names(w_biloxi$alpha_diversity)
names(w_pascagoula$alpha_diversity)


alpha_combined = trans_alpha$new(dataset = w_combined, group = "Location")
alpha_biloxi = trans_alpha$new(dataset = w_biloxi, group = "Group")
alpha_pascagoula = trans_alpha$new(dataset = w_pascagoula, group = "Group")

df_alpha = alpha_combined$data_alpha %>%
  filter(Measure %in% c("Observed", "Shannon", "Chao1")) %>%
  mutate(Group = factor(Group, levels = group_order),
         Location = factor(Location, levels = names(location_colors)),
         Salinity_PSU = as.numeric(as.character(Salinity_PSU)),
         Measure = factor(Measure, levels = c("Observed", "Shannon", "Chao1")))

write_csv(df_alpha, "results/tables/alpha_beta/Alpha_diversity_all.csv")

############################################################################################
###### ALPHA DIVERSITY — OBSERVED RICHNESS ONLY
############################################################################################

# Filter Observed richness
df_observed = df_alpha %>% filter(Measure == "Observed")

# 1. Combined: Biloxi vs Pascagoula
alpha_location = ggplot(df_observed, aes(Location, Value, fill = Location)) +
  geom_boxplot(width = 0.55, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.10, size = 2.2, alpha = 0.8) +
  scale_fill_manual(values = location_colors) +
  labs(x = "Bay", y = "Observed ASV richness") +
  theme_classic(base_size = 13)  +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA)
  )

alpha_location

ggsave("results/figures/alpha_diversity/Observed_combined_bays.pdf", alpha_location, width = 6, height = 5, bg = "white")


# 2. Salinity groups within each bay
alpha_salinity = df_observed %>%
  filter(!is.na(Group)) %>%
  ggplot(aes(Group, Value, fill = Group)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.10, size = 2.5, alpha = 0.85) +
  facet_wrap(~Location, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = group_colors, breaks = group_order, drop = TRUE) +
  scale_x_discrete(drop = TRUE) +
  labs(x = "Salinity regime", y = "Observed ASV richness", fill = "Salinity regime") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size =12, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA)
  )

alpha_salinity

ggsave("results/figures/alpha_diversity/Observed_salinity_within_bays.pdf", 
       alpha_salinity, width = 10, height = 5, bg = "white")


# 3. Continuous salinity vs Observed richness — faceted by bay
# Calculate regression equation, R2, and P for each bay
reg_stats = df_observed %>%
  group_by(Location) %>%
  group_modify(~{
    fit = lm(Value ~ Salinity_PSU, data = .x)
    s = summary(fit)
    tibble(
      intercept = unname(coef(fit)[1]),
      slope = unname(coef(fit)[2]),
      R2 = s$r.squared,
      p = s$coefficients[2, 4]
    )
  }) %>%
  mutate(
    equation = paste0("y = ", round(intercept, 1),
                      if_else(slope >= 0, " + ", " - "),
                      round(abs(slope), 2), "x"),
    r2_label = paste0("R² = ", sprintf("%.2f", R2)),
    p_label = paste0("P ", if_else(p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", p)))),
    label = paste(equation, r2_label, p_label, sep = "\n")
  )

# Plot
alpha_regression = ggplot(df_observed, aes(Salinity_PSU, Value, color = Location)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9) +
  geom_text(data = reg_stats, aes(x = -Inf, y = Inf, label = label, color = Location),
            inherit.aes = FALSE, hjust = -0.12, vjust = 1.2,
            fontface = "bold", size = 3.8, lineheight = 1.3) +
  facet_wrap(~Location, nrow = 1, scales = "fixed") +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_color_manual(values = location_colors) +
  labs(x = "Salinity (PSU)", y = "Observed ASV richness") +
  theme_bw(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA)
  )

alpha_regression


alpha_regression

ggsave("results/figures/alpha_diversity/Observed_salinity_regression.pdf", alpha_regression, width = 7, height = 5, bg = "white")

# 4. Separate regression plots with equation, R², and P-value
plot_bay_observed = function(bay) {
  df_bay = df_observed %>% filter(Location == bay)
  stats_bay = reg_stats %>% filter(Location == bay)
  
  ggplot(df_bay, aes(Salinity_PSU, Value)) +
    geom_point(color = location_colors[[bay]], size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                color = location_colors[[bay]], linewidth = 0.9) +
    geom_text(data = stats_bay, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.12, vjust = 1.2,
              color = location_colors[[bay]], fontface = "bold",
              size = 4.2, lineheight = 1.4) +
    scale_x_continuous(breaks = seq(0, 30, by = 5)) +
    labs(title = bay, x = "Salinity (PSU)", y = "Observed ASV richness") +
    theme_bw(base_size = 13) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      axis.text = element_text(color = "black", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

alpha_biloxi_plot = plot_bay_observed("Biloxi Bay")
alpha_pascagoula_plot = plot_bay_observed("Pascagoula Bay")

alpha_biloxi_plot
alpha_pascagoula_plot

ggsave("results/figures/alpha_diversity/Observed_Biloxi.pdf", alpha_biloxi_plot, width = 6, height = 5, bg = "white")
ggsave("results/figures/alpha_diversity/Observed_Pascagoula.pdf", alpha_pascagoula_plot, width = 6, height = 5, bg = "white")

##----------------------------------------------------------------------------------------------------------------

beta_combined = trans_beta$new(dataset = w_combined, group = "Group", measure = "bray")
beta_biloxi = trans_beta$new(dataset = w_biloxi, group = "Group", measure = "bray")
beta_pascagoula = trans_beta$new(dataset = w_pascagoula, group = "Group", measure = "bray")

set.seed(123)
beta_combined$cal_ordination(method = "NMDS")

set.seed(123)
beta_biloxi$cal_ordination(method = "NMDS")

set.seed(123)
beta_pascagoula$cal_ordination(method = "NMDS")

theme_nmds = theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )


nmds_combined = beta_combined$plot_ordination(
  plot_color = "Group",
  plot_shape = "Location",
  plot_group_order = group_order,
  plot_type = "point",
  point_size = 3.2,
  point_alpha = 0.9,
  NMDS_stress_text_prefix = "Stress = "
) +
  scale_color_manual(values = group_colors, breaks = group_order, drop = TRUE) +
  scale_shape_manual(values = c("Biloxi Bay" = 16, "Pascagoula Bay" = 17)) +
  labs(title = "Ordination between two estuaries",
       x = "NMDS1", y = "NMDS2", color = "Salinity group", shape = "Bay")+
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.3),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

nmds_combined

ggsave("results/figures/beta_diversity/NMDS_combined.pdf", nmds_combined, width = 7.5, height = 5.5, bg = "white")

# Extract stress value from the NMDS result
# Try this first:
stress_biloxi = 0.012
stress_pascagoula = 0.011

# If the above does not work, inspect:
# str(beta_biloxi$res_ordination, max.level = 1)

nmds_biloxi = beta_biloxi$plot_ordination(
  plot_color = "Group",
  plot_group_order = group_order,
  plot_type = c("point", "centroid"),
  point_size = 3.5,
  point_alpha = 0.9,
  centroid_segment_alpha = 0.5,
  centroid_segment_size = 0.6,
  centroid_segment_linetype = "solid"
) +
  scale_color_manual(values = group_colors, breaks = group_order, drop = TRUE) +
  labs(title = "Biloxi Bay", x = "NMDS1", y = "NMDS2", color = "Salinity group") +
  annotate("text",
           x = -Inf, y = Inf,
           label = paste0("Stress = ", stress_biloxi),
           hjust = -0.1, vjust = 1.3,
           size = 4.2, fontface = "bold") +
  theme_bw(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.3),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )


nmds_biloxi = nmds_biloxi +
  scale_x_continuous(expand = expansion(mult = 0.30)) +
  scale_y_continuous(expand = expansion(mult = 0.50))

nmds_biloxi

ggsave("results/figures/beta_diversity/NMDS_Biloxi.pdf", nmds_biloxi, width = 7, height = 5.5, dpi = 1000)

nmds_pascagoula = beta_pascagoula$plot_ordination(
  plot_color = "Group",
  plot_group_order = group_order,
  plot_type = c("point", "centroid"),
  point_size = 3.5,
  point_alpha = 0.9,
  centroid_segment_alpha = 0.5,
  centroid_segment_size = 0.6,
  centroid_segment_linetype = "solid"
) +
  scale_color_manual(values = group_colors, breaks = group_order, drop = TRUE) +
  labs(title = "Pascagoula Bay", x = "NMDS1", y = "NMDS2", color = "Salinity group") +
  annotate("text",
           x = -Inf, y = Inf,
           label = paste0("Stress = ", stress_pascagoula),
           hjust = -0.1, vjust = 1.3,
           size = 4.2, fontface = "bold") +
  theme_bw(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.3),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )


nmds_pascagoula = nmds_pascagoula +
  scale_x_continuous(expand = expansion(mult = 0.30)) +
  scale_y_continuous(expand = expansion(mult = 0.50))

nmds_pascagoula

ggsave("results/figures/beta_diversity/NMDS_Pascagoula.pdf", nmds_pascagoula, width = 7, height = 5.5, dpi = 1000)

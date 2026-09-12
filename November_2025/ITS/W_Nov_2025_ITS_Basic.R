# =================================================================================================
# November 2025 water fungal ITS analysis with microeco
# Adapted from the November 2025 16S workflow
# Input files:
#   Data files/w_nov_2025_its.biom
#   Data files/w_nov25_its_repseqs.fasta
#   Data files/w_nov25_its_rooted-tree.nwk
#   Data files/w_nov25_its_unrooted-tree.nwk
#   Data files/metadata.txt
#   Data files/nov_2025_env_data.xlsx
# =================================================================================================

# ---- Core data handling / plotting ----
library(tidyverse)
library(magrittr)
library(readxl)
library(RColorBrewer)
library(patchwork)
library(aplot)

# ---- Microbiome data processing ----
library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(ape)
library(picante)
library(Biostrings)

# ---- microeco plotting / environmental associations ----
library(ggtree)
library(ggnested)
library(linkET)

# ---- Null-model / classifier dependencies ----
library(NST)
library(Boruta)
library(randomForest)
library(caret)
library(multiROC)
library(rfPermute)

# Install missing packages only once if needed; do not reinstall them every analysis run.
# remotes::install_github("gmteunisse/ggnested")
# remotes::install_github("WandeRum/multiROC")

set.seed(123)

getwd()

out_base <- "November_2025/ITS"
dir.create(file.path(out_base, "Diversity_Metrics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "RDA-Analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "Abundance_results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "Output figures"), recursive = TRUE, showWarnings = FALSE)

location_colors <- c("#7B3294", "#C2A5CF")

salinity_colors <- c("#FDAE61", "#F46D43", "#A50026")


##### Importing files ########

biom = import_biom("Data files/w_nov_2025_its.biom")

metadata = import_qiime_sample_data("Data files/metadata.txt")

tree = read_tree("Data files/w_nov25_its_rooted-tree.nwk")

rep_fasta = readDNAStringSet("Data files/w_nov25_its_repseqs.fasta", format = "fasta")

sample_names(biom) <- sub("_[A-Z]{2}$", "", sample_names(biom))

biom_ids <- sample_names(biom)
metadata <- metadata[biom_ids, , drop = FALSE]
stopifnot(identical(sample_names(biom), sample_names(metadata)))

w_biom = merge_phyloseq(biom, rep_fasta, metadata, tree)


#rename columns in taxonomy table
tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if(ncol(tax_table(w_biom)) != length(tax_ranks)) stop("Expected 7 ITS taxonomy columns, but found ", ncol(tax_table(w_biom)), ". Check the BIOM taxonomy table.")
colnames(tax_table(w_biom)) <- tax_ranks

w_nov_2025_its <- phyloseq2meco(w_biom)

# removing zero-abundance taxa and make sure sample and taxonomy tables are consistent.
w_nov_2025_its$tidy_dataset()

# check sequencing depth 

sample_depth <- colSums(w_nov_2025_its$otu_table)

summary(sample_depth)
min(sample_depth)
quantile(sample_depth, probs = c(0.05, 0.10, 0.25, 0.50, 0.75))

hist(sample_depth, breaks = 30,
     main = "Sequencing depth per sample",
     xlab = "Reads per sample")

tmp <- trans_norm$new(dataset = w_nov_2025_its)

# Rarefaction depth retained from the 16S workflow.
# Review sample_depth above before running; change this value if ITS depth requires a different cutoff.
rarefaction_depth <- 40000
low_depth_samples <- names(sample_depth[sample_depth < rarefaction_depth])
if(length(low_depth_samples) > 0) message(length(low_depth_samples), " ITS sample(s) have fewer than ", rarefaction_depth, " reads and microeco will remove them during rarefaction: ", paste(low_depth_samples, collapse = ", "))

w_nov_2025_its_rarefied <- w_nov_2025_its

w_nov_2025_its_rarefied$sample_table$Salinity <- factor(
  w_nov_2025_its_rarefied$sample_table$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)

w_nov_2025_its_rarefied$tidy_dataset()

# add_data is used to add the environmental data

env <- readxl::read_excel ("Data files/nov_2025_env_data.xlsx") %>% as.data.frame()

rownames(env) <- env[, 1]

env = env[ ,-1]


env_nov25_its <- trans_env$new(dataset = w_nov_2025_its_rarefied, add_data = env)

w_nov_2025_its_rarefied$cal_abund()
w_nov_2025_its_rarefied$cal_alphadiv()
w_nov_2025_its_rarefied$cal_betadiv()

############################################################################################
######  Alpha diversity
############################################################################################


loc_alpha <- trans_alpha$new(dataset = w_nov_2025_its_rarefied, group = "Location")
loc_alpha$cal_diff(method = "anova", formula = "Location+Salinity")

loc_alpha

loc_alpha <- trans_alpha$new(dataset = w_nov_2025_its_rarefied, group = "Location")
loc_alpha$cal_diff(method = "t.test")
loc_alpha$plot_alpha(measure = "Shannon")

alpha_dat <- w_nov_2025_its_rarefied$alpha_diversity %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::left_join(
    w_nov_2025_its_rarefied$sample_table %>% as.data.frame() %>% tibble::rownames_to_column("Sample"),
    by = "Sample"
  )

write.csv(alpha_dat, file = "Diversity_Metrics/Alpha.csv", row.names = FALSE)

alpha_dat$Salinity <- factor(alpha_dat$Salinity,
                             levels = c("Freshwater", "Moderate Salinity", "High Salinity"))

w_nov_2025_its_rarefied$sample_table$Salinity <- factor(
  w_nov_2025_its_rarefied$sample_table$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)

alpha_shannon <- ggplot(alpha_dat, 
                                  aes(x = factor(Salinity), 
                                      y = Shannon, 
                                      fill = factor(Salinity))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  scale_fill_manual(values = salinity_colors)+
  facet_wrap(~ Location, scales = "free_y") +
  labs(x = "Salinity",
       y = "Shannon",
       fill = "Salinity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

alpha_shannon_loc <- ggplot(alpha_dat, 
                        aes(x = Location, 
                            y = Shannon, 
                            fill = factor(Location))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  scale_fill_manual(values = location_colors)+
  labs(x = "Location",
       y = "Shannon",
       fill = "Location") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

alpha_chao1 <- ggplot(alpha_dat, 
                        aes(x = factor(Salinity), 
                            y = Chao1, 
                            fill = factor(Salinity))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  scale_fill_manual(values = salinity_colors)+
  facet_wrap(~ Location, scales = "free_y") +
  labs(x = "Salinity",
       y = "Chao1",
       fill = "Salinity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

alpha_chao1_loc <- ggplot(alpha_dat, 
                      aes(x = Location, 
                          y = Chao1, 
                          fill = factor(Location))) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2)  +
  scale_fill_manual(values = location_colors)+
  labs(x = "Location",
       y = "Chao1",
       fill = "Location") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

alpha_chao1_loc

ggsave("Diversity_Metrics/alpha_shannon.pdf", plot = alpha_shannon, width = 6, height = 4, dpi = 1000)
ggsave("Diversity_Metrics/alpha_chao1.pdf", plot = alpha_chao1, width = 6, height = 4, dpi = 1000)
ggsave("Diversity_Metrics/alpha_chao1_loc.pdf", plot = alpha_chao1_loc, width = 6, height = 4, dpi = 1000)
ggsave("Diversity_Metrics/alpha_shannon_loc.pdf", plot = alpha_shannon_loc, width = 6, height = 4, dpi = 1000)

alpha_chao1_loc

############################################################################################
######  Beta diversity
############################################################################################

beta_loc <- trans_beta$new(dataset = w_nov_2025_its_rarefied, group = "Location", measure = "bray")

beta_loc$cal_ordination(method = "NMDS")

nov_2025_loc = beta_loc$plot_ordination(plot_color = "Salinity", plot_shape = "Location", plot_type = c("point", "ellipse"), point_size = 2, color_values = salinity_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), # Increase x-axis text size
        axis.text.y = element_text(size = 15), # Increase y-axis text size
        axis.title.x = element_text(size = 15), # Increase x-axis label size
        axis.title.y = element_text(size = 15), # Increase y-axis label size
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

nov_2025_loc

nov_2025_loc_elips <- beta_loc$plot_ordination(plot_color = "Location", plot_shape = "Salinity", plot_type = "point", point_size = 2, color_values = location_colors) +
  stat_ellipse(
    aes(group = Location),
    type = "t",
    linewidth = 1,
    linetype = 1
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

nov_2025_loc_elips

ggsave("Diversity_Metrics/nov_2025_loc.pdf", plot = nov_2025_loc, device = "pdf", width = 7, 
       height = 5, units = "in", dpi = 1000)

ggsave("Diversity_Metrics/nov_2025_loc_elips.pdf", plot = nov_2025_loc_elips, device = "pdf", width = 7, 
       height = 5, units = "in", dpi = 1000)

manova_loc <- trans_beta$new(dataset = w_nov_2025_its_rarefied, group = "Location", measure = "bray")
manova_sal <- trans_beta$new(dataset = w_nov_2025_its_rarefied, group = "Salinity", measure = "bray")
# manova for all groups when manova_all = TRUE
manova_loc$cal_manova(manova_all = TRUE)
manova_sal$cal_manova(manova_all = TRUE)

manova_loc$res_manova
manova_sal$res_manova

# manova for specified group set: such as "Group + Type"
manova <- trans_beta$new(dataset = w_nov_2025_its_rarefied, group = "Location", measure = "bray")
manova$cal_manova(manova_set = "Location + Salinity")
manova$res_manova

write.csv(manova$res_manova, "Diversity_Metrics/permanova_all.csv")
write.csv(manova_loc$res_manova, "Diversity_Metrics/permanova_location.csv")
write.csv(manova_sal$res_manova, "Diversity_Metrics/permanova_salinity.csv")

############################################################################################
######  Water variables
############################################################################################

env_data <- readxl::read_excel("Data files/nov_2025_env_data.xlsx")
env_nov25_its <- trans_env$new(dataset = w_nov_2025_its_rarefied, add_data = env[,1:6])

env_plot <- env_data %>% mutate(
    Location = case_when( str_detect(`Sample ID`, "_BIL_") ~ "Biloxi Bay", str_detect(`Sample ID`, "_PAS_") ~ "Pascagoula Bay"),
    Station = as.numeric(str_extract(`Sample ID`, "(?<=_ST)\\d+"))) %>%
  # Environmental measurements are identical for _01 and _02
  # so retain one observation per station
  distinct(Location, Station, .keep_all = TRUE)

env_long <- env_plot %>% dplyr::select(Location, Station, `Salinity (PSU)`, `Temperature (°C)`, pH, `ODO (mg/L)`) %>%
  tidyr::pivot_longer(cols = -c(Location, Station), names_to = "Variable", values_to = "Value") %>%
  dplyr::mutate(Variable = factor(Variable, levels = c("Salinity (PSU)", "Temperature (°C)", "pH", "ODO (mg/L)"),
                                  labels = c("Salinity (PSU)", "Temperature (°C)", "pH", "Dissolved oxygen (mg/L)")))

env_var_plot = ggplot(env_long, aes(x = Station, y = Value, color = Location, group = Location)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 1:10, labels = paste0("ST", 1:10)) +
  scale_color_manual(values = location_colors) +
  labs(x = "Station", y = NULL, color = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

ggsave("RDA-Analysis/env_varible_plot.pdf", env_var_plot, w = 10, h = 7, dpi = 1000)

# use bray-curtis distance for dbRDA
env_nov25_its$cal_ordination(method = "dbRDA", use_measure = "bray")

# show the orginal results
env_nov25_its$trans_ordination()
env_nov25_its$plot_ordination(plot_color = "Location")

# the main results of RDA are related with the projection and angles between arrows
# adjust the length of the arrows to show them better
env_nov25_its$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 0.5)

# t1$res_rda_trans is the transformed result for plotting
dbrda_plot= env_nov25_its$plot_ordination(plot_color = "Location", plot_shape = "Salinity")

ggsave("RDA-Analysis/db-rda_plot.pdf", dbrda_plot, w= 7, height = 7, dpi = 1000)

env_nov25_its$cal_ordination_anova()
env_nov25_its$cal_ordination_envfit()

env_nov25_its$res_ordination_envfit
write.csv(env_nov25_its$res_ordination_envfit, "RDA-Analysis/dbRDA_envfit.csv", row.names = FALSE)

# use Genus
env_nov25_its$cal_ordination(method = "RDA", taxa_level = "Genus")
# select 10 features and adjust the arrow length
env_nov25_its$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
# t1$res_rda_trans is the transformed result for plot
rda_plot = env_nov25_its$plot_ordination(plot_color = "Location", plot_shape = "Salinity")

ggsave("RDA-Analysis/rda_plot.pdf", rda_plot, w= 7, height = 7, dpi = 1000)


## Mantel tests

env_nov25_its$cal_mantel(use_measure = "bray")

# return t1$res_mantel
head(env_nov25_its$res_mantel)
write.csv(env_nov25_its$res_mantel, "November_2025/ITS/RDA-Analysis/mantel_overall.csv", row.names = FALSE)

env_use <- env_nov25_its$data_env %>% dplyr::select(`Salinity (PSU)`, `Temperature (°C)`, pH, `ODO (mg/L)`)

# Overall trans_env object for environmental correlation matrix
t_env <- trans_env$new(dataset = w_nov_2025_its_rarefied, add_data = env_use, standardize = TRUE)

# Identify top 10 abundant phyla
top_phyla <- tibble(
  ASV = rownames(w_nov_2025_its_rarefied$tax_table),
  Phylum = w_nov_2025_its_rarefied$tax_table$Phylum,
  Abundance = rowSums(w_nov_2025_its_rarefied$otu_table[rownames(w_nov_2025_its_rarefied$tax_table), , drop = FALSE])
) %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  dplyr::filter(!is.na(Phylum), Phylum != "", Phylum != "p__") %>%
  dplyr::arrange(desc(Abundance))

phyla <- head(top_phyla$Phylum, 10)
phyla

# Mantel test separately for each phylum
plot_table <- lapply(phyla, function(ph){
  d <- clone(w_nov_2025_its_rarefied); taxa <- rownames(d$tax_table)[d$tax_table$Phylum == ph]
  if(length(taxa) == 0) return(NULL)
  d$tax_table <- d$tax_table[taxa, , drop = FALSE]; d$otu_table <- d$otu_table[taxa, , drop = FALSE]; d$phylo_tree <- NULL
  d$tidy_dataset(); d$cal_betadiv()
  te <- trans_env$new(dataset = d, add_data = env_use, standardize = TRUE); te$cal_mantel(use_measure = "bray", partial_mantel = TRUE)
  te$res_mantel %>%
    dplyr::transmute(
      spec = gsub("^p__", "", ph),
      env = Variables,
      r = `Correlation coefficient`,
      p = p.adjusted
    )
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
                spec = factor(spec, levels = rev(gsub("^p__", "", phyla))))

# Check Mantel results
plot_table
table(plot_table$spec)
write.csv(plot_table, "November_2025/ITS/RDA-Analysis/mantel_top10_phyla.csv", row.names = FALSE)

# Mantel + environmental Pearson correlation plot
g1 <- qcorrplot(correlate(t_env$data_env), type = "upper", diag = FALSE) +
  geom_square() +
  geom_mark(sig_thres = 0.05, sig_level = c(0.05, 0.01, 0.001), size = 4) +
  geom_couple(aes(colour = pd, size = rd), data = plot_table, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-1, 1)) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#BDBDBD")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 2), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

g1
ggsave("November_2025/ITS/RDA-Analysis/mantel_phylum_env.pdf", plot = g1, width = 10, height = 7, dpi = 1000)

############################################################################################
######  Abundance 
############################################################################################

w_nov_2025_its$cal_abund()
colSums(w_nov_2025_its$taxa_abund$Phylum) %>% summary()
rowMeans(w_nov_2025_its$taxa_abund$Phylum) %>% sort(decreasing = TRUE) %>% head(20)
rownames(w_nov_2025_its$taxa_abund$Phylum)[str_detect(rownames(w_nov_2025_its$taxa_abund$Phylum), regex("unclassified|unknown|NA", ignore_case = TRUE))]

w_nov_2025_its$taxa_abund$Phylum["k__Fungi|p__", ] %>% sort(decreasing = TRUE)


top15_palette <- c(
  "#7F3C8D", "#11A579", "#3969AC", "#F2B701", "#E73F74",
  "#80BA5A", "#E68310", "#008695", "#CF1C90", "#f97b72",
  "#4B4B8F", "#A5AA99", "#8C564B", "#BCBD22", "#17BECF",
  "#C49A6C", "#6B8E23", "#D2691E", "#9C6ADE", "#5F9EA0"
)

Abundance_L <- trans_abund$new(dataset = w_nov_2025_its, taxrank = "Phylum", ntaxa = 15)


box_abundance_loc = Abundance_L$plot_box(group = "Location",  xtext_angle = 30)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

box_abundance_sal = Abundance_L$plot_box(group = "Salinity",  xtext_angle = 30)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

ggsave("Abundance_results/box_abundance_loc.pdf", plot = box_abundance_loc, device = "pdf", width = 9, 
       height = 7, units = "in", dpi = 1000)

ggsave("Abundance_results/box_abundance_sal.pdf", plot = box_abundance_sal, device = "pdf", width = 9, 
       height = 7, units = "in", dpi = 1000)

stack_loc <- trans_abund$new(dataset = w_nov_2025_its_rarefied, taxrank = "Phylum", ntaxa = 15, groupmean = "Location")

stack_sal <- trans_abund$new(dataset = w_nov_2025_its_rarefied, taxrank = "Phylum", ntaxa = 15, groupmean = "Salinity")

stack_loc_plot <- stack_loc$plot_bar(others_color = "grey70", legend_text_italic = FALSE, color_values = top15_palette)  +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border



stack_sal_plot <- stack_sal$plot_bar(others_color = "grey70", legend_text_italic = FALSE, color_values = top15_palette)  +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

ggsave("Abundance_results/stack_loc_plot.pdf", plot = stack_loc_plot, device = "pdf", width = 8, 
       height = 7, units = "in", dpi = 1000)

ggsave("Abundance_results/stack_sal_plot.pdf", plot = stack_sal_plot, device = "pdf", width = 9, 
       height = 7, units = "in", dpi = 1000)



make_stamp_plot <- function(dataset,
                            tax_rank = "Phylum",
                            group_var = "Location",
                            groups = NULL,
                            top_n = 10) {
  
  otu <- as.data.frame(dataset$otu_table)
  tax <- as.data.frame(dataset$tax_table)
  meta <- as.data.frame(dataset$sample_table)
  
  meta$SampleID <- rownames(meta)
  meta[[group_var]] <- trimws(as.character(meta[[group_var]]))
  
  if (!tax_rank %in% colnames(tax)) {
    stop("Taxonomy rank not found. Available columns: ",
         paste(colnames(tax), collapse = ", "))
  }
  
  if (!group_var %in% colnames(meta)) {
    stop("Grouping variable not found. Available metadata columns: ",
         paste(colnames(meta), collapse = ", "))
  }
  
  if (is.null(groups)) {
    groups <- unique(meta[[group_var]])
  }
  
  groups <- trimws(as.character(groups))
  
  if (length(groups) != 2) {
    stop("This plot requires exactly two groups. Current groups are: ",
         paste(groups, collapse = ", "))
  }
  
  meta_sub <- meta[meta[[group_var]] %in% groups, , drop = FALSE]
  
  cat("\nSamples per location:\n")
  print(table(meta_sub[[group_var]]))
  
  if (length(unique(meta_sub[[group_var]])) != 2) {
    stop("Only one Location was retained after filtering. Check exact spelling in groups.")
  }
  
  # Make sure OTU table has taxa as rows and samples as columns
  if (all(meta_sub$SampleID %in% colnames(otu))) {
    otu_sub <- otu[, meta_sub$SampleID, drop = FALSE]
  } else if (all(meta_sub$SampleID %in% rownames(otu))) {
    otu <- as.data.frame(t(otu))
    otu_sub <- otu[, meta_sub$SampleID, drop = FALSE]
  } else {
    stop("Sample IDs do not match between otu_table and sample_table.")
  }
  
  # Relative abundance percentage
  otu_rel <- sweep(otu_sub, 2, colSums(otu_sub), FUN = "/") * 100
  
  # Match taxonomy
  tax_sub <- tax[rownames(otu_rel), , drop = FALSE]
  
  taxon_vec <- as.character(tax_sub[[tax_rank]])
  taxon_vec[is.na(taxon_vec) | taxon_vec == "" | taxon_vec == " "] <- "Unclassified"
  
  # Aggregate ASVs/OTUs to selected taxonomic rank
  tax_abund <- rowsum(as.matrix(otu_rel), group = taxon_vec)
  
  # Create long table using base R
  long_df <- expand.grid(
    Taxon = rownames(tax_abund),
    SampleID = colnames(tax_abund),
    stringsAsFactors = FALSE
  )
  
  long_df$Abundance <- as.vector(tax_abund)
  
  long_df <- merge(
    long_df,
    meta_sub[, c("SampleID", group_var)],
    by = "SampleID",
    all.x = TRUE
  )
  
  colnames(long_df)[colnames(long_df) == group_var] <- "Group"
  long_df$Group <- factor(long_df$Group, levels = groups)
  
  long_df <- long_df[!is.na(long_df$Group), ]
  
  cat("\nLong table columns:\n")
  print(colnames(long_df))
  
  cat("\nSamples in final long table:\n")
  print(table(long_df$Group))
  
  # Select top taxa
  mean_taxa <- aggregate(
    Abundance ~ Taxon,
    data = long_df,
    FUN = mean
  )
  
  mean_taxa <- mean_taxa[order(mean_taxa$Abundance, decreasing = TRUE), ]
  top_taxa <- head(mean_taxa$Taxon, top_n)
  
  long_df <- long_df[long_df$Taxon %in% top_taxa, ]
  
  # Mean abundance by group
  mean_df <- aggregate(
    Abundance ~ Taxon + Group,
    data = long_df,
    FUN = mean
  )
  
  colnames(mean_df)[colnames(mean_df) == "Abundance"] <- "mean_abund"
  
  # Statistics
  stat_list <- lapply(unique(long_df$Taxon), function(tx) {
    
    dat <- long_df[long_df$Taxon == tx, ]
    
    g1 <- dat$Abundance[dat$Group == groups[1]]
    g2 <- dat$Abundance[dat$Group == groups[2]]
    
    pval <- wilcox.test(g1, g2)$p.value
    
    data.frame(
      Taxon = tx,
      mean_group1 = mean(g1, na.rm = TRUE),
      mean_group2 = mean(g2, na.rm = TRUE),
      diff = mean(g1, na.rm = TRUE) - mean(g2, na.rm = TRUE),
      p_value = pval
    )
  })
  
  stat_df <- do.call(rbind, stat_list)
  
  stat_df$sig <- ifelse(stat_df$p_value < 0.001, "***",
                        ifelse(stat_df$p_value < 0.01, "**",
                               ifelse(stat_df$p_value < 0.05, "*", "")))
  
  # Order taxa
  tax_order <- aggregate(
    mean_abund ~ Taxon,
    data = mean_df,
    FUN = sum
  )
  
  tax_order <- tax_order[order(tax_order$mean_abund), "Taxon"]
  
  mean_df$Taxon <- factor(mean_df$Taxon, levels = tax_order)
  stat_df$Taxon <- factor(stat_df$Taxon, levels = tax_order)
  
  # Plot A: mean relative abundance
  p1 <- ggplot(mean_df, aes(x = mean_abund, y = Taxon, fill = Group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = location_colors)+
    labs(x = "Proportions (%)", y = NULL, fill = "Location") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, colour = "black"),
      axis.title.x = element_text(size = 14),
      legend.position = "top",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)
    )
  
  # Plot B: difference
  p2 <- ggplot(stat_df, aes(x = diff, y = Taxon)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2.8) +
    geom_text(aes(label = sig), nudge_x = 0.3, color = "red", size = 5) +
    labs(
      x = paste0("Difference between proportions (%)"),
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 14, colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)
    )
  
  # Plot C: p-values
  p3 <- ggplot(stat_df, aes(y = Taxon, x = 1, label = signif(p_value, 3))) +
    geom_text(size = 3.5) +
    labs(x = NULL, y = NULL, title = "P-value") +
    theme_void() +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  final_plot <- p1 + p2 + p3 +
    plot_layout(widths = c(1.5, 1.1, 0.45))
  
  return(list(
    plot = final_plot,
    stats = stat_df,
    abundance_table = mean_df,
    long_table = long_df
  ))
}

unique(w_nov_2025_its_rarefied$sample_table$Location)
table(w_nov_2025_its_rarefied$sample_table$Location)

phylum_location <- make_stamp_plot( w_nov_2025_its_rarefied,
  tax_rank = "Phylum",
  group_var = "Location",
  groups = c("Biloxi Bay", "Pascagoula Bay"),
  top_n = 15
)

loc_stat_abund <- phylum_location$plot
write.csv(phylum_location$stats, "November_2025/ITS/Abundance_results/loc_stamp_stats.csv", row.names = FALSE)
write.csv(phylum_location$abundance_table, "November_2025/ITS/Abundance_results/loc_stamp_abundance.csv", row.names = FALSE)

ggsave("November_2025/ITS/Abundance_results/loc_stat_abund.pdf", plot = loc_stat_abund, device = "pdf", width = 9, 
       height = 7, units = "in", dpi = 1000)


# show 40 taxa at Genus level
w_nov_2025_its_rarefied$sample_table$Salinity <- factor(
  w_nov_2025_its_rarefied$sample_table$Salinity,
  levels = c(
    "Freshwater",
    "Moderate Salinity",
    "High Salinity"
  )
)

abun_genus <- trans_abund$new(dataset = w_nov_2025_its_rarefied, taxrank = "Genus", ntaxa = 40)
plot_abun_genus_loc <- abun_genus$plot_heatmap(facet = "Location", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
plot_abun_genus_sal <- abun_genus$plot_heatmap(facet = "Salinity", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))

plot_abun_genus_loc
plot_abun_genus_sal
p1 = plot_abun_genus_loc + theme(axis.text.y = element_text(face = 'italic'))

plot_abun_genus_loc_sal <- abun_genus$plot_heatmap(
  xtext_keep = FALSE,
  withmargin = FALSE,
  plot_breaks = c(0.01, 0.1, 1, 10)
) +
  facet_grid(
    ~ Location + Salinity,
    scales = "free_x",
    space = "free_x"
  ) +
  theme(
    axis.text.y = element_text(face = "italic")
  )

plot_abun_genus_loc_sal

ggsave("November_2025/ITS/Abundance_results/plot_abun_genus_loc.pdf", plot = plot_abun_genus_loc, device = "pdf", width = 10, 
       height = 10, units = "in", dpi = 1000)

ggsave("November_2025/ITS/Abundance_results/plot_abund_loc_sal.pdf", plot = plot_abun_genus_loc_sal, device = "pdf", width = 12, 
       height = 10, units = "in", dpi = 1000)


# Top 10 genera, grouped/nested within Phylum
genus_nested <- trans_abund$new(
  dataset = w_nov_2025_its_rarefied,
  taxrank = "Genus",
  ntaxa = 10,
  high_level = "Phylum",
  delete_taxonomy_prefix = FALSE
)

genus_nested_plot <- genus_nested$plot_bar(
  ggnested = TRUE,
  facet = c("Location", "Salinity"),
  xtext_angle = 30
)

ggsave("November_2025/ITS/Abundance_results/genus_nested_location_salinity.pdf", plot = genus_nested_plot,
       width = 12, height = 8, dpi = 1000)


## Differential abundance tests ####

lefse_nov_2025_loc <- trans_diff$new(dataset = w_nov_2025_its_rarefied, method = "lefse", group = "Location", alpha = 0.01, lefse_subgroup = NULL)
lefse_nov_2025_sal <- trans_diff$new(dataset = w_nov_2025_its_rarefied, method = "lefse", group = "Salinity", alpha = 0.05, lefse_subgroup = NULL)

loc_lefse= lefse_nov_2025_loc$plot_diff_bar(threshold = 3.5)
sal_lefse= lefse_nov_2025_sal$plot_diff_bar(threshold = 3)
sal_lefse

ggsave("November_2025/ITS/Output figures/loc_lefse.pdf", plot = loc_lefse, device = "pdf", width =8, 
       height = 10, units = "in", dpi = 1000)

ggsave("November_2025/ITS/Output figures/sal_lefse.pdf", plot = sal_lefse, device = "pdf", width =8, 
       height = 10, units = "in", dpi = 1000)

# LEfSe cladogram for fungal ITS taxa.
# If labels overlap, first inspect the plot and then use select_show_labels with fungal taxa of interest.
lefse_clado_loc <- lefse_nov_2025_loc$plot_diff_cladogram(
  use_taxa_num = 200,
  use_feature_num = 50,
  clade_label_level = 5,
  group_order = c("Biloxi Bay", "Pascagoula Bay")
)

ggsave("November_2025/ITS/Output figures/lefse_clado_loc.pdf", plot = lefse_clado_loc, device = "pdf", width =18, 
       height = 10, units = "in", dpi = 1000)

## Random Forest + Differential abundance test

# use Genus level for parameter taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
rf_loc <- trans_diff$new(dataset = w_nov_2025_its_rarefied, method = "rf", group = "Location", taxa_level = "Genus")
rf_sal <- trans_diff$new(dataset = w_nov_2025_its_rarefied, method = "rf", group = "Salinity", taxa_level = "Genus")

# plot the MeanDecreaseGini bar
# group_order is designed to sort the groups

g1 <- rf_loc$plot_diff_bar(use_number = 1:30, group_order = c("Biloxi Bay", "Pascagoula Bay"))
r1 <- rf_sal$plot_diff_bar(use_number = 1:30, group_order = c("Freshwater", "Moderate Salinity", "High Salinity"))

# plot the abundance using same taxa in g1
g2 <- rf_loc$plot_diff_abund(group_order = c("Biloxi Bay", "Pascagoula Bay"), select_taxa = rf_loc$plot_diff_bar_taxa, plot_type = "barerrorbar", add_sig = F, errorbar_addpoint = FALSE, errorbar_color_black = TRUE)
r2 <- rf_sal$plot_diff_abund(group_order = c("Freshwater", "Moderate Salinity", "High Salinity"), select_taxa = rf_sal$plot_diff_bar_taxa, plot_type = "barerrorbar", add_sig = F, errorbar_addpoint = FALSE, errorbar_color_black = TRUE)

# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())
p <- g1 %>% aplot::insert_right(g2)
p

r1 <- r1 + theme(legend.position = "none")
r2 <- r2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())
q <- r1 %>% aplot::insert_right(r2)
q

ggsave("November_2025/ITS/Output figures/rf_loc.pdf", plot = p, device = "pdf", width =9, 
       height = 10, units = "in", dpi = 1000)


ggsave("November_2025/ITS/Output figures/rf_sal.pdf", plot = q, device = "pdf", width =9, 
       height = 10, units = "in", dpi = 1000)

############ Null Model analysis #########################
# IMPORTANT FOR ITS:
# betaNRI/betaNTI and the betaNTI + RCbray process partition depend on phylogenetic distances.
# ITS can be difficult to align reliably across deep fungal lineages. Interpret these metrics only if
# the supplied ITS tree is biologically defensible for the taxa retained in this dataset.
###------ beta NRI -------#######

# generate trans_nullmodel object
# as an example, we only use high abundance OTU with mean relative abundance > 0.0005
# Set salinity order BEFORE trans_beta$new()
w_nov_2025_its_rarefied$sample_table$Salinity <- factor(
  w_nov_2025_its_rarefied$sample_table$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)

unique(w_nov_2025_its_rarefied$sample_table$Salinity)

nm<- trans_nullmodel$new(w_nov_2025_its_rarefied, filter_thres = 0.0005)

# see null.model parameter for other null models
# null model run 500 times for the example
nm$cal_ses_betampd(runs = 500, abundance.weighted = TRUE)
# return t1$res_ses_betampd

# add betaNRI matrix to beta_diversity list
w_nov_2025_its_rarefied$beta_diversity[["betaNRI"]] <- nm$res_ses_betampd

# create trans_beta class, use measure "betaNRI"
nm_beta_loc <- trans_beta$new(dataset = w_nov_2025_its_rarefied, group = "Location", measure = "betaNRI")
nm_beta_sal <- trans_beta$new(dataset = w_nov_2025_its_rarefied, group = "Salinity", measure = "betaNRI")

# transform the distance for each group
nm_beta_loc$cal_group_distance()
nm_beta_sal$cal_group_distance()

# see the help document for more methods, e.g. "anova" and "KW_dunn"
nm_beta_loc$cal_group_distance_diff(method = "wilcox")

nm_beta_sal$cal_group_distance_diff(method = "wilcox")

# Make sure salinity group is ordered correctly
nm_beta_sal$sample_table$Salinity <- factor(
  nm_beta_sal$sample_table$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)

# plot the results
g1 <- nm_beta_loc$plot_group_distance(add = "mean")
g1 + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)

q1 <- nm_beta_sal$plot_group_distance(group = "Salinity", add = "mean") + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
q1 

# Extract result table
df <- nm_beta_sal$res_group_distance

df$Salinity <- factor(
  df$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)

beta_NRI_sal <- ggplot(df, aes(x = Salinity, y = Value, fill = Salinity)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4) +
  geom_hline(yintercept = -2, linetype = 2) +
  geom_hline(yintercept = 2, linetype = 2) +
  labs(x = "Salinity", y = "betaNRI") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 14)
  )

q1

ggsave("November_2025/ITS/Output figures/betaNRI_sal.pdf", plot = beta_NRI_sal , device = "pdf", width =7, 
       height = 7, units = "in", dpi = 1000)

###------ beta NTI -------#######

# null model run 500 times
nm$cal_ses_betamntd(runs = 500, abundance.weighted = TRUE, null.model = "taxa.labels")
# return t1$res_ses_betamntd

# add betaNTI matrix to beta_diversity list
w_nov_2025_its_rarefied$beta_diversity[["betaNTI"]] <- nm$res_ses_betamntd

nm_NTI_sal <- trans_beta$new(dataset = w_nov_2025_its_rarefied, group = "Salinity", measure = "betaNTI")
nm_NTI_loc <- trans_beta$new(dataset = w_nov_2025_its_rarefied, group = "Location", measure = "betaNTI")

nm_NTI_sal$cal_group_distance()
nm_NTI_loc$cal_group_distance()

nm_NTI_sal$cal_group_distance_diff(method = "wilcox")
nm_NTI_loc$cal_group_distance_diff(method = "wilcox")

# Extract result table
df_nti <- nm_NTI_sal$res_group_distance
df_nti_loc<- nm_NTI_loc$res_group_distance

df_nti$Salinity <- factor(
  df_nti$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)

beta_NTI_sal <- ggplot(df_nti, aes(x = Salinity, y = Value, fill = Salinity)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4) +
  geom_hline(yintercept = -2, linetype = 2) +
  geom_hline(yintercept = 2, linetype = 2) +
  labs(x = "Salinity", y = "betaNTI") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 14)
  )

beta_NTI_loc <- ggplot(df_nti_loc, aes(x = Location, y = Value, fill = Location)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4) +
  geom_hline(yintercept = -2, linetype = 2) +
  geom_hline(yintercept = 2, linetype = 2) +
  labs(x = "Location", y = "betaNTI") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 14)
  )

beta_NTI_sal
beta_NTI_loc

ggsave("November_2025/ITS/Output figures/beta_NTI_sal.pdf", plot = beta_NTI_sal , device = "pdf", width =7, 
       height = 7, units = "in", dpi = 1000)

ggsave("November_2025/ITS/Output figures/beta_NTI_loc.pdf", plot = beta_NTI_loc , device = "pdf", width =7, 
       height = 7, units = "in", dpi = 1000)

## RC Bray

# result stored in t1$res_rcbray
nm$cal_rcbray(runs = 1000)
# return t1$res_rcbray

# use betaNTI and rcbray to evaluate processes
rcbray_loc = nm$cal_process(use_betamntd = TRUE, group = "Location")
rcbray_sal = nm$cal_process(use_betamntd = TRUE, group = "Salinity")

rcbray_loc$res_process
rcbray_sal$res_process


df_loc_rcbray <- rcbray_loc$res_process
df_sal_rcbray <- rcbray_sal$res_process
# order processes
df_loc_rcbray$process <- factor(
  df_loc_rcbray$process,
  levels = c(
    "variable selection",
    "homogeneous selection",
    "dispersal limitation",
    "homogeneous dispersal",
    "drift"
  )
)

# order salinity groups
df_sal_rcbray$Salinity <- factor(
  df_sal_rcbray$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)

# order locations if needed

p_loc <- ggplot(df_loc_rcbray, aes(x = Location, y = percentage, fill = process)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(percentage > 3, paste0(round(percentage, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "variable selection" = "#D55E00",
    "homogeneous selection" = "#0072B2",
    "dispersal limitation" = "#009E73",
    "homogeneous dispersal" = "#CC79A7",
    "drift" = "#F0E442"
  )) +
  labs(
    x = "Location",
    y = "Percentage (%)",
    fill = "Process"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank()
  )

p_sal <- ggplot(df_sal_rcbray, aes(x = Salinity, y = percentage, fill = process)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(percentage > 3, paste0(round(percentage, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "variable selection" = "#D55E00",
    "homogeneous selection" = "#0072B2",
    "dispersal limitation" = "#009E73",
    "homogeneous dispersal" = "#CC79A7",
    "drift" = "#F0E442"
  )) +
  labs(
    x = "Salinity",
    y = "Percentage (%)",
    fill = "Process"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank()
  )

p_loc
p_sal

ggsave("November_2025/ITS/Output figures/rcbray_sal.pdf", plot = p_sal , device = "pdf", width =7, 
       height = 8, units = "in", dpi = 1000)

ggsave("November_2025/ITS/Output figures/rcbray_loc.pdf", plot = p_loc , device = "pdf", width =7, 
       height = 8, units = "in", dpi = 1000)

## NST

# require NST package to be installed
nst_loc = nm$cal_NST(method = "tNST", group = "Location", dist.method = "bray", abundance.weighted = TRUE, output.rand = TRUE, SES = TRUE)
nst_sal = nm$cal_NST(method = "tNST", group = "Salinity", dist.method = "bray", abundance.weighted = TRUE, output.rand = TRUE, SES = TRUE)

nst_loc$res_NST$index.grp
nst_sal$res_NST$index.grp
write.csv(nst_loc$res_NST$index.grp, "November_2025/ITS/Output figures/NST_location.csv", row.names = FALSE)
write.csv(nst_sal$res_NST$index.grp, "November_2025/ITS/Output figures/NST_salinity.csv", row.names = FALSE)

###### Classifier based analysis ############

rf_loc_1 <- trans_classifier$new(dataset = w_nov_2025_its_rarefied, y.response = "Location", x.predictors = "All")

# generate train and test set
rf_loc_1$cal_split(prop.train = 3/4)

# require caret package
rf_loc_1$set_trainControl()

# use default parameter method = "rf"
rf_loc_1$cal_train(method = "rf")

rf_loc_1$cal_predict()
# plot the confusionMatrix to check out the performance
conf_mat_plot = rf_loc_1$plot_confusionMatrix()

ggsave("November_2025/ITS/Output figures/conf_mat_plot.pdf", plot = conf_mat_plot , device = "pdf", width =7, 
       height = 6, units = "in", dpi = 1000)

rf_loc_1$cal_ROC()
# select one group to plot ROC
rf_loc_1$plot_ROC(plot_group = "Biloxi.Bay")
rf_loc_1$plot_ROC(plot_group = "Biloxi.Bay", color_values = "black")
# default all groups
rf_loc_1$plot_ROC(size = 0.5, alpha = 0.7)


# require Boruta package
rf_loc_1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.01)

rf_loc_2 <- trans_classifier$new(dataset = w_nov_2025_its_rarefied, y.response = "Location", x.predictors = "All")

rf_loc_2$cal_split(prop.train = 3/4)

rf_loc_2$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.01)

rf_loc_2$set_trainControl()

rf_loc_2$cal_train()

rf_loc_2$cal_predict()

rf_loc_2$plot_confusionMatrix()

rf_loc_2$cal_ROC()

rf_loc_2$plot_ROC(size = 0.5, alpha = 0.7)


# default method in caret package without significance
rf_loc_2$cal_feature_imp()
rf_loc_2$plot_feature_imp(colour = "red", fill = "red", width = 0.6)

# generate significance with rfPermute package
rf_loc_2$cal_feature_imp(rf_feature_sig = TRUE, num.rep = 1000)

# add_sig = TRUE: add significance label
rf_loc_2$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6, add_sig = TRUE)

# show_sig_group = TRUE: show different colors in groups with different significance labels
rf_loc_2$plot_feature_imp(show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE)

rf_loc_2$plot_feature_imp(show_sig_group = TRUE, coord_flip = TRUE, width = 0.6, add_sig = TRUE)

# rf_sig_show = "MeanDecreaseGini": switch to MeanDecreaseGini
rf_loc_2$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE)

# group_aggre = FALSE: donot aggregate features for each group
rf_loc_2$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)

#### ----------------------------------------------------------------------------------------------------------------------------------------
# beta clustering
# Location
w_nov2025_its_beta <- clone(w_nov_2025_its_rarefied)

w_nov2025_its_beta$sample_table %<>% subset(Location %in% c("Biloxi Bay", "Pascagoula Bay"))

w_nov2025_its_beta$tidy_dataset()

# calculate beta diversity again after subsetting
w_nov2025_its_beta$cal_betadiv(method = "bray", unifrac = FALSE)

loc_beta_cluster <- trans_beta$new(dataset = w_nov2025_its_beta, group = "Location", measure = 'bray')

# use replace_name to set the label name, group parameter used to set the color
loc_clus = loc_beta_cluster$plot_clustering(group = "Location",replace_name = c("Location"))

# Salt

w_nov2025_its_beta <- clone(w_nov_2025_its_rarefied)

w_nov2025_its_beta$sample_table %<>% subset(Salinity %in% c("Freshwater", "Moderate Salinity", "High Salinity"))

w_nov2025_its_beta$tidy_dataset()

# calculate beta diversity again after subsetting
w_nov2025_its_beta$cal_betadiv(method = "bray", unifrac = FALSE)

sal_beta_cluster <- trans_beta$new(dataset = w_nov2025_its_beta, group = "Salinity", measure = 'bray')

# use replace_name to set the label name, group parameter used to set the color
sal_clus = sal_beta_cluster$plot_clustering(group = "Salinity", replace_name = c("Salinity"))

ggsave("November_2025/ITS/Output figures/loc_clus.pdf", loc_clus , width =7, height = 8, units = "in", dpi = 1000)
ggsave("November_2025/ITS/Output figures/sal_clus.pdf", sal_clus , width =7, height = 8, units = "in", dpi = 1000)


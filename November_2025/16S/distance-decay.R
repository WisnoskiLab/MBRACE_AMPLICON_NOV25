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

# ---- Null-model / classifier dependencies ----
library(tidyverse) 
library(vegan) 
library(geosphere) 
library(mgcv)
library(ecodist)
library(patchwork)
library(ape)
library(NST)
library(spaa)

set.seed(123)

out_base <- "November_2025/16S"

dir.create(file.path(out_base, "distance-decay"), recursive = TRUE, showWarnings = FALSE)

location_colors <- c("#0072B2","#E69F00")

salinity_colors <- c( "#009E73", "#56B4E9", "#0072B2")
subcommunity_colors <- c("Abundant" = "#F8766D", "CRT" = "#00BA38", "Rare" = "#619CFF")

##### Importing files ########

biom = import_biom("Data files/w_nov_2025_table.biom")

metadata = import_qiime_sample_data("Data files/metadata.txt")

tree = read_tree("Data files/w_nov_2025-unrooted_tree.nwk")

rep_fasta = readDNAStringSet("Data files/w_nov_2025_repseqs.fasta", format = "fasta")

w_biom = merge_phyloseq(biom, rep_fasta, metadata, tree)

#rename columns in taxonomy table
colnames(tax_table(w_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

w_nov_2025 <- phyloseq2meco(w_biom)

# removing zero-abundance taxa and make sure sample and taxonomy tables are consistent.
w_nov_2025$tidy_dataset()

# check sequencing depth 

sample_depth <- colSums(w_nov_2025$otu_table)

summary(sample_depth)
min(sample_depth)
quantile(sample_depth, probs = c(0.05, 0.10, 0.25, 0.50, 0.75))

hist(sample_depth, breaks = 30,
     main = "Sequencing depth per sample",
     xlab = "Reads per sample")

tmp <- trans_norm$new(dataset = w_nov_2025)

# rarefaction
w_nov_2025_rarefied <- tmp$norm(method = "rarefy", sample.size = 40000)

w_nov_2025_rarefied$sample_table$Salinity <- factor(
  w_nov_2025_rarefied$sample_table$Salinity,
  levels = c("Freshwater", "Moderate Salinity", "High Salinity")
)
# Preparing the data
w_nov_2025_rarefied$tidy_dataset()

otu_raw <- as.matrix(w_nov_2025$otu_table)
meta <- w_nov_2025$sample_table %>% as.data.frame() %>% rownames_to_column("SampleID")
tax <- w_nov_2025$tax_table %>% as.data.frame() %>% rownames_to_column("ASV")

# Safety check for OTU-table orientation
if(!all(meta$SampleID %in% colnames(otu_raw)) && all(meta$SampleID %in% rownames(otu_raw))) otu_raw <- t(otu_raw)
stopifnot(all(meta$SampleID %in% colnames(otu_raw)))

# Derive bay/station directly from sample IDs
meta <- meta %>% mutate(Location = case_when(str_detect(SampleID, "_BIL_") ~ "Biloxi Bay", str_detect(SampleID, "_PAS_") ~ "Pascagoula Bay", TRUE ~ as.character(Location)),
                        Station = as.integer(str_match(SampleID, "_ST(\\d+)_")[,2]),
                        StationID = paste0(if_else(Location == "Biloxi Bay", "BIL", "PAS"), "_ST", Station))

# Check sequencing depths
tibble(SampleID = colnames(otu_raw), Reads = colSums(otu_raw)) %>% summary()

# ============================================================
# 2. REMOVE ABSOLUTE SINGLETONS
# Paper removed ASVs with total abundance = 1
# ============================================================

singleton_summary <- tibble(ASV = rownames(otu_raw), Total_reads = rowSums(otu_raw), Occurrence_n = rowSums(otu_raw > 0)) %>%
  summarise(ASVs_before = n(), Singletons = sum(Total_reads == 1), Doubletons = sum(Total_reads == 2),
            Reads_1_5 = sum(Total_reads <= 5), Reads_1_10 = sum(Total_reads <= 10), Occurs_once = sum(Occurrence_n == 1))

singleton_summary

otu <- otu_raw[rowSums(otu_raw) > 1, , drop = FALSE]
dim(otu)

# ============================================================
# 3. CONVERT FILTERED COUNTS TO RELATIVE ABUNDANCE
# ============================================================

rel <- sweep(otu, 2, colSums(otu), "/")

# Confirm each sample sums to 1
summary(colSums(rel))

# ============================================================
# 4. CLASSIFY ABUNDANT, CRT AND RARE ASVs
#
# Abundant = >=1% in at least one sample
# CRT      = >=0.1% in at least one sample, but always <1%
# Rare     = always <0.1%
#
# Occurrence:
# Broad        >=75%
# Intermediate >10% and <75%
# Narrow       <=10%
# ============================================================

asv_class <- tibble(ASV = rownames(rel), Max_RA = apply(rel, 1, max), Mean_RA = rowMeans(rel),
                    Occurrence = rowMeans(rel > 0), Occurrence_n = rowSums(rel > 0)) %>%
  mutate(Subcommunity = case_when(Max_RA >= 0.01 ~ "Abundant", Max_RA >= 0.001 ~ "CRT", TRUE ~ "Rare"),
         Distribution = case_when(Occurrence >= 0.75 ~ "Broad", Occurrence > 0.10 ~ "Intermediate", TRUE ~ "Narrow"),
         Subcommunity = factor(Subcommunity, levels = c("Abundant", "CRT", "Rare")),
         Distribution = factor(Distribution, levels = c("Broad", "Intermediate", "Narrow")))

# Number and percentage of ASVs
class_summary <- asv_class %>% count(Subcommunity, name = "ASVs") %>% mutate(Percent_ASVs = 100 * ASVs / sum(ASVs))
class_summary

# Relative contribution of each subcommunity to total sequence abundance
read_summary <- asv_class %>% mutate(Reads = rowSums(otu)[match(ASV, rownames(otu))]) %>%
  group_by(Subcommunity) %>% summarise(ASVs = n(), Reads = sum(Reads), .groups = "drop") %>%
  mutate(Relative_reads = 100 * Reads / sum(Reads))
read_summary

# ============================================================
# 5. CHECK LOW-ABUNDANCE / OCCURRENCE STRUCTURE AFTER FILTERING
# ============================================================

rare_check <- asv_class %>% mutate(Total_reads = rowSums(otu)[match(ASV, rownames(otu))])

rare_check_summary <- rare_check %>% group_by(Subcommunity) %>%
  summarise(ASVs = n(), Median_reads = median(Total_reads), Q1_reads = quantile(Total_reads, 0.25),
            Q3_reads = quantile(Total_reads, 0.75), Median_samples = median(Occurrence_n),
            One_sample = sum(Occurrence_n == 1), One_sample_pct = 100 * mean(Occurrence_n == 1),
            Reads_1_5 = sum(Total_reads <= 5), Reads_1_5_pct = 100 * mean(Total_reads <= 5),
            Reads_1_10 = sum(Total_reads <= 10), Reads_1_10_pct = 100 * mean(Total_reads <= 10))
rare_check_summary

# ============================================================
# 6. RANK-ABUNDANCE CURVE
# ============================================================

rank_ab <- asv_class %>% arrange(desc(Max_RA)) %>% mutate(Rank = row_number())

p_rank <- ggplot(rank_ab, aes(Rank, Max_RA, color = Subcommunity)) +
  geom_point(size = 1.2, alpha = 0.75) +
  geom_hline(yintercept = c(0.001, 0.01), linetype = 2) +
  scale_y_log10() +
  scale_color_manual(values = subcommunity_colors) +
  labs(x = "ASV rank", y = "Maximum relative abundance", color = NULL) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA))

p_rank

# ============================================================
# 7. OCCURRENCE STRUCTURE WITHIN EACH SUBCOMMUNITY
# ============================================================

occ_summary <- asv_class %>% count(Subcommunity, Distribution) %>% group_by(Subcommunity) %>%
  mutate(Percent_ASVs = 100 * n / sum(n))

p_occ <- ggplot(occ_summary, aes(Subcommunity, Percent_ASVs, fill = Distribution)) +
  geom_col(width = 0.7) +
  labs(x = NULL, y = "ASVs (%)", fill = "Occurrence") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA))

p_occ

# Relative read contribution by abundance class × occurrence class
occ_reads <- asv_class %>% mutate(Reads = rowSums(otu)[match(ASV, rownames(otu))]) %>%
  group_by(Subcommunity, Distribution) %>% summarise(ASVs = n(), Reads = sum(Reads), .groups = "drop") %>%
  mutate(Relative_reads = 100 * Reads / sum(Reads))

occ_reads

p_occ_reads <- ggplot(occ_reads, aes(Subcommunity, Relative_reads, fill = Distribution)) +
  geom_col(width = 0.7) +
  labs(x = NULL, y = "Relative sequence abundance (%)", fill = "Occurrence") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA))

p_occ_reads

# ============================================================
# 8. CONVERT SAMPLE-LEVEL RELATIVE ABUNDANCE TO LONG FORMAT
# ============================================================

ra_long <- rel %>% as.data.frame(check.names = FALSE) %>% rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "SampleID", values_to = "RA") %>%
  left_join(meta %>% select(SampleID, Location, Station, StationID), by = "SampleID") %>%
  left_join(asv_class %>% select(ASV, Subcommunity, Distribution), by = "ASV")

# ============================================================
# 9. AVERAGE BIOLOGICAL REPLICATES AT EACH STATION
#
# Classification above used individual samples.
# Spatial distance-decay below uses one community per station.
# ============================================================

station_ra <- ra_long %>% group_by(Location, Station, StationID, ASV, Subcommunity, Distribution) %>%
  summarise(RA = mean(RA), .groups = "drop")

# Relative contribution of abundant / CRT / rare communities along transects
station_subcommunity <- station_ra %>% group_by(Location, Station, Subcommunity) %>%
  summarise(Relative_abundance = sum(RA), .groups = "drop")

p_subcommunity <- ggplot(station_subcommunity, aes(factor(Station), Relative_abundance * 100, fill = Subcommunity)) +
  geom_col() +
  facet_wrap(~Location, scales = "free_x") +
  scale_fill_manual(values = subcommunity_colors) +
  labs(x = "Station", y = "Relative abundance (%)", fill = NULL) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA))

p_subcommunity

# ============================================================
# 10. ENVIRONMENTAL DATA
#
# Assumes your existing env object has:
# Salinity_PSU, Temperature_C, pH, ODO_mg_L,
# Conductivity_uS_cm, TDS_mg_L
#
# For multivariate environmental distance we intentionally use
# salinity, temperature, pH and ODO because conductivity/TDS
# are strongly redundant with salinity.
# ============================================================

env <- readxl::read_excel ("Data files/nov_2025_env_data copy.xlsx") %>% as.data.frame()

rownames(env) <- env[, 1]

env = env[ ,-1]

env_df <- as.data.frame(env)
if(!"SampleID" %in% names(env_df)) env_df <- env_df %>% rownames_to_column("SampleID")

env_vars <- c("Salinity_PSU", "Temperature_C", "pH", "ODO_mg_L", "Conductivity_µS_cm", "TDS_mg_L")
stopifnot(all(env_vars %in% names(env_df)))

env_station <- env_df %>% left_join(meta %>% select(SampleID, Location, Station, StationID), by = "SampleID") %>%
  group_by(Location, Station, StationID) %>%
  summarise(across(all_of(env_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")

env_station

# ============================================================
# 11. STATION COORDINATES
#
# Required file columns:
# Location, Station, Latitude, Longitude
#
# Optional preferred column:
# Transect_km
#
# Transect_km = distance along the estuarine transect from a
# common origin. If present, the analysis uses it instead of
# straight-line geographic distance.
# ============================================================

station_coords <- read_csv("Data files/coordinates.csv", show_col_types = FALSE) %>%
  mutate(Location = case_when(str_detect(Sample_ID, "_BIL_") ~ "Biloxi Bay",
                              str_detect(Sample_ID, "_PAS_") ~ "Pascagoula Bay",
                              TRUE ~ Location),
         Station = as.integer(str_match(Sample_ID, "_ST(\\d+)_")[,2]),
         StationID = paste0(if_else(Location == "Biloxi Bay", "BIL", "PAS"), "_ST", Station)) %>%
  select(Location, Station, StationID, Latitude, Longitude) %>%
  distinct()

station_coords

# ============================================================
# 12. FUNCTION TO CREATE COMMUNITY MATRIX
# ============================================================

make_comm <- function(location, subcommunity = "Total", rank = "ASV"){
  d <- station_ra %>% filter(Location == location)
  if(subcommunity != "Total") d <- d %>% filter(Subcommunity == subcommunity)
  
  if(rank == "ASV"){
    d <- d %>% transmute(StationID, Taxon = ASV, RA)
  } else {
    d <- d %>% left_join(tax %>% select(ASV, Taxon = all_of(rank)), by = "ASV") %>%
      filter(!is.na(Taxon), Taxon != "", !str_detect(Taxon, "__$|unassigned$|Unassigned$")) %>%
      group_by(StationID, Taxon) %>% summarise(RA = sum(RA), .groups = "drop")
  }
  
  x <- d %>% group_by(StationID, Taxon) %>% summarise(RA = sum(RA), .groups = "drop") %>%
    pivot_wider(names_from = Taxon, values_from = RA, values_fill = 0)
  
  comm <- x %>% column_to_rownames("StationID") %>% as.matrix()
  comm[, colSums(comm) > 0, drop = FALSE]
}

# Example
make_comm("Biloxi Bay", "Abundant")[1:5, 1:min(5, ncol(make_comm("Biloxi Bay", "Abundant")))]

# ============================================================
# 13. FUNCTION: BRAY + GEOGRAPHIC + ENVIRONMENTAL DISTANCE
# ============================================================

get_spatial <- function(location, subcommunity = "Total", rank = "ASV"){
  comm <- make_comm(location, subcommunity, rank)
  ids <- rownames(comm)
  
  crd <- station_coords %>% filter(Location == location, StationID %in% ids) %>% arrange(match(StationID, ids))
  ev <- env_station %>% filter(Location == location, StationID %in% ids) %>% arrange(match(StationID, ids))
  
  if(nrow(crd) != length(ids)) stop(paste("Missing coordinates for", location))
  if(nrow(ev) != length(ids)) stop(paste("Missing environmental data for", location))
  if(anyNA(ev[, env_vars])) stop(paste("Missing environmental values in", location))
  
  bray <- vegdist(comm, method = "bray")
  
  if("Transect_km" %in% names(crd)){
    geo <- dist(matrix(crd$Transect_km, ncol = 1))
  } else {
    geo <- as.dist(distm(as.matrix(crd[, c("Longitude", "Latitude")]), fun = distHaversine) / 1000)
  }
  
  env_mat <- ev %>% select(all_of(env_vars))
  keep_env <- map_lgl(env_mat, ~sd(.x, na.rm = TRUE) > 0)
  env_dist <- dist(scale(env_mat[, keep_env, drop = FALSE]))
  sal_dist <- dist(matrix(ev$Salinity_PSU, ncol = 1))
  
  list(bray = bray, geo = geo, env = env_dist, salinity = sal_dist, ids = ids)
}

# ============================================================
# 14. CONVERT DISTANCE MATRICES TO PAIRWISE DATA FOR PLOTS
# ============================================================

pair_data <- function(location, subcommunity = "Total", rank = "ASV"){
  x <- get_spatial(location, subcommunity, rank)
  b <- as.matrix(x$bray); g <- as.matrix(x$geo); e <- as.matrix(x$env); s <- as.matrix(x$salinity)
  z <- which(lower.tri(b), arr.ind = TRUE)
  
  tibble(Location = location, Subcommunity = subcommunity, TaxRank = rank,
         Site1 = rownames(b)[z[,1]], Site2 = colnames(b)[z[,2]],
         Bray = b[z], Geo_km = g[z], Env_distance = e[z], Salinity_difference = s[z])
}

locations <- sort(unique(meta$Location[meta$Location %in% c("Biloxi Bay", "Pascagoula Bay")]))
subcommunities <- c("Total", "Abundant", "CRT", "Rare")

dd_pairs <- crossing(Location = locations, Subcommunity = subcommunities) %>%
  pmap_dfr(~pair_data(..1, ..2))

dd_pairs %>% group_by(Location, Subcommunity) %>%
  summarise(Pairs = n(), Min_km = min(Geo_km), Max_km = max(Geo_km), Mean_Bray = mean(Bray), .groups = "drop")

# ============================================================
# 15. COMPARE OVERALL BRAY-CURTIS DISSIMILARITY
# Analogous to Fig. 2A of reference paper
# ============================================================

p_bray <- ggplot(dd_pairs, aes(Subcommunity, Bray, fill = Subcommunity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  facet_wrap(~Location) +
  scale_fill_manual(values = c("Total" = "grey60", subcommunity_colors)) +
  labs(x = NULL, y = "Bray–Curtis dissimilarity") +
  theme_minimal() +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA))

p_bray

# ============================================================
# 16. MAIN SPATIAL DISTANCE-DECAY FIGURE
# Analogous to time-decay Fig. 2B, but distance in km
# ============================================================

p_dd <- ggplot(dd_pairs, aes(Geo_km, Bray, color = Location)) +
  geom_point(alpha = 0.45, size = 1.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), method.args = list(method = "REML"), se = TRUE, linewidth = 1) +
  facet_wrap(~Subcommunity, ncol = 2) +
  scale_color_manual(values = location_colors) +
  labs(x = "Geographic distance (km)", y = "Bray–Curtis dissimilarity", color = NULL) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(size = 11), axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))

p_dd

# ============================================================
# 17. DESCRIPTIVE SPEARMAN + LINEAR DECAY SLOPES
#
# Spearman parallels the reference paper.
# LM slope is useful for comparing decay strength.
# Inferential conclusions should rely more heavily on Mantel/MRM.
# ============================================================

dd_stats <- dd_pairs %>% group_by(Location, Subcommunity) %>%
  group_modify(~{
    ct <- cor.test(.x$Geo_km, .x$Bray, method = "spearman", exact = FALSE)
    lm1 <- lm(Bray ~ Geo_km, data = .x)
    tibble(Spearman_rho = unname(ct$estimate), Spearman_p = ct$p.value,
           Slope = unname(coef(lm1)[2]), Intercept = unname(coef(lm1)[1]), R2 = summary(lm1)$r.squared)
  }) %>% ungroup() %>% mutate(Spearman_q = p.adjust(Spearman_p, method = "BH"))

dd_stats

# ============================================================
# 18. MANTEL TESTS
#
# Geo_r       = spatial distance-decay
# Env_r       = community turnover with environmental distance
# Salinity_r  = community turnover with salinity difference
# Geo_partial = geographic effect after controlling environment
# ============================================================

run_mantel <- function(location, subcommunity){
  x <- get_spatial(location, subcommunity)
  mg <- vegan::mantel(x$bray, x$geo, method = "spearman", permutations = 9999)
  me <- vegan::mantel(x$bray, x$env, method = "spearman", permutations = 9999)
  ms <- vegan::mantel(x$bray, x$salinity, method = "spearman", permutations = 9999)
  mp <- vegan::mantel.partial(x$bray, x$geo, x$env, method = "spearman", permutations = 9999)
  
  tibble(Location = location, Subcommunity = subcommunity,
         Geo_r = unname(mg$statistic), Geo_p = mg$signif,
         Env_r = unname(me$statistic), Env_p = me$signif,
         Salinity_r = unname(ms$statistic), Salinity_p = ms$signif,
         Geo_partial_r = unname(mp$statistic), Geo_partial_p = mp$signif)
}

mantel_results <- crossing(Location = locations, Subcommunity = subcommunities) %>%
  pmap_dfr(~run_mantel(..1, ..2)) %>%
  mutate(Geo_q = p.adjust(Geo_p, "BH"), Env_q = p.adjust(Env_p, "BH"),
         Salinity_q = p.adjust(Salinity_p, "BH"), Geo_partial_q = p.adjust(Geo_partial_p, "BH"))

mantel_results

# ============================================================
# 19. ENVIRONMENTAL DISTANCE VS COMMUNITY DISSIMILARITY
# Analogous to environmental heterogeneity analysis in paper
# ============================================================

p_env <- ggplot(dd_pairs, aes(Env_distance, Bray, color = Location)) +
  geom_point(alpha = 0.45, size = 1.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  facet_wrap(~Subcommunity, ncol = 2) +
  scale_color_manual(values = location_colors) +
  labs(x = "Environmental distance", y = "Bray–Curtis dissimilarity", color = NULL) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA))

p_env

# ============================================================
# 20. SALINITY DIFFERENCE VS COMMUNITY DISSIMILARITY
# Particularly important for your estuarine gradient
# ============================================================

p_salinity <- ggplot(dd_pairs, aes(Salinity_difference, Bray, color = Location)) +
  geom_point(alpha = 0.45, size = 1.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), method.args = list(method = "REML"), se = TRUE, linewidth = 1) +
  facet_wrap(~Subcommunity, ncol = 2) +
  scale_color_manual(values = location_colors) +
  labs(x = expression(Delta*" salinity (PSU)"), y = "Bray–Curtis dissimilarity", color = NULL) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA))

p_salinity

# ============================================================
# 21. MULTIPLE REGRESSION ON DISTANCE MATRICES (MRM)
#
# Model A: geographic + multivariate environmental distance
# Model B: geographic + salinity difference
#
# MRM is especially useful because geographic and environmental
# distance covary strongly along an estuary.
# ============================================================

run_mrm_env <- function(location, subcommunity){
  x <- get_spatial(location, subcommunity)
  bray <- x$bray; geo <- x$geo; envd <- x$env
  ecodist::MRM(bray ~ geo + envd, nperm = 9999, mrank = TRUE)
}

run_mrm_salinity <- function(location, subcommunity){
  x <- get_spatial(location, subcommunity)
  bray <- x$bray; geo <- x$geo; sal <- x$salinity
  ecodist::MRM(bray ~ geo + sal, nperm = 9999, mrank = TRUE)
}

mrm_env_models <- crossing(Location = locations, Subcommunity = subcommunities) %>%
  mutate(Model = map2(Location, Subcommunity, run_mrm_env))

mrm_salinity_models <- crossing(Location = locations, Subcommunity = subcommunities) %>%
  mutate(Model = map2(Location, Subcommunity, run_mrm_salinity))

# Extract coefficient tables
mrm_env_coef <- mrm_env_models %>% mutate(Result = map(Model, ~as.data.frame(.x$coef) %>% rownames_to_column("Term"))) %>%
  select(-Model) %>% unnest(Result)

mrm_salinity_coef <- mrm_salinity_models %>% mutate(Result = map(Model, ~as.data.frame(.x$coef) %>% rownames_to_column("Term"))) %>%
  select(-Model) %>% unnest(Result)

mrm_env_coef
mrm_salinity_coef

# View model R2 and permutation P
mrm_env_models %>% transmute(Location, Subcommunity, R2 = map(Model, "r.squared")) %>% unnest_wider(R2)
mrm_salinity_models %>% transmute(Location, Subcommunity, R2 = map(Model, "r.squared")) %>% unnest_wider(R2)

# ============================================================
# 22. DISTANCE-DECAY ACROSS TAXONOMIC RESOLUTION
#
# Spatial analogue of reference-paper Fig. 4
# ASV -> Genus -> Family -> Order -> Class -> Phylum
# ============================================================

ranks <- c("ASV", c("Genus", "Family", "Order", "Class", "Phylum")[c("Genus", "Family", "Order", "Class", "Phylum") %in% names(tax)])

rank_stats <- crossing(Location = locations, Subcommunity = c("Abundant", "CRT", "Rare"), TaxRank = ranks) %>%
  pmap_dfr(function(Location, Subcommunity, TaxRank){
    x <- get_spatial(Location, Subcommunity, TaxRank)
    d <- tibble(Bray = as.vector(x$bray), Geo_km = as.vector(x$geo))
    lm1 <- lm(Bray ~ Geo_km, data = d)
    ct <- cor.test(d$Geo_km, d$Bray, method = "spearman", exact = FALSE)
    mt <- vegan::mantel(x$bray, x$geo, method = "spearman", permutations = 9999)
    
    tibble(Location = Location, Subcommunity = Subcommunity, TaxRank = TaxRank,
           Slope = unname(coef(lm1)[2]), R2 = summary(lm1)$r.squared,
           Spearman_rho = unname(ct$estimate), Spearman_p = ct$p.value,
           Mantel_r = unname(mt$statistic), Mantel_p = mt$signif)
  }) %>%
  mutate(TaxRank = factor(TaxRank, levels = ranks),
         Mantel_q = p.adjust(Mantel_p, method = "BH"))

rank_stats

p_rank_decay <- ggplot(rank_stats, aes(TaxRank, Slope, color = Location, group = Location)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = Mantel_q < 0.05), size = 3) +
  facet_wrap(~Subcommunity) +
  scale_color_manual(values = location_colors) +
  labs(x = "Taxonomic resolution", y = "Distance–decay slope", color = NULL, shape = "FDR < 0.05") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA))

p_rank_decay

# ============================================================
# 23. NICHE WIDTH
# Spatial analogue of reference-paper niche-width analysis
# ============================================================

calc_niche <- function(location, subcommunity){
  d <- station_ra %>% filter(Location == location, Subcommunity == subcommunity) %>%
    select(StationID, ASV, RA) %>% pivot_wider(names_from = ASV, values_from = RA, values_fill = 0)
  
  mat <- d %>% column_to_rownames("StationID") %>% as.matrix()
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  nw <- spaa::niche.width(mat, method = "levins")
  
  tibble(ASV = names(nw), Niche_width = as.numeric(nw), Location = location, Subcommunity = subcommunity)
}

niche_width <- crossing(Location = locations, Subcommunity = c("Abundant", "CRT", "Rare")) %>%
  pmap_dfr(~calc_niche(..1, ..2))

p_niche <- ggplot(niche_width, aes(Subcommunity, log10(Niche_width), fill = Subcommunity)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~Location) +
  scale_fill_manual(values = subcommunity_colors) +
  labs(x = NULL, y = expression(log[10]*" Levins' niche width")) +
  theme_minimal() +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA))

p_niche

# Statistical comparison of niche width
niche_tests <- niche_width %>% group_by(Location) %>%
  group_modify(~{
    kw <- kruskal.test(Niche_width ~ Subcommunity, data = .x)
    tibble(Kruskal_chisq = unname(kw$statistic), df = unname(kw$parameter), P = kw$p.value)
  })

niche_tests

# Pairwise Wilcoxon if Kruskal-Wallis is significant
niche_pairwise <- niche_width %>% group_split(Location) %>%
  set_names(unique(niche_width$Location)) %>%
  map(~pairwise.wilcox.test(.x$Niche_width, .x$Subcommunity, p.adjust.method = "BH"))

niche_pairwise

# ============================================================
# 24. OPTIONAL: pNST FOR ASSEMBLY PROCESSES
#
# Uses individual samples rather than station-averaged data.
# Runs each bay separately so each estuary is evaluated as its
# own spatial metacommunity.
# ============================================================

tree <- w_nov_2025$phylo_tree

run_pnst <- function(location, subcommunity){
  samples <- meta %>% filter(Location == location) %>% pull(SampleID)
  taxa_keep <- if(subcommunity == "Total") asv_class$ASV else asv_class %>% filter(Subcommunity == subcommunity) %>% pull(ASV)
  
  comm <- otu[intersect(taxa_keep, rownames(otu)), intersect(samples, colnames(otu)), drop = FALSE] %>% t()
  taxa_common <- intersect(colnames(comm), tree$tip.label)
  comm <- comm[, taxa_common, drop = FALSE]
  tr <- ape::keep.tip(tree, taxa_common)
  
  grp <- data.frame(Group = rep(location, nrow(comm)), row.names = rownames(comm))
  
  pd_dir <- tempfile(pattern = paste0("pNST_", gsub(" ", "_", location), "_", subcommunity, "_"))
  dir.create(pd_dir, recursive = TRUE)
  on.exit(unlink(pd_dir, recursive = TRUE, force = TRUE), add = TRUE)
  
  NST::pNST(comm = comm, tree = tr, group = grp, pd.wd = pd_dir,
            abundance.weighted = TRUE, rand = 1000, phylo.shuffle = TRUE,
            nworker = 4, between.group = FALSE)
}

pnst_models <- crossing(Location = locations, Subcommunity = subcommunities) %>%
  mutate(pNST = map2(Location, Subcommunity, run_pnst))

# Inspect group-level pNST
pnst_models %>% mutate(Group_result = map(pNST, "index.grp")) %>% select(-pNST) %>% unnest(Group_result)

# ============================================================
# 25. SAVE MAIN RESULTS
# ============================================================

dir.create("Results/Distance_decay", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures/Distance_decay", recursive = TRUE, showWarnings = FALSE)

write_csv(class_summary, "distance-decay/subcommunity_ASV_summary.csv")
write_csv(read_summary, "distance-decay/subcommunity_read_summary.csv")
write_csv(rare_check_summary, "distance-decay/rare_ASV_QC_summary.csv")
write_csv(dd_pairs, "distance-decay/pairwise_distance_decay_data.csv")
write_csv(dd_stats, "distance-decay/distance_decay_statistics.csv")
write_csv(mantel_results, "distance-decay/mantel_results.csv")
write_csv(rank_stats, "distance-decay/taxonomic_resolution_distance_decay.csv")
write_csv(niche_width, "distance-decay/niche_width.csv")

ggsave("distance-decay/rank_abundance.png", p_rank, width = 8, height = 6, dpi = 600)
ggsave("distance-decay/subcommunity_spatial_abundance.png", p_subcommunity, width = 10, height = 6, dpi = 600)
ggsave("distance-decay/bray_subcommunities.png", p_bray, width = 9, height = 6, dpi = 600)
ggsave("distance-decay/geographic_distance_decay.png", p_dd, width = 10, height = 8, dpi = 600)
ggsave("distance-decay/environmental_distance_decay.png", p_env, width = 10, height = 8, dpi = 600)
ggsave("distance-decay/salinity_distance_decay.png", p_salinity, width = 10, height = 8, dpi = 600)
ggsave("distance-decay/taxonomic_resolution_distance_decay.png", p_rank_decay, width = 10, height = 6, dpi = 600)
ggsave("distance-decay/niche_width.png", p_niche, width = 9, height = 6, dpi = 600)

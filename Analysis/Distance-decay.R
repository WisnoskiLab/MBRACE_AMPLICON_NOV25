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
library(geosphere)

# ---- microeco plotting / environmental associations ----
library(ggtree)
library(ggplot2)
library(ggnested)
library(linkET)

# ---- Null-model / classifier dependencies ----
library(NST)
library(Boruta)
library(randomForest)
library(caret)
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
library(rstatix)

set.seed(123)

location_colors <- c("#0072B2","#E69F00")

subcommunity_colors <- c("Abundant" = "#F8766D", "CRT" = "#00BA38", "Rare" = "#619CFF")

##### Importing files ########

biom = import_biom("Data files 200_230/w_nov_2025_200-230_feture.biom")
metadata = import_qiime_sample_data("Data files 200_230/metadata.txt")
tree = read_tree("Data files 200_230/w_nov_2025_200_230_unrooted-tree.nwk")
rep_fasta = readDNAStringSet("Data files 200_230/dna-sequences_230_200.fasta", format = "fasta")

w_biom = merge_phyloseq(biom, rep_fasta, metadata, tree)

colnames(tax_table(w_biom)) = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

w_nov_2025 = phyloseq2meco(w_biom)
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

# Retain ASVs present in at least 5% of samples
prev_threshold = 0.05
min_samples = max(1, ceiling(ncol(otu) * prev_threshold))

otu = otu[rowSums(otu > 0) >= min_samples, , drop = FALSE]

# Check filtering results
filter_summary = tibble(
  Step = c("Original", "After singleton removal", "After prevalence filtering"),
  ASVs = c(nrow(otu_raw), sum(rowSums(otu_raw) > 1), nrow(otu))
)

filter_summary

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

sal_levels = c("Limnetic (0-0.5)", "Oligohaline (0.5-5)", "Mesohaline (5-18)", "Polyhaline (18-30)")
sub_levels = c("Total", "Abundant", "CRT", "Rare")
sub_colors = c("Total" = "#333333", "Abundant" = "#E64B35", "CRT" = "#00A087", "Rare" = "#3C78D8")
bay_colors = c("Biloxi Bay" = "#0072B2", "Pascagoula Bay" = "#E69F00")
env_vars = c("Salinity_PSU", "Temperature_C", "pH", "ODO_mg_L")

# 04. ASSIGN ABUNDANCE AND OCCURRENCE CLASSES PER ASV --------------------------
# Abundant: >=1% relative abundance in >=1 sample.
# CRT: reaches >=0.1% but NEVER >=1% in any sample.
# Rare: <0.1% in EVERY sample.
# Distribution describes sample prevalence, not species niche width.
asv_class = tibble(ASV = rownames(rel), Max_RA = apply(rel, 1, max), Mean_RA = rowMeans(rel),
                   Occurrence = rowMeans(rel > 0), Occurrence_n = rowSums(rel > 0)) %>%
  mutate(Subcommunity = case_when(Max_RA >= 0.01 ~ "Abundant", Max_RA >= 0.001 ~ "CRT", TRUE ~ "Rare"),
         Distribution = case_when(Occurrence >= 0.75 ~ "Broad", Occurrence > 0.10 ~ "Intermediate", TRUE ~ "Narrow"),
         Subcommunity = factor(Subcommunity, levels = sub_levels[-1]),
         Distribution = factor(Distribution, levels = c("Broad", "Intermediate", "Narrow")))

class_summary = asv_class %>% count(Subcommunity, name = "ASVs") %>% mutate(Percent_ASVs = 100 * ASVs / sum(ASVs))
read_summary = asv_class %>% mutate(Reads = rowSums(otu)[match(ASV, rownames(otu))]) %>%
  group_by(Subcommunity) %>% summarise(ASVs = n(), Reads = sum(Reads), .groups = "drop") %>% mutate(Read_percent = 100 * Reads / sum(Reads))
occ_summary = asv_class %>% count(Subcommunity, Distribution, name = "ASVs") %>% group_by(Subcommunity) %>% mutate(Percent_ASVs = 100 * ASVs / sum(ASVs)) %>% ungroup()

p_rank = asv_class %>% arrange(desc(Max_RA)) %>% mutate(Rank = row_number()) %>%
  ggplot(aes(Rank, Max_RA, color = Subcommunity)) + geom_point(alpha = 0.65, size = 0.7) +
  geom_hline(yintercept = c(.001, .01), linetype = 2) + scale_y_log10() +
  scale_color_manual(values = sub_colors) + labs(x = "ASV rank", y = "Maximum sample relative abundance", color = "Subcommunity") + 
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.6, "cm")
  )

p_rank

ggsave("rank_RT_Abund.pdf", p_rank, w = 6, h = 6, dpi = 1000)


# 05. PREPARE INDIVIDUAL-SAMPLE AND STATION-AVERAGED DATA -----------------------
# Classification uses individual samples, whereas distance-decay analyses
# use the mean ASV relative abundance across replicates at each station.
ra_long = as.data.frame(rel, check.names = FALSE) %>% rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "SampleID", values_to = "RA") %>%
  left_join(meta %>% select(SampleID, Location, Station, StationID, Group), by = "SampleID") %>%
  left_join(asv_class %>% select(ASV, Subcommunity, Distribution), by = "ASV")

station_ra = ra_long %>% group_by(Location, Station, StationID, ASV, Subcommunity, Distribution) %>%
  summarise(RA = mean(RA), .groups = "drop")

stopifnot(all(abs(station_ra %>% group_by(StationID) %>% summarise(SumRA = sum(RA), .groups = "drop") %>% pull(SumRA) - 1) < 1e-8))

# 06. DESCRIPTIVE FIGURES: ASV CLASSES ALONG STATIONS AND SALINITY -------------
# Percentages show contribution to ALL retained reads per station or sample.
station_subcommunity = station_ra %>% group_by(Location, Station, Subcommunity) %>%
  summarise(Relative_abundance = sum(RA), .groups = "drop") %>%
  mutate(Station = factor(paste0("ST", Station), levels = paste0("ST", 1:10)))

p_station = ggplot(station_subcommunity, aes(Station, 100 * Relative_abundance, fill = Subcommunity)) +
  geom_col() + facet_wrap(~Location) + scale_fill_manual(values = sub_colors) +
  labs(x = "Station", y = "Relative contribution to all reads (%)", fill = "Subcommunity") +
  theme_bw(base_size = 12) + 
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.6, "cm")
  )

p_station

ggsave("station-wise-contrb.pdf", p_station, w =10, h = 6, dpi = 1000)

salinity_subcommunity = ra_long %>% group_by(SampleID, Location, Group, Subcommunity) %>%
  summarise(RA = sum(RA), .groups = "drop") %>% group_by(Location, Group, Subcommunity) %>%
  summarise(Mean_RA = mean(RA), .groups = "drop")

sal_levels = c("Limnetic (0-0.5)", "Oligohaline (0.5-5)", "Mesohaline (5-18)", "Polyhaline (18-30)")

salinity_subcommunity = ra_long %>%
  group_by(SampleID, Location, Group, Subcommunity) %>%
  summarise(RA = sum(RA), .groups = "drop") %>%
  group_by(Location, Group, Subcommunity) %>%
  summarise(Mean_RA = mean(RA), .groups = "drop") %>%
  mutate(Group = factor(Group, levels = sal_levels))


salinity_subcommunity

p_salinity_composition = ggplot(salinity_subcommunity, aes(Group, 100 * Mean_RA, fill = Subcommunity)) +
  geom_col() + facet_wrap(~Location, scales = "free_x") + scale_fill_manual(values = sub_colors) +
  labs(x = "Salinity category", y = "Mean relative contribution (%)", fill = "Subcommunity") +
  theme_bw(base_size = 12) + theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.6, "cm")
  )

p_salinity_composition

ggsave("salinity-wise-contribution.pdf", p_salinity_composition, w =7, h = 6, dpi = 1000)

# Station-level phylum composition PER SUBCOMMUNITY. All bars retain original
# whole-community RA, so Abundant, CRT, Rare panels are not forced to 100%.
station_tax = station_ra %>% left_join(tax %>% select(ASV, Phylum), by = "ASV") %>%
  mutate(Phylum = str_remove(as.character(Phylum), "^p__"),
         Phylum = if_else(is.na(Phylum) | Phylum == "" | str_detect(Phylum, "__$"), "Unclassified", Phylum))

top_phyla = station_tax %>% filter(Phylum != "Unclassified") %>% group_by(Phylum) %>%
  summarise(Total_RA = sum(RA), .groups = "drop") %>% slice_max(Total_RA, n = 25, with_ties = FALSE) %>% pull(Phylum)

phylum_plot = station_tax %>% mutate(Phylum = case_when(Phylum == "Unclassified" ~ "Unclassified", Phylum %in% top_phyla ~ Phylum, TRUE ~ "Others")) %>%
  group_by(Location, Station, Subcommunity, Phylum) %>% summarise(RA = sum(RA), .groups = "drop") %>%
  mutate(Station = factor(paste0("ST", Station), levels = paste0("ST", 1:10)))

phylum_plot

p_phylum = ggplot(phylum_plot, aes(Station, 100 * RA, fill = Phylum)) + geom_col() +
  facet_grid(Subcommunity ~ Location, scales = "free_y") +
  labs(x = "Station", y = "Relative contribution to total reads (%)", fill = "Phylum") +
  theme_bw(base_size = 11) +
  theme(
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.6, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "cm")
  )

p_phylum

ggsave("p_phylum.pdf", p_phylum, width = 15, height = 8)


# 07. ENVIRONMENTAL METADATA, ALIGNMENT, AND PREPROCESSING ---------------------
# Use 4 variables consistently in ALL environmental distance calculations.
# Conductivity/TDS excluded here because of redundancy with salinity.
# Modify name aliases only if the workbook uses different column names.

env <- readxl::read_excel ("Data files/nov_2025_env_data copy.xlsx") %>% as.data.frame()

rownames(env) <- env[, 1]

env = env[ ,-1]

env_df <- as.data.frame(env)
if(!"SampleID" %in% names(env_df)) env_df <- env_df %>% rownames_to_column("SampleID")
# Support common alternate spellings / units without silently altering values.

for(target in env_vars){
  if(!target %in% names(env_df)){
    possibilities = switch(target,
                           Salinity_PSU = c("Salinity (PSU)", "Salinity"),
                           Temperature_C = c("Temperature (°C)", "Temperature"),
                           pH = c("PH"),
                           ODO_mg_L = c("ODO (mg/L)", "Dissolved oxygen", "DO_mg_L"))
    matched = intersect(possibilities, names(env_df))
    if(length(matched) == 1) names(env_df)[match(matched, names(env_df))] = target
  }
}

stopifnot(all(env_vars %in% names(env_df)), !anyDuplicated(env_df$SampleID))

env_station = env_df %>% select(SampleID, all_of(env_vars)) %>%
  mutate(across(all_of(env_vars), as.numeric)) %>%
  inner_join(meta %>% select(SampleID, Location, Station, Group, StationID), by = "SampleID") %>%
  group_by(Location, Station, StationID) %>%
  summarise(across(all_of(env_vars), ~if(all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)), .groups = "drop")

stopifnot(!anyNA(env_station[, env_vars]))

# 08. IMPORT NEW CUMULATIVE TRANSECT DISTANCES -------------------------------
# CRITICAL: this replaces the old pairwise Haversine geographic matrix.
# The supplied CSV includes StationID, Location, Segment_km and Transect_km.
# Transect_km is cumulative distance from ST1 using station-ordered segments.
# It does NOT yet represent bends between neighboring stations: update later
# using a map-traced channel polyline if this matters to final inference.
station_coords = read_csv("Data files/station_transect_distances.csv", show_col_types = FALSE) %>%
  select(Location, Station, StationID, Latitude, Longitude, Segment_km, Transect_km) %>%
  distinct() %>% arrange(Location, Station)

stopifnot(!anyDuplicated(station_coords$StationID), !anyNA(station_coords$Transect_km), all(station_coords$Segment_km >= 0))

stopifnot(setequal(unique(meta$StationID), station_coords$StationID))

transect_check = station_coords %>% group_by(Location) %>%
  summarise(N_stations = n(), Transect_length_km = max(Transect_km) - min(Transect_km),
            Sum_segments_km = sum(Segment_km), .groups = "drop")

transect_check

station_coords

# 09. PREPARE ONE COMPLETE STATION x ASV COMMUNITY MATRIX ----------------------
# make_comm() keeps every station within its bay, including stations where a
# particular subset or taxon is absent (a zero row). This prevents mismatches
# across dissimilarity and geographic/environmental matrices.
make_comm = function(location, subcommunity = "Total", rank = "ASV", normalize_subset = FALSE){
  ids = station_coords %>% filter(Location == location) %>% arrange(Station) %>% pull(StationID)
  d = station_ra %>% filter(Location == location)
  if(subcommunity != "Total") d = d %>% filter(Subcommunity == subcommunity)
  if(rank == "ASV") d = d %>% transmute(StationID, Taxon = ASV, RA)
  else d = d %>% left_join(tax %>% select(ASV, Taxon = all_of(rank)), by = "ASV") %>%
    filter(!is.na(Taxon), Taxon != "", !str_detect(Taxon, "__$|unassigned$|Unassigned$")) %>%
    transmute(StationID, Taxon, RA)
  d = d %>% group_by(StationID, Taxon) %>% summarise(RA = sum(RA), .groups = "drop") %>%
    pivot_wider(names_from = Taxon, values_from = RA, values_fill = 0)
  x = d %>% column_to_rownames("StationID") %>% as.matrix()
  x = x[match(ids, rownames(x)), , drop = FALSE]
  rownames(x) = ids
  x[is.na(x)] = 0
  x = x[, colSums(x) > 0, drop = FALSE]
  if(normalize_subset){
    totals = rowSums(x)
    if(any(totals <= 0)) stop("Zero-sum station; cannot normalize this subset: ", location, " / ", subcommunity, " / ", rank)
    x = sweep(x, 1, totals, "/")
  }
  x
}

# Sanity check: 9 Biloxi rows, 10 Pascagoula rows if all stations represented.
print(dim(make_comm("Biloxi Bay", "Abundant")))
print(dim(make_comm("Pascagoula Bay", "Abundant")))


# 10. BRAY + TRANSECT + ENVIRONMENTAL + SALINITY DISTANCE MATRICES -------------
# Pairwise distances must share EXACTLY the same ordered station IDs.
# PRIMARY Bray: normalize each subcommunity to 100% at each station to assess
# compositional differences within Abundant, CRT, and Rare subsets.
# SENSITIVITY Bray: original whole-community RA (different total subset fractions
# can influence dissimilarity). For Total both approaches are identical.
# Environmental distance = Euclidean distance after z-scoring 4 station means.
get_spatial = function(location, subcommunity = "Total", rank = "ASV", mode = bray_mode){
  stopifnot(mode %in% c("within_subset", "original"))
  comm = make_comm(location, subcommunity, rank, normalize_subset = (mode == "within_subset"))
  ids = rownames(comm)
  crd = station_coords %>% filter(Location == location) %>% arrange(match(StationID, ids))
  ev = env_station %>% filter(Location == location) %>% arrange(match(StationID, ids))
  stopifnot(identical(crd$StationID, ids), identical(ev$StationID, ids))
  if(any(rowSums(comm) == 0)) stop("Zero abundance at a station: ", location, " / ", subcommunity)
  bray = vegan::vegdist(comm, method = "bray")
  geo = dist(matrix(crd$Transect_km, ncol = 1))
  env_mat = ev %>% select(all_of(env_vars))
  keep_env = vapply(env_mat, sd, numeric(1), na.rm = TRUE) > 0
  if(!any(keep_env)) stop("No varying environmental covariates")
  envd = dist(scale(as.matrix(env_mat[, keep_env, drop = FALSE])))
  sal = dist(matrix(ev$Salinity_PSU, ncol = 1))
  list(bray = bray, geo = geo, env = envd, salinity = sal, ids = ids)
}

# 11. CONVERT DISTANCE MATRICES INTO ONE ROW PER STATION PAIR -------------------
# For n stations there are n(n-1)/2 UNIQUE unordered pairs. In Biloxi n=9,
# so 36 pairs; in Pascagoula n=10, so 45 pairs (if all stations retained).
# 11. BRAY-CURTIS DISSIMILARITY FOR STATION PAIRS ----------------------------

pair_data <- function(location, subcommunity = "Total", rank = "ASV"){
  x <- get_spatial(location, subcommunity, rank, mode = "original")
  b <- as.matrix(x$bray); g <- as.matrix(x$geo)
  e <- as.matrix(x$env); s <- as.matrix(x$salinity)
  z <- which(lower.tri(b), arr.ind = TRUE)
  
  tibble(Location = location, Subcommunity = subcommunity, TaxRank = rank,
         Site1 = rownames(b)[z[, 1]], Site2 = colnames(b)[z[, 2]],
         Bray = b[z], Geo_km = g[z],
         Env_distance = e[z], Salinity_difference = s[z])
}

locations <- c("Biloxi Bay", "Pascagoula Bay")
subcommunities <- sub_levels

dd_pairs <- crossing(Location = locations, Subcommunity = subcommunities) %>%
  pmap_dfr(~pair_data(..1, ..2))

p_bray <- dd_pairs %>%
  ggplot(aes(x = Subcommunity, y = Bray, fill = Subcommunity)) +
  geom_boxplot(alpha = 0.75, width = 0.6, outlier.shape = NA) +
  facet_wrap(~Location, nrow = 1) +
  labs(x = "Subcommunity", y = "Bray–Curtis dissimilarity") +
  scale_fill_manual(values = sub_colors) +
  theme_bw(base_size = 12) + theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.6, "cm")
  )

p_bray

ggsave("bray-curtis-differences.pdf", p_bray, width = 7, height = 6)


# 12. MAIN SPATIAL DISTANCE-DECAY: CURVES AND BAY COMPARISON --------------------
# Each point is a station PAIR. GAM (k=4) permits a curved association.
# Descriptive only: pairwise distances are NOT independent observations.
# Calculate Mantel r and permutation p-value for each subcommunity and bay
dd_plot = dd_pairs %>%
  mutate(Subcommunity = factor(Subcommunity, levels = sub_levels))

bray_mode = "within_subset"

geo_stats = crossing(Location = locations, Subcommunity = sub_levels) %>%
  pmap_dfr(function(Location, Subcommunity){
    x = get_spatial(Location, Subcommunity, mode = bray_mode)
    mt = vegan::mantel(x$bray, x$geo, method = "spearman", permutations = 9999)
    tibble(Location = Location, Subcommunity = Subcommunity,
           Mantel_r = unname(mt$statistic), Mantel_p = mt$signif)
  }) %>%
  mutate(Subcommunity = factor(Subcommunity, levels = sub_levels),
         Label = paste0("r = ", sprintf("%.3f", Mantel_r),
                        ", p ", if_else(Mantel_p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", Mantel_p)))))

# Display statistics separately for each bay and subcommunity
geo_labels = geo_stats %>%
  group_by(Location) %>%
  arrange(Subcommunity, .by_group = TRUE) %>%
  mutate(x = -Inf, y = seq(0.98, 0.77, length.out = n())) %>%
  ungroup()

# Plot
p_geo = ggplot(dd_plot, aes(Geo_km, Bray, color = Subcommunity)) +
  geom_point(alpha = 0.45, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              method.args = list(method = "REML"), se = FALSE, linewidth = 1.2) +
  geom_text(data = geo_labels, aes(x = x, y = y, label = Label, color = Subcommunity),
            inherit.aes = FALSE, hjust = -0.05, size = 4, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~Location, nrow = 1) +
  scale_color_manual(values = sub_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Cumulative along-transect separation (km)",
       y = "Bray–Curtis dissimilarity", color = "Subcommunity") +
  theme_bw(base_size = 13) + theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.6, "cm")
  )

p_geo

ggsave("geographic_distance_decay_Mantel.pdf", p_geo, width = 9, height = 5)

# 14. ENVIRONMENTAL AND SALINITY DISTANCE-DECAY PLOTS ---------------------------
# Environmental distance combines standardized salinity, temperature, pH, ODO.
# Salinity distance alone is absolute PSU difference. These are descriptive
# linear trend plots; assess actual associations by Mantel/MRM below.

sub_colors = c(
  "Abundant" = "#E64B35",
  "CRT" = "#00A087",
  "Rare" = "#3C78D8",
  "Total" = "#333333"
)

dd_plot = dd_pairs |>
  mutate(
    Location = factor(Location, levels = c("Biloxi Bay", "Pascagoula Bay")),
    Subcommunity = factor(Subcommunity, levels = c("Abundant", "CRT", "Rare", "Total"))
  )


env_labels = dd_plot |>
  group_by(Location, Subcommunity) |>
  group_modify(~{
    fit = lm(Bray ~ Env_distance, data = .x)
    tibble(R2 = summary(fit)$r.squared)
  }) |>
  ungroup() |>
  mutate(
    Label = paste0(as.character(Subcommunity), ": R² = ", sprintf("%.3f", R2)),
    x = -Inf,
    y = c(0.98, 0.91, 0.84, 0.77)[match(Subcommunity, levels(dd_plot$Subcommunity))]
  )

p_env2 = ggplot(dd_plot, aes(Env_distance, Bray, color = Subcommunity)) +
  geom_point(alpha = 0.45, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  geom_text(
    data = env_labels,
    aes(x = x, y = y, label = Label, color = Subcommunity),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~Location, nrow = 1) +
  scale_color_manual(values = sub_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Environmental distance (standardized Euclidean)",
    y = "Bray–Curtis dissimilarity",
    color = "Subcommunity"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.6, "cm")
  )

p_env2

# Calculate linear regression R²
sal_labels = dd_plot |>
  group_by(Location, Subcommunity) |>
  group_modify(~{
    fit = lm(Bray ~ Salinity_difference, data = .x)
    tibble(R2 = summary(fit)$r.squared)
  }) |>
  ungroup() |>
  mutate(
    Label = paste0(as.character(Subcommunity), ": R² = ", sprintf("%.3f", R2)),
    x = -Inf,
    y = c(0.98, 0.91, 0.84, 0.77)[match(Subcommunity, levels(dd_plot$Subcommunity))]
  )

# Salinity distance-decay plot
p_sal2 = ggplot(dd_plot, aes(Salinity_difference, Bray, color = Subcommunity)) +
  geom_point(alpha = 0.45, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  geom_text(data = sal_labels, aes(x = x, y = y, label = Label, color = Subcommunity),
            inherit.aes = FALSE, hjust = -0.05, size = 4, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~Location, nrow = 1) +
  scale_color_manual(values = sub_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Absolute salinity difference (PSU)", y = "Bray–Curtis dissimilarity", color = "Subcommunity") +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.6, "cm")
  )

p_sal2

ggsave("salinity_distance_combined_LM.pdf",
       p_sal2, width = 10, height = 5.5, dpi= 1000)

ggsave("env_distance_combined_LM.pdf",
       p_env2, width = 10, height = 5.5, dpi= 1000)

# ============================================================
# GEOGRAPHIC–ENVIRONMENTAL DISTANCE RELATIONSHIP
# ============================================================

# STEP 1: Extract each unique station pair only once.
# Environmental and geographic distances are identical across
# subcommunities, so use Total to avoid counting pairs four times.

geo_env_pairs = dd_pairs %>%
  filter(Subcommunity == "Total") %>%
  select(Location, Site1, Site2, Geo_km, Env_distance,
         Salinity_difference) %>%
  distinct()

# STEP 2: Fit descriptive regressions separately for each bay.
# R2 measures how much pairwise environmental-distance variation
# is described by geographic separation in an ordinary linear fit.
# These regressions are descriptive, not independent-pair tests.

geo_env_stats = geo_env_pairs %>%
  group_by(Location) %>%
  summarise(Slope = unname(coef(lm(Env_distance ~ Geo_km))[2]),
            R2 = summary(lm(Env_distance ~ Geo_km))$r.squared,
            Pearson_r = cor(Geo_km, Env_distance),
            Salinity_r = cor(Geo_km, Salinity_difference),
            .groups = "drop")

print(geo_env_stats)

write_csv(geo_env_stats, "geography_environment_relationship.csv")

# STEP 3: Plot geographic separation against environmental distance.
# A strong positive relationship indicates that these predictors
# share substantial variation along the estuarine transect.

p_geo_env = ggplot(geo_env_pairs, aes(Geo_km, Env_distance)) +
  geom_point(color = "#0072B2", size = 3, alpha = .7) +
  geom_smooth(method = "lm", color = "#D55E00", se = TRUE, linewidth = 1.1) +
  facet_wrap(~Location, nrow = 1) +
  labs(x = "Cumulative along-transect separation (km)",
       y = "Environmental distance (standardized Euclidean)") +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold"))

p_geo_env

ggsave("geographic_vs_environmental_distance.pdf",
       p_geo_env, width = 9, height = 6, device = cairo_pdf)

# ============================================================
# SPATIALLY DETRENDED ENVIRONMENTAL DISTANCE
# ============================================================

# STEP 1: Within each bay, regress environmental distance on
# geographic separation and retain the regression residuals.
#
# Positive residual = stations are more environmentally different
# than expected from their geographic separation.
#
# Negative residual = stations are more environmentally similar
# than expected from their geographic separation.

geo_env_residuals = geo_env_pairs %>%
  group_by(Location) %>%
  mutate(Env_predicted = fitted(lm(Env_distance ~ Geo_km)),
         Env_residual = resid(lm(Env_distance ~ Geo_km))) %>%
  ungroup()

# STEP 2: Check that residuals are centered near zero and
# uncorrelated with geographic distance under the linear model.
# This property is mathematical, not evidence that all spatial
# autocorrelation has been eliminated.

residual_check = geo_env_residuals %>%
  group_by(Location) %>%
  summarise(Mean_residual = mean(Env_residual),
            Cor_geo_residual = cor(Geo_km, Env_residual),
            .groups = "drop")

print(residual_check)

# STEP 3: Join the residual environmental distances back to
# community dissimilarities for all four abundance classes.

dd_detrended = dd_pairs %>%
  left_join(geo_env_residuals %>%
              select(Location, Site1, Site2, Env_predicted, Env_residual),
            by = c("Location", "Site1", "Site2"))

stopifnot(!anyNA(dd_detrended$Env_residual))

write_csv(dd_detrended, "spatially_detrended_distance_pairs.csv")

dd_detrended

# ============================================================
#  ENVIRONMENTAL DISTANCE-DECAY AFTER SPATIAL DETRENDING
# ============================================================

# STEP 1: Calculate descriptive slopes and R2 values.
# Because these are station pairs, lm p-values are not used.

detrended_stats = dd_detrended %>%
  group_by(Location, Subcommunity) %>%
  summarise(Residual_slope = unname(coef(lm(Bray ~ Env_residual))[2]),
            Residual_R2 = summary(lm(Bray ~ Env_residual))$r.squared,
            Original_slope = unname(coef(lm(Bray ~ Env_distance))[2]),
            Original_R2 = summary(lm(Bray ~ Env_distance))$r.squared,
            .groups = "drop")

print(detrended_stats, n = Inf)

write_csv(detrended_stats, "environmental_detrending_summary.csv")

# STEP 2: Prepare figure labels.

detrended_labels = detrended_stats %>%
  mutate(Subcommunity = factor(Subcommunity, levels = sub_levels),
         Label = paste0("Slope = ", sprintf("%.3f", Residual_slope),
                        " | R² = ", sprintf("%.3f", Residual_R2))) %>%
  group_by(Location) %>%
  arrange(Subcommunity, .by_group = TRUE) %>%
  mutate(x = -Inf, y = seq(.98, .77, length.out = n())) %>%
  ungroup()

# STEP 3: Plot microbial community dissimilarity against
# environmental distance after accounting for the geographic trend.

p_env_detrended = dd_detrended %>%
  mutate(Subcommunity = factor(Subcommunity, levels = sub_levels)) %>%
  ggplot(aes(Env_residual, Bray, color = Subcommunity)) +
  geom_point(alpha = .45, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = detrended_labels,
            aes(x = x, y = y, label = Label, color = Subcommunity),
            inherit.aes = FALSE, hjust = -0.05, size = 3.7,
            fontface = "bold", show.legend = FALSE) +
  facet_wrap(~Location, nrow = 1) +
  scale_color_manual(values = sub_colors) +
  scale_y_continuous(breaks = seq(0, 1, .2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Spatially detrended environmental distance (residual)",
       y = "Bray–Curtis dissimilarity") +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")

p_env_detrended

ggsave("environmental_distance_spatially_detrended.pdf",
       p_env_detrended, width = 9, height = 6, device = cairo_pdf)


# ============================================================
# COMBINED RDA FIGURE
# Top row: Biloxi Bay
# Bottom row: Pascagoula Bay
# Columns: Total, Abundant, CRT, Rare
# ============================================================

sub_levels = c("Total", "Abundant", "CRT", "Rare")
bay_levels = c("Biloxi Bay", "Pascagoula Bay")
env_vars = c("Salinity_PSU", "Temperature_C", "pH", "ODO_mg_L")

# 18.1. Prepare environmental variables -----------------------

stopifnot(all(c("Location", "StationID", env_vars) %in% names(env_station)))

env_rda = env_station %>%
  select(Location, StationID, all_of(env_vars)) %>%
  mutate(across(all_of(env_vars), as.numeric))

stopifnot(!anyDuplicated(env_rda$StationID))

# 18.2. Run Hellinger-RDA -------------------------------------

run_rda = function(location, subcommunity){
  
  # Station-level ASV community matrix
  comm = make_comm(location, subcommunity, "ASV", normalize_subset = TRUE)
  comm = as.matrix(comm)
  storage.mode(comm) = "numeric"
  
  # Match station-level environmental data
  env_sub = env_rda %>%
    filter(Location == location, StationID %in% rownames(comm)) %>%
    arrange(match(StationID, rownames(comm)))
  
  stopifnot(identical(env_sub$StationID, rownames(comm)))
  stopifnot(!anyNA(env_sub[, env_vars]))
  stopifnot(all(is.finite(comm)), all(comm >= 0))
  
  # Remove absent ASVs
  comm = comm[, colSums(comm) > 0, drop = FALSE]
  stopifnot(all(rowSums(comm) > 0))
  
  # Hellinger transformation: sqrt(relative abundance)
  comm_hel = vegan::decostand(comm, method = "hellinger")
  
  # Standardize environmental variables
  env_scaled = env_sub %>%
    select(all_of(env_vars)) %>%
    mutate(across(everything(), ~as.numeric(scale(.x))))
  
  stopifnot(all(is.finite(as.matrix(env_scaled))))
  
  # RDA model
  set.seed(123)
  
  model = vegan::rda(comm_hel ~ ., data = env_scaled)
  
  # Permutation test and adjusted R2
  test = vegan::anova.cca(model, permutations = 9999)
  adj_r2 = vegan::RsquareAdj(model)$adj.r.squared
  p_value = test$`Pr(>F)`[1]
  
  # Extract site scores
  site_scores = as.data.frame(vegan::scores(model, display = "sites",
                                            choices = 1:2, scaling = 2)) %>%
    rownames_to_column("StationID") %>%
    left_join(env_sub %>% select(StationID, Salinity_PSU), by = "StationID")
  
  # Extract environmental arrows
  arrows = as.data.frame(vegan::scores(model, display = "bp",
                                       choices = 1:2, scaling = 2)) %>%
    rownames_to_column("Variable") %>%
    mutate(Label = recode(Variable,
                          Salinity_PSU = "Salinity",
                          Temperature_C = "Temperature",
                          ODO_mg_L = "DO",
                          pH = "pH"))
  
  # Axis percentages relative to total community inertia
  eig = model$CCA$eig
  axis_pct = 100 * eig[1:2] / model$tot.chi
  
  list(model = model, sites = site_scores, arrows = arrows,
       adj_r2 = adj_r2, p_value = p_value, axis_pct = axis_pct,
       n_sites = nrow(comm), n_asvs = ncol(comm))
}

# 18.3. Run all eight RDA models -------------------------------

rda_results = crossing(Location = bay_levels, Subcommunity = sub_levels) %>%
  mutate(Location = factor(Location, levels = bay_levels),
         Subcommunity = factor(Subcommunity, levels = sub_levels)) %>%
  arrange(Location, Subcommunity) %>%
  mutate(Result = map2(as.character(Location), as.character(Subcommunity), run_rda))

saveRDS(rda_results, "distance-decay2/RDA/RDA_models.rds")

# 18.4. Extract RDA statistics --------------------------------

rda_stats = rda_results %>%
  transmute(Location, Subcommunity,
            Adj_R2 = map_dbl(Result, "adj_r2"),
            P_value = map_dbl(Result, "p_value"),
            RDA1_percent = map_dbl(Result, ~.x$axis_pct[1]),
            RDA2_percent = map_dbl(Result, ~.x$axis_pct[2]),
            N_stations = map_int(Result, "n_sites"),
            N_ASVs = map_int(Result, "n_asvs"))

print(rda_stats, n = Inf)

write_csv(rda_stats, "RDA_statistics.csv")

# ============================================================
# 19. COMBINED RDA FIGURE
# ============================================================

# 19.1. Consistent salinity scale ------------------------------

sal_range = range(env_rda$Salinity_PSU, na.rm = TRUE)

# 19.2. RDA plotting function ---------------------------------

plot_rda = function(result){
  
  ggplot(result$sites, aes(RDA1, RDA2)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.4) +
    geom_point(aes(color = Salinity_PSU), size = 2, alpha = 0.9) +
    geom_segment(data = result$arrows,
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 inherit.aes = FALSE, color = "#253DB8", linewidth = 0.8,
                 arrow = grid::arrow(length = unit(0.18, "cm"))) +
    geom_text_repel(data = result$arrows,
                    aes(RDA1, RDA2, label = Label),
                    inherit.aes = FALSE, color = "#253DB8",
                    size = 2, fontface = "bold",
                    max.overlaps = Inf, min.segment.length = 0) +
    scale_color_gradient(low = "#0072B2", high = "#D55E00",
                         limits = sal_range, name = "Salinity (PSU)") +
    labs(x = paste0("RDA1 (", sprintf("%.1f", result$axis_pct[1]), "%)"),
         y = paste0("RDA2 (", sprintf("%.1f", result$axis_pct[2]), "%)"),
         subtitle = paste0("Adj. R² = ", sprintf("%.3f", result$adj_r2),
                           " | p ", ifelse(result$p_value < 0.001, "< 0.001",
                                           paste0("= ", sprintf("%.3f", result$p_value))))) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11, color = "black"),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 11))
}

# 19.3. Generate all eight panels ------------------------------

rda_plots = rda_results %>%
  arrange(Location, Subcommunity) %>%
  mutate(Plot = pmap(list(Result, Location, Subcommunity), function(Result, Location, Subcommunity){
    plot_rda(Result) +
      labs(title = paste0(Location, " | ", Subcommunity)) +
      theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
  }))

p_rda_combined = wrap_plots(rda_plots$Plot, ncol = 4, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right",
        plot.tag = element_text(size = 11, face = "bold"))

p_rda_combined

ggsave("RDA_Biloxi_Pascagoula_combined.pdf", p_rda_combined, width = 10, height = 7, device = cairo_pdf)


#  LEVINS' NICHE WIDTH BY ASV AND SUBCOMMUNITY -----------------------------
# Width describes how evenly ASV abundances are distributed across stations.
# This is an inferred occupancy breadth, not an experimentally measured niche.
# 1. Prepare log-transformed niche widths

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


# Pairwise Wilcoxon if Kruskal-Wallis is significant
niche_pairwise <- niche_width %>% group_split(Location) %>%
  set_names(unique(niche_width$Location)) %>%
  map(~pairwise.wilcox.test(.x$Niche_width, .x$Subcommunity, p.adjust.method = "BH"))

niche_pairwise


niche_plot = niche_width |> mutate(
    Log_Niche = log10(Niche_width),
    Subcommunity = factor(Subcommunity, levels = c("Abundant", "CRT", "Rare")))

# 2. Overall Kruskal-Wallis test
niche_tests = niche_plot |>
  group_by(Location) |>
  kruskal_test(Log_Niche ~ Subcommunity)

niche_tests

# 3. Pairwise Wilcoxon tests with BH correction
niche_pairwise = niche_plot |>
  group_by(Location) |>
  wilcox_test(Log_Niche ~ Subcommunity, p.adjust.method = "BH") |>
  adjust_pvalue(method = "BH") |>
  add_significance("p.adj") |>
  add_xy_position(x = "Subcommunity", fun = "max", step.increase = 0.08)

niche_pairwise |>
  select(Location, group1, group2, p, p.adj, p.adj.signif, y.position)

# 4. Plot with correctly positioned brackets
p_niche = ggplot(niche_plot, aes(Subcommunity, Log_Niche, fill = Subcommunity)) +
  geom_violin(alpha = 0.7, trim = TRUE, width = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  stat_pvalue_manual(niche_pairwise, label = "p.adj.signif",
                     tip.length = 0.01, bracket.size = 0.5, size = 5,
                     hide.ns = FALSE) +
  facet_wrap(~Location) +
  scale_fill_manual(values = sub_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = expression(log[10]*" Levins' niche width")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_niche

ggsave("niche_width_significance_corrected.pdf", p_niche, width = 8, height = 6, device = cairo_pdf)

# 22.2. Run pNST for each bay and subcommunity -----------------

run_pnst = function(location, subcommunity){
  
  samples = meta %>% filter(Location == location) %>% pull(SampleID)
  
  taxa_keep = if(subcommunity == "Total") asv_class$ASV else
    asv_class %>% filter(Subcommunity == subcommunity) %>% pull(ASV)
  
  comm = otu[intersect(taxa_keep, rownames(otu)),
             intersect(samples, colnames(otu)), drop = FALSE] %>% t()
  
  # Remove taxa absent from the selected subset and samples with no reads
  comm = comm[, colSums(comm) > 0, drop = FALSE]
  comm = comm[rowSums(comm) > 0, , drop = FALSE]
  
  # Match community taxa with phylogenetic tree
  taxa_common = intersect(colnames(comm), tree$tip.label)
  comm = comm[, taxa_common, drop = FALSE]
  
  if(nrow(comm) < 3 || ncol(comm) < 2)
    stop("Insufficient samples or taxa: ", location, " / ", subcommunity)
  
  if(any(rowSums(comm) == 0))
    stop("Zero-read samples after matching tree: ", location, " / ", subcommunity)
  
  tr = ape::keep.tip(tree, taxa_common)
  comm = comm[, tr$tip.label, drop = FALSE]
  storage.mode(comm) = "numeric"
  
  grp = data.frame(Group = rep(location, nrow(comm)), row.names = rownames(comm))
  
  pd_dir = tempfile(pattern = paste0("pNST_", gsub(" ", "_", location), "_", subcommunity, "_"))
  dir.create(pd_dir, recursive = TRUE)
  on.exit(unlink(pd_dir, recursive = TRUE, force = TRUE), add = TRUE)
  
  message("Running pNST: ", location, " | ", subcommunity,
          " | Samples: ", nrow(comm), " | ASVs: ", ncol(comm))
  
  NST::pNST(comm = comm, tree = tr, group = grp, pd.wd = pd_dir,
            abundance.weighted = TRUE, rand = 1000, phylo.shuffle = TRUE,
            nworker = 4, between.group = FALSE)
}

# 22.3. Run or load saved models -------------------------------

model_file = "pNST_models_corrected.rds"
rerun_pnst = TRUE

safe_run_pnst = purrr::possibly(run_pnst, otherwise = NULL, quiet = FALSE)

if(file.exists(model_file) && !rerun_pnst) {
  pnst_models = readRDS(model_file)
} else {
  pnst_models = tidyr::crossing(Location = locations, Subcommunity = subcommunities) %>%
    mutate(pNST = purrr::map2(Location, Subcommunity, safe_run_pnst)) %>%
    filter(!purrr::map_lgl(pNST, is.null))
  
  if(nrow(pnst_models) == 0) stop("No valid pNST models. Check run_pnst() before saving.")
  
  saveRDS(pnst_models, model_file)
}

print(pnst_models %>% select(Location, Subcommunity))

# 22.4. Extract group-level normalized pNST -------------------

pnst_summary = pnst_models %>%
  mutate(Group_result = map(pNST, ~as_tibble(.x$index.grp))) %>%
  select(Location, Subcommunity, Group_result) %>%
  unnest(Group_result)

stopifnot("NST.i.bMNTD" %in% names(pnst_summary))

pnst_plot = pnst_summary %>%
  mutate(pNST_percent = `NST.i.bMNTD` * 100,
         MST_percent = `MST.i.bMNTD` * 100,
         Subcommunity = factor(Subcommunity, levels = subcommunities),
         Location = factor(Location, levels = locations),
         Process = case_when(pNST_percent > 50 ~ "Stochastic tendency",
                             pNST_percent < 50 ~ "Deterministic tendency",
                             TRUE ~ "Threshold"))

print(pnst_plot, n = Inf)

write_csv(pnst_plot, "pNST_group_summary.csv")

# 22.5. Group-level pNST bar plot ------------------------------

p_pnst = ggplot(pnst_plot, aes(Subcommunity, pNST_percent, fill = Subcommunity)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pNST_percent)),
            vjust = -0.5, size = 5, fontface = "bold") +
  facet_wrap(~Location) +
  scale_fill_manual(values = sub_colors) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, 105)) +
  labs(x = "Subcommunity", y = "Phylogenetic normalized stochasticity ratio (pNST, %)") +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.6, "cm"))

p_pnst

ggsave("pNST_subcommunity.pdf", p_pnst, width = 8, height = 6, device = cairo_pdf)

# 22.6. Extract valid pairwise MST and ST values --------------

pnst_pair = pnst_models %>%
  mutate(Pair_result = map(pNST, ~as_tibble(.x$index.pair))) %>%
  select(Location, Subcommunity, Pair_result) %>%
  unnest(Pair_result, names_repair = "unique")

# Identify pairwise index columns across NST package versions
mst_col = grep("^MST\\.ij", names(pnst_pair), value = TRUE)
st_col = grep("^ST\\.ij", names(pnst_pair), value = TRUE)

stopifnot(length(mst_col) == 1, length(st_col) == 1)

pnst_pair_plot = pnst_pair %>%
  mutate(MST_percent = .data[[mst_col]] * 100,
         ST_percent = .data[[st_col]] * 100,
         Subcommunity = factor(Subcommunity, levels = subcommunities),
         Location = factor(Location, levels = locations)) %>%
  filter(is.finite(MST_percent))

write_csv(pnst_pair_plot, "pNST_pairwise_values.csv")

# 22.7. Pairwise modified stochasticity ratio (MST) violin ----

p_mst_violin = ggplot(pnst_pair_plot, aes(Subcommunity, MST_percent, fill = Subcommunity)) +
  geom_violin(trim = TRUE, alpha = 0.75, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = NA, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.8, color = "black") +
  facet_wrap(~Location) +
  scale_fill_manual(values = sub_colors) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(mult = c(0.03, 0.05))) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Subcommunity", y = "Pairwise modified stochasticity ratio (MST, %)") +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.6, "cm"))

p_mst_violin

ggsave("dMST_pairwise_violin.pdf", p_mst_violin, width = 8, height = 6, device = cairo_pdf)

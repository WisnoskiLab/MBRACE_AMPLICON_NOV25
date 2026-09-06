Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
install.packages(Packages)
lapply(Packages, library, character.only = TRUE)

Sys.unsetenv("GITHUB_TOKEN")
Sys.unsetenv("GITHUB_PAT")


remotes::install_github("cozygene/FEAST", auth_token = NULL)

library(FEAST)
library(dplyr)

metadata_path <- "~/Documents/R-Git/MBRACE Analysis/November_2025/FEAST-data/metadata-feast.txt"
otu_path <- "~/Documents/R-Git/MBRACE Analysis/November_2025/FEAST-data/merged_otutable.txt"

# -----------------------------
# 1. Read metadata
# -----------------------------
metadata <- read.delim(
  metadata_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(metadata) <- trimws(colnames(metadata))

# Make all source and sink samples part of same FEAST source-tracking group
metadata$id <- 1

# Define soil locations as separate source environments
metadata <- metadata %>%
  mutate(
    Env = case_when(
      grepl("^S_EPP", SampleID) ~ "Soil_EPP",
      grepl("^S_GCP", SampleID) ~ "Soil_GCP",
      grepl("^S_HLP", SampleID) ~ "Soil_HLP",
      grepl("^S_PCP", SampleID) ~ "Soil_PCP",
      grepl("^S_RSP", SampleID) ~ "Soil_RSP",
      grepl("^S_RKM", SampleID) ~ "Soil_RKM",
      grepl("^W_", SampleID) ~ "Water",
      TRUE ~ Env
    )
  )

# Set sample IDs as metadata row names
rownames(metadata) <- metadata$SampleID

# -----------------------------
# 2. Read QIIME2 exported OTU table
# -----------------------------

first_line <- readLines(otu_path, n = 1)
skip_n <- ifelse(grepl("^# Constructed", first_line), 1, 0)

otu_raw <- read.delim(
  otu_path,
  header = TRUE,
  sep = "\t",
  skip = skip_n,
  check.names = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE
)

# First column is #OTUID / ASV ID
colnames(otu_raw)[1] <- "FeatureID"

rownames(otu_raw) <- otu_raw$FeatureID
otu_raw$FeatureID <- NULL

otu_raw <- as.matrix(otu_raw)
storage.mode(otu_raw) <- "numeric"

# -----------------------------
# 3. IMPORTANT: transpose for FEAST
# -----------------------------
# QIIME table: ASVs rows, samples columns
# FEAST C object: samples rows, ASVs columns

C <- t(otu_raw)

# -----------------------------
# 4. Match OTU table and metadata
# -----------------------------

common_samples <- intersect(rownames(C), metadata$SampleID)

cat("Common samples:", length(common_samples), "\n")

# Check missing samples
cat("Metadata samples missing from OTU table:\n")
print(setdiff(metadata$SampleID, rownames(C)))

cat("OTU table samples missing from metadata:\n")
print(setdiff(rownames(C), metadata$SampleID))

# Keep only common samples
C <- C[common_samples, ]

metadata2 <- metadata[common_samples, ]

# Final check
stopifnot(all(rownames(C) == metadata2$SampleID))

# -----------------------------
# 5. Check source-sink structure
# -----------------------------

table(metadata2$SourceSink)
table(metadata2$Env, metadata2$SourceSink)
table(metadata2$id, metadata2$SourceSink)

# -----------------------------
# 6. Run FEAST
# -----------------------------

dir.create("~/FEAST/Data_files/", recursive = TRUE, showWarnings = FALSE)

# Replace NA with zero
C[is.na(C)] <- 0

# Make sure counts are numeric integers
storage.mode(C) <- "numeric"

# Remove ASVs with zero counts across all matched samples
C <- C[, colSums(C) > 0]

# Remove samples with zero total reads, if any
C <- C[rowSums(C) > 0, ]

# Match metadata again after filtering
metadata2 <- metadata2[rownames(C), ]

stopifnot(all(rownames(C) == metadata2$SampleID))

source_samples <- metadata2$SampleID[metadata2$SourceSink == "Source"]
sink_samples   <- metadata2$SampleID[metadata2$SourceSink == "Sink"]

length(source_samples)
length(sink_samples)

dir.create("~/FEAST/Data_files/", recursive = TRUE, showWarnings = FALSE)

feast_list <- list()

for (sink in sink_samples) {
  
  message("Running FEAST for sink: ", sink)
  
  # Keep all soil sources + one water sink
  keep_samples <- c(source_samples, sink)
  
  C_sub <- C[keep_samples, , drop = FALSE]
  meta_sub <- metadata2[keep_samples, , drop = FALSE]
  
  # Remove ASVs absent in this subset
  C_sub <- C_sub[, colSums(C_sub) > 0, drop = FALSE]
  
  # FEAST metadata rownames must match C rownames
  rownames(meta_sub) <- meta_sub$SampleID
  
  # Make sure source and sink are in same id group
  meta_sub$id <- 1
  
  # Final checks
  stopifnot(all(rownames(C_sub) == meta_sub$SampleID))
  stopifnot(sum(meta_sub$SourceSink == "Sink") == 1)
  stopifnot(sum(meta_sub$SourceSink == "Source") >= 1)
  
  # Run FEAST
  feast_result <- FEAST(
    C = C_sub,
    metadata = meta_sub,
    different_sources_flag = 0,
    dir_path = "~/FEAST/Data_files/",
    outfile = paste0("FEAST_", sink)
  )
  
  feast_list[[sink]] <- feast_result
}




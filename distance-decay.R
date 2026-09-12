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

out_base <- "November_2025/16S"
# GEO Data read in  -----------
#
#   Jennifer Franks - March 2023 
#
# -----------------------

library(dplyr)
library(monocle3)
library(Matrix)
library(tidyr)
library(SingleCellExperiment)

source("h5_into_monocle.R")
setwd("")

# ---- Load in Data & Process ------ 

folder_path <- getwd()


matrix_files <- list.files(folder_path, pattern = "\\_matrix.mtx.gz$", full.names = TRUE)

sample_names <- sub("\\_matrix.mtx.gz$",  "", 
                    list.files(folder_path, pattern = "\\_matrix.mtx.gz$", full.names = F))



cds.list<- list()

for (i in 1:length(sample_names)){
  cds <- load_mm_data(mat_path = paste0(sample_names[i], "_matrix.mtx.gz"), 
                      feature_anno_path = paste0(sample_names[i], "_features.tsv.gz"), # may be "genes" not "features"
                      cell_anno_path = paste0(sample_names[i], "_barcodes.tsv.gz"))
  colnames(fData(cds)) <- c("gene_short_name")
  fData(cds)$id <- rownames(fData(cds))
  samplenameparts <- strsplit(sample_names[i], "_")[[1]]
  pData(cds)$GSM <- samplenameparts[1]
  pData(cds)$Treatment <- samplenameparts[2]
  pData(cds)$Timepoint <- samplenameparts[3]
  pData(cds)$ID <- paste(samplenameparts[2], samplenameparts[3], sep="_")
  
  counts_matrix = t(counts(cds))
  genes = as.array(rownames(fData(cds)))
  
  
  library(reticulate)
  use_condaenv("myenv")
  py_run_string("import scrublet as scr")
  py_run_string("import scipy.io")
  py_run_string("import matplotlib.pyplot as plt")
  py_run_string("import numpy as np")
  py_run_string("import os")
  py_run_string("scrub = scr.Scrublet(r.counts_matrix, expected_doublet_rate=0.06)")
  py_run_string("scrubresults = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)")
  
  scrubresults <- py$scrubresults
  pData(cds)$scrublet_score <-  as.numeric(scrubresults[[1]])
  pData(cds)$scrublet_call <- scrubresults[[2]]
  
  
  mt <- c(grep("^MT-",fData(cds)$gene_short_name,value = TRUE), grep("^mt-",fData(cds)$gene_short_name,value = TRUE))
  mt_cds <- cds[which(fData(cds)$gene_short_name %in% mt),]
  pData(cds)$perc_mitochondrial_umis <- Matrix::colSums(exprs(mt_cds))/Matrix::colSums(exprs(cds)) * 100
  
  rb <- grep("^RP[LS]",fData(cds)$gene_short_name,value = TRUE)
  rb_cds <- cds[which(fData(cds)$gene_short_name %in% rb),]
  pData(cds)$perc_ribosomal_umis <- Matrix::colSums(exprs(rb_cds))/Matrix::colSums(exprs(cds)) * 100
  
  cds.list <- c(cds.list, cds)
  remove(cds, scrubresults, counts_matrix, genes, rb, rb_cds, mt, mt_cds)
  
}

raw.cds <- combine_cds(cds.list)


# --- parameters -----
PCs_to_use <- 100 # 
UMImax <- 15000
UMImin <- 100
mitoMax <- 10.
set.seed(8248959)
cluster.k <- 20
scrublet.score.max <- 0.2
# --------------------

# Filtering QC ----

# Doublet filtering using scrublet calls/scores
pData(raw.cds)$scrublet_call[is.na(pData(raw.cds)$scrublet_call)] <- "Unk"
raw.cds <- raw.cds[,which(pData(raw.cds)$scrublet_call != "Doublet")]
#raw.cds <- raw.cds[,which(pData(raw.cds)$scrublet_score < 0.2)]



# UMI filtering 
#raw.cds <- raw.cds[,which(pData(raw.cds)$n.umi > UMImin)]


# Mitochondrial reads filtering 
raw.cds <- raw.cds[,pData(raw.cds)$perc_mitochondrial_umis < mitoMax]


# Basic Preprocessing --- 

cds <- raw.cds %>% 
  estimate_size_factors() %>%
  preprocess_cds(num_dim=PCs_to_use) %>% 
  detect_genes() %>%
  reduce_dimension() %>%
  cluster_cells(k=50, random_seed = 3752)


pData(cds)$original_cluster <- clusters(cds)


# Save Monocle Object ----

save_monocle_objects(cds = cds, 
                     directory_path = "")




# --------------------------------------------------------------
# Spatial visualization of Hallmark and custom gene signatures
# 
# Jennifer Franks (Dec 2025)
#
# Overview:
#   - Pull gene sets from MSigDB
#   - Compute per-spot signature scores using monocle3
#   - Store scores as numeric columns in colData(cds)
#   - Visualize signatures by overlaying scores on Visium H&E images
#     using plot_signature_spatial_sample() and
#     plot_all_samples_signature_spatial() from visium_monocle_helpers.R
#
# Requirements:
#   - cds:
#       * monocle3 cell_data_set built with build_visium_cds()
#       * attach_spatial_to_cds() has been run so colData(cds) contains
#           x_highres, y_highres, x_lowres, y_lowres, sample, sample_label
#   - imgs:
#       * list of image data.frames returned by build_visium_images(cds, ...)
#           imgs$highres : high-resolution H&E pixels
#           imgs$lowres  : low-resolution H&E pixels
#       * each image df has columns:
#           x, y, rgb.val, sample, sample_label, sample_group, resolution,
#           contains_spots ("yes"/"no")
#   - Signature scores:
#       * for each gene set (e.g. "HALLMARK_INFLAMMATORY_RESPONSE"),
#         a numeric column is added to colData(cds):
#           colData(cds)[, "HALLMARK_INFLAMMATORY_RESPONSE"] <- <numeric vector>
# --------------------------------------------------------------

library(msigdbr)
library(monocle3)
library(ggplot2)
library(dplyr)

# --------------------------------------------------------------
# 1) Get MSigDB gene sets of interest
# --------------------------------------------------------------

all_gene_sets <- msigdbr(species = "Homo sapiens")

# Hallmark gene sets
msigdbr_hallmark <- all_gene_sets %>%
  dplyr::filter(grepl("HALLMARK", gs_name)) %>%
  distinct()

# Osteoclast-related GO BP gene sets
msigdbr_osteoclast <- all_gene_sets %>%
  dplyr::filter(grepl("GOBP", gs_name),
                grepl("OSTEOCLAST", gs_name)) %>%
  distinct()

# Combine all of these into one data.frame
msigdbr_df <- bind_rows(
  msigdbr_hallmark,
  msigdbr_osteoclast,
)

# Optional: table of gene set -> gene mappings (if needed elsewhere)
# msigdbr_t2g <- msigdbr_df %>%
#   distinct(gs_name, gene_symbol) %>%
#   as.data.frame()

# Split into a named list: each element is a character vector of genes for that gene set
msigdbr_list <- split(msigdbr_df$gene_symbol, msigdbr_df$gs_name)

# --------------------------------------------------------------
# 2) Optionally add a custom gene signature
# --------------------------------------------------------------

custom_gene_list       <- c("SPP1", "ACP5", "TGFB1")
custom_gene_list_name  <- "CUSTOM_SPP1_ACP5_TGFB1"

# Append to the list of gene sets
msigdbr_list[[custom_gene_list_name]] <- custom_gene_list

# Inspect gene set names
names(msigdbr_list)

# --------------------------------------------------------------
# 3) Compute signature scores per spot using monocle3
#    aggregate_gene_expression expects a gene grouping matrix:
#       rownames = genes in cds
#       column 1 = gene names
#       column 2 = group labels (0 or 1 here, in vs. out of gene set)
# --------------------------------------------------------------

# Initialize results with per-spot metadata
results <- pData(cds)  # or as.data.frame(colData(cds)) if you are fully in monocle3

# 3) Compute signature scores per spot and store in colData(cds)
#    aggregate_gene_expression expects a "gene.group" matrix with:
#      rownames = gene ids
#      column 1 = gene_short_name
#      column 2 = group label (0 or 1 for "in set")
#
#    We then pull the row "1" (genes in set) and attach scores
#    as new columns in colData(cds).

# For convenience, also keep a plain data.frame copy for fast inspection
results <- as.data.frame(colData(cds))

for (gs_name in names(msigdbr_list)) {
  message("Computing signature: ", gs_name)
  
  gene_set_gene_list <- msigdbr_list[[gs_name]]
  
  # membership matrix
  gene.group <- matrix(
    nrow = nrow(cds),
    ncol = 2
  )
  gene.group[, 1] <- as.matrix(rowData(cds)$gene_short_name)
  gene.group[, 2] <- as.numeric(gene.group[, 1] %in% gene_set_gene_list)
  rownames(gene.group) <- rownames(rowData(cds))
  
  # aggregate expression
  agg_expr <- aggregate_gene_expression(
    cds,
    gene.group,
    scale_agg_values = TRUE,
    norm_method      = "log",
    gene_agg_fun     = "sum"
  )
  
  # "1" row = genes in set
  sig_vec <- as.numeric(agg_expr["1", colnames(cds)])
  
  # attach to colData(cds) and results data.frame
  colData(cds)[, gs_name] <- sig_vec
  results[[gs_name]]      <- sig_vec
  
  rm(gene.group, agg_expr, sig_vec)
}

# Quick check
colnames(results)[grepl("HALLMARK|GOBP|CUSTOM", colnames(results))]
head(results[, grepl("HALLMARK|GOBP|CUSTOM", colnames(results)), drop = FALSE])


# --------------------------------------------------------------
# 4) Plot spatial signature scores over high-res images
# --------------------------------------------------------------
#one sample only-
plot_signature_spatial_sample(
  cds           = cds,
  images        = imgs,
  signature_col = "HALLMARK_INFLAMMATORY_RESPONSE",
  sample_id     = "Sample_9874-WS-2",
  resolution    = "highres",
  palette       = "diverging",
  diverging_midpoint = 0
) 

#iterate over all samples-
plots_hallmark_inflam <- plot_all_samples_signature_spatial(
  cds           = cds,
  images        = imgs,
  signature_col = "HALLMARK_INFLAMMATORY_RESPONSE",
  resolution    = "highres",
  output_dir    = "plots/signatures_highres",
  file_prefix   = "HALLMARK_INFLAMMATORY_RESPONSE_",
  palette       = "diverging",
  diverging_midpoint = 0
)

plots_hallmark_inflam[["Sample_9874-WS-2"]]



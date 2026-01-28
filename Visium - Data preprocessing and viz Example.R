# --------------------------------------------------------------
# Build Visium monocle3 cds for Silica-CWP project
#
# Pipeline overview:
#   1. Discover Visium samples and paths on disk
#   2. Build and preprocess a combined monocle3 cds
#   3. Attach spatial coordinates and a "spatial tSNE" embedding
#   4. Build H&E image data frames (highres and lowres)
#   5. Run a small set of QC and visualization plots
#   6. Save cds for downstream analyses
#
# Note:
#   All reusable logic lives in visium_monocle_helpers.R
#   This script is intended to be a readable outline of the workflow.
# --------------------------------------------------------------

project_dir <- "/data/hps/assoc/private/franks_lab/user/jfra11/projects/Silica/spatial-prelim/visium"
setwd(project_dir)

source(file.path(project_dir, "visium_monocle_helpers.R"))

# --------------------------------------------------------------
# 1. Discover Visium samples and build sample_table
# --------------------------------------------------------------
# Convention:
#   sample_id    = full sequencing sample name (for example, "Sample_9874-WS-2")
#   sample_label = short colloquial label (for example, "CWP", "Control")
#   sample_group = grouping variable (for example, "CWP", "Control", "Case", "Timepoint1")
#
# Here there are only two samples:
#   Sample_9874-WS-2  -> CWP
#   Sample_9874-WS-4  -> Control

sample_table <- discover_visium_sample_table(
  root_dir      = file.path(project_dir, "data"),
  sample_dirs   = c("Sample_9874-WS-2", "Sample_9874-WS-4"),
  sample_ids    = c("Sample_9874-WS-2", "Sample_9874-WS-4"),
  sample_labels = c("CWP", "Control"),
  sample_groups = c("CWP", "Control")
)

sample_table

# --------------------------------------------------------------
# 2. Build and preprocess monocle3 cds
# --------------------------------------------------------------
# Steps inside build_visium_cds:
#   * Read each filtered_feature_bc_matrix.h5 into a cell_data_set
#   * combine_cds to merge per sample cds
#   * Map numeric sample IDs to sample_id / sample_label / sample_group
#   * Basic QC:
#       - Filter on UMI (min_umi)
#       - detect_genes and estimate_size_factors
#       - Filter on num_genes_expressed (min_genes)
#   * preprocess_cds, align_cds by "sample", reduce_dimension (UMAP), cluster_cells

cds <- build_visium_cds(
  sample_table = sample_table,
  min_umi      = 100,
  min_genes    = 10,
  num_dim      = 25,
  cluster_res  = 1e-3
)

cds

# --------------------------------------------------------------
# 3. Attach spatial coordinates and "spatial tSNE"
# --------------------------------------------------------------
# attach_spatial_to_cds:
#   * Reads tissue_positions*.csv and scalefactors_json.json
#   * Computes x/y_highres and x/y_lowres coordinates
#   * Joins spatial metadata into colData(cds) by barcode
#
# add_spatial_tsne:
#   * Writes full resolution pixel coordinates into the "tSNE" slot
#     so that plot_cells(..., reduction_method = "tSNE") shows the
#     physical spatial layout of spots on the slide

cds <- attach_spatial_to_cds(cds, sample_table)
cds <- add_spatial_tsne(cds)

# Quick sanity check: spatial columns should exist now
colnames(colData(cds))[grepl("pxl_|_lowres$|_highres$", colnames(colData(cds)))]

# --------------------------------------------------------------
# 4. Build H&E image data frames
# --------------------------------------------------------------
# build_visium_images:
#   * Loads tissue_hires_image.png and tissue_lowres_image.png
#   * Converts them to per pixel data frames with columns:
#       x, y, rgb.val, sample, sample_label, sample_group,
#       resolution (highres or lowres), contains_spots
#   * contains_spots is TRUE inside a buffered bounding box around
#     the spot coordinates

imgs <- build_visium_images(cds, sample_table, img_buffer = 5)
images_highres <- imgs$highres
images_lowres  <- imgs$lowres

# --------------------------------------------------------------
# 5. Quick QC plots
# --------------------------------------------------------------
# These do not use images, they just sanity check cluster and gene
# level signals in lowres and highres coordinate systems.

# Clusters on lowres coordinates
plot_visium_clusters_lowres(cds)

# Clusters on highres coordinates
plot_visium_clusters_highres(cds)

# Example gene, lowres scatter, no H&E overlay
plot_visium_gene_lowres(cds, gene.name = "CTSK")

# Highres H&E + single gene for CWP sample
plot_visium_gene_spatial_sample(
  cds       = cds,
  images    = imgs,
  gene.name = "CTSK",
  sample_id = "Sample_9874-WS-2",
  resolution = "highres"
) # + ggtitle("tester") you can also amend ggplot properties like this

# Optional: use "spatial tSNE" embedding (fullres pixels) with monocle
plot_cells(
  cds,
  reduction_method  = "tSNE",
  color_cells_by    = "clusters",
  label_cell_groups = FALSE,
  cell_size         = 1
) + facet_wrap(~sample_label)

# --------------------------------------------------------------
# 6. Image based spatial plots (H&E plus spots)
# --------------------------------------------------------------
# The helper plot_all_samples_spatial wraps:
#   * plot_visium_spatial_sample for each sample
#   * Optionally saves one file per sample
#
# resolution can be "highres" or "lowres"
# color_by can be any column in colData(cds)

# 6a. Clusters over highres H&E, one plot per sample, also save to disk
plots_clusters_highres <- plot_all_samples_spatial(
  cds        = cds,
  images     = imgs,
  resolution = "highres",
  color_by   = "clusters",
  output_dir = "plots/spatial_clusters_highres",  # set to NULL to skip saving
  file_prefix = "clusters_highres_"
)

# Inspect first sample
plots_clusters_highres[[1]]
# Or by sample_id
plots_clusters_highres[["Sample_9874-WS-2"]]

# 6b. Example: sample_group over highres images
plots_group_highres <- plot_all_samples_spatial(
  cds        = cds,
  images     = imgs,
  resolution = "highres",
  color_by   = "sample_group",
  output_dir = "plots/spatial_samplegroup_highres",
  file_prefix = "samplegroup_highres_"
)

# 6c. Example: SLC40A1 levels over lowres H&E
# First attach the gene as a column in colData if desired, otherwise color_by can be "SLC40A1" directly
# Here we demonstrate using color_by with an existing metadata column:
#   For custom genes, better to build a separate helper or reuse plot_visium_gene_lowres.

# --------------------------------------------------------------
# 7. Save cds object for downstream analyses
# --------------------------------------------------------------
# Note:
#   writes cds plus metadata to a directory.

save_monocle_objects(
  cds,
  directory_path = "cds objects/2024-12-11_CDS",
  comment        = "Data read in from filtered h5 files, preprocessed and aligned by sample using default Monocle3 parameters."
)

# --------------------------------------------------------------
# 8. Optional: manual ggplot overlay example
# --------------------------------------------------------------
# Keeping a version here as a reference for custom tweaks.

ggplot() + 
  geom_raster(
    data = images_highres[images_highres$contains_spots == "yes", ],
    aes(x = x, y = y, fill = rgb.val)
  ) +
  scale_fill_identity() + 
  scale_y_reverse() +
  geom_point(
    data  = as.data.frame(colData(cds)),
    aes(x = x_highres, y = y_highres, colour = clusters),
    size  = 0.7,
    alpha = 0.8
  ) + 
  facet_wrap(~sample_label, scales = "free") + 
  theme_void() + 
  theme(
    text         = element_text(size = 24, color = "black"), 
    legend.text  = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    aspect.ratio = 1
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

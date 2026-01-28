# --------------------------------------------------------------
# Visium -> monocle3 reusable helpers
#
#   - Auto discover Visium samples on disk
#   - Load CellRanger h5 into monocle3 cell_data_set objects
#   - Attach spatial coordinates and H&E images
#   - Provide basic spatial QC and overlay plotting helpers
#
# Conventions:
#   colData(cds)$sample        = full sequencing sample ID
#                                (for example, "Sample_9874-WS-2")
#   colData(cds)$sample_label  = colloquial label
#                                (for example, "CWP", "Control"; may be non-unique)
#   colData(cds)$sample_group  = grouping label
#                                (for example, "Case", "Control", "Timepoint1")
#
# Jennifer Franks
# Started: 2024-12-11
# Refactored: 2025-12-08
# --------------------------------------------------------------

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(hdf5r)
  library(Matrix)
  library(jsonlite)
  library(tibble)
  library(imager)
})

# --------------------------------------------------------------
# 0) Auto discover Visium sample metadata
# --------------------------------------------------------------
# Looks for, under each sample directory:
#   - filtered_feature_bc_matrix.h5
#   - tissue_positions*.csv
#   - scalefactors_json.json
#   - tissue_hires_image.png
#   - tissue_lowres_image.png
#
# Arguments:
#   root_dir      : directory containing per sample subdirectories
#   sample_dirs   : optional vector of subdirectory names or full paths;
#                   if NULL, all first level subdirectories under root_dir
#                   are treated as samples
#   sample_ids    : full sample IDs (go to colData(cds)$sample);
#                   if NULL, defaults to basename(sample_dirs)
#   sample_labels : colloquial labels (go to colData(cds)$sample_label);
#                   if NULL, defaults to sample_ids
#   sample_groups : grouping labels (go to colData(cds)$sample_group);
#                   if NULL, defaults to sample_labels
#   barcode_suffix: suffix appended per sample to barcodes (for example, "_1");
#                   if NULL, defaults to "_1", "_2", ...
#
# Returns:
#   tibble with columns:
#     sample_id, sample_label, sample_group,
#     h5_file, positions_csv, scalefactors_json,
#     hires_image, lowres_image, barcode_suffix

discover_visium_sample_table <- function(root_dir,
                                         sample_dirs    = NULL,
                                         sample_ids     = NULL,
                                         sample_labels  = NULL,
                                         sample_groups  = NULL,
                                         barcode_suffix = NULL) {
  if (!dir.exists(root_dir)) {
    stop("root_dir does not exist: ", root_dir)
  }
  
  # Determine sample directories
  if (is.null(sample_dirs)) {
    sample_paths <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
  } else {
    sample_paths <- ifelse(
      dir.exists(sample_dirs),
      sample_dirs,
      file.path(root_dir, sample_dirs)
    )
  }
  
  if (length(sample_paths) == 0) {
    stop("No sample directories found under: ", root_dir)
  }
  
  n <- length(sample_paths)
  
  # Sample IDs
  if (is.null(sample_ids)) {
    sample_ids <- basename(sample_paths)
  } else if (length(sample_ids) != n) {
    stop("sample_ids must have length equal to number of sample_paths")
  }
  
  # Sample labels and groups
  if (is.null(sample_labels)) {
    sample_labels <- sample_ids
  } else if (length(sample_labels) != n) {
    stop("sample_labels must have length equal to number of sample_paths")
  }
  
  if (is.null(sample_groups)) {
    sample_groups <- sample_labels
  } else if (length(sample_groups) != n) {
    stop("sample_groups must have length equal to number of sample_paths")
  }
  
  # Barcode suffixes
  if (is.null(barcode_suffix)) {
    barcode_suffix <- paste0("_", seq_len(n))
  } else if (length(barcode_suffix) != n) {
    stop("barcode_suffix must have length equal to number of sample_paths")
  }
  
  # Discover per sample files
  h5_files          <- character(n)
  positions_files   <- character(n)
  scalefactors_json <- character(n)
  hires_imgs        <- character(n)
  lowres_imgs       <- character(n)
  
  for (i in seq_along(sample_paths)) {
    sp <- sample_paths[i]
    
    if (!dir.exists(sp)) {
      warning("Sample directory does not exist: ", sp)
      next
    }
    
    # filtered_feature_bc_matrix.h5
    h5_candidate <- list.files(
      sp,
      pattern = "filtered_feature_bc_matrix\\.h5$",
      full.names = TRUE,
      recursive = TRUE
    )
    h5_files[i] <- if (length(h5_candidate) > 0) h5_candidate[1] else NA_character_
    if (is.na(h5_files[i])) warning("No filtered_feature_bc_matrix.h5 found under ", sp)
    
    # tissue_positions*.csv
    pos_candidate <- list.files(
      sp,
      pattern = "tissue_positions.*\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    )
    positions_files[i] <- if (length(pos_candidate) > 0) pos_candidate[1] else NA_character_
    if (is.na(positions_files[i])) warning("No tissue_positions*.csv found under ", sp)
    
    # scalefactors_json.json
    sf_candidate <- list.files(
      sp,
      pattern = "scalefactors_json\\.json$",
      full.names = TRUE,
      recursive = TRUE
    )
    scalefactors_json[i] <- if (length(sf_candidate) > 0) sf_candidate[1] else NA_character_
    if (is.na(scalefactors_json[i])) warning("No scalefactors_json.json found under ", sp)
    
    # tissue_hires_image.png
    hi_candidate <- list.files(
      sp,
      pattern = "tissue_hires_image\\.png$",
      full.names = TRUE,
      recursive = TRUE
    )
    hires_imgs[i] <- if (length(hi_candidate) > 0) hi_candidate[1] else NA_character_
    if (is.na(hires_imgs[i])) warning("No tissue_hires_image.png found under ", sp)
    
    # tissue_lowres_image.png
    lo_candidate <- list.files(
      sp,
      pattern = "tissue_lowres_image\\.png$",
      full.names = TRUE,
      recursive = TRUE
    )
    lowres_imgs[i] <- if (length(lo_candidate) > 0) lo_candidate[1] else NA_character_
    if (is.na(lowres_imgs[i])) warning("No tissue_lowres_image.png found under ", sp)
  }
  
  tibble(
    sample_id         = sample_ids,
    sample_label      = sample_labels,
    sample_group      = sample_groups,
    h5_file           = h5_files,
    positions_csv     = positions_files,
    scalefactors_json = scalefactors_json,
    hires_image       = hires_imgs,
    lowres_image      = lowres_imgs,
    barcode_suffix    = barcode_suffix
  )
}

# --------------------------------------------------------------
# 1) Read CellRanger h5 into monocle3 cell_data_set
# --------------------------------------------------------------

read.cds.cellranger.h5.file <- function(h5.file) {
  if (!file.exists(h5.file)) {
    stop("File ", h5.file, " not found")
  }
  
  s <- H5File$new(h5.file, mode = "r")
  
  barcodes    <- s[["matrix/barcodes"]][]
  gene_ids    <- s[["matrix/features/id"]][]
  gene_names  <- s[["matrix/features/name"]][]
  featuretype <- s[["matrix/features/feature_type"]][]
  data        <- s[["matrix/data"]][]
  indices     <- s[["matrix/indices"]][] + 1
  indptr      <- s[["matrix/indptr"]][]
  shape       <- s[["matrix/shape"]][]
  
  h5close(s)
  
  if (length(levels(factor(featuretype))) > 1) {
    warning("Multiple feature types found in h5 (for example, gene expression plus antibody capture)")
  }
  
  gbm <- sparseMatrix(
    x    = data,
    i    = indices,
    p    = indptr,
    dims = shape
  )
  
  pData.df <- data.frame(
    barcode = barcodes,
    stringsAsFactors = FALSE
  )
  
  fData.df <- data.frame(
    id             = gene_ids,
    gene_short_name = gene_names,
    feature_type   = featuretype,
    stringsAsFactors = FALSE
  )
  
  rownames(pData.df) <- pData.df$barcode
  rownames(fData.df) <- fData.df$id
  rownames(gbm)      <- rownames(fData.df)
  colnames(gbm)      <- rownames(pData.df)
  
  suppressWarnings({
    cds <- new_cell_data_set(
      gbm,
      cell_metadata = pData.df,
      gene_metadata = fData.df
    )
  })
  
  cds
}

# --------------------------------------------------------------
# 2) Read Visium spatial coordinates and scale factors
# --------------------------------------------------------------

read_visium_positions <- function(positions_csv,
                                  scalefactors_json,
                                  barcode_suffix = NULL,
                                  sample_id = NULL) {
  dt <- fread(positions_csv)
  
  if (!is.null(barcode_suffix)) {
    dt$barcode <- paste0(dt$barcode, barcode_suffix)
  }
  
  scale.factors <- jsonlite::fromJSON(scalefactors_json)
  
  dt$y_highres <- dt$pxl_row_in_fullres * scale.factors$tissue_hires_scalef
  dt$x_highres <- dt$pxl_col_in_fullres * scale.factors$tissue_hires_scalef
  dt$y_lowres  <- dt$pxl_row_in_fullres * scale.factors$tissue_lowres_scalef
  dt$x_lowres  <- dt$pxl_col_in_fullres * scale.factors$tissue_lowres_scalef
  
  if (!is.null(sample_id)) {
    dt$sample <- sample_id
  }
  
  dt
}

# --------------------------------------------------------------
# 3) Build and preprocess a combined Visium cds
# --------------------------------------------------------------
# sample_table must contain:
#   sample_id, sample_label, sample_group,
#   h5_file, positions_csv, scalefactors_json, barcode_suffix

build_visium_cds <- function(sample_table,
                             min_umi     = 100,
                             min_genes   = 10,
                             num_dim     = 25,
                             cluster_res = 1e-3) {
  n <- nrow(sample_table)
  cds.list <- vector("list", length = n)
  
  # Load each h5 into its own cds
  for (i in seq_len(n)) {
    message("Reading h5: ", sample_table$h5_file[i])
    cds_i <- read.cds.cellranger.h5.file(sample_table$h5_file[i])
    
    # If needed, enforce unique barcodes by suffixing:
    # colnames(exprs(cds_i)) <- paste0(colnames(exprs(cds_i)), sample_table$barcode_suffix[i])
    
    cds.list[[i]] <- cds_i
  }
  
  # Combine cds objects
  cds <- combine_cds(cds.list)
  
  # Map numeric sample labels to sample_table
  orig_sample    <- as.character(colData(cds)$sample)
  numeric_levels <- sort(unique(orig_sample))
  
  if (length(numeric_levels) != n) {
    stop("Number of unique cds$sample levels (", length(numeric_levels),
         ") does not match number of rows in sample_table (", n, ").")
  }
  
  map_sample_id    <- stats::setNames(sample_table$sample_id,    numeric_levels)
  map_sample_label <- stats::setNames(sample_table$sample_label, numeric_levels)
  map_sample_group <- stats::setNames(sample_table$sample_group, numeric_levels)
  
  colData(cds)$sample_numeric <- orig_sample
  colData(cds)$sample         <- unname(map_sample_id[orig_sample])
  colData(cds)$sample_label   <- unname(map_sample_label[orig_sample])
  colData(cds)$sample_group   <- unname(map_sample_group[orig_sample])
  
  # Basic QC
  cds <- cds[, Matrix::colSums(counts(cds)) > min_umi]
  cds <- cds %>%
    detect_genes() %>%
    estimate_size_factors()
  cds <- cds[, colData(cds)$num_genes_expressed > min_genes]
  
  # Preprocess, align by sample, reduce dimension, cluster
  cds <- cds %>%
    preprocess_cds(num_dim = num_dim) %>%
    align_cds(alignment_group = "sample") %>%
    reduce_dimension(umap.fast_sgd = TRUE) %>%
    cluster_cells(resolution = cluster_res)
  
  colData(cds)$clusters <- clusters(cds)
  
  cds
}

# --------------------------------------------------------------
# 4) Attach spatial coordinates to cds
# --------------------------------------------------------------

attach_spatial_to_cds <- function(cds, sample_table) {
  
  dt_list <- vector("list", length = nrow(sample_table))
  
  for (i in seq_len(nrow(sample_table))) {
    message("Reading spatial for: ", sample_table$positions_csv[i])
    dt_list[[i]] <- read_visium_positions(
      positions_csv     = sample_table$positions_csv[i],
      scalefactors_json = sample_table$scalefactors_json[i],
      barcode_suffix    = sample_table$barcode_suffix[i],
      sample_id         = sample_table$sample_id[i]
    )
  }
  
  dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  
  colData(cds)$barcode <- rownames(colData(cds))
  
  join_df <- as.data.frame(colData(cds)[, "barcode", drop = FALSE]) %>%
    dplyr::left_join(dt, by = "barcode")
  
  if (!all(rownames(colData(cds)) == join_df$barcode)) {
    warning("Row ordering mismatch after join; check barcode suffixing and combine_cds.")
  }
  
  existing_cols <- colnames(colData(cds))
  join_df_clean <- join_df %>%
    dplyr::select(-barcode, -sample) %>%
    dplyr::select(dplyr::any_of(setdiff(colnames(.), existing_cols)))
  
  colData(cds)$barcode <- NULL
  
  colData(cds) <- cbind(colData(cds), join_df_clean)
  
  cds
}

# --------------------------------------------------------------
# 5) Store unscaled spatial coordinates in tSNE slot
# --------------------------------------------------------------
# Uses full resolution pixel coordinates as the "tSNE" embedding so
# plot_cells(cds, reduction_method = "tSNE") shows spatial layout.

add_spatial_tsne <- function(cds) {
  if (!all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(colData(cds)))) {
    stop("pxl_col_in_fullres / pxl_row_in_fullres not present in colData(cds). ",
         "Did you run attach_spatial_to_cds()?")
  }
  
  reducedDim(x = cds, type = "tSNE") <-
    matrix(
      cbind(colData(cds)$pxl_col_in_fullres,
            -colData(cds)$pxl_row_in_fullres),
      ncol = 2
    )
  
  cds
}

# --------------------------------------------------------------
# 6) Quick QC plotting helpers (without H&E)
# --------------------------------------------------------------

plot_visium_clusters_lowres <- function(cds) {
  df <- as.data.frame(colData(cds))
  
  ggplot(df) +
    geom_point(aes(x = x_lowres, y = -y_lowres, color = clusters), size = 0.5) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme_void() +
    facet_wrap(~sample_label)
}

plot_visium_clusters_highres <- function(cds) {
  df <- as.data.frame(colData(cds))
  
  ggplot(df) +
    geom_point(aes(x = x_highres, y = -y_highres, color = clusters), size = 0.5) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme_void() +
    facet_wrap(~sample_label)
}

plot_visium_gene_lowres <- function(cds,
                                    gene.name,
                                    facet_by = c("sample_label", "sample")) {
  facet_by <- match.arg(facet_by)
  
  df <- as.data.frame(colData(cds))
  
  # get expression vector for this gene
  idx <- fData(cds)$gene_short_name == gene.name
  if (!any(idx)) {
    stop("Gene not found in fData(cds)$gene_short_name: ", gene.name)
  }
  expr_vec <- exprs(cds)[idx, , drop = TRUE]
  
  df[[gene.name]] <- log2(as.numeric(expr_vec) + 0.01)
  
  facet_formula <- as.formula(paste("~", facet_by))
  
  ggplot(df) +
    geom_point(
      aes(
        x = x_lowres,
        y = -y_lowres,               # flip y so it matches image orientation
        color = .data[[gene.name]]
      ),
      size = 0.5
    ) +
    scale_color_gradient(low = "white", high = "red", na.value = "grey90") +
    coord_fixed() +                 # keep spatial aspect ratio
    facet_wrap(facet_formula) +     # fixed scales to play nicely with coord_fixed
    theme_void() +
    labs(
      color = paste0("log2(", gene.name, " + 0.01)"),
      title = paste0("Lowres spatial expression: ", gene.name)
    ) +
    theme(
      plot.title  = element_text(hjust = 0.5),
      legend.text = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 10, color = "black")
    )
}


# --------------------------------------------------------------
# 7) Read H&E images and align to spatial coordinates
# --------------------------------------------------------------
# Returns:
#   list(highres = df, lowres = df)
#   where each df has columns:
#     x, y, rgb.val, sample, sample_label, sample_group,
#     resolution, contains_spots

build_visium_images <- function(cds, sample_table, img_buffer = 5) {
  required_cols <- c("x_highres", "y_highres", "x_lowres", "y_lowres", "sample")
  if (!all(required_cols %in% colnames(colData(cds)))) {
    stop("cds is missing one of: ", paste(required_cols, collapse = ", "),
         ". Did you run attach_spatial_to_cds() and build_visium_cds()?")
  }
  
  df_spots <- as.data.frame(colData(cds))
  df_spots$sample <- as.character(df_spots$sample)
  
  sample_table$sample_id    <- as.character(sample_table$sample_id)
  sample_table$sample_label <- as.character(sample_table$sample_label)
  sample_table$sample_group <- as.character(sample_table$sample_group)
  
  cds_samples    <- sort(unique(df_spots$sample))
  table_samples  <- sort(unique(sample_table$sample_id))
  common_samples <- intersect(cds_samples, table_samples)
  
  if (length(common_samples) == 0) {
    stop("No overlapping sample IDs between cds$sample and sample_table$sample_id.\n",
         "cds samples: ", paste(cds_samples, collapse = ", "), "\n",
         "sample_table$sample_id: ", paste(table_samples, collapse = ", "))
  }
  
  hi_list <- list()
  lo_list <- list()
  
  for (sample_id in common_samples) {
    df_samp <- df_spots[df_spots$sample == sample_id, , drop = FALSE]
    if (nrow(df_samp) == 0) next
    
    row_idx <- which(sample_table$sample_id == sample_id)[1]
    if (is.na(row_idx)) next
    
    hires_image_file  <- sample_table$hires_image[row_idx]
    lowres_image_file <- sample_table$lowres_image[row_idx]
    
    samp_label <- sample_table$sample_label[row_idx]
    samp_group <- sample_table$sample_group[row_idx]
    
    xh_min <- min(df_samp$x_highres, na.rm = TRUE)
    xh_max <- max(df_samp$x_highres, na.rm = TRUE)
    yh_min <- min(df_samp$y_highres, na.rm = TRUE)
    yh_max <- max(df_samp$y_highres, na.rm = TRUE)
    
    xl_min <- min(df_samp$x_lowres, na.rm = TRUE)
    xl_max <- max(df_samp$x_lowres, na.rm = TRUE)
    yl_min <- min(df_samp$y_lowres, na.rm = TRUE)
    yl_max <- max(df_samp$y_lowres, na.rm = TRUE)
    
    if (!is.na(hires_image_file) && file.exists(hires_image_file)) {
      img_hi <- imager::load.image(hires_image_file)
      df_hi <- as.data.frame(img_hi, wide = "c") %>%
        dplyr::mutate(rgb.val = rgb(c.1, c.2, c.3))
      
      df_hi$sample       <- sample_id
      df_hi$sample_label <- samp_label
      df_hi$sample_group <- samp_group
      df_hi$resolution   <- "highres"
      df_hi$contains_spots <- ifelse(
        df_hi$x >= (xh_min - img_buffer) & df_hi$x <= (xh_max + img_buffer) &
          df_hi$y >= (yh_min - img_buffer) & df_hi$y <= (yh_max + img_buffer),
        "yes", "no"
      )
      
      hi_list[[length(hi_list) + 1]] <- df_hi
    }
    
    if (!is.na(lowres_image_file) && file.exists(lowres_image_file)) {
      img_lo <- imager::load.image(lowres_image_file)
      df_lo <- as.data.frame(img_lo, wide = "c") %>%
        dplyr::mutate(rgb.val = rgb(c.1, c.2, c.3))
      
      df_lo$sample       <- sample_id
      df_lo$sample_label <- samp_label
      df_lo$sample_group <- samp_group
      df_lo$resolution   <- "lowres"
      df_lo$contains_spots <- ifelse(
        df_lo$x >= (xl_min - img_buffer) & df_lo$x <= (xl_max + img_buffer) &
          df_lo$y >= (yl_min - img_buffer) & df_lo$y <= (yl_max + img_buffer),
        "yes", "no"
      )
      
      lo_list[[length(lo_list) + 1]] <- df_lo
    }
  }
  
  images_highres <- if (length(hi_list) > 0) dplyr::bind_rows(hi_list) else NULL
  images_lowres  <- if (length(lo_list) > 0) dplyr::bind_rows(lo_list) else NULL
  
  list(
    highres = images_highres,
    lowres  = images_lowres
  )
}

# --------------------------------------------------------------
# 8) Aspect ratio helper based on fullres pixels
# --------------------------------------------------------------

compute_spatial_xy_ratio <- function(cds,
                                     sample_id,
                                     use_fullres = TRUE) {
  df <- as.data.frame(colData(cds))
  df <- df[df$sample == sample_id, , drop = FALSE]
  
  if (nrow(df) == 0) {
    stop("No spots found in cds for sample_id: ", sample_id)
  }
  
  if (use_fullres) {
    if (!all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(df))) {
      stop("pxl_col_in_fullres / pxl_row_in_fullres not found in colData(cds).")
    }
    x_range <- range(df$pxl_col_in_fullres, na.rm = TRUE)
    y_range <- range(df$pxl_row_in_fullres, na.rm = TRUE)
  } else {
    if (!all(c("x_highres", "y_highres") %in% colnames(df))) {
      stop("x_highres / y_highres not found in colData(cds).")
    }
    x_range <- range(df$x_highres, na.rm = TRUE)
    y_range <- range(df$y_highres, na.rm = TRUE)
  }
  
  (x_range[2] - x_range[1]) / (y_range[2] - y_range[1])
}

# --------------------------------------------------------------
# 9) Single sample H&E overlay
# --------------------------------------------------------------
# Generic per sample plot that works with either highres or lowres
# image data frames, by controlling which spot coordinates to use.

plot_visium_spatial_sample <- function(cds,
                                       image_df,
                                       sample_id,
                                       color_by        = "clusters",
                                       coord_type      = c("highres", "lowres"),
                                       point_size      = 0.5,
                                       point_alpha     = 0.8,
                                       img_buffer_only_spots = TRUE) {
  coord_type <- match.arg(coord_type)
  
  df_spots <- as.data.frame(colData(cds))
  df_spots <- df_spots[df_spots$sample == sample_id, , drop = FALSE]
  
  if (nrow(df_spots) == 0) {
    stop("No spots found in cds for sample_id: ", sample_id)
  }
  
  df_img <- image_df[image_df$sample == sample_id, , drop = FALSE]
  if (nrow(df_img) == 0) {
    stop("No image pixels found for sample_id: ", sample_id)
  }
  
  if (img_buffer_only_spots && "contains_spots" %in% colnames(df_img)) {
    df_img <- df_img[df_img$contains_spots == "yes", , drop = FALSE]
  }
  
  # Which spot coordinates to use
  if (coord_type == "highres") {
    x_col <- "x_highres"
    y_col <- "y_highres"
  } else {
    x_col <- "x_lowres"
    y_col <- "y_lowres"
  }
  
  if (!all(c(x_col, y_col) %in% colnames(df_spots))) {
    stop("Required coordinate columns ", x_col, ", ", y_col, " not found in colData(cds).")
  }
  
  xy_ratio <- compute_spatial_xy_ratio(cds, sample_id, use_fullres = TRUE)
  sample_label <- unique(df_spots$sample_label)
  
  ggplot() +
    geom_raster(
      data = df_img,
      aes(x = x, y = y, fill = rgb.val)
    ) +
    scale_fill_identity() +
    scale_y_reverse() +
    geom_point(
      data = df_spots,
      aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        colour = .data[[color_by]]
      ),
      size  = point_size,
      alpha = point_alpha
    ) +
    theme_void() +
    theme(
      aspect.ratio = 1 / xy_ratio,
      text         = element_text(size = 16, color = "black"),
      legend.text  = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 10, color = "black")
    ) +
    labs(
      title  = paste0("Sample: ", sample_label, " (", sample_id, ")"),
      colour = color_by
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4)))
}

# --------------------------------------------------------------
# Gene expression over H&E for a single sample
# --------------------------------------------------------------
# Arguments:
#   cds        : monocle3 cell_data_set
#   images     : list from build_visium_images() (images$highres / images$lowres)
#   gene.name  : gene_short_name to visualize
#   sample_id  : full sample ID (matches colData(cds)$sample)
#   resolution : "highres" or "lowres" image
#   point_size : point size for spots
#   point_alpha: transparency of the spots
#
# Returns:
#   ggplot object with image + gene expression overlay

plot_visium_gene_spatial_sample <- function(cds,
                                            images,
                                            gene.name,
                                            sample_id,
                                            resolution  = c("highres", "lowres"),
                                            point_size  = 0.5,
                                            point_alpha = 0.8,
                                            img_buffer_only_spots = TRUE) {
  resolution <- match.arg(resolution)
  
  # choose the appropriate image df
  image_df <- images[[resolution]]
  if (is.null(image_df)) {
    stop("No images found for resolution: ", resolution)
  }
  
  # subset spots for this sample
  df_spots <- as.data.frame(colData(cds))
  df_spots <- df_spots[df_spots$sample == sample_id, , drop = FALSE]
  if (nrow(df_spots) == 0) {
    stop("No spots found in cds for sample_id: ", sample_id)
  }
  
  # get expression vector for this gene across ALL spots
  idx <- fData(cds)$gene_short_name == gene.name
  if (!any(idx)) {
    stop("Gene not found in fData(cds)$gene_short_name: ", gene.name)
  }
  expr_vec <- exprs(cds)[idx, , drop = TRUE]  # named by colnames(cds)
  
  # align gene expression to just this sample's spots
  # rownames(df_spots) are barcodes, which match colnames(cds)
  expr_sample <- expr_vec[rownames(df_spots)]
  
  if (any(is.na(expr_sample))) {
    warning("NA values in expression for gene ", gene.name,
            " after subsetting to sample ", sample_id,
            ". Proceeding but check barcodes if this looks odd.")
  }
  
  df_spots[[gene.name]] <- log2(as.numeric(expr_sample) + 0.01)
  
  # pick coordinates based on resolution
  if (resolution == "highres") {
    x_col <- "x_highres"
    y_col <- "y_highres"
  } else {
    x_col <- "x_lowres"
    y_col <- "y_lowres"
  }
  
  if (!all(c(x_col, y_col) %in% colnames(df_spots))) {
    stop("Required coordinate columns ", x_col, ", ", y_col, " not found in colData(cds).")
  }
  
  # subset image pixels for this sample
  df_img <- image_df[image_df$sample == sample_id, , drop = FALSE]
  if (nrow(df_img) == 0) {
    stop("No image pixels found for sample_id: ", sample_id)
  }
  
  if (img_buffer_only_spots && "contains_spots" %in% colnames(df_img)) {
    df_img <- df_img[df_img$contains_spots == "yes", , drop = FALSE]
  }
  
  # aspect ratio from fullres pixels
  xy_ratio     <- compute_spatial_xy_ratio(cds, sample_id, use_fullres = TRUE)
  sample_label <- unique(df_spots$sample_label)
  
  ggplot() +
    geom_raster(
      data = df_img,
      aes(x = x, y = y, fill = rgb.val)
    ) +
    scale_fill_identity() +
    scale_y_reverse() +
    geom_point(
      data = df_spots,
      aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        colour = .data[[gene.name]]
      ),
      size  = point_size,
      alpha = point_alpha
    ) +
    scale_colour_gradient(low = "white", high = "red", na.value = "grey90") +
    theme_void() +
    theme(
      aspect.ratio = 1 / xy_ratio,
      text         = element_text(size = 16, color = "black"),
      legend.text  = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 10, color = "black")
    ) +
    labs(
      title  = paste0("Sample: ", sample_label, " (", sample_id, ") – ", gene.name),
      colour = paste0("log2(", gene.name, " + 0.01)")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4)))
}

# --------------------------------------------------------------
# 10) Multi sample convenience wrapper
# --------------------------------------------------------------
# Plot all samples on highres or lowres images and optionally save.

plot_all_samples_spatial <- function(cds,
                                     images,
                                     resolution = c("highres", "lowres"),
                                     sample_ids = NULL,
                                     color_by   = "clusters",
                                     output_dir = NULL,
                                     file_prefix = NULL,
                                     width = 6,
                                     height = 6) {
  resolution <- match.arg(resolution)
  
  image_df <- images[[resolution]]
  if (is.null(image_df)) {
    stop("No images found for resolution: ", resolution)
  }
  
  if (is.null(sample_ids)) {
    sample_ids <- sort(unique(as.character(colData(cds)$sample)))
  }
  
  if (is.null(file_prefix)) {
    file_prefix <- paste0("spatial_", resolution, "_")
  }
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  coord_type <- if (resolution == "highres") "highres" else "lowres"
  
  plot_list <- lapply(sample_ids, function(sid) {
    p <- plot_visium_spatial_sample(
      cds       = cds,
      image_df  = image_df,
      sample_id = sid,
      color_by  = color_by,
      coord_type = coord_type,
      point_size = ifelse(resolution == "highres", 0.5, 1),
      point_alpha = 0.8
    )
    
    if (!is.null(output_dir)) {
      out_file <- file.path(output_dir, paste0(file_prefix, sid, ".pdf"))
      ggsave(out_file, p, width = width, height = height)
    }
    
    p
  })
  
  names(plot_list) <- sample_ids
  plot_list
}

# --------------------------------------------------------------
# 11) Signature over H&E for a single sample
# --------------------------------------------------------------

# Purpose:
#   Generate a single spatial overlay plot for one sample, showing a
#   numeric gene-signature score projected onto a Visium highres or
#   lowres tissue image.
#
# Inputs:
#   cds              monocle3 cell_data_set with signature columns.
#   images           list returned by build_visium_images():
#                      list(highres = df, lowres = df)
#   signature_col    name of numeric signature column in colData(cds).
#   sample_id        sample to plot (matches colData(cds)$sample).
#   resolution       "highres" or "lowres".
#   palette          "diverging", "sequential", "viridis", "custom".
#   diverging_midpoint   midpoint for diverging scale (default 0).
#   sequential_low       lower color for sequential scale.
#   sequential_high      upper color for sequential scale.
#   viridis_opt          viridis palette option ("magma", "inferno", ...).
#   custom_cols          character vector of custom colors.
#   point_size           spot size.
#   point_alpha          spot transparency.
#
# Output:
#   A ggplot object.
#
# Example:
#   p <- plot_signature_spatial_sample(
#          cds           = cds,
#          images        = imgs,
#          signature_col = "HALLMARK_INFLAMMATORY_RESPONSE",
#          sample_id     = "Sample_9874-WS-2",
#          resolution    = "highres",
#          palette       = "diverging",
#          diverging_midpoint = 0
#        )
#   p
#
#   # Sequential palette example:
#   plot_signature_spatial_sample(
#       cds           = cds,
#       images        = imgs,
#       signature_col = "GOBP_FIBROBLAST_ACTIVATION",
#       sample_id     = "Sample_9874-WS-4",
#       palette       = "sequential",
#       sequential_low  = "white",
#       sequential_high = "red"
#   )
#

plot_signature_spatial_sample <- function(cds,
                                          images,
                                          signature_col,
                                          sample_id,
                                          resolution  = c("highres", "lowres"),
                                          point_size  = 0.5,
                                          point_alpha = 0.8,
                                          img_buffer_only_spots = TRUE,
                                          palette      = c("diverging", "sequential", "viridis", "custom"),
                                          custom_cols  = NULL,
                                          viridis_opt  = "magma",
                                          diverging_midpoint = 0) {
  resolution <- match.arg(resolution)
  palette    <- match.arg(palette)
  
  # pick image data frame for this resolution
  image_df <- images[[resolution]]
  if (is.null(image_df)) {
    stop("No images available for resolution: ", resolution)
  }
  
  # subset spots for this sample
  df_spots <- as.data.frame(colData(cds))
  df_spots <- df_spots[df_spots$sample == sample_id, , drop = FALSE]
  if (nrow(df_spots) == 0) {
    stop("No spots found in cds for sample_id: ", sample_id)
  }
  
  # signature column must exist and be numeric
  if (!signature_col %in% colnames(df_spots)) {
    stop("signature_col '", signature_col, "' not found in colData(cds).")
  }
  if (!is.numeric(df_spots[[signature_col]])) {
    stop("signature_col '", signature_col, "' is not numeric. It should be a numeric signature score.")
  }
  
  # spatial coordinates based on resolution
  if (resolution == "highres") {
    x_col <- "x_highres"
    y_col <- "y_highres"
  } else {
    x_col <- "x_lowres"
    y_col <- "y_lowres"
  }
  if (!all(c(x_col, y_col) %in% colnames(df_spots))) {
    stop("Required coordinate columns ", x_col, ", ", y_col, " not found in colData(cds).")
  }
  
  # subset image pixels for this sample
  df_img <- image_df[image_df$sample == sample_id, , drop = FALSE]
  if (nrow(df_img) == 0) {
    stop("No image pixels found for sample_id: ", sample_id)
  }
  
  if (img_buffer_only_spots && "contains_spots" %in% colnames(df_img)) {
    df_img <- df_img[df_img$contains_spots == "yes", , drop = FALSE]
  }
  
  # aspect ratio using full resolution pixel ranges
  xy_ratio     <- compute_spatial_xy_ratio(cds, sample_id, use_fullres = TRUE)
  sample_label <- unique(df_spots$sample_label)
  
  # choose continuous color scale for numeric signature
  color_scale <- switch(
    palette,
    "sequential" = scale_colour_gradient(
      low      = "white",
      high     = "red",
      na.value = "grey80"
    ),
    "diverging" = scale_colour_gradient2(
      low       = "blue",
      mid       = "grey90",
      high      = "red",
      midpoint  = diverging_midpoint,
      na.value  = "grey80"
    ),
    "viridis" = scale_colour_viridis_c(
      option   = viridis_opt,
      na.value = "grey80"
    ),
    "custom" = {
      if (is.null(custom_cols) || length(custom_cols) < 2) {
        stop("For palette = 'custom', supply at least 2 colors in custom_cols.")
      }
      scale_colour_gradientn(
        colors   = custom_cols,
        na.value = "grey80"
      )
    }
  )
  
  ggplot() +
    geom_raster(
      data = df_img,
      aes(x = x, y = y, fill = rgb.val)
    ) +
    scale_fill_identity() +
    scale_y_reverse() +
    geom_point(
      data = df_spots,
      aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        colour = .data[[signature_col]]
      ),
      size  = point_size,
      alpha = point_alpha
    ) +
    color_scale +
    guides(colour = guide_colorbar()) +
    theme_void() +
    theme(
      aspect.ratio = 1 / xy_ratio,
      text         = element_text(size = 16, color = "black"),
      legend.text  = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 10, color = "black")
    ) +
    labs(
      title  = paste0("Sample: ", sample_label, " (", sample_id, ")"),
      colour = signature_col
    )
}

# --------------------------------------------------------------
# 12) Signature over H&E for all samples
# --------------------------------------------------------------
# Purpose:
#   Wrapper around plot_signature_spatial_sample().
#   Produces one spatial signature plot per sample and optionally saves
#   them to disk.
#
# Inputs:
#   cds              monocle3 cell_data_set.
#   images           list(highres = df, lowres = df).
#   signature_col    signature column name in colData(cds).
#   resolution       "highres" or "lowres".
#   palette          "diverging", "sequential", "viridis", "custom".
#   diverging_midpoint   midpoint for diverging palette.
#   sequential_low       lower color.
#   sequential_high      upper color.
#   viridis_opt          viridis palette choice.
#   custom_cols          custom color vector.
#   point_size           spot size.
#   point_alpha          transparency.
#   output_dir           directory to save plots (NULL = no saving).
#   file_prefix          filename prefix when saving.
#
# Output:
#   Named list of ggplot objects:
#      $Sample_9874-WS-2
#      $Sample_9874-WS-4
#      ...
#
# Example:
#   plots_hallmark_inflam <- plot_all_samples_signature_spatial(
#       cds            = cds,
#       images         = imgs,
#       signature_col  = "HALLMARK_INFLAMMATORY_RESPONSE",
#       resolution     = "highres",
#       palette        = "diverging",
#       output_dir     = "plots/signatures_highres",
#       file_prefix    = "HALLMARK_INFLAMMATORY_"
#   )
#
#   # Inspect one result:
#   plots_hallmark_inflam[["Sample_9874-WS-2"]]
#
plot_all_samples_signature_spatial <- function(cds,
                                               images,
                                               signature_col,
                                               resolution  = c("highres", "lowres"),
                                               sample_ids  = NULL,
                                               output_dir  = NULL,
                                               file_prefix = NULL,
                                               width       = 6,
                                               height      = 6,
                                               point_size  = 0.5,
                                               point_alpha = 0.8,
                                               img_buffer_only_spots = TRUE,
                                               palette      = c("diverging", "sequential", "viridis", "custom"),
                                               custom_cols  = NULL,
                                               viridis_opt  = "magma",
                                               diverging_midpoint = 0) {
  resolution <- match.arg(resolution)
  palette    <- match.arg(palette)
  
  # basic check for images at this resolution
  image_df <- images[[resolution]]
  if (is.null(image_df)) {
    stop("No images available for resolution: ", resolution)
  }
  
  # default: use all sample IDs from cds
  if (is.null(sample_ids)) {
    sample_ids <- sort(unique(as.character(colData(cds)$sample)))
  }
  
  # default file prefix if saving
  if (is.null(file_prefix)) {
    file_prefix <- paste0("signature_", signature_col, "_", resolution, "_")
  }
  
  # optional output directory
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  plot_list <- lapply(sample_ids, function(sid) {
    p <- plot_signature_spatial_sample(
      cds           = cds,
      images        = images,
      signature_col = signature_col,
      sample_id     = sid,
      resolution    = resolution,
      point_size    = point_size,
      point_alpha   = point_alpha,
      img_buffer_only_spots = img_buffer_only_spots,
      palette       = palette,
      custom_cols   = custom_cols,
      viridis_opt   = viridis_opt,
      diverging_midpoint = diverging_midpoint
    )
    
    if (!is.null(output_dir)) {
      out_file <- file.path(output_dir, paste0(file_prefix, sid, ".pdf"))
      ggsave(out_file, p, width = width, height = height)
    }
    
    p
  })
  
  names(plot_list) <- sample_ids
  plot_list
}


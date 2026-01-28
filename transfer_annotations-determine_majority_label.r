# --------------------------------------------------------------
# Transfer cell-level annotations between monocle3 CDS objects
# --------------------------------------------------------------
# Input: reference CDS object with celltype annotations
#          target CDS object
#
# Output: target CDS object with transferred celltype annotations
#
# Usage:
#.  cds2 <- transfer_cell_annotations(
#.    cds_ref    = cds1,
#.    cds_target = cds2,
#.    annotation_col = "cell_type",
#.    cell_id_col = "cell"   # set to NULL to use rownames
#.  )


transfer_cell_annotations <- function(
    cds_ref,  # reference (old) dataset
    cds_target,   # target (new) dataset
    annotation_col,    #column name where original annotations are written
    new_annotation_col = NULL,    #column name to be added to the target/new dataset
    cell_id_col = NULL,
    verbose = TRUE
) {
  # Extract pData
  pd_ref <- as.data.frame(pData(cds_ref))
  pd_target <- as.data.frame(pData(cds_target))
  
  # Determine output column name
  if (is.null(new_annotation_col)) {
    new_annotation_col <- annotation_col
  }
  
  # Determine cell identifiers
  if (!is.null(cell_id_col)) {
    if (!cell_id_col %in% colnames(pd_ref) ||
        !cell_id_col %in% colnames(pd_target)) {
      stop("cell_id_col not found in both cds_ref and cds_target pData")
    }
    ref_ids <- pd_ref[[cell_id_col]]
    target_ids <- pd_target[[cell_id_col]]
  } else {
    ref_ids <- rownames(pd_ref)
    target_ids <- rownames(pd_target)
  }
  
  # Check annotation column exists in reference
  if (!annotation_col %in% colnames(pd_ref)) {
    stop("annotation_col not found in cds_ref pData")
  }
  
  # Warn if overwriting existing column
  if (new_annotation_col %in% colnames(pd_target)) {
    warning(
      "Column '", new_annotation_col,
      "' already exists in cds_target pData and will be overwritten"
    )
  }
  
  # Create mapping table
  mapping <- data.frame(
    cell_id = ref_ids,
    annotation = pd_ref[[annotation_col]],
    stringsAsFactors = FALSE
  )
  
  # Match target cells to reference
  match_idx <- match(target_ids, mapping$cell_id)
  transferred <- mapping$annotation[match_idx]
  
  # Assign annotation
  pData(cds_target)[[new_annotation_col]] <- transferred
  
  # Diagnostics
  if (verbose) {
    n_total <- length(target_ids)
    n_matched <- sum(!is.na(match_idx))
    n_unmatched <- sum(is.na(match_idx))
    
    message("Annotation transfer summary:")
    message("  Reference annotation: ", annotation_col)
    message("  Target column name:   ", new_annotation_col)
    message("  Total target cells:   ", n_total)
    message("  Matched cells:        ", n_matched)
    message("  Unmatched cells:      ", n_unmatched)
    
    if (n_unmatched > 0) {
      message("  Unmatched cells assigned NA")
    }
  }
  
  return(cds_target)
}



# --------------------------------------------------------------
# Determine majority annotation per cluster
# --------------------------------------------------------------
# Input:  CDS object with transfered celltype annotations
#         must specify which columns contain cluster, celltype label
#
# Output: CDS object with column majority cluster label
#
#Usage:
# cluster_labels <- majority_label_by_cluster(
# cds,
# annotation_col = "cell_type",
# cluster_col = "cluster",
# min_prop = 0.6
# )
# 
# # Or directly annotate the CDS
# cluster_labels <- majority_label_by_cluster(
#   cds,
#   annotation_col = "cell_type",
#   annotate_cds = TRUE
# )

majority_label_by_cluster <- function(
    cds,
    annotation_col,
    cluster_col = "cluster",
    min_prop = 0,
    annotate_cds = FALSE
) {
  pd <- as.data.frame(pData(cds))
  
  if (!annotation_col %in% colnames(pd)) {
    stop("annotation_col not found in pData(cds)")
  }
  if (!cluster_col %in% colnames(pd)) {
    stop("cluster_col not found in pData(cds)")
  }
  
  # Remove NA annotations
  pd <- pd[!is.na(pd[[annotation_col]]), ]
  
  # Tabulate
  tab <- table(
    cluster = pd[[cluster_col]],
    annotation = pd[[annotation_col]]
  )
  
  # Compute majority label
  majority_df <- do.call(rbind, lapply(rownames(tab), function(cl) {
    counts <- tab[cl, ]
    prop <- counts / sum(counts)
    top_label <- names(which.max(counts))
    top_prop <- max(prop)
    
    if (top_prop < min_prop) {
      top_label <- NA
    }
    
    data.frame(
      cluster = cl,
      majority_label = top_label,
      majority_prop = top_prop,
      stringsAsFactors = FALSE
    )
  }))
  
  # Optionally annotate CDS
  if (annotate_cds) {
    cluster_to_label <- setNames(
      majority_df$majority_label,
      majority_df$cluster
    )
    pData(cds)$majority_label <- cluster_to_label[
      as.character(pData(cds)[[cluster_col]])
    ]
  }
  
  return(majority_df)
}
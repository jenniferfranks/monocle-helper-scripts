facet_gene_dotplot_group_by_Timepoint <- function(cds,
                                markers,
                                group_cells_by="clusters",
                                reduction_method = "UMAP",
                                norm_method = c("log", "size_only"),
                                lower_threshold = 0,
                                max.size = 10,
                                axis_order = c('group_marker', 'marker_group'),
                                flip_percentage_mean = FALSE,
                                pseudocount = 1,
                                scale_max = 3,
                                scale_min = -3) {
  gene_ids = as.data.frame(fData(cds)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname %in% markers | gene_short_name %in% markers) %>%
    dplyr::pull(rowname)
  
  major_axis <- 1
  minor_axis <- 2
  
  
  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  exprs_mat$Timepoint <- pData(cds)$Timepoint
  colnames(exprs_mat) <- c('Cell', 'Gene',  'Expression', 'Timepoint')
  exprs_mat$Gene <- as.character(exprs_mat$Gene)
  cell_group <- colData(cds)[,group_cells_by]
  
  names(cell_group) = colnames(cds)
  
  
  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene, Timepoint) %>%
    dplyr::summarize(mean = mean(log(Expression + pseudocount)),
                     percentage = sum(Expression > lower_threshold) /
                       length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, ExpVal$mean)
  
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, 'gene_short_name']
  
  
  
  g <- ggplot(ExpVal, aes(y = Group,  x = Timepoint)) +
    geom_point(aes(colour = mean,  size = percentage)) +
    viridis::scale_color_viridis(name = 'log(mean + 0.1)') +
    scale_size(name = 'percentage', range = c(0, max.size)) + facet_wrap(~Gene) + 
    theme_bw() + xlab("Timepoint") + ylab(group_cells_by)
  
  g
  
  
}



 
facet_gene_dotplot_group_by_Timepoint2 <- function(cds,
                                                  markers,
                                                  group_cells_by="clusters",
                                                  reduction_method = "UMAP",
                                                  norm_method = c("log", "size_only"),
                                                  lower_threshold = 0,
                                                  max.size = 10,
                                                  axis_order = c('group_marker', 'marker_group'),
                                                  flip_percentage_mean = FALSE,
                                                  pseudocount = 1,
                                                  scale_max = 3,
                                                  scale_min = -3) {
  gene_ids = as.data.frame(fData(cds)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname %in% markers | gene_short_name %in% markers) %>%
    dplyr::pull(rowname)
  
  major_axis <- 1
  minor_axis <- 2
  
  
  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  exprs_mat$Timepoint <- pData(cds)$Timepoint
  colnames(exprs_mat) <- c('Cell', 'Gene',  'Expression', 'Timepoint')
  exprs_mat$Gene <- as.character(exprs_mat$Gene)
  cell_group <- colData(cds)[,group_cells_by]
  
  names(cell_group) = colnames(cds)
  
  
  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene, Timepoint) %>%
    dplyr::summarize(mean = mean(log(Expression + pseudocount)),
                     percentage = sum(Expression > lower_threshold) /
                       length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, ExpVal$mean)
  
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, 'gene_short_name']
  
  
  
  g <- ggplot(ExpVal, aes(y = Gene,  x = Timepoint)) +
    geom_point(aes(colour = mean,  size = percentage)) +
    viridis::scale_color_viridis(name = 'log(mean + 0.1)') +
    scale_size(name = 'percentage', range = c(0, max.size)) + #facet_wrap(~Gene) + 
    theme_bw() + xlab("Timepoint") + ylab(group_cells_by)
  
  g
  
  
}

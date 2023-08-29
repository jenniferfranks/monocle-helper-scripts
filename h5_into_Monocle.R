library(hdf5r)
library(Matrix)


read.cds.cellranger.h5.file = function(h5.file) {
  if(!file.exists(h5.file)){stop(paste0("File ", h5.file," not found"))}
  #s<-h5file(h5.file, mode="r")
  s <- H5File$new(h5.file, mode="r")
  
  #s$ls(recursive=T)
  
  barcodes = s[["matrix/barcodes"]][]
  gene_ids = s[["matrix/features/id"]][]
  gene_names =s[["matrix/features/name"]][]
  featuretype = s[["matrix/features/feature_type"]][]
  data = s[["matrix/data"]][]
  indices = s[["matrix/indices"]][]+1
  indptr = s[["matrix/indptr"]][]
  shape = s[["matrix/shape"]][]
  #browser()
  h5close(s)
  # gbm = new(
  #   "dgCMatrix",
  #   x = data, i = indices, p = indptr,
  #   Dim = shape)
  if(length(levels(factor(featuretype)))>1){
    warning("Warning - Multiple feature types found")
  }
  
  gbm = sparseMatrix(x = data, i = indices, p = indptr,
                     dims = shape)
  
  
  pData.df = data.frame(
    barcode = barcodes,
    stringsAsFactors = F
  )
  
  fData.df = data.frame(
    id = gene_ids,
    gene_short_name = gene_names,
    feature_type = featuretype,
    stringsAsFactors = F
  )
  
  rownames(pData.df) = pData.df$barcode
  rownames(fData.df) = fData.df$id
  
  rownames(gbm) = rownames(fData.df)
  colnames(gbm) = rownames(pData.df)
  
  suppressWarnings({res = new_cell_data_set(
    gbm,
    cell_metadata = pData.df,
    gene_metadata = fData.df
  )})
  
  return(res)
}
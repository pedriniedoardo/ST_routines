# -------------------------------------------------------------------------
# read in the h5ad file and convert it to a Seurat object

h5ad_to_seurat = function(file_in) {
  
  stopifnot(require(rhdf5))
  stopifnot(require(Matrix))
  stopifnot(require(Seurat))
  
  cell_counts = h5read(file = file_in, 'X')
  cell_obs = h5read(file = file_in, 'obs')
  cell_var = h5read(file = file_in, 'var')
  
  n_rows <- length(cell_counts$indptr) - 1  # Number of rows in the matrix
  n_cols <- max(cell_counts$indices) + 1    # Number of columns (assuming 0-based indices)
  
  i = rep(1:n_rows, times = diff(cell_counts$indptr))
  j = cell_counts$indices + 1
  x = as.numeric(cell_counts$data)
  
  dcg_mat = sparseMatrix(
    i = i, j = j, x = as.numeric(x), dims = c(n_rows, n_cols),
    dimnames = list(
      as.character(cell_obs$`_index`),
      as.character(cell_var$`_index`)
    ))
  
  CreateSeuratObject(CreateAssayObject(counts = t(dcg_mat)), meta.data = cell_obs)
  
}

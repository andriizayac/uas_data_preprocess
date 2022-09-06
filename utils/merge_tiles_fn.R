# helper function
merge_tiles <- function(file_path, output_path) {
  # merge multiple laz/las tiles into a single large file
  # input 1: a path to where the tiles are, will use all files in the dir
  # input 2: a path, including the name without extension, where the merged file will go 
  # output: a single .las file
  
  pkgs <- library("lidR", logical.return = TRUE, quietly = TRUE)
  if( pkgs == FALSE ) {
    stop(cat("Package `lidR` is not available"))
  }
  ctg <- readLAScatalog(file_path)
  
  opt_output_files(ctg) <- output_path
  opt_chunk_buffer(ctg) <- 0
  opt_chunk_size(ctg) <- 1e6
  singlefile_ctg <- catalog_retile(ctg)
  
  rm(ctg, singlefile_ctg)
  
}

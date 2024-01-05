#'Read counts from multiple files 
#'
#'date 06/01/2023
#'@param dir Directory of count files
#'
#'@return Data Frame
#'@export

generate_counts <- function(count_files_dir, file_ext_to_remove, remove_last_rows=FALSE,rows_to_remove=4, sep = '\t') {
  
  countfiles <- list.files(path = count_files_dir,
                           pattern = '.txt')
  raw_counts_list <- list()
  #Read your expression data files and make a list
  for (f in countfiles ) {
    x <- as.data.frame( read.table(paste0(count_files_dir,'/',f), row.names=1, sep = sep) )
    raw_counts_list[[ length( raw_counts_list ) + 1]] <- x
    row.names( raw_counts_list[[ length(raw_counts_list) ]] ) <-  row.names(x)
  }
  
  ##Convert datalist to 1 table
  data <- as.data.frame(raw_counts_list)
  
  if (!is.null(file_ext_to_remove)) {
    names(data) <- gsub(file_ext_to_remove,'',countfiles)
  } else {
    names(data) <- countfiles
  }
  
  if (remove_last_rows){
    data <- data[-c((dim(data)[1]-rows_to_remove):dim(data)[1]),] 
  }
 data 
}

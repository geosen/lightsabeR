#'Isolate a subset of genes from a gene expression table
#'
#'date 08/02/2023
#'@param ensg_ids Vector of Ensembl IDs with version number
#'
#'
#'@return Vector of Ensembl IDs without version number
#'@export
#'
remove_version <- function(ensg_ids) {
  gsub('\\.[[:digit:]]{1+}$','',ensg_ids)
}

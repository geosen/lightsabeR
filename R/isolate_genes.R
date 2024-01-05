#'Isolate a subset of genes from a gene expression table
#'
#'date 08/02/2023
#'@param object Gene expression table
#'@param gene_list Vector with gene names to be isolated
#'@param gene_column Column in the gene expression table that contains the gene names
#'
#'
#'@return Same class as object but subsetted
#'@export
isolate_genes <- function(object,
                          gene_list,
                          gene_column = "rownames")
{
  

  if(gene_column == 'rownames') {
    object[rownames(object) %in% gene_list,]
  } else if(gene_column != 'rownames' & gene_column %in% names(object)) {
    object[object[,gene_column] %in% gene_list,]
  }

}

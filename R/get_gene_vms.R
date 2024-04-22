#' Get gene Variance, Mean and Sum.
#' 
#' George Sentis
#'
#'date 31/01/2024
#'
#' @param df Data frame containing gene expression with gene names as rownames.
#' @param genes The list of genes to gather Variance, Mean and Sum.
#'
#' @return A Data frame containing gene Variance, Mean and Sum.
#' @export
#'
#'
#'
#'
#'

get_gene_vms <- function(df,genes) {
  
  geneinfo <- data.frame( vars = unlist(sapply(as.data.frame(t(df[genes,])),var)),
                          means = unlist(sapply(as.data.frame(t(df[genes,])),mean)),
                          sums = unlist(sapply(as.data.frame(t(df[genes,])),sum)))
  
}


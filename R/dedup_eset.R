#'Deduplicate probes of Expression Set Microarray object
#'
#'This function takes an Expression Set object, finds probes that correspond to the same gene 
#'according to fData, and deduplicates.
#'
#'Deduplication is performed based on higher expression and if that is a tie then 
#'it chooses higher variation.
#'If only 1 probe has a mean expression over 1, removes all others. 
#'If 2 or more probes have expression over 1, or 2 or all probes have expression 
#'lower than 1, then it keeps the one with the highest variance.
#'
#'Dependencies are:
#' library(Biobase)
#'
#'date 10/07/2023
#'@param eset Input ExpressionSet object 
#'@param fData_gene_symbol_column_name The name of the column containing the gene symbols in feature Data of the ExpressionSet object 
#'
#'@return An ExpressionSet object with 1 probe per gene
#'@export

dedup_eset <- function(eset, fData_gene_symbol_column_name = 'Gene Symbol') {
  
  table <-  as.data.frame(exprs(eset))
  table$genes <- fData(eset)[,fData_gene_symbol_column_name]
  double_genes <- fData(eset)[,fData_gene_symbol_column_name][duplicated(fData(eset)[,fData_gene_symbol_column_name])]
  genes_to_remove <- c()
  for (i in 1:length(double_genes)) {
    dg <- double_genes[i]
    possible_indices <- which(table$genes == dg)
    dedup <- table[table$genes == dg, colnames(table) != "genes"]
    tdedup <- as.data.frame(t(dedup))
    dmeans <- data.frame(colMeans(tdedup, na.rm = TRUE))
    dvars <- data.frame(lapply(tdedup, var, na.rm = TRUE))
    
    #choosing which genes to remove
    if (sum(dmeans > 1) == 1) {
      genes_to_remove <- c(genes_to_remove, possible_indices[dmeans < 
                                                               1])
    } else if (sum(!is.na(dvars)) > 0 & length(possible_indices[dvars == 
                                                                max(dvars)]) == 1) {
      genes_to_remove <- c(genes_to_remove, possible_indices[dvars != 
                                                               max(dvars)])
    } else if ((sum(dmeans) != 0) & (sum(dmeans == max(dmeans)) == 
                                     1)) {
      genes_to_remove <- c(genes_to_remove, possible_indices[dmeans != 
                                                               max(dmeans)])
    } else {
      genes_to_remove <- c(genes_to_remove, possible_indices[2:length(possible_indices)])
    }
  }
  
  #removing genes
  eset2 <- eset[-genes_to_remove]
  eset2
  
}
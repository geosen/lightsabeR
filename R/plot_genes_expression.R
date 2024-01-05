#'Plot histograms of expression of specific genes
#'
#'Dependencies are:
#' library(ggplot2)
#' library(tidyverse)
#'
#'
#'date 10/02/2023
#'@param expression Object with expression values. The row names of the object must be the gene identifiers.
#'@param genes The genes that will be plotted. Identifiers must match some row names of the expression object.
#'
#'@return Data Frame
#'@export

plot_genes_expression <-  function(expression,
                                        genes){
#loading libraries  
library(ggplot2)
library(tidyverse)
    
  
gene <- rownames(expression)[match(genes, rownames(expression))]  
    
if(is.null(gene)){
      errorCondition('The genes you provided cannot be found in the expression table. 
                     Make sure the identifiers in the expression table and the gene list are in a common format and same organism')
    } else {
      
      #isolating expression
      gene_expression <- isolate_genes(expression, gene_list = gene)
      
      #reformating table for plotting
      gene_expression <- as.data.frame(t(gene_expression)) %>%
        rownames_to_column() %>%
        gather(key = rowname)
      
      #plotting expression
      ggplot(gene_expression, aes(value)) + 
        geom_histogram() + 
        labs(x = 'Expression',y = 'Number of samples')+
        theme_minimal() + 
        facet_wrap(~rowname)
    } 
    
}
  





  
  

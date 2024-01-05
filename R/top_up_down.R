#'Choose a number of up and down regulated genes based on logFC from a certain comparison.
#'
#'Dependencies are:
#' library(dplyr)
#'
#'date 21/07/2023
#'@param lrt_table
#'@param order_factor 
#'@param ntop 
#'@param ntop_up
#'@param ntop_down
#'
#'@return A vector or data.frame(when both ntop_up and ntop_down are used)
#'@export
#'
#'

top_up_down <- function(lrt_table,
                                order_factor = 'logFC',
                                ntop = 100,
                                ntop_up = NULL,
                                ntop_down = NULL,
                                gene_id_column = NULL) {
  if(!is.null(gene_id_column)){
    rownames(lrt_table) <- NULL
    lrt_table <- lrt_table %>%
      column_to_rownames(gene_id_column)
  }
  
  if(is.null(ntop_up) & is.null(ntop_down)) {
    #selecting ntop up and downregulated genes
    upgenes <- lrt_table %>%
        slice_max(order_by = get(order_factor), n = ntop)
    topgenes <- lrt_table %>%
        slice_min(order_by = get(order_factor), n = ntop) %>% 
        bind_rows(upgenes)
    
    #returning gene list
      genes <- rownames(topgenes)
      genes
      
  }  else if(!is.null(ntop_up) & is.null(ntop_down)) {
    ntop_up <- as.numeric(ntop_up)
    #selecting upregulated genes
    upgenes <- lrt_table %>%
      slice_max(order_by = get(order_factor), n = ntop_up)
    
    #returning upgene list
    genes <- rownames(upenes)
    genes
  } else if(is.null(ntop_up) & !is.null(ntop_down)) {
    ntop_down <- as.numeric(ntop_down)
    
    #selecting downregulated genes
    downgenes <- lrt_table %>%
      slice_min(order_by = get(order_factor), n = ntop_down)
    
    genes <- rownames(downgenes)
    genes
  } else {
    ntop_up <- as.numeric(ntop_up)
    ntop_down <- as.numeric(ntop_down)
    
    #selecting ntop_up upregulated genes and ntop_down downregulated genes
    upgenes <- lrt_table %>%
      slice_max(order_by = get(order_factor), n = ntop_up) %>%
      mutate(DE = 'upregulated')
    topgenes <- lrt_table %>%
      slice_min(order_by = get(order_factor), n = ntop_down) %>% 
      mutate(DE = 'downregulated') %>%
      bind_rows(upgenes) %>%
      rownames_to_column('Genes') %>%
      dplyr::select(Genes, DE)
    
    #returning gene list
    topgenes
    
  }
    

  
  }



#'
#'
#'Dependencies are:
#' none
#' 
#'
#'date 21/07/2023
#'@param lrt_table
#'@param logFC_col_name 
#'@param sig_measure_col_name 
#'@param remove_version
#'@param gene_id_column 
#'@param output_table_filename
#'@param return_data
#'
#'@return Creates an rnk file in the directory (current or specified)
#'@export
#'
#'
create_rnk <- function(lrt_table, 
                       logFC_col_name = 'logFC' , 
                       sig_measure_col_name = 'fdr', 
                       remove_version = FALSE,
                       gene_id_column = NULL,
                       output_table_filename = NULL,
                       return_data = FALSE) {

  if(is.null(gene_id_column)) {
  DEGlist <- lrt_table 
  DEGlist$genes <- rownames(DEGlist)
  } else {
    DEGlist <- lrt_table
  }
  
  #calculate logFc*-log10FDR (or pvalue)
  DEGlist$log_FC_p<-DEGlist[,logFC_col_name]*-log10(DEGlist[,sig_measure_col_name])
  
  #isolate columns of interest and sort based on product
  if(!is.null(gene_id_column)) {
    nDEGlist<-DEGlist[c(as.character(gene_id_column), ("log_FC_p"))]
    sorted_nDEGlist<-nDEGlist[with(nDEGlist, order(nDEGlist$log_FC_p, decreasing = TRUE)),]  
  } else {
    nDEGlist<-DEGlist[c('genes', ("log_FC_p"))]
    sorted_nDEGlist<-nDEGlist[with(nDEGlist, order(nDEGlist$log_FC_p, decreasing = TRUE)),]
  }
  
  if(remove_version) {
    sorted_nDEGlist$genes <- gsub('\\.[[:digit:]]+$','',sorted_nDEGlist$genes)
  }
  
  
  #output
  if(return_data) {
    sorted_nDEGlist
  } else {
    
    if(is.null(output_table_filename)) {
      output_table_filename = paste0(as.numeric(Sys.time()),'_R_generated.rnk')
    }   
    
    write.table(sorted_nDEGlist,
              file =  output_table_filename, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE ) #output filepath here
    }
  }

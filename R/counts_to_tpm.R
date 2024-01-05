#'Convert ENSG_ID rownames of a table to HGNC_ID IDs
#'
#'This function takes a table with Ensembl Gene IDs or HGNC symbolsn and converts the raw expression values to HGNC symbols TPM 
#'This function does not take into account effective gene length. Maybe you want to look more into this for scientific and publication purposes.
#'Dependencies are:
#' library(biomaRt)
#'
#'date 30/04/2023
#'@param df Input data frame with ENSG IDs or HGNC symbols as rownames
#'@param ensembl_version The version of ensembl to be used. If none is provided the function uses the default (latest) version
#'@param organism The organism to search for: human or mouse
#'
#'@return A dataframe with TPM expression values
#'@export




counts_to_tpm <- function(df,
                          identifier_type = c('hgnc_symbol','ensembl_gene_id','mgi_symbol'),
                          ensembl_version = NULL,
                          gene_lengths = NULL, 
                          organism = c('human','mouse')){
  
  
  
  library(biomaRt)
  
  if(organism == 'human'){
    
    if(is.null(ensembl_version)) {
      mart = useEnsembl(biomart = 'ensembl',dataset = 'hsapiens_gene_ensembl')
      gene_list <- getBM(mart = mart, attributes=c('ensembl_gene_id','hgnc_symbol','entrezgene','start_position','end_position','description'))
      gene_list$gene_length = abs(gene_list$start_position-gene_list$end_position)
      
      if(identifier_type == 'hgnc_symbol'){
        df_rpk <- df/gene_list$gene_length[match(rownames(df),gene_list$hgnc_symbol)]
      } else if (identifier_type == 'ensembl_gene_id') {
        df_rpk <- df/gene_list$gene_length[match(rownames(df),gene_list$ensembl_gene_id)]
      }
      
    } else {
      mart = useEnsembl(biomart = 'ensembl',dataset = 'hsapiens_gene_ensembl', version = ensembl_version)
      gene_list <- getBM(mart = mart, attributes=c('ensembl_gene_id_version','ensembl_gene_id','hgnc_symbol','entrezgene','start_position','end_position','description'))
      gene_list$gene_length = abs(gene_list$start_position-gene_list$end_position)
     
      if(identifier_type == 'hgnc_symbol'){
        df_rpk <- df/gene_list$gene_length[match(rownames(df),gene_list$hgnc_symbol)]
      } else if (identifier_type == 'ensembl_gene_id') {
        df_rpk <- df/gene_list$gene_length[match(rownames(df),gene_list$ensembl_gene_id_version)]
      } 
      
    }
  }
  
  if(organism == 'mouse'){
    if(is.null(ensembl_version)) {
      mart = useEnsembl(biomart = 'ensembl',dataset = 'mmusculus_gene_ensembl')
      gene_list <- getBM(mart = mart, attributes=c('ensembl_gene_id','mgi_symbol','entrezgene','start_position','end_position','description'))
      gene_list$gene_length = abs(gene_list$start_position-gene_list$end_position)
      
      if(identifier_type == 'mgi_symbol'){
        df_rpk <- df/gene_list$gene_length[match(rownames(df),gene_list$mgi_symbol)]
      } else if (identifier_type == 'ensembl_gene_id') {
        df_rpk <- df/gene_list$gene_length[match(rownames(df),gene_list$ensembl_gene_id)]
      }
      
    } else {
      mart = useEnsembl(biomart = 'ensembl',dataset = 'mmusculus_gene_ensembl', version = ensembl_version)
      gene_list <- getBM(mart = mart, attributes=c('ensembl_gene_id_version','mgi_symbol','entrezgene','start_position','end_position','description'))
      gene_list$gene_length = abs(gene_list$start_position-gene_list$end_position)
      
      if(identifier_type == 'mgi_symbol'){
        df_rpk <- df/gene_list$gene_length[match(rownames(df),gene_list$mgi_symbol)] 
      } else if (identifier_type == 'ensembl_gene_id') {
        df_rpk <- df/gene_list$gene_length[match(rownames(df),gene_list$ensembl_gene_id_version)]
      }
      
    }
    
  }
  
  #divide every column of the rpk data frame with a "per milliion" normalizing factor
  df_tpm <- sweep(df_rpk,2,colSums(df_rpk, na.rm = TRUE)/10^6,'/')
  #return the data frame
  df_tpm
}
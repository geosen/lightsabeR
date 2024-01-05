#'Convert ENSG_ID rownames of a table to HGNC_ID IDs
#'
#'This function takes a table with Ensembl Gene IDs, converts them to HGNC symbols, 
#'leaves out all non-converted IDs (Blanks) and deduplicates 
#'
#'Deduplication is performed based on higher expression and if that is a tie then 
#'it chooses higher variation.
#'If only 1 gene has a mean expression over 1, removes all others. 
#'If 2 or more genes have expression over 1, or 2 or all genes have expression 
#'lower than 1, then keeps the one with the highest variance.
#'
#'Dependencies are:
#' library(biomaRt)
#'
#'date 04/07/2022
#'@param table Input table with ENSG IDs as rownames
#'@param ensembl_version The version of ensembl to be used. If none is provided the function uses the default (latest) version
#'@param organism The organism to search for: hsapiens or mmusculus
#'@param version TRUE/FALSE. Whether the identifiers of the original dataframe have an ensembl gene id version number
#'
#'@return A dataframe with deduplicated HGNC IDs as rownames
#'@export


ensg_to_hgnc <- function(table, ensembl_version = 0, organism = 'hsapiens',
                         version = FALSE) {
  
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package \"biomaRt\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  table <- as.data.frame(table)
  
  
  if(organism == 'hsapiens')
  {
    dataset = 'hsapiens_gene_ensembl'
    
    
    
    if (ensembl_version != 0){
      #if an ensembl version is provided, the mart is defined by that version number
      mart = useEnsembl(biomart = 'ensembl', dataset = dataset, version = ensembl_version)  
      
      #the process of converting and deduplicating is different if the original identifiers have a version number or not.
      if(version){
        conv_table <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id_version'), filters = 'ensembl_gene_id_version', values = rownames(table), mart = mart)
        conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id_version),]
        table$genes <- conv_table$hgnc_symbol[match(rownames(table),conv_table$ensembl_gene_id_version)]  
      } else {
        conv_table <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id', values = rownames(table), mart = mart)
        conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id),]
        table$genes <- conv_table$hgnc_symbol[match(rownames(table),conv_table$ensembl_gene_id)]  
      }
      
      #if an ensembl version if not provided then the mart will use the latest ensembl version.
    } else {
      #rownames(table) <- gsub('\\.[[:digit:]]{1+}$','',rownames(table)) This removed the version number for safety. Now it completely dependend upon the user.
      mart = useEnsembl(biomart = 'ensembl', dataset = dataset)
      
      if(version){
        conv_table <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id_version'), filters = 'ensembl_gene_id_version', values = rownames(table), mart = mart)
        conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id_version),]
        table$genes <- conv_table$hgnc_symbol[match(rownames(table),conv_table$ensembl_gene_id_version)]  
      } else {
      conv_table <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id', values = rownames(table), mart = mart)
      conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id),]
      table$genes <- conv_table$hgnc_symbol[match(rownames(table),conv_table$ensembl_gene_id)]
      }
    }
    
    
    
  } else if (organism == 'mmusculus') {
    
    
    
    dataset = 'mmusculus_gene_ensembl'
    
    
    if (ensembl_version != 0){
      ##If an ensembl version is provided, the mart is defined by it.
      mart = useEnsembl(biomart = 'ensembl', dataset = dataset, version = ensembl_version)  
      
      if(version){
      conv_table <- getBM(attributes = c('mgi_symbol','ensembl_gene_id_version'), filters = 'ensembl_gene_id_version', values = rownames(table), mart = mart)
      conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id_version),]
      table$genes <- conv_table$mgi_symbol[match(rownames(table),conv_table$ensembl_gene_id_version)]
      } else {
      conv_table <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id', values = rownames(table), mart = mart)
      conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id),]
      table$genes <- conv_table$mgi_symbol[match(rownames(table),conv_table$ensembl_gene_id)] 
      }
    
      
      } else {
      ##if an ensembl version is not provided, the mart will use the latest available ensembl version
        #rownames(table) <- gsub('\\.[[:digit:]]{1+}$','',rownames(table))
      mart = useEnsembl(biomart = 'ensembl', dataset = dataset)
      
      if(version){
      conv_table <- getBM(attributes = c('mgi_symbol','ensembl_gene_id_version'), filters = 'ensembl_gene_id_version', values = rownames(table), mart = mart)
      conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id_version),]
      table$genes <- conv_table$mgi_symbol[match(rownames(table),conv_table$ensembl_gene_id_version)]
      } else {
      conv_table <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id', values = rownames(table), mart = mart)
      conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id),]
      table$genes <- conv_table$mgi_symbol[match(rownames(table),conv_table$ensembl_gene_id)]
      }
    }
    
  } else {
    errorCondition(message = 'Invalid organism, please select hsapiens or mmusculus')
  }
  
  #Discarding blank and NA hgnc symbols
  table = table[table$genes != '',]
  table = table[!is.na(table$genes),]
  
  double_genes <- table$genes[which(duplicated(table$genes))]  
  table_no_genes <- data.frame(table[,-length(table)])
  rownames(table_no_genes) <- rownames(table)
  colnames(table_no_genes) <- colnames(table)[-length(table)]
  genes_to_remove  <- c()
  
  if(length(double_genes)>0) {
      
    for (i in 1:length(double_genes)) {
      
      dg <- double_genes[i]
      possible_indices <- which(table$genes == dg)
      
      #subset the data frame with only the duplicate genes
      dedup <- table[table$genes == dg,colnames(table)!='genes']
      tdedup <- as.data.frame(t(dedup))
      
      #calculate mean expression and variance of duplicate genes
      dmeans <- data.frame(colMeans(tdedup, na.rm=TRUE))
      dvars <- data.frame(lapply(tdedup, var, na.rm = TRUE))
      
      
      #check which genes to keep
      if(sum(dmeans > 1) == 1) {
        #if only 1 gene has a mean expression over 1, remove all others 
        genes_to_remove <-  c(genes_to_remove,possible_indices[dmeans <1])
      } else if(sum(!is.na(dvars)) > 0 & length(possible_indices[dvars == max(dvars)]) == 1) {
        #if 2 or more genes have expression over 1 or 2 or all genes have expression lower than 1 then keep the one with the highest variance
        genes_to_remove <- c(genes_to_remove, possible_indices[dvars != max(dvars)])
      } else if((sum(dmeans) != 0) & (sum(dmeans == max(dmeans)) == 1 ))  {
        genes_to_remove <- c(genes_to_remove, possible_indices[dmeans != max(dmeans)])
      } else {
        genes_to_remove <- c(genes_to_remove, possible_indices[2:length(possible_indices)])
      }
    
    } ## end of duplicate gene indexes decisive loop
    if(organism == 'hsapiens') {
      
      if(version) {
          final_table <- data.frame(table_no_genes[-genes_to_remove,])
          rownames(final_table) <- rownames(table_no_genes)[-genes_to_remove]
          rownames(final_table) <- conv_table$hgnc_symbol[match(rownames(final_table),conv_table$ensembl_gene_id_version)]
          } else {
          final_table <- data.frame(table_no_genes[-genes_to_remove,])
          rownames(final_table) <- rownames(table_no_genes)[-genes_to_remove]
          rownames(final_table) <- conv_table$hgnc_symbol[match(rownames(final_table),conv_table$ensembl_gene_id)]      
          }
    } else if (organism == 'mmusculus') 
      {
      if(version) {
          final_table <- data.frame(table_no_genes[-genes_to_remove,])
          rownames(final_table) <- rownames(table_no_genes)[-genes_to_remove]
          rownames(final_table) <- conv_table$mgi_symbol[match(rownames(final_table),conv_table$ensembl_gene_id_version)]
          } else {
          final_table <- data.frame(table_no_genes[-genes_to_remove,])
          rownames(final_table) <- rownames(table_no_genes)[-genes_to_remove]
          rownames(final_table) <- conv_table$mgi_symbol[match(rownames(final_table),conv_table$ensembl_gene_id)]
          }
        }
    ####end of if length_double>0 bracket
  } else 
    {
    
    if(organism == 'hsapiens') {
        if(version) {
        final_table <- data.frame(table_no_genes)
        rownames(final_table) <- conv_table$hgnc_symbol[match(rownames(table_no_genes),conv_table$ensembl_gene_id_version)]
        } else {
        final_table <- data.frame(table_no_genes)
        rownames(final_table) <- conv_table$hgnc_symbol[match(rownames(table_no_genes),conv_table$ensembl_gene_id)]      
        }
      } else if(organism == 'mmusculus') {
      
      if(version){
        final_table <- data.frame(table_no_genes)
        rownames(final_table) <- conv_table$mgi_symbol[match(rownames(table_no_genes),conv_table$ensembl_gene_id_version)]
      } else {
        final_table <- data.frame(table_no_genes)
        rownames(final_table) <- conv_table$mgi_symbol[match(rownames(table_no_genes),conv_table$ensembl_gene_id)]
      }
    
      }
    
      } ##end of if length_double>0 bracket else version
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    colnames(final_table) <- colnames(table_no_genes)
    return(final_table)
    
    }  ##end of function bracket
  
  


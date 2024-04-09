#'Convert raw counts to tpm counts
#'
#'This function takes a table with raw expression counts and converts them to TPM 
#'This function does not take into account effective gene length.
#'Dependencies are:
#'
#'
#'date 09/04/2024
#'@param df Input data frame with raw counts and gene ids as rownames
#'@param annot A dataframe with annotation containing Gene start (bp) or Gene end (bp) columns to determine gene length. This should also contain a column with the gene ids (as found on the df input) to match the dfs.
#'@param check.names Whether to use 'Gene.start.(bp)' or 'Gene start (bp)'
#'@param gene_col_id The name of the column that corresponds to the gene ids as they are found on the rownames of the df input.
#'
#'@return A dataframe with TPM expression values
#'@export


counts_to_tpm_94 <- function (df, annot, check.names = FALSE, gene_col_id = 'Gene stable ID version') 
{
  if (check.names) {
    annot$gene_length = abs(annot$`Gene.start.(bp)` - annot$`Gene.end.(bp)`)  
  } else {
    annot$gene_length = abs(annot$`Gene start (bp)` - annot$`Gene end (bp)`)
  }
  
  #Length normalization
  df_rpk <- df/annot$gene_length[match(rownames(df), annot[,gene_col_id])]
  
  #Library size normalization
  df_tpm <- sweep(df_rpk, 2, colSums(df_rpk, na.rm = TRUE)/10^6, 
                  "/")
  #return the data frame
  df_tpm
}
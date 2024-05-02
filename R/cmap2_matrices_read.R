#'Parse Cmap2 results into a table
#'
#'Dependencies are:
#' library(cmapR)
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#'
#'date 19/04/2024
#'@param clue_results_dir Directory path of the downloaded CLUE Cmap2.0 results
#'@param keep_only_morpheus Whether to keep only rows corresponding to the online CLUE Morpheus visualization (essentially the same table as the query_result.gct in the ARFS folder, with the addition of the fdr column)
#'
#'
#'@return A dataframe
#'@export

cmap2_matrices_read <- function(clue_results_dir, keep_only_morpheus = FALSE) {
  if(!requireNamespace("cmapR")) {
    stop('Please load cmapR to continue')
  }
  if(!requireNamespace("dplyr")) {
    stop('Please load dplyr to continue')
  }
  if(!requireNamespace("tibble")) {
    stop('Please load tibble to continue')
  }
  if(!requireNamespace("tidyr")) {
    stop('Please load tidyr to continue')
  }
  
  
  
  #identifying files
  query_dir <- paste0(clue_results_dir,'/matrices/query/')
  filelist <- list.files(path = query_dir)
  
  ncs_file <- paste0(query_dir,filelist[grepl('ncs.gct',filelist)])
  fdr_file <- paste0(query_dir,filelist[grepl('fdr_qvalue.gct',filelist)])
  raw_cs_file <- paste0(query_dir,filelist[grepl('^cs.gct',filelist)])
  
  #reading normalized connectivity score
  clue_res_query <- parse_gctx(ncs_file)
  
  #melting
  melt_norm_cs <- melt_gct(clue_res_query) %>%
    mutate(norm_cs = round(value,4), 
           id = id.x) %>% 
    dplyr::select(-id.y, -id.x, -value)
  
  #reading fdr values
  clue_fdrs <- parse_gctx(fdr_file)
  
  #meltin
  melt_fdrs <- melt_gct(clue_fdrs) %>%
    dplyr::select(fdr = value, id = id.x, - id.y) %>%
    mutate(fdr_q_nlog10 = round(-log10(fdr + 2.22024039534507e-16),4)) #2.22024039534507e-16 was identified as the lowest fdr q value (besides 0) and 
  #it was determined that in order to match
  #the fdr_q_nlog10 in the ARFS query_result, this value should be added to 
  #the fdr values of the matrices/query/fdr.gct values
  #to find this I examined the histogram and the tables of fdr_q_nlog10 > 5 values
  #of each table, and then I converted the same values back to fdr to check their difference.
  
  #reading raw_cs
  clue_raw_cs <- parse_gctx(raw_cs_file)
  
  #melting
  melt_raw_cs <- melt_gct(clue_raw_cs) %>%
    dplyr::select(raw_cs = value, id = id.x, - id.y) %>%
    mutate(raw_cs = round(raw_cs,4))
  
  
  #combining into 1 dataframe
  melt_complete <- melt_norm_cs %>%
    left_join(melt_fdrs, by = 'id') %>%
    left_join(melt_raw_cs, by = 'id') %>%
    mutate(id.y = 'norm_cs') %>%
    mutate(ss_ngene = as.integer(ss_ngene),
           cc_q75 = as.numeric(cc_q75),
           tas = as.numeric(tas),
           qc_pass = as.integer(qc_pass))
  
  #replacing -666 with NAs
    #from connectopedia: Note that any missing metadata value is represented by -666, 
    #which indicates that the information is not available or not applicable.
  melt_complete[melt_complete=='-666'] <- NA
  
  if(keep_only_morpheus) {
    #keeping only CLUE MORPHEUS (ARFS) rows with is_exemplar_sig
    melt_final <- melt_complete %>%
      filter(is_exemplar_sig == 1) %>%
      arrange(id)
    rownames(melt_final) <- 1:nrow(melt_final)
  } else {
    melt_final <- melt_complete
  }
  
  melt_final
}


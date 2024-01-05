#'Read multiple gProfiler output csv files
#'
#'Dependencies are:
#' none
#' 
#'
#'date 12/07/2023
#'@param enrich_folder The directory containing the gprofiler files to be read.
#'@param no_of_top_results If you need filtering, how many top results to keep from each file. 10000 usually results in all terms being collected.
#'@param pvalue_thres The p value significance threshold for the enriched terms
#'@param file_labels The labels of the files. This is going to be converted into a column of the final table. This will make it easier to identify which file the enriched term originated from.
#'
#'
#'
#'
#'@return A dataframe
#'@export


read_gprofiler_multi <- function(enrich_folder,
                                  no_of_top_results = 10000,
                                  pvalue_thres = 0.05,
                                  file_labels = NULL) {
  
  
  #read files list
  enrich_files <- list.files(enrich_folder, pattern = '.csv')
  #the loop below reads CSVs exported from gProfiler, from a specific folder and produces a table with the top N term from each file
  
  
  
  for (i in 1:length(enrich_files)) {
    file = enrich_files[i]
    
    #read file
    ids <- read.csv(paste0(enrich_folder, enrich_files[i]),header = TRUE)
    
    #order according to increasing p-value - most significants first
    ids1 <- ids[order(ids$adjusted_p_value),]
    #deduplicate terms
    ids2 <- ids1[!duplicated(ids1$term_name),]
    
    #set name - Warning! Name is inferred from enrichment file if not provided.
    
    if(is.null(file_labels)){
      #Name should be the first word and separated by -> _ <- from the all other characters in the filename.
      file_label <- data.frame(c(rep(strsplit(paste0(file),'_')[[1]][1],no_of_top_results)))
      names(file_label) <- 'File_label'  
    } else if (!is.null(file_labels)) {
      file_label <- data.frame(c(rep(file_labels[i],no_of_top_results)))
      names(file_label) <- 'File_label'  
    }
    
    #creating temporary table with top results
    ids3 <- ids2[1:no_of_top_results,c('source','term_name','intersection_size','adjusted_p_value')]
    ids4 <- cbind(ids3,file_label)
    assign(paste0('table',i),ids4)
    
    #assigning to the final enrichment table
    if(i>1){
      enrichment_all1 <- rbind(enrichment_all1,get(paste0('table',i)))   
    } else {
      enrichment_all1 <- rbind(get(paste0('table',i)))  
    }
    
    #clearing unneeded variables
    rm(file, ids, ids1, ids2, ids3, ids4, file_label, list = paste0('table',i))
    
    #filtering for p_value
    if(i == length(enrich_files)){
      enrichment_all <- enrichment_all1[enrichment_all1$adjusted_p_value < pvalue_thres,]
      rm(enrichment_all1)
    }
    
  } # end of for loop
  
  #clearing out NA rows
  enrichment_all <- enrichment_all[!rowSums(is.na(enrichment_all)) == ncol(enrichment_all),]
  enrichment_all
  
} ##end of read multi function
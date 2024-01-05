#' Perform enrichment analysis on WGCNA gene lists using the latest version of gprofiler2
#' 
#' 
#' 
#'Dependencies are:
#' library(gprofiler2)
#' library(openxlsx)
#'
#'date 02/08/2023
#'@param wgcna_gene_info Table containing gene names and the module they belong
#'@param gene_id_cols The name of the column containing gene ids
#'@param organism The organism to search for: hsapiens or mmusculus
#'@param modules A vector containing the modules to be tested. If this is given, any other module parameter is overriden. 
#'@param module_Results A data frame with module names as rownames containing columns with correlation p.values between module eigengenes and phenotypical traits.
#'@param significant_only TRUE/FALSE. Whether to check only modules that are significantly correlated with at least 1 phenotypical trait or check every module detected.
#'@param export_to_xlsx TRUE/FALSE. Whether to export results for each module to an excel file.
#'
#'
#'@return A list with enrichment results as elements.
#'@export
#'
#'
#'
#'
#'
#'

wgcna_enrich_modules <- function(wgcna_gene_info,
                                 gene_id_cols = 'EnsemblGeneID',
                                 module_color_column = 'moduleColor',
                                 organism = 'hsapiens',
                                 modules = NULL,
                                 module_Results = NULL,
                                 significant_only = TRUE,
                                 export_to_xlsx = FALSE) {
  
  if(!requireNamespace("gprofiler2")) {
    stop('Please load gprofiler2 to continue')
  }
  
  gost_res <- list()
  
  if(significant_only == TRUE & is.null(modules)) {
    stopifnot('To check for significance you have to provide module_Results data frame' =  !is.null(module_Results))
    
    #check for modules significant in at least 1 characteristic
    mods_sig <- module_Results[grepl('^p\\.',colnames(module_Results))] < 0.05
    #isolate names
    mods <- rownames(mods_sig)[rowSums(mods_sig)>0]
    #fix names
    modules <- gsub('^ME','',mods)
      
  } else if(significant_only == FALSE & is.null(modules)) {
    stopifnot('To check all modules you have to provide module_Results data frame' =  !is.null(module_Results))
    #Check all names
    modules = gsub('^ME','', rownames(module_Results))
    }
  
  if(!is.null(modules)) {
    cat('Proceeding with selected modules:\n')
    cat(paste0(modules))
  }
  
  
  if(export_to_xlsx) {
    wb <- createWorkbook()
  }
  
  
  
  ##enrichment check loop
  for (i in 1:length(modules)) {
    #selecting module
    module <- modules[i]
    #isolating gene list
    genes <- geneInfo %>%
      filter(get(module_color_column) == module) %>%
      mutate(ensg = remove_version(get(gene_id_cols))) %>%
      dplyr::select(ensg, paste0('MM.',module), paste0('p.MM.',module))
    #running enrichment
    gostMod <- gost(genes$ensg, organism = organism)
    
    gost_res[[module]] <- gostMod
    
    if(export_to_xlsx) {
      
      if(!requireNamespace("openxlsx")) {
        stop('Please load openxlsx to export to excel file')
      }
      
      addWorksheet(wb = wb, sheetName = module)
      
      writeData(wb = wb, sheet = module, x = gostMod$result)
                 
    }
    
    }
  
  
  if(export_to_xlsx){
    saveWorkbook(wb, 
                 paste0('gprofiler_enrichment_',
                        gsub('\\.[[:digit:]]+','',
                             as.numeric(Sys.time())),
                        '.xlsx')
                 )
    
  }
  gost_res
  
  
}



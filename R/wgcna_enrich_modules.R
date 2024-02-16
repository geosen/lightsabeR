#' Perform enrichment analysis on WGCNA gene lists using the latest version of clusterProfiler
#' 
#' 
#' 
#'Dependencies are:
#' library(clusterProfiler)
#' library(openxlsx)
#' library(dplyr)
#' library(stringr)
#' library(msigdbr)
#'
#'date 02/08/2023
#'@param wgcna_gene_info Table containing gene names and the module they belong
#'@param gene_id_cols The name of the column containing gene ids (HGNC_Symbols)
#'@param organism The organism to search for: Homo sapiens or Mus musculus
#'@param modules A vector containing the modules to be tested. If this is given, any other module parameter is overriden. 
#'@param module_Results A data frame with module names as rownames containing columns with correlation p.values between module eigengenes and phenotypical traits.
#'@param significant_only TRUE/FALSE. Whether to check only modules that are significantly correlated with at least 1 phenotypical trait or check every module detected.
#'@param export_to_xlsx TRUE/FALSE. Whether to export results for each module to an excel file.
#'@param genesets The collections of MSigDB R package to check enrichment of in the format of Category;Subcategory.
#'@param fdr_threshold The FDR threshold for results exported in Excel file
#'@param pvalueCutoff ClusterProfiler's adjusted pvalue cutoff on enrichment tests to report
#'@param pAdjustMethod ClusterProfiler's: one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#'@param export_fdr_thres Used when export to xlsx is TRUE to filter the final list of geneset results to be exported. By default, it equals the pvalueCutoff.
#'@param minGSSize ClusterProfiler's minimal size of genes annotated for testing
#'@param maxGSSize ClusterProfiler's maximal size of genes annotated for testing
#'@param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
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
                                 gene_id_cols = 'geneSymbol',
                                 module_color_column = 'moduleColor',
                                 organism = 'Homo sapiens',
                                 modules = NULL,
                                 module_Results = NULL,
                                 significant_only = TRUE,
                                 export_to_xlsx = FALSE,
                                 collections = c('H','C2;CP','C5;GO:BP'),
                                 pvalueCutoff = 0.05,
                                 export_fdr_thres = 0.05,
                                 pAdjustMethod = 'BH',
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 qvalueCutoff = 1) {
  
  if(!requireNamespace("clusterProfiler")) {
    stop('Please load clusterProfiler to continue')
  }
  
  if(export_to_xlsx & !requireNamespace("openxlsx")) {
    stop('Please load openxlsx to continue')
  }
  if(!requireNamespace("dplyr")) {
    stop('Please load dplyr to continue')
  }
  if(!requireNamespace("stringr")) {
    stop('Please load stringr to continue')
  }
  if(!requireNamespace("msigdbr")) {
    stop('Please load msigdbr to continue')
  }
  
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
  
  ##Defining genesets 
  categories <- as.data.frame(stringr::str_split(collections,';',n = 2,simplify = TRUE))
  
  #Checking categories
  stopifnot('Invalid category entered. Please check spelling!' = all(pull(categories,1) %in% unique(pull(msigdbr_collections(),'gs_cat'))))
  #checking subcategories
  stopifnot('Invalid subcategory entered. Please check spelling!' = all(pull(categories,2) %in% unique(pull(msigdbr_collections(),'gs_subcat'))))
  
  #Gathering genesets
  for (i in 1:dim(categories)[1]) {
    
    if(i == 1) {
      genesets <- msigdbr(species = organism,
                          category= categories[i,1],
                          subcategory = categories[i,2])
    } else {
      genesets <- genesets %>%
        bind_rows(msigdbr(species = organism,
                          category= categories[i,1],
                          subcategory = categories[i,2]))
    }
    
  }
  
  ##creating background
  background <- geneInfo %>%
    filter(!duplicated(geneSymbol),!is.na(geneSymbol), geneSymbol != '') %>%
    pull(geneSymbol)
    
  
  ##Creating Term2Gene
  t2g_cp <- genesets %>%
    dplyr::distinct(gs_name,
                    gene_symbol) %>%
    as.data.frame()
  
    #creating list for storing results
    cp_res <- list()
  
    ##enrichment check loop
  for (i in 1:length(modules)) {
    #selecting module
    module <- modules[i]
    #isolating gene list
    genes <- geneInfo %>%
      filter(get(module_color_column) == module) %>%
      mutate(geneSymb = lightsabeR::remove_version(get(gene_id_cols))) %>%
      dplyr::select(geneSymb, paste0('MM.',module), paste0('p.MM.',module))
    
    
    #running enrichment
    cpMod <- enricher(genes$geneSymb,
                        pvalueCutoff = pvalueCutoff,
                        pAdjustMethod = pAdjustMethod,
                        universe = background,
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        qvalueCutoff = qvalueCutoff,
                        TERM2GENE = t2g_cp)
    
    
    cp_res[[module]] <- cpMod
    
    if(export_to_xlsx) {
      
      if(!requireNamespace("openxlsx")) {
        stop('Please load openxlsx to export to excel file')
      }
      
      addWorksheet(wb = wb, sheetName = module)
      
      writeDataTable(wb = wb, sheet = module, x = cpMod@result %>%
                       filter(p.adjust < export_fdr_thres))
                 
    }
    
    }
  
  
  if(export_to_xlsx){
    saveWorkbook(wb, 
                 paste0('clusterProfiler_enrichment_',
                        gsub('\\.[[:digit:]]+','',
                             as.numeric(Sys.time())),
                        '.xlsx')
                 )
    
  }
  cp_res
  
  
}



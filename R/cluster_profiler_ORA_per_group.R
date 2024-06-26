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
#'date 02/05/2024
#'@param df Table containing gene names and the group they belong
#'@param gene_id_cols The name of the column containing gene ids (HGNC_Symbols)
#'@param group_column The name of the column defining the gene groups
#'@param organism The organism to search for: Homo sapiens or Mus musculus
#'@param groups A vector containing the modules to be tested. If this is given, any other module parameter is overriden. 
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

cluster_profiler_ORA_per_group <- function(df,
                                 gene_id_cols = 'hgnc_symbol',
                                 group_column = 'gene_groups',
                                 organism = 'Homo sapiens',
                                 groups = NULL,
                                 ids_type = c('hgnc_symbol','ensembl_gene_id'),
                                 export_to_xlsx = FALSE,
                                 collections = c('H','C2;CP','C5;GO:BP'),
                                 pvalueCutoff = 0.05,
                                 export_fdr_thres = 0.05,
                                 pAdjustMethod = 'BH',
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 qvalueCutoff = 1,
                                 useBackground = T) {
  
  if(!requireNamespace("clusterProfiler")) {
    stop('Please load clusterProfiler to continue')
  }
  
  if(!requireNamespace("openxlsx")) {
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
  
  #handle ids_type
  ids_type <- match.arg(ids_type)
  
  #removing ensembl version in case there is one
  if(ids_type == 'ensembl_gene_id') {
    df[,gene_id_cols] <- lightsabeR::remove_version(df[,gene_id_cols])
  }
  
  if(!is.null(groups)) {
    cat('Proceeding with selected groups:\n')
    cat(paste0(groups))
  } else {
      groups = levels(as.factor(df[,group_column]))
      }
  
  
  if(export_to_xlsx) {
    wb <- createWorkbook()
  }
  
  ##Defining genesets 
  categories <- as.data.frame(stringr::str_split(collections,';',n = 2,simplify = TRUE))
  
  #Checking categories
  stopifnot('Invalid category entered. Please check spelling!' = all(pull(categories,1) %in% unique(pull(msigdbr_collections(),'gs_cat'))))
  
  if(!all(pull(categories,2) == '')) { #this performs subcategory check only if subcategories are present
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
  } else { #if subcategories are not present
    for (i in 1:dim(categories)[1]) {
      
      if(i == 1) {
        genesets <- msigdbr(species = organism,
                            category= categories[i,1])
      } else {
        genesets <- genesets %>%
          bind_rows(msigdbr(species = organism,
                            category= categories[i,1]))
      }
      
    }
    
  }
  
  
  if(useBackground){
    
  ##creating background
  background <- df %>%
    filter(!duplicated(get(gene_id_cols)),!is.na(get(gene_id_cols)), get(gene_id_cols) != '') %>%
    pull(get(gene_id_cols))
  }
  
  ##Creating Term2Gene
  if(ids_type == 'hgnc_symbol'){
  t2g_cp <- genesets %>%
    dplyr::distinct(gs_name,
                    gene_symbol) %>%
    as.data.frame()
  } else if(ids_type == 'ensembl_gene_id') {
    t2g_cp <- genesets %>%
      dplyr::distinct(gs_name,
                      ensembl_gene) %>%
      dplyr::select(gs_name, gene_symbol = ensembl_gene) %>%
      as.data.frame()
  } else {
    stop('Please provide a valid ids_type value')
  }
  
  #creating list for storing results
  cp_res <- list()
  
  ##enrichment check loop
  for (i in 1:length(groups)) {
    
    #selecting group
    group <- groups[i]
    
    #isolating gene list
    genes <- df %>%
      filter(get(group_column) == group) %>%
      mutate(geneSymb = get(gene_id_cols)) %>%
      dplyr::select(geneSymb)
    
    if(length(intersect(genes$geneSymb,t2g_cp$gene_symbol)) != 0){
      
    #running enrichment
      if(useBackground){
        
    cpMod <- enricher(genes$geneSymb,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      universe = background,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      qvalueCutoff = qvalueCutoff,
                      TERM2GENE = t2g_cp)
      } else {
      cpMod <- enricher(genes$geneSymb,
                        pvalueCutoff = pvalueCutoff,
                        pAdjustMethod = pAdjustMethod,
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        qvalueCutoff = qvalueCutoff,
                        TERM2GENE = t2g_cp)
      }
    
    cp_res[[group]] <- cpMod
    
    if(export_to_xlsx) {
      
      if(!requireNamespace("openxlsx")) {
        stop('Please load openxlsx to export to excel file')
      }
      
      addWorksheet(wb = wb, sheetName = group)
      
      writeDataTable(wb = wb, sheet = group, x = cpMod@result %>%
                       filter(p.adjust < export_fdr_thres))
      
    }
    } else {
      warning(paste0('No genes mapped to pathways for group: ', group,'.\nDo not expect an entry for this group in the relevant object and/or xlsx file.'))
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



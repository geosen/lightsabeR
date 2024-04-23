#'Create a STRINGdb PPI network based on genes in selected WGCNA modules
#'
#'Dependencies are:
#' library(tidyverse)
#' library(STRINGdb)
#' library(igraph)
#'
#'date 23/04/2024
#'@param wgcna_gene_info Data frame containing the results of gene-module membership of WGCNA results
#'@param moduleColors Names of the modules to be used for the network
#'@param module_color_col Name of the column containing the gene module color
#'@param string_version Version of the STRINGdb to be queried
#'@param species Species to be queried - 9606 is H.sapiens
#'@param score_threshold Interaction score threshold
#'@param network_type Physical of full
#'@param geneNameCol Name of the column containing the gene names
#'
#'
#'@return An igraph network object
#'@export




wgcna_module_stringdb_network_create <- function(wgcna_gene_info,
                                                 moduleColors, 
                                                 module_color_col = 'moduleColor',
                                                 string_version = '12.0',
                                                 species = 9606,
                                                 score_threshold = 500,
                                                 network_type = c('physical','full'),
                                                 geneNameCol = 'geneSymbol') {
  
  if(!requireNamespace("STRINGdb")) {
    stop('Please load STRINGdb to continue')
  }
  if(!requireNamespace("tidyverse")) {
    stop('Please load tidyverse to continue')
  }
  if(!requireNamespace("igraph")) {
    stop('Please load igraph to continue')
  }

#isolating modules gene lists
isolMod <- wgcna_gene_info %>%
  filter(get(module_color_col) %in% moduleColors)

#creating String DB instance
string <- STRINGdb(species = species,
                   score_threshold = score_threshold, #>=600 is kind of strict
                   network_type = network_type,
                   version = string_version)



#first step - mapping gene names
mapped <- string$map(isolMod, 
                     geneNameCol, 
                     removeUnmappedRows = TRUE)

#second step - get interactions
interactions <- string$get_interactions(mapped$STRING_id)

#third step = changing String IDs to Gene Symbols and creating final network df
ann <- interactions %>%
  left_join(mapped %>% 
              dplyr::select(paste0(geneNameCol),STRING_id, paste0(module_color_col)),
            by = c('from' = 'STRING_id')) %>%
  mutate(from_geneSymbol = get(geneNameCol),
         from_moduleColor = get(module_color_col)) %>%
  dplyr::select(-paste0(geneNameCol),-paste0(module_color_col)) %>% ##end of 'from' genes
  left_join(mapped %>% 
              dplyr::select(paste0(geneNameCol),STRING_id, paste0(module_color_col)),
            by = c('to' = 'STRING_id')) %>%
  mutate(to_geneSymbol = get(geneNameCol),
         to_moduleColor = get(module_color_col)) %>%
  dplyr::select(-paste0(geneNameCol),-paste0(module_color_col)) %>% ##end of 'to' genes
  mutate(interaction_pair = paste0(from_geneSymbol,'_',to_geneSymbol)) %>%
  filter(!duplicated(interaction_pair))


vertices1 <- data.frame(vertices = unique(c(ann$from_geneSymbol,
                                            ann$to_geneSymbol)))

vertices <- vertices1 %>%
  left_join(wgcna_gene_info, by = c('vertices' = geneNameCol)) %>%
  filter(!duplicated(vertices)) %>%
  dplyr::select(vertices,paste0(module_color_col))

#fourth step - creating network objects
net <- graph_from_data_frame(d = ann[,c('from_geneSymbol','to_geneSymbol')],
                             vertices = vertices,
                             directed = FALSE)
net
}


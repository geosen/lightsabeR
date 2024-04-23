#'Create a STRINGdb PPI network based on genes in selected WGCNA modules
#'
#'Dependencies are:
#' library(tidyverse)
#' library(ggraph)
#' library(igraph)
#'
#'date 23/04/2024
#'@param net The network object to be plotted
#'@param layout layout package
#'@param algorithm igraph layout algorithm to be used
#'@param seed Seed to be used
#'@param edge_alpha Alpha value for edge transparency
#'@param node_size_vec A vector specifying node sizes. By default it is a function of their degree
#'@param node_shape Shape of the nodes
#'@param node_stroke Node circle line width
#'@param node_alpha Alpha value for node transparency
#'@param node_text_size_vec A vector specifying node text sizes. By default it is a function of their degree
#'@param pdf_width If print_to_pdf = T, specifies pdf width
#'@param pdf_height If print_to_pdf = T, specifies pdf height
#'@param pdf_filename If print_to_pdf = T, specifies pdf filename
#'@param print_to_pdf Whether to export to pdf file
#'
#'
#'@return A ggplot object
#'@export
#'
wgcna_module_stringdb_network_plot <- function(net,layout = 'igraph',algorithm = 'kk',seed = 7,
                                               edge_alpha = 0.05, node_size_vec = NULL,node_shape = 21,
                                               node_stroke = 1, node_alpha = 0.7, node_text_size_vec = NULL,
                                               pdf_width = 14, pdf_height = 14, pdf_filename = NULL, print_to_pdf = F) {
  
  if(!requireNamespace("ggraph")) {
    stop('Please load ggraph to continue')
  }
  if(!requireNamespace("tidyverse")) {
    stop('Please load tidyverse to continue')
  }
  if(!requireNamespace("igraph")) {
    stop('Please load igraph to continue')
  }
  
  if(is.null(node_size_vec)) {
    node_size_vec <- (1.5*degree(net)/(mean(degree(net)))+4)
  }
  
  if(is.null(node_text_size_vec)) {
    node_text_size_vec <- (2*degree(net)/(mean(degree(net))**2)+1)
  }
  #fifth step - visualizing network
  set.seed(7)
  #layout
  layout <- create_layout(net, layout = "igraph", algorithm = igraph_layout_algorithm)
  #layout2 <- create_layout(net, layout = "igraph", algorithm = 'dh')
  
  g <- ggraph(net, layout = layout2) +
    geom_edge_link(alpha = edge_alpha) + 
    geom_node_point(size = node_size_vec, 
                    shape = node_shape, 
                    stroke = node_stroke, 
                    fill = vertices$moduleColor,
                    alpha = node_alpha) +
    geom_node_text(aes(label = name), 
                   size = node_text_size_vec) + 
    theme_void() 
  
  if(print_to_pdf) {
    if(is.null(pdf_filename)) {
      pdf(paste0(paste(moduleColors, collapse = '_'),'_network_String_v',string_version,'.pdf'), width = pdf_width, height=pdf_height)
    } else {
      pdf(pdf_filename, width = pdf_width, height = pdf_height)
    }
    print(g)
    dev.off()
  } else {
    g 
    }
  
}

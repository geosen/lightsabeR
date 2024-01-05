#'Plot Trait Gene Significance over Module Membership using the output of WGCNA
#'
#'Dependencies are:
#' library(ggplot2)
#' library(ggrepel)
#' library(stringr)
#'
#'date 11/07/2023
#'@param wgcna_res Input wgcna results data.frame 
#'@param module The name of the module you wish to plot. It has to match the MM.module column in your wgnca_res object
#'@param trait The name of the trait you wish to plot. It has to match the GS.trait column in your wgnca_res object
#'@param module_column The name of the column containing the information of which gene belongs to which column
#'@param MM.threshold_low The lower threshold of Module membership to limit the X axis
#'@param MM.threshold_high The upper threshold of Module membership to limit the X axis
#'@param GS.threshold_low The lower threshold of Gene significance to limit the Y axis
#'@param GS.threshold_high The upper threshold of Gene significance to limit the Y axis
#'@param label_col The name of the column containing the labels you wish to annotate the plot with (i.e. Gene Symbols)
#'@param max.overlaps Parameter for geom_label_repel
#'@param min.segment.length Parameter for geom_label_repel
#'@param add_labels Whether to use geom_label_repel to add labels to the points
#'@param color_by_module Whether to color the points according to the module color
#'
#'
#'@return A ggplot object
#'@export


wgcna_hub_gene_plot <- function(wgcna_res,
                                module,
                                trait,
                                module_column = 'moduleColor',
                                MM.threshold_low = 0.7,
                                MM.threshold_high = 1,
                                GS.threshold_low = -0.5,
                                GS.threshold_high = 0.5,
                                label_col = 'geneSymbol',
                                max.overlaps = 10,
                                min.segment.length = 0.01,
                                add_labels = TRUE,
                                color_by_module = TRUE
                                )
{
  library(ggplot2)
  library(ggrepel)
  
g <- ggplot(wgcna_res[wgcna_res[,module_column] == module,], 
       aes(x = get(paste0('MM.',module)),
           y = get(paste0('GS.', trait))
           )
       ) + 
  scale_x_continuous(limits = c(MM.threshold_low, MM.threshold_high)) + 
  scale_y_continuous(limits = c(GS.threshold_low, GS.threshold_high)) + 
  geom_point() + 
  theme_minimal() + 
  labs(x = paste0(stringr::str_to_title(module), '  module membership'),
       y = paste0(trait, ' gene significance')
       )

if(add_labels) {
  
g <- g + geom_label_repel(aes(label = get(label_col)),
                   max.overlaps = max.overlaps,
                  min.segment.length = min.segment.length) 
}

if (color_by_module) {
g <- g + geom_point(aes(color = moduleColor)) + 
    scale_color_manual(values = module) + 
    guides(color = 'none')
}
g
}
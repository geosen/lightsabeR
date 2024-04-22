#'Plot a heatmap based on a Cmap2 results table
#'
#'Dependencies are:
#' library(cmapR)
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(ComplexHeatmap)
#' library(circlize)
#' library(RColorBrewer)
#'
#'date 22/04/2024
#'@param cmap_melted A melted data.frame containing Connectivity Map 2 results
#'@param n_top Number of top results to keep
#'@param direction What direction of signatures should be kept. Normal keeps mimicking signatures, reverse keeps reversing signatures, both keeps both.
#'@param CS Heatmap color scale function
#'@param brewer_palette_qual Qualitative brewer palette for mapping qualitative variables. Check out brewer_pal_info for more info.
#'@param row_label_size The size of row labels
#'@param col_label_size The size of column labels
#'@param show_row_dend Whether to show row dendrogram
#'@param heatmap_title Heatmap title
#'@param legend_title How to name the heatmap color scale
#'@param legend_breaks Which legend breaks to display
#'@param label_column Which column to use for row labels
#'@param row_annot_cols Which columns to keep as row annotation
#'@param border_color Color of the heatmap border
#'@param border_linetype Linetype of the heatmap border
#'@param border_linewidth Linewidth of the heatmap border
#'@param cell_border_color Color of the cell border
#'@param cell_border_linewidth Linewidth of the cell border
#'
#'
#'@return A heatmap
#'@export

cmap2_heatmap_create <- function(cmap_melted,
                                 n_top = 20,
                                 direction = c('both','normal','reverse'),
                                 CS = NULL,
                                 brewer_palette_qual = 'Paired',
                                 row_label_size = 7,
                                 col_label_size = 9,
                                 show_row_dend = T,
                                 heatmap_title = 'NULL',
                                 legend_title = 'Normalized\nConnectivity\nScore',
                                 legend_breaks = c(-3,-1.5,0,1.5,3),
                                 label_column = c('pert_iname','pert_id'),
                                 row_annot_cols = c("moa","pert_idose",  "pert_itime",  "pert_type", 
                                                    "target_name",  "tas", "cc_q75", "fdr_q_nlog10",
                                                    "raw_cs","cell_iname"),
                                 border_color = 'black',
                                 border_linetype = 1 ,
                                 border_linewidth = 1,
                                 cell_border_color = 'grey90',
                                 cell_border_linewidth = 1)
{
  
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
  if(!requireNamespace("circlize")) {
    stop('Please load circlize to continue')
  }
  if(!requireNamespace("ComplexHeatmap")) {
    stop('Please load ComplexHeatmap to continue')
  }
  if(!requireNamespace("RColorBrewer")) {
    stop('Please load RColorBrewer to continue')
  }
  
  
  
#Force data frame
cmap_melted <- as.data.frame(cmap_melted)

#handle direction and label column parameter
direction <- match.arg(direction)
label_column <- match.arg(label_column)

#determine top entries
if(direction == 'both') {
  cmap_top <- cmap_melted %>%
    slice_max(abs(norm_cs), n = n_top)
  } else if(direction == 'normal') {
    cmap_top <- cmap_melted %>%
      slice_max(norm_cs, n = n_top)
  } else if(direction == 'reverse') {
    cmap_top <- cmap_melted %>%
      slice_min(norm_cs, n = n_top)
  } else {
  stop('direction should be one of: both, normal, reverse')
    }

if(is.null(CS)) {
  #defining color scale
  CS <- circlize::colorRamp2(breaks = c(-3,-1,0,1,3),
                             colors = c('steelblue3','steelblue1','white','firebrick1','firebrick3'))
}


#Correcting column labels####
row_annot_labels <- data.frame(old = c("moa","pert_idose",  "pert_itime",  "pert_type", "target_name",  "tas", "cc_q75", "fdr_q_nlog10","raw_cs","cell_iname"),
                               new = c("Mech. of action","Pert. dose",  "Pert. time",  "Pert. type", "Target",  "TAS", "CC_q75", "FDR_q_nlog10","Raw_cs","Cell line"))
#determining columns to keep
row_annot_kept <- row_annot_labels %>% filter(old %in% row_annot_cols)
#renaming
colnames(cmap_top)[match(row_annot_kept$old, colnames(cmap_top))] <- row_annot_kept$new

#Creating annotation colors named vectors list
annot_list <- as.list(cmap_top[,row_annot_kept$new])

annot_colors <- list()

for (i in 1:length(annot_list)) {
    x <- annot_list[[i]]
  if(class(x) == 'character'){
      no_of_colors <- length(unique(x))
      if(no_of_colors <= RColorBrewer::brewer.pal.info[brewer_palette_qual,"maxcolors"]){
        colors1 <- RColorBrewer::brewer.pal(no_of_colors, brewer_palette_qual)
        names(colors1) <- unique(x)
      
      } else {
        set.seed(7)
        colors1 <- colors()[sample(1:657,no_of_colors)]
        names(colors1) <- unique(x)
      
      }
      annot_colors[[i]] <- colors1[!is.na(names(colors1))] #return any values without NA colnames
      
  } else if(class(x) == 'numeric') {
    if(0 > min(x, na.rm=T) & 0 < max(x, na.rm=T)) {
      colors1 <- circlize::colorRamp2(breaks = c(min(x, na.rm=T),0,max(x, na.rm=T)),
                                        colors = c('darkolivegreen3','white','mediumorchid3'))
    } else if(0 > max(x, na.rm=T)) {
      colors1 <- circlize::colorRamp2(breaks = c(min(x, na.rm=T),0),
                                        colors = c('darkolivegreen3','white'))
    } else if(0 < min(x, na.rm=T)) {
      colors1 <- circlize::colorRamp2(breaks = c(0,max(x, na.rm=T)),
                                        colors = c('white','mediumorchid3'))
    }
    annot_colors[[i]] <- colors1 #return the function
  }
                            
                      }
names(annot_colors) <- names(annot_list)

#Creating heatmap####
ht <- ComplexHeatmap::Heatmap(as.matrix(cmap_top[,c('norm_cs')]),
                        width = unit(15, "mm"),
                        row_labels = cmap_top[,label_column],
                        column_labels = 'Normalized\nConnectivity\nScore',
                        col = CS, 
                        border_gp = gpar(col = border_color, 
                                         lty = border_linetype, 
                                         lwd = border_linewidth),
                        rect_gp = gpar(col = cell_border_color, 
                                       lwd = cell_border_linewidth),
                        row_names_gp = gpar(fontsize = row_label_size),
                        column_names_gp = gpar(fontsize = col_label_size),
                        show_row_dend = show_row_dend,
                        column_title = heatmap_title,
                        heatmap_legend_param = list(title = legend_title, 
                                                    at = legend_breaks),
                        right_annotation = HeatmapAnnotation(df = cmap_top[,row_annot_kept$new],
                                                             which = 'row',
                                                             col = annot_colors )
                                                            )

draw(ht, merge_legend = TRUE)

}





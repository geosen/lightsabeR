#'Plot Module - Trait correlation using the output of WGCNA
#'
#'Dependencies are:
#' library(ComplexHeatmap)
#' library(circlize)
#' library(stringr)
#'
#'date 11/07/2023
#'@param wgcna_trait_res Input wgcna results data.frame 
#'@param row_labels A vector to name the row labels
#'@param col_labels A vector to name the column labels
#'@param row_label_size The size of the row labels
#'@param col_label_size The size of the column labels
#'@param colf A function for the heatmap colors
#'@param row_angle The angle of the x axis labels
#'@param column_km How many clusters should the columns be divided into
#'@param row_km How many clusters should the rows be divided into
#'@param show_row_dend Whether to show the row dendrogram
#'@param show_column_dend Whether to show the column dendrogram
#'@param title The title of the plot
#'@param legend_title The title of the legend
#'@param cluster_cols Whether to cluster columns or not
#'@param cluster_rows Whether to cluster rows or not
#'@param border_color Outer rectangle line color
#'@param border_linetype Outer rectangle linetype
#'@param border_linewidth Outer rectangle linewidth
#'@param cell_color Inner rectangle line color
#'@param cell_linetype Inner rectangle linetype
#'@param cell_linewidth Inner rectangle linewidth
#'
#'
#'@return A ComplexHeatmap
#'@export

wgcna_module_heatmap_plot <- function(wgcna_trait_res, 
                                      row_labels = NULL,
                                      col_labels = NULL,
                                      row_label_size = 14,
                                      col_label_size = 7,
                                      colf = NULL,
                                      row_angle = 0,
                                      column_km = 1,
                                      row_km = 1,
                                      show_row_dend = FALSE,
                                      show_column_dend = TRUE,
                                      title = 'WGCNA Module - Trait Correlation',
                                      legend_title = 'Correlation',
                                      cluster_rows = FALSE,
                                      cluster_cols = FALSE,
                                      border_color = 'black',
                                      border_linetype = 1,
                                      border_linewidth = 1,
                                      cell_color = 'white',
                                      cell_linetype = 1,
                                      cell_linewidth = 1)

{
    #filtering out P-value columns and converting to a matrix
  wgcna_mat <- as.matrix(wgcna_trait_res[,!grepl('p\\.', colnames(wgcna_trait_res))])
  
  #keeping P-values seprately for significance annotation
  wgcna_pvalues <- as.matrix(ifelse(wgcna_trait_res[,grepl('p\\.', colnames(wgcna_trait_res))] < 0.05,
                                    TRUE,FALSE))
  
  #setting row labels
  if(is.null(row_labels)) {
    row_labels = rownames(wgcna_mat)
    row_labels <- stringr::str_to_title(gsub('^ME','',row_labels[grepl('^ME', row_labels)]))
    
  }
  
  
  #setting column labels
  if(is.null(col_labels)) {
    col_labels = colnames(wgcna_mat)
  }
  
  #plotting the heatmap
  
  if(!is.null(colf)) {
    ComplexHeatmap::Heatmap(wgcna_mat,
            col = colf, 
            row_labels = row_labels, 
            column_labels = col_labels,
            row_names_gp = gpar(fontsize = row_label_size), 
            column_names_gp = gpar(fontsize = col_label_size), 
            row_names_rot = row_angle,
            column_km = column_km,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_cols,
            row_km = row_km,
            border_gp = gpar(col = border_color, 
                             lty = border_linetype,
                             lwd = border_linewidth),
            rect_gp = gpar(col = cell_color,
                           lty = cell_linetype,
                           lwd = cell_linewidth),
            show_row_dend = show_row_dend, 
            show_column_dend = show_column_dend, 
            column_title = gt_render(paste0("<span style='color:black'>**",title,"**</span>"),
                                     r = unit(2, "pt"),
                                     padding = unit(c(2, 2, 2, 2), "pt")), 
            column_title_gp = gpar(box_fill = "#EEEEEE"),
            heatmap_legend_param = list(title = legend_title),
            cell_fun = function(j, i, x, y, w, h, fill){
              if(wgcna_pvalues[i, j]) {
                grid.text('*', x, y)
              }
            }
    )
    
  } else if (is.null(colf)) {
    ComplexHeatmap::Heatmap(wgcna_mat,
            #col = colf, 
            row_labels = row_labels, 
            column_labels = col_labels,
            row_names_gp = gpar(fontsize = row_label_size), 
            column_names_gp = gpar(fontsize = col_label_size), 
            row_names_rot = row_angle,
            column_km = column_km,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_cols,
            row_km = row_km,
            border_gp = gpar(col = border_color, 
                             lty = border_linetype,
                             lwd = border_linewidth),
            rect_gp = gpar(col = cell_color,
                           lty = cell_linetype,
                           lwd = cell_linewidth),
            show_row_dend = show_row_dend, 
            show_column_dend = show_column_dend, 
            column_title = gt_render(paste0("<span style='color:black'>**",title,"**</span>"),
                                     r = unit(2, "pt"),
                                     padding = unit(c(2, 2, 2, 2), "pt")), 
            column_title_gp = gpar(box_fill = "#EEEEEE"),
            heatmap_legend_param = list(title = legend_title),
            cell_fun = function(j, i, x, y, w, h, fill){
              if(wgcna_pvalues[i, j]) {
                grid.text('*', x, y)
              }
            }
    )
  }
  
}




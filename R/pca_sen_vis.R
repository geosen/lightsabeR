#'Create a pca plot from an existing pca_object
#'
#'Creates a pca plot based on an existing pca object
#'
#'
#'Dependencies are:
#' library(ggplot2)
#' library(FactoMineR)
#'
#'
#'date 07/11/2022
#'@param pca_object PCA list created using pca_sen with return_data=TRUE
#'@param color_factor A vector of a parameter to color the samples by
#'@param color_factor_name A title for the legend of the color parameter
#'@param colors_list A list of colors if you wish to manually specify group colors
#'@param title Title of the plot
#'@param PCx The x-axis Principal Component
#'@param PCy The y-axis Principal Component
#'@param point_alpha Set the transparency of points
#'@param point_size Set the size of points
#'@param add_ellipse If you wish to add an ellipse around the points
#'@param save_pdf Whether to save the plot to a pdf file
#'@param pdf_height The height of the output pdf
#'@param pdf_width The width of the output pdf
#'@param pdf_file The name or path of the output pdf
#'
#'
#'
#'@return ggplot object
#'@export

pca_sen_vis <- function(pca_object, color_factor,
                    color_factor_name = "Color",
                    colors_list = NULL,
                    shape_vector = NULL,
                    stroke =2,
                    return_data = FALSE,
                    title = "PCA",
                    PCx=1,
                    PCy=2,
                    point_alpha = 0.7,
                    point_size = 7,
                    add_ellipse = FALSE,
                    save_pdf = 'FALSE', 
                    pdf_height = 10, 
                    pdf_width = 10,
                    pdf_file = "PCA_sen_output.pdf"
) {
  ##load libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    stop(
      "Package \"FactoMineR\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  pcacpm <- pca_object
  
  ##choose colors
  if (is.null(colors_list)) {
    color <-  grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    color <- color[!grepl('white',color)]
    color_choice <- seq(1,433,443/length(levels(factor(color_factor)))) +70
    colors_plot <- color[color_choice]
  } else {
    colors_plot <- colors_list
  }
  
  if (save_pdf) {
    pdf(pdf_file, height = pdf_height, width = pdf_width)  
  }
  
  if (length(color_factor) != dim(pcacpm$ind$coord)[1]) {
    print("Color factor length is not equal to sample size")
  }
  
  data_pca <- as.data.frame(pcacpm$ind$coord)
  ####Plotting####
  if(!is.null(shape_vector)){
    g <- ggplot(data_pca, 
                aes(data_pca[,PCx],data_pca[,PCy], color = factor(color_factor))) +
      geom_point(alpha = point_alpha, size= point_size, key_glyph='point', shape = shape_vector, stroke = stroke) + labs(title = title, x=paste0('PC',PCx, '(',round(pcacpm$eig[PCx,2],2),')'), y= paste0('PC',PCy, '(',round(pcacpm$eig[PCy,2],2),')')) +
      scale_color_manual(name = color_factor_name, labels = levels(factor(color_factor)),values = colors_plot)+
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    g <- ggplot(data_pca, 
                aes(data_pca[,PCx],data_pca[,PCy], color = factor(color_factor))) +
      geom_point(alpha = point_alpha, size= point_size, key_glyph='point') + labs(title = title, x=paste0('PC',PCx, '(',round(pcacpm$eig[PCx,2],2),')'), y= paste0('PC',PCy, '(',round(pcacpm$eig[PCy,2],2),')')) +
      scale_color_manual(name = color_factor_name, labels = levels(factor(color_factor)),values = colors_plot)+
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
  }
  
  if (add_ellipse) {
    g <- g + stat_ellipse()
  }
  
  if (save_pdf){
    print(g)  
    dev.off()  
  } else {
    g    
  }
  
}

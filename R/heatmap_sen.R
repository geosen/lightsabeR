#'
#'
#'Dependencies are:
#' library(ComplexHeatmap)
#' library(gridtext)
#'
#'date 21/07/2023
#'@param norm_expr 
#'@param lrt_table
#'@param col_fun 
#'@param z_scaling 
#'@param samples
#'@param ntop 
#'@param order_factor
#'@param title
#'@param legend_title 
#'@param col_annotation 
#'@param annotcolors
#'@param col_name_size 
#'@param row_name_size
#'@param convert_to_hgnc
#'@param ensembl_version
#'@param identifier_versions
#'@param conversion_organism
#'@param cluster_rows
#'@param cluster_columns
#'@param annotation_legend
#'@param title_render
#'@param box_fill
#'
#'@return A ComplexHeatmap
#'@export
#'
#'

heatmap_sen <- function(norm_expr,
                        lrt_table = NULL,
                        col_fun = NULL,
                        z_scaling = TRUE,
                        samples = NULL,
                        ntop = NULL,
                        order_factor = 'logFC',
                        title = NULL,
                        legend_title = NULL,
                        col_annotation = NULL,
                        annotcolors = NULL,
                        col_name_size = 8,
                        row_name_size = 8,
                        convert_to_hgnc=FALSE,
                        ensembl_version = 0,
                        identifier_versions = FALSE,
                        conversion_organism = 'hsapiens',
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        annotation_legend = NULL,
                        title_render = TRUE,
                        box_fill = '#EEEEEE') {


##choosing top n genes up and down

  if(is.null(samples)) {
    samples = colnames(norm_expr)
  }
  
  if(!is.null(ntop)){
    
    if(is.null(lrt_table)){
      stop('lrt_table is missing, cannot select top n genes')
    } else {
      upgenes <- lrt_table %>%
        slice_max(order_by = get(order_factor), n = ntop)

      topgenes <- lrt_table %>%
        slice_min(order_by = get(order_factor), n = ntop) %>% 
        bind_rows(upgenes)
      
      genes <- rownames(topgenes)
      
      ##isolating expression
      norm_genes <- norm_expr[genes,samples]
    }
    
    } else {
      norm_genes <- norm_expr[,samples]
      }



##deduplicating if needed
if(convert_to_hgnc) {
  norm_genes <- ensg_to_hgnc(norm_genes, 
                             ensembl_version = ensembl_version, 
                             organism = conversion_organism, 
                             version = identifier_versions)
} 


##scaling expression
if (z_scaling){
z_genes <- t(scale(t(norm_genes)))
} else {
  z_genes <- norm_genes
}

if(title_render){
  column_title = gt_render(paste0("<span style='color:black'>**",title,"**</span>"),
                           r = unit(2, "pt"),
                           padding = unit(c(2, 2, 2, 2), "pt"))
} else {
  column_title = title
}




#plotting heatmap
if (!is.null(col_annotation)){

  ##defining annotation colors - converting to named vector
  if(!is.null(annotcolors)) {
    ##if provided, check if it is the same number as the categories
    if(length(annotcolors) != length(unique(col_annotation))){
      stop('Number of annotation colors does not match the number of annotation categories')
    
      } else {
        #if it is the same number
      names(annotcolors) <- unique(col_annotation)
    }
      
  } else {
    annotcolors <- sample(colors()[!grepl('grey',colors()) & !grepl('gray',colors())],
                          length(unique(col_annotation))
                          )
    names(annotcolors) <- unique(col_annotation)
  }
  
  
  
  #converting to data.frame
  col_annotation <- as.data.frame(col_annotation)
  
  ##fixing annotation legend if it is provided
  if(!is.null(annotation_legend)){
    names(col_annotation) <- annotation_legend
  }
  
  #converting colors to named list with name matching the col_annotation
  annotlist <- list(annotcolors)
  names(annotlist) <- colnames(col_annotation)
  
  if(!is.null(col_fun)) {
    Heatmap(z_genes,
            col = col_fun,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_columns,
            row_names_gp = gpar(fontsize = row_name_size), 
            column_names_gp = gpar(fontsize = col_name_size),
            column_title = column_title, 
            column_title_gp = gpar(box_fill = box_fill),
            heatmap_legend_param = list(title = legend_title),
            top_annotation = HeatmapAnnotation(df = col_annotation,
                                               which = 'column',
                                               col = annotlist)
    )
  } else {
    
  Heatmap(z_genes,
          cluster_rows = cluster_rows,
          cluster_columns = cluster_columns,
          row_names_gp = gpar(fontsize = row_name_size), 
          column_names_gp = gpar(fontsize = col_name_size),
          column_title = column_title, 
          column_title_gp = gpar(box_fill = box_fill),
          heatmap_legend_param = list(title = legend_title),
          top_annotation = HeatmapAnnotation(df = col_annotation,
                                              which = 'column',
                                             col = annotlist)
  )
  }

  } else {

    if(!is.null(col_fun)){
      Heatmap(z_genes,
              col = col_fun,
              cluster_rows = cluster_rows,
              cluster_columns = cluster_columns,
              row_names_gp = gpar(fontsize = row_name_size), 
              column_names_gp = gpar(fontsize = col_name_size),
        column_title = column_title, 
        column_title_gp = gpar(box_fill = box_fill),
        heatmap_legend_param = list(title = legend_title)
        )
    } else {
      Heatmap(z_genes,
              cluster_rows = cluster_rows,
              cluster_columns = cluster_columns,
              row_names_gp = gpar(fontsize = row_name_size), 
              column_names_gp = gpar(fontsize = col_name_size),
              column_title = column_title, 
              column_title_gp = gpar(box_fill = box_fill),
              heatmap_legend_param = list(title = legend_title)
      )
    }
}
}

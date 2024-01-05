#'Create a volcano plot
#'
#'Creates a volcano plot based on a gene list with logFC and (fdr or PValue) columns
#'Works best with edgeR results. 
#'
#'Dependencies are:
#' library(ggplot2)
#' library(dplyr)
#' library(ggrepel)
#'
#'
#'date 24/10/2022
#'@param table Input table with logFC and (fdr or PValue) columns
#'@param title Title of the plot
#'@param significanceMeasure the significance Measure to be used for statistical significance. Must match the corresponding column' name
#'@param thresSig The statistical Measure significance value
#'@param thresLogFC The Log Fold change significance threshold
#'@param colorUp The upregulated genes colors
#'@param colorDown The downregulated genes colors
#'@param colorNeutral The color of non-statistically significant or non-logFC significant genes
#'@param xlim The limits of the x-axis
#'@param ylim The limits of the y-axis
#'@param breaks_x_default If you wish to use the default x-axis breaks
#'@param breaks_y_default If you wish to use the default y-axis breaks
#'@param genes_to_label The list of genes to label
#'@param label_column_names The columns of the gene labels
#'@param ids If you wish to have a mixed type of ids (i.e. Ensembl and HGNC symbols) or a single type
#'@param label_size The size of the gene labels
#'@param label_alpha The transparency of the gene labels
#'@param max_overlaps The number of geom_label_repel max.overlaps parameter
#'
#'@return ggplot object
#'@export


volcano_sen <- function(table, title = "Cond_Up vs Cond_Down", significanceMeasure = c("fdr","PValue"), thresSig = 0.05, thresLogFC = 0.58, 
                        colorUp="darkgreen",colorDown = "steelblue", colorNeutral = "lightgrey" , 
                        xlim = c(-10,10), ylim = c(0,15), breaks_x_default = FALSE, breaks_y_default = FALSE,
                        genes_to_label = NULL, label_columns_names = NULL, ids=c("single","mixed"),
                        label_size = 4, label_alpha = 0.7, max_overlaps = 40)
  
  {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Package \"ggrepel\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  
  if (length(significanceMeasure) >1) {
    significanceMeasure <- 'fdr'
  }
  ##Volcano code
  upregulated <- table %>%
    filter(get(significanceMeasure) < thresSig, logFC > thresLogFC) %>%
    mutate(
      change = 'Upregulated'
    )
  
  downregulated <- table %>%
    filter(get(significanceMeasure) < thresSig, logFC < -thresLogFC) %>%
    mutate(
      change = 'Downregulated'
    )
  
  insignificant <- table %>%
    filter(get(significanceMeasure) > thresSig | get(significanceMeasure) < thresSig & logFC <= thresLogFC & logFC >= -thresLogFC ) %>%
    mutate(
      change = 'Insignificant'
    )
  
  table <- rbind(upregulated,downregulated,insignificant)
  
  ##creating base plot
g <-ggplot(table, aes(logFC, -log10(get(significanceMeasure)),color = as.factor(change))) + 
    geom_point() +
    geom_hline(yintercept= -log10(thresSig))+
    geom_vline(xintercept = -thresLogFC)+
    geom_vline(xintercept = thresLogFC)+
    scale_color_manual(values = c('Upregulated' = colorUp, 'Downregulated' = colorDown, 'Insignificant' = colorNeutral))+
    labs(title = title, x = 'LogFC', y = paste0('-Log10(',toupper(significanceMeasure),')'), color = 'Differential \nExpression \nResult')+
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5))


##adding x scale
if (!breaks_x_default) {
  g <- g + scale_x_continuous(breaks = c(seq(xlim[1],xlim[2],2),-thresLogFC,0,thresLogFC), labels = c(seq(xlim[1],xlim[2],2),-thresLogFC,0,thresLogFC), limits = xlim)
} else {
  g <- g + scale_x_continuous(limits = xlim) 
}

##adding y scale
if (!breaks_y_default) {
  g <- g + scale_y_continuous(breaks = c(seq(ylim[1],ylim[2],2),-log10(thresSig)), labels = c(seq(ylim[1],ylim[2],2),paste0(toupper(significanceMeasure)," = ", thresSig)), limits = ylim)
} else {
  g <- g + scale_y_continuous(limits = ylim)
}

##adding genes with labels
if (!is.null(genes_to_label)) {
  
  if(length(ids) >1) {
    ids <- "single"
  }
  
  if(is.null(label_columns_names)){
    if(sum(genes_to_label %in% rownames(table))==0){
      stop("Genes not found in rownames. Please provide the column name of the label column.")
    } else {
    g <- g + geom_label_repel(
      data = table[rownames(table) %in% genes_to_label,],
      aes(label= rownames(table[rownames(table) %in% genes_to_label,]),color = as.factor(change)),
      min.segment.length = unit(0.01,"cm"),
      size = label_size,
      alpha = label_alpha,
      max.overlaps = max_overlaps
      )
    }
    
  }  else if (!is.null(label_columns_names) & !(label_columns_names %in% colnames(table))) {
    stop("Please provide valid column names")
    
  } else {
    
  if (ids == "mixed"){
    if (length(label_columns_names) >2){
      stop("Please choose only two columns for mixed IDs")
    } else if (length(label_columns_names <2)) {
      stop("Please choose two columns for IDs")
    } else {
      
      if(sum(is.null(table[,label_columns_names[1]]))< sum(is.null(table[,label_columns_names[2]]))){
        
      table$mixed <- table[,label_columns_names[2]]
      table$mixed(is.null(table$mixed)) <- table[is.null(table$mixed),label_columns_names[1]]

      } else {
        table$mixed <- table[,label_columns_names[1]]
        table$mixed(is.null(table$mixed)) <- table[is.null(table$mixed),label_columns_names[2]]
      }
      
      g <- g + geom_label_repel(
        data = table[table$mixed %in% genes_to_label,],
        aes(label= mixed,color = as.factor(change)),
        min.segment.length = unit(0.01,"cm"),
        size = label_size,
        alpha = label_alpha,
        max.overlaps = max_overlaps
      )
      
    }
    
  } else if (ids =="single") {
    if (length(label_columns_names) !=1){
      stop("Please choose only one column for hgnc/ensembl IDs")
    } else {
      
      keep <- match(genes_to_label,table[,label_columns_names])
      
      g <- g + geom_label_repel(
        data = table[keep,],
        aes(label= table[keep,label_columns_names],color = as.factor(change)),
        min.segment.length = unit(0.01,"cm"),
        size = label_size,
        alpha = label_alpha,
        max.overlaps = max_overlaps
      )
    }
  } else {
      stop("ids must be one of mixed or single")
    }
    } 
  
}

##return value
g

}

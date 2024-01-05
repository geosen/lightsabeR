#'Run an interactive ComplexHeatmap session
#'
#'Dependencies are:
#' library(shiny)
#' library(ComplexHeatmap)
#' library(circlize)
#' library(dplyr)
#' library(openxlsx)
#' library(edgeR)
#' library(stringr)
#' 
#'date 14/02/2023
#'
#'
#'@return Opens a shiny session
#'@export

shinyHeatmap_sen <- function(){
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(openxlsx)
library(edgeR)
library(stringr)
library(shiny)

  options(shiny.maxRequestSize = 200*1024^2)
ui <- fluidPage(
  splitLayout(cellWidths = c('45%','65%'),
    mainPanel(
      splitLayout(cellWidths = 300,
        mainPanel(
          actionButton('click','Generate Heatmap'),
          fileInput('expression_table','Expression file\n(xlsx format)', accept = '.xlsx'),
          checkboxInput('filter',label='Remove lowly-expressed genes',value = TRUE),
      selectInput('normalization','Normalize the values with', 
                  choices = c('None','CPM','LogCPM','Z-score','Log only'), selected = 'Z-score'),
      selectInput('colors_scale_low', label = 'Color low', choices = colors(), selected = 'green'),
      selectInput('colors_scale_mid', label = 'Color mid', choices = c('white','black'), selected = 'white'),
      selectInput('colors_scale_high', label = 'Color high', choices = colors(), selected = 'purple')
      ,
      fileInput('metadata','Metadata file (xlsx format)', accept = '.xlsx'),
      uiOutput('annot_sel'),
      selectInput('annot_row_or_col','Annotate rows or columns?', choices = c('column','row')),
      selectInput('annotcolors', label = 'Annotation colors', choices = colors(), multiple = TRUE)
      ),
      mainPanel(
      textInput('heatmap_title','Heatmap title', value = ''),
      textInput('legend_title','Legend title', value = NULL),
      textInput('legend_annotation_title','Annotation title', value = 'Annotation'),
      
      sliderInput('colsize','Column labels size', min=0, max=20, step=1, value=8 ),
      sliderInput('colkm','Number of Column clusters', min=1, max=10, step=1, value=1 ),
      checkboxInput('coldend','Show column dendrogram'),

      sliderInput('rowsize','Row labels size', min=0, max=20, step=1, value=8 ),
      sliderInput('rowkm','Number of Row clusters', min=1, max=10, step=1, value=1 ),
      checkboxInput('rowdend','Show row dendrogram'),
      textAreaInput('gene_list','List of genes to plot(1 gene per line)', value = NULL),
      textAreaInput('highlight_gene','List of genes to highlight(1 gene per line)', value = NULL)
    )
    )
    ),
    
    mainPanel(
     plotOutput('heatmap', height= '600px')
    )
  )
  
)

server <- function(input, output, session) {
    
    ####expression file####
    rval_counts <-   reactive({
      
      counts_read <- read.xlsx(input$expression_table$datapath, rowNames = TRUE, colNames=TRUE)
      
      if(input$filter) {
        counts_read %>%
          filter(
            (rowSums(. > 10)) >= ceiling(0.10*dim(counts_read)[2])
          )
      } else {
        counts_read
      }
      
      })
    
    ####normalization type####
    rval_norm_counts <- reactive({
      
    if(input$normalization == 'CPM'){
      d <- DGEList(counts = rval_counts())
      d <- calcNormFactors(d)
      as.data.frame(cpm(d))
    } else if(input$normalization == 'LogCPM'){
      d <- DGEList(counts = rval_counts())
      d <- calcNormFactors(d)
      as.data.frame(log(cpm(d)+1))
    } else if(input$normalization == 'Z-score'){
      as.data.frame(t(scale(t(rval_counts()))))
    } else if(input$normalization == 'CPM Z-score'){
      d <- DGEList(counts = rval_counts())
      d <- calcNormFactors(d)
      as.data.frame(t(scale(t(cpm(d)))))
    } else if(input$normalization == 'None'){
      rval_counts()
    } else if(input$normalization == 'Log only'){
      as.data.frame(log(rval_counts() + 1))
    }
    
    })
    
    
    ####legend title####
    legend_title <- reactive ({
      
    if(input$normalization == 'CPM' & input$legend_title == ''){
      'CPM\nvalues'
    } else if(input$normalization == 'LogCPM' & input$legend_title == ''){
      'Log Normalized\nCPM values'
    } else if(input$normalization == 'Z-score' & input$legend_title == ''){
      'Z-score\nvalues'
    } else if(input$normalization == 'CPM Z-score' & input$legend_title == ''){
      'CPM\nZ-score\nvalues'
    } else if(input$normalization == 'None' & input$legend_title == ''){
      'Raw\nvalues'
    } else if(input$normalization == 'Log only' & input$legend_title == ''){
      'Log\nNormalized\nValues'
    } else {
      input$legend_title
    }
    
    })
    
    ####metadata file####
    metadata <- reactive({
      
    if(!is.null(input$metadata)){
      metadata <- read.xlsx(input$metadata$datapath, rowNames = TRUE, colNames=TRUE)
    } 
    
    })
    
    ##render selector
    output$annot_sel <- renderUI({
      validate(need(input$metadata,message = 'Metadata must be provided to choose annotation variable', label = 'Metadata'))
      tagList(selectInput('annotcol', 'Annotation column(s)\nof metadata file', choices = colnames(metadata())),
              checkboxInput('annotate','Add annotation column'))
      })
    
    
    
    ####Gene list####
    expression_isolated <- reactive({
      
    if(input$gene_list!=''){
      genes_list <- unlist(strsplit(input$gene_list, split = '\n'))
      rval_norm_counts() %>% dplyr::filter(rownames(.) %in% genes_list)
    } else {
      rval_norm_counts()
    }
    
    })
    
    ####Expression color scale####
    CS <- reactive({
      CS <-colorRamp2(
                        c(min(expression_isolated(), na.rm = TRUE),
                        0,
                        max(expression_isolated(), na.rm = TRUE)
                        ),
                      c(input$colors_scale_low,
                        input$colors_scale_mid,
                        input$colors_scale_high)
                      )
    })
    
    
    ##Generate heatmap
    observeEvent(input$click, {
      
    output$heatmap <-  renderPlot( expr = {
      
      ####Highlight gene list####
      if(input$highlight_gene!=''){
        highlight_gene <- unlist(strsplit(input$highlight_gene, split = '\n'))
        new_rownames <- case_when(rownames(expression_isolated()) %in% highlight_gene ~ rownames(expression_isolated()),
                                  !rownames(expression_isolated()) %in% highlight_gene ~ '')
        row_labels = structure(paste0(new_rownames), names = rownames(expression_isolated()))
      } else {
        highlight_gene <- NULL
      }
      
      ####Plot heatmap####
    if(!is.null(input$annotcol) & input$annotate){
      
        if(length(unlist(input$annotcolors))==length(levels(as.factor(metadata()[,input$annotcol])))){
        
        annotcolors <- unlist(input$annotcolors)
        ####Annotation colors####
        names(annotcolors)<-levels(as.factor(metadata()[,input$annotcol]))
        ####Annotation list####
        annotlist <- list(annotcolors)
        names(annotlist) <- input$legend_annotation_title
        
        ####Annotation data.frame####
        annotdf <- data.frame(metadata()[,input$annotcol]) 
        names(annotdf) <- input$legend_annotation_title
        } else {
        annotcolors <- sample(colors(),
                              size = length(levels(as.factor(metadata()[,input$annotcol]))),
                              replace = FALSE)
        ####Annotation colors####
        names(annotcolors)<-levels(as.factor(metadata()[,input$annotcol]))
        ####Annotation list####
        annotlist <- list(annotcolors)
        names(annotlist) <- input$legend_annotation_title
          
        ####Annotation data.frame####
        annotdf <- data.frame(metadata()[,input$annotcol]) 
        names(annotdf) <- input$legend_annotation_title
        }
        
        if (is.null(highlight_gene)) {
        Heatmap(as.matrix(expression_isolated()),
            col = CS(),
            row_names_gp = gpar(fontsize = input$rowsize),
            column_names_gp = gpar(fontsize = input$colsize),
            column_km = input$colkm,
            row_km= input$rowkm,
            show_row_dend = input$rowdend,
            show_column_dend = input$coldend,
            column_title = input$heatmap_title,
            heatmap_legend_param = list(title = legend_title()),
            top_annotation = HeatmapAnnotation(df = annotdf, which = input$annot_row_or_col,
                                               col = annotlist )
        )
        } else {
          Heatmap(as.matrix(expression_isolated()),
                  col = CS(),
                  row_labels = row_labels[rownames(expression_isolated())], 
                  row_names_gp = gpar(fontsize = input$rowsize),
                  column_names_gp = gpar(fontsize = input$colsize),
                  column_km = input$colkm,
                  row_km= input$rowkm,
                  show_row_dend = input$rowdend,
                  show_column_dend = input$coldend,
                  column_title = input$heatmap_title,
                  heatmap_legend_param = list(title = legend_title()),
                  top_annotation = HeatmapAnnotation(df = annotdf, which = input$annot_row_or_col,
                                                     col = annotlist ))
        }
          
        
       
      
    } else { ####DONT TOUCH BELOW HERE
      
      
      if (is.null(highlight_gene)) {
      Heatmap(as.matrix(expression_isolated()), 
              col = CS(), 
              row_names_gp = gpar(fontsize = input$rowsize),
              column_names_gp = gpar(fontsize = input$colsize),
              column_km = input$colkm,
              row_km= input$rowkm,
              show_row_dend = input$rowdend,
              show_column_dend = input$coldend,
              column_title = input$heatmap_title,
              heatmap_legend_param = list(title = legend_title())
              )
      } else  {
        Heatmap(as.matrix(expression_isolated()), 
                col = CS(),
                row_labels = row_labels[rownames(expression_isolated())], 
                row_names_gp = gpar(fontsize = input$rowsize),
                column_names_gp = gpar(fontsize = input$colsize),
                column_km = input$colkm,
                row_km= input$rowkm,
                show_row_dend = input$rowdend,
                show_column_dend = input$coldend,
                column_title = input$heatmap_title,
                heatmap_legend_param = list(title = legend_title()))
      }
    }
    
  })
  })
    
}

shinyApp(ui, server)
}
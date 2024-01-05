#'Run an interactive PCA plot session
#'
#'Dependencies are:
#' library(shiny)
#' library(DT)
#' library(plotly)
#' library(openxlsx)
#' library(tidyverse)
#' library(sentisR)
#' library(FactoMineR)
#' library(edgeR)
#' 
#'date 10/12/2023
#'
#'
#'@return Opens a shiny session
#'@export




shinyPCA_sen <- function() {
  library(shiny)
  library(DT)
  library(plotly)
  library(openxlsx)
  library(tidyverse)
  library(sentisR)
  library(FactoMineR)
  library(edgeR)
  
  
  options(shiny.maxRequestSize = 200*1024^2)
  
  ui <- fluidPage(
    sidebarPanel(
      fileInput('counts',label= 'Input you count file (.xlsx)', accept = '.xlsx'),
      checkboxInput('filter',label='Remove lowly-expressed genes',value = TRUE),
      uiOutput('selectNormalization'),
      fileInput('metadata',label = 'Input your metadata table (.xlsx)',accept = '.xlsx'),
      uiOutput('selectColor'),
      checkboxInput('PCA',label='Create PCA plot',value = FALSE),
      uiOutput('PCA_variables')
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(title='Counts',
                 DTOutput('Counts')
        ),
        tabPanel(title='Metadata',
                 DTOutput('Metadata')
        ),
        tabPanel(title='PCA Plot',
                 plotlyOutput('PCA_plot')
        ),
        tabPanel(title='PCA Data',
                 DTOutput('DTPCA')
        )
      )
    )
    
    
    
  )
  
  server <- function(input, output, session) {
    #Counts
    rval_counts <- reactive({
      counts_read <- read.xlsx(xlsxFile = input$counts$datapath, rowNames = TRUE, colNames = TRUE)
      
      #filtering - genes with more than 10 raw counts in at least 10% of the samples (rounded up)
      if(input$filter) {
        counts_read %>%
          filter(
            (rowSums(. > 10)) >= ceiling(0.10*dim(counts_read)[2])
          )
      } else {
        counts_read
      }
      
    })
    
    output$selectNormalization <- renderUI({validate(need(input$counts, message = ''));
      selectInput('norm_method', label='Select normalization method \n(logCPM recommened for PCAs)', 
                  choices = c('None','CPM (EdgeR)','LogCPM (EdgeR)','Z-Score'),
                  selected = 'None')
    })
    
    
    rval_norm_counts <- reactive({
      if(input$norm_method == 'None')
      {
        rval_counts()
      } else if(input$norm_method == 'CPM (EdgeR)') {
        
        d <- DGEList(counts=rval_counts())
        d <- calcNormFactors(d)
        cpm(d)
        
      } else if(input$norm_method == 'LogCPM (EdgeR)') {
        d <- DGEList(counts=rval_counts())
        d <- calcNormFactors(d)
        log(cpm(d)+1)      
        
      } else if(input$norm_method == 'Z-Score') {
        t(scale(t(rval_counts())))
      }
    })
    
    output$Counts <- renderDT({ validate(need(input$norm_method, 
                                              message = 'Counts must be provided. If providing a large file, it may take some time.'));
      rval_norm_counts()
    })
    
    
    
    #Metadata
    rval_metadata <- reactive({
      read.xlsx(xlsxFile = input$metadata$datapath, rowNames = TRUE, colNames = TRUE) 
    })
    
    output$Metadata <- renderDT({ validate(need(input$metadata, message = 'Metadata must be provided. If providing a large file, it may take some time.'));
      rval_metadata()
    })
    
    output$selectColor <- renderUI({validate(need(input$metadata, message = ''));
      selectInput('color_factor', label='Select PCA color variable', choices = colnames(rval_metadata()))
    })
    
    #PCA
    pca_data <- reactive({
      
      PCA(t(rval_norm_counts()), scale.unit = T, graph = F)
      
    })
    
    
    output$DTPCA <- renderDT({
      validate(need(input$color_factor, 
                    message = 'Counts and metadata must be provided. If providing a large file, it may take some time.'));
      if(input$PCA){
        pca_data()$eig
      } else {
        data.frame('Column' = c('Press create PCA plot to analyze'))
      }
      
    })
    
    #generating PCA options menu
    output$PCA_variables <- renderUI({
      if(input$PCA){
        tagList(
          numericInput('PCx',label = 'PCx', value = 1, min = 1, max = dim(pca_data()$ind$coord)[2]),
          numericInput('PCy',label = 'PCy', value = 2, min = 1, max = dim(pca_data()$ind$coord)[2]),
          textInput('color_factor_name', label = 'Color legend name', value = 'Color'),
          selectInput('colors_list', label = 'Group colors', choices = colors(),multiple = TRUE),
          selectInput('shape_vector', label = 'Group shapes', choices = c(1:20),multiple = TRUE, selected = 16),
          numericInput('stroke', label = 'Point width', value = 2, min = 1, max = 20),
          textInput('title',label = 'Plot title', value = 'PCA'),
          numericInput('point_alpha',label = 'Point transparency', value = 0.7, min = 0, max = 1),
          numericInput('point_size',label = 'Point size', value = 7, min = 0, max = 20),
          checkboxInput('add_ellipse',label = 'Add ellipse', value = FALSE)
          
        )
      } 
      
    })
    
    output$PCA_plot <- renderPlotly({validate(need(input$color_factor, 
                                                   message = 'Counts, metadata and a color factor must be provided. If providing a large file, it may take some time.'));
      if(input$PCA){
        #### pca_sen_vis -- adjusted
        pcacpm <- pca_data()
        
        #creating color list
        if (is.null(input$colors_list)) {
          color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), 
                                            invert = T)]
          color <- color[!grepl("white", color)]
          color_choice <- seq(1, 433, 443/length(levels(factor(rval_metadata()[,input$color_factor])))) + 
            70
          colors_plot <- color[color_choice]
        }
        else {
          colors_plot <- input$colors_list
        }
        
        #checking color factor length
        if (length(rval_metadata()[,input$color_factor]) != dim(pcacpm$ind$coord)[1]) {
          print("Color factor length is not equal to sample size")
        }
        
        #plotting
        data_pca <- as.data.frame(pcacpm$ind$coord) %>%
          mutate(`Point name` = rownames(.))
        
        if (!is.null(input$shape_vector)) {
          g <- ggplot(data_pca, aes(data_pca[, input$PCx], 
                                    data_pca[, input$PCy], 
                                    color = factor(rval_metadata()[,input$color_factor]),
                                    label = `Point name`)
          ) +
            geom_point(alpha = input$point_alpha, 
                       size = input$point_size, 
                       #key_glyph = "point", 
                       shape = input$shape_vector, 
                       stroke = input$stroke) + 
            labs(title = input$title, 
                 x = paste0("PC", input$PCx, "(", round(pcacpm$eig[input$PCx, 2], 2), ")"), 
                 y = paste0("PC", input$PCy, "(", round(pcacpm$eig[input$PCy, 2], 2), ")")) + 
            scale_color_manual(name = input$color_factor_name,
                               labels = levels(factor(rval_metadata()[,input$color_factor])), 
                               values = colors_plot) + 
            theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
        }
        else {
          g <- ggplot(data_pca, aes(data_pca[, input$PCx], 
                                    data_pca[, input$PCy], 
                                    color = factor(rval_metadata()[,input$color_factor]),
                                    label = `Point name`)
          ) +
            geom_point(alpha = input$point_alpha, 
                       size = input$point_size, 
                       #key_glyph = "point", 
                       shape = input$shape_vector) + 
            labs(title = input$title, 
                 x = paste0("PC", input$PCx, "(", round(pcacpm$eig[input$PCx, 2], 2), ")"), 
                 y = paste0("PC", input$PCy, "(", round(pcacpm$eig[input$PCy, 2], 2), ")")) + 
            scale_color_manual(name = input$color_factor_name,
                               labels = levels(factor(rval_metadata()[,input$color_factor])), 
                               values = colors_plot) + 
            theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
        }
        if (input$add_ellipse) {
          g <- g + stat_ellipse()
        }
        ####
        
        ggplotly(g, tooltip = 'label')
      }
      
      
    })
    
  }
  
  shinyApp(ui, server)
}
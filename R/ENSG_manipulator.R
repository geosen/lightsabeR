#'Run an interactive session to remove ENSG ID version number
#'
#'Dependencies are:
#' library(shiny)
#' library(biomaRt)
#' 
#'date 24/02/2023
#'
#'to specify repositories for deployment use options(repos = BiocManager::repositories())
#'
#'
#'@return Opens a shiny session
#'@export
ENSG_manipulator_shiny <- function() {
  library(shiny)
  options(shiny.maxRequestSize = 200*1024^2)
  ui <- fluidPage(
    splitLayout(cellWidths = c('45%','65%'),
                mainPanel(textAreaInput('gene_list',
                                        'List of genes (1 gene per line)', 
                                        value = NULL),
                          selectInput('what','What to do',choices = c('Remove version','Convert to Gene Symbol','Convert to Ensembl ID'),
                                      selected = 'Remove version'),
                          sliderInput('ensembl_version','Choose your ensembl version:',min = 0, max = 106, value = 105, step = 1),
                          selectInput('organism','Organism',choices = c('Homo sapiens','Mus musculus'),
                                      selected = 'Homo sapiens'),
                          actionButton('click','Run')
                ),
                mainPanel(
                  verbatimTextOutput('deversioned')
                )
    )
    
  )
  
  server <- function(input, output, session) {
    #output
    library(biomaRt)
    
    observeEvent(input$click, {
      
      if(input$what == 'Remove version') {
        output$deversioned <-  renderText( expr = {
          genes <- unlist(strsplit(input$gene_list, split = '\n'))
          paste(gsub('\\.[[:digit:]]{1+}$','',genes), collapse = '\n')
        })
      } else if (input$what == 'Convert to Gene Symbol') {
        
        output$deversioned <-  renderText( expr = {
          
          withProgress(message = 'Converting...', value = 0, {
            
            genes <- unlist(strsplit(input$gene_list, split = '\n'))
            
            
            if(input$organism == 'Homo sapiens') {
              dataset = 'hsapiens_gene_ensembl'
              
              if (input$ensembl_version != 0){
                mart = useEnsembl(biomart = 'ensembl', dataset = dataset, version = input$ensembl_version)  
                incProgress(1/3, detail = "Mart ready, getting gene symbols...")
                conv_table <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id_version'), filters = 'ensembl_gene_id_version', values = genes, mart = mart)
                incProgress(1/3, detail = "Symbols ready, preparing output")
                conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id_version),]
                genes2 <- conv_table$hgnc_symbol[match(genes,conv_table$ensembl_gene_id_version)]
                incProgress(1/3, detail = "All done!")
                paste(genes2, collapse = '\n')
              } else {
                genes <- gsub('\\.[[:digit:]]{1+}$','',genes)
                mart = useEnsembl(biomart = 'ensembl', dataset = dataset)
                incProgress(1/3, detail = "Mart ready, getting gene symbols...")
                conv_table <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id', values = genes, mart = mart)
                incProgress(1/3, detail = "Symbols ready, preparing output")
                conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id),]
                genes2 <- conv_table$hgnc_symbol[match(genes,conv_table$ensembl_gene_id)]
                incProgress(1/3, detail = "All done!")
                paste(genes2, collapse = '\n')
              }
            } else if (organism == 'Mus musculus') {
              
              dataset = 'mmusculus_gene_ensembl'
              if (input$ensembl_version != 0){
                mart = useEnsembl(biomart = 'ensembl', dataset = dataset, version = input$ensembl_version)  
                incProgress(1/3, detail = "Mart ready, getting gene symbols...")
                conv_table <- getBM(attributes = c('mgi_symbol','ensembl_gene_id_version'), filters = 'ensembl_gene_id_version', values = genes, mart = mart)
                incProgress(1/3, detail = "Symbols ready, preparing output")
                conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id_version),]
                genes2 <- conv_table$mgi_symbol[match(genes,conv_table$ensembl_gene_id_version)]
                incProgress(1/3, detail = "All done!")
                paste(genes2, collapse = '\n')
              } else {
                genes <- gsub('\\.[[:digit:]]{1+}$','',genes)
                mart = useEnsembl(biomart = 'ensembl', dataset = dataset)
                incProgress(1/3, detail = "Mart ready, getting gene symbols...")
                conv_table <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id', values = genes, mart = mart)
                incProgress(1/3, detail = "Symbols ready, preparing output")
                conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id),]
                genes2 <- conv_table$mgi_symbol[match(genes,conv_table$ensembl_gene_id)]
                incProgress(1/3, detail = "All done!")
                paste(genes2, collapse = '\n')
              }
            }
            
            
            
            
            
          })
        })
      } else if (input$what == 'Convert to Ensembl ID') {
        
        output$deversioned <-  renderText( expr = {
          
          withProgress(message = 'Converting...', value = 0, {
            
            genes <- unlist(strsplit(input$gene_list, split = '\n'))
            
            
            if(input$organism == 'Homo sapiens') {
              dataset = 'hsapiens_gene_ensembl'
              
              if (input$ensembl_version != 0){
                mart = useEnsembl(biomart = 'ensembl', dataset = dataset, version = input$ensembl_version)  
                incProgress(1/3, detail = "Mart ready, getting ensembl ids...")
                conv_table <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id_version'), filters = 'hgnc_symbol', values = genes, mart = mart)
                incProgress(1/3, detail = "IDs ready, preparing output")
                conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id_version),]
                genes2 <- conv_table$ensembl_gene_id_version[match(genes,conv_table$hgnc_symbol)]
                incProgress(1/3, detail = "All done!")
                paste(genes2, collapse = '\n')
              } else {
                genes <- gsub('\\.[[:digit:]]{1+}$','',genes)
                mart = useEnsembl(biomart = 'ensembl', dataset = dataset)
                incProgress(1/3, detail = "Mart ready, getting ensembl ids...")
                conv_table <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'), filters = 'hgnc_symbol', values = genes, mart = mart)
                incProgress(1/3, detail = "IDs ready, preparing output")
                conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id),]
                genes2 <- conv_table$ensembl_gene_id[match(genes,conv_table$hgnc_symbol)]
                incProgress(1/3, detail = "All done!")
                paste(genes2, collapse = '\n')
              }
            } else if (organism == 'Mus musculus') {
              
              dataset = 'mmusculus_gene_ensembl'
              if (input$ensembl_version != 0){
                mart = useEnsembl(biomart = 'ensembl', dataset = dataset, version = input$ensembl_version)  
                incProgress(1/3, detail = "Mart ready, getting ensembl ids...")
                conv_table <- getBM(attributes = c('mgi_symbol','ensembl_gene_id_version'), filters = 'mgi_symbol', values = genes, mart = mart)
                incProgress(1/3, detail = "IDs ready, preparing output")
                conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id_version),]
                genes2 <- conv_table$ensembl_gene_id_version[match(genes,conv_table$mgi_symbol)]
                incProgress(1/3, detail = "All done!")
                paste(genes2, collapse = '\n')
              } else {
                genes <- gsub('\\.[[:digit:]]{1+}$','',genes)
                mart = useEnsembl(biomart = 'ensembl', dataset = dataset)
                incProgress(1/3, detail = "Mart ready, getting ensembl ids...")
                conv_table <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'), filters = 'mgi_symbol', values = genes, mart = mart)
                incProgress(1/3, detail = "IDs ready, preparing output")
                conv_table <- conv_table[!duplicated(conv_table$ensembl_gene_id),]
                genes2 <- conv_table$ensembl_gene_idd[match(genes,conv_table$mgi_symbol)]
                incProgress(1/3, detail = "All done!")
                paste(genes2, collapse = '\n')
              }
            }
            
            
            
            
            
          })
        })
      }
      
    })
    
  }
  
  
  shinyApp(ui, server)
  
}

library(shiny)
library(magrittr)
library(devtools)
library(tidyr)
library(dplyr)


setwd("/Users/noam/Desktop/Technion/bioinformatics project/xCell2")
options(shiny.maxRequestSize = 100*1024^2) # 100MB

ui <- fluidPage(
  titlePanel("xCell2 Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("ref_file", "Upload reference file"),
      fileInput("mix_file", "Upload mixture file"),
      actionButton("run_button", "Run Analysis")
    ),
    mainPanel(
      tableOutput("output_table")
    )
  )
)

server <- function(input, output) {
  
  # Load reference file
  ref_data <- reactive({
    readRDS(input$ref_file$datapath)
  })
  
  # Load mixture file
  mix_data <- reactive({
    if (is.null(input$mix_file))
      return(NULL)
    bulk <- read.table(input$mix_file$datapath, check.names = FALSE, row.names = 1, header = TRUE)
    return(bulk)
  })
  
  # Run xCell2Analysis when run button is clicked
  xcell_output <- eventReactive(input$run_button, {
    if (is.null(ref_data()) || is.null(mix_data()))
      return(NULL)
    xCell2Analysis(mix = mix_data(), xcell2ref = ref_data())
  })
  
  # Display xCell2Analysis output as a table
  output$output_table <- renderTable({
    if (is.null(xcell_output()))
      return(NULL)
    xcell_output()
  })
}

shinyApp(ui, server)

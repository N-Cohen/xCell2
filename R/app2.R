library(shiny)
library(dplyr)
library(readr)
library(devtools)
library(shinyjs)
library(xCell2)
options(shiny.maxRequestSize = 100*1024^2) # 100MB

setwd("/Users/noam/Desktop/Technion/bioinformatics project/xCell2")


# Define UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("xCell2Train App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("ref", "Reference file"),
      fileInput("labels", "Labels file"),
      actionButton("run", "Run xCell2Train"),
      br(),
      helpText("Note: Both files must be in csv format."),
      br(),
      tags$div(
        style = "margin-bottom: 10px",
        HTML("Instructions for the labels file:"),
        HTML("<ol>
              <li>The first column 'ont' for cell type
                (i.e., 'CL:0000775', 'CL:0000576') <b>**</b></li>
              <li>The second column 'label' can be any string that describes the cell type
                (i.e., 'Neutrophils', 'Monocytes') <b>***</b></li>
              <li>The third column 'sample' is the sample name in the reference - it must be the same as
                the reference column names
                all(bp_labels$sample == colnames(bp_ref))</li>
              <li>The fourth column 'dataset' is the name of the reference/dataset the sample originate from</li>
            </ol>
            <b>**</b> Use NA if you don't know the cell type ontology<br>
            <b>***</b> You can use any label you like but it is recommended to use the same labels as in onyology<br>
            <a href='https://www.ebi.ac.uk/ols/index' target='_blank'>an example screenshot</a>"),
        tags$img(src = "/Users/noam/Desktop/Technion/bioinformatics project/xCell2/xCell2/xCell2screenshot.png", style = "max-width: 100%; height: auto;")
      ),
      downloadButton("download", "Download results")
    ),
    mainPanel(
      tableOutput("output_table")
    )
  )
)

# Define server
server <- function(input, output) {

  # Load ref file
  ref_data <- reactive({
    if (is.null(input$ref_file))
      return(NULL)
    ref_matrix <- read.table(input$ref_file$datapath, check.names = FALSE, row.names = 1, header = TRUE)
    return(ref_matrix)
  })

  # Load labels file
  labels_data <- reactive({
    if (is.null(input$labels))
      return(NULL)
    labels_file <- read.table(input$labels$datapath, check.names = FALSE, row.names = 1, header = TRUE)
    return(ref_matrix)
  })

  # Run xCell2Train when run button is clicked
  xcell_output <- eventReactive(input$run_button, {
    if (is.null(ref_data()) || is.null(labels_data()))
      return(NULL)
    xCell2Train(ref = ref_matrix, labels = labels_file, data_type = "rnaseq", lineage_file = NULL, mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02)),
                probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                min_genes = 5, max_genes = 500)
  })

  # Display xCell2Analysis output as a table
  output$output_table <- renderTable({
    if (is.null(xcell_output()))
      return(NULL)
    xcell_output()
  })
}

# Run the app
shinyApp(ui, server)

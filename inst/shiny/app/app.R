library(shiny)
library(visNetwork)
library(getNC)


ui <- fluidPage(
  tags$head(tags$style(HTML("
      /* Sidebar styling */
      #sidebar {
        position: fixed;
        top: 70px;        /* below the title */
        left: 0;
        width: 260px;
        bottom: 0;
        padding: 15px;
        background-color: #f8f9fa;
        border-right: 1px solid #ccc;
        overflow-y: auto;
        z-index: 1000;
      }

      /* Main content area */
      #main_content {
        margin-left: 280px; /* match sidebar width + padding */
        padding: 10px;
      }
  "))),
  
  titlePanel("Interactive Gene Network (GLgraph)"),
  
  # Sidebar (fixed)
  div(
    id = "sidebar",
    h4("Settings"),
    
    numericInput(
      "top_k", 
      "Top K genes:", 
      value = 20,
      min = 5, 
      max = 200, 
      step = 1
    ),
    
    fileInput(
      "data_file", 
      "Upload raw matrix / Seurat:",
      accept = c(".csv", ".rds")
    ),
    
    helpText("Leave empty to use the default dataset.")
  ),
  
  # Main content
  div(
    id = "main_content",
    fluidRow(
      visNetworkOutput("network", height = "700px")
    )
  )
)

server <- function(input, output, session) {

  # -------------------------
  # 1. Reactive data loading
  data_reactive <- reactive({

    req(input$top_k)

    # -------------------------
    # Case 1: No file uploaded → use default data
    # -------------------------
    if (is.null(input$data_file)) {
      result <- auto_fit_glasso(x = NULL)

    } else {

      ext <- tools::file_ext(input$data_file$name)

      # -------------------------
      # Case 2: CSV
      # -------------------------
      if (ext == "csv") {
        mat <- as.matrix(read.csv(input$data_file$datapath, row.names = 1))
        result <- auto_fit_glasso(x = mat)

      # -------------------------
      # Case 3: RDS
      # -------------------------
      } else if (ext == "rds") {
        obj <- readRDS(input$data_file$datapath)

        if (inherits(obj, "Seurat")) {
          mat <- preprocess_seurat_data(obj)
          result <- auto_fit_glasso(x = mat)

        } else if (is.matrix(obj)) {
          result <- auto_fit_glasso(x = obj)

        } else {
          stop("RDS must contain a matrix or Seurat object.")
        }

      # -------------------------
      # Case 4: Unsupported file
      # -------------------------
      } else {
        stop("File type not supported.")
      }
    }

    # -------------------------
    # Clamp top_k to valid range for this dataset
    # -------------------------
    k <- input$top_k
    gene_count <- length(result$features)

    if (k > gene_count) {
      updateNumericInput(session, "top_k", value = gene_count)
    }

    if (k < 1) {
      updateNumericInput(session, "top_k", value = 1)
    }

    return(result)
  })


  # -------------------------
  # 2. Build GLgraph
  # -------------------------
  G_reactive <- reactiveVal(NULL)

  observeEvent(data_reactive(), {
    req(data_reactive()) 
    fit <- data_reactive()
    k   <- input$top_k
    G0  <- new_GLgraph(fit, k = k, target = 1)
    G_reactive(G0)
  })

  observeEvent(input$top_k, {
    req(data_reactive()) 
    fit <- data_reactive()
    k   <- input$top_k
    G0  <- new_GLgraph(fit, k = k, target = 1)
    G_reactive(G0)
  })

  # -------------------------
  # 3. Render network
  # -------------------------
  output$network <- renderVisNetwork({
    req(G_reactive())

    visnetwork_from_GLgraph(G_reactive()) %>%
      visNetwork::visEvents(
        click = "function(params) {
          if (params.nodes.length > 0) {
            Shiny.setInputValue('node_click',
              params.nodes[0], {priority:'event'});
          }
        }"
      )
  })

  # -------------------------
  # 4. Node-click knockout toggle
  # -------------------------
  observeEvent(input$node_click, {
    clicked <- input$node_click

    G_reactive(
      visnetwork_toggle_knock(G_reactive(), clicked)
    )

    output$network <- renderVisNetwork({
      visnetwork_from_GLgraph(G_reactive()) %>%
        visNetwork::visEvents(
          click = "function(params) {
            if (params.nodes.length > 0) {
              Shiny.setInputValue('node_click',
                params.nodes[0], {priority:'event'});
            }
          }"
        )
    })
  })
}


shinyApp(ui, server)
library(shiny)
library(visNetwork)
library(getNC)

ui <- fluidPage(
  tags$head(tags$style(HTML("
      /* Sidebar styling */
      #sidebar {
        position: fixed;
        top: 70px;
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
        margin-left: 280px;
        padding: 10px;
      }
  "))),
  
  titlePanel("Interactive Gene Network (GLgraph)"),
  
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
  
  div(
    id = "main_content",
    fluidRow(
      visNetworkOutput("network", height = "700px")
    )
  )
)


server <- function(input, output, session) {

  # =====================================================
  # 1. Load data (reactive)
  # =====================================================
  data_reactive <- reactive({

    req(input$top_k)

    # Default case: no file uploaded
    if (is.null(input$data_file)) {
      result <- auto_fit_glasso(x = NULL)

    } else {

      ext <- tools::file_ext(input$data_file$name)

      if (ext == "csv") {

        mat <- as.matrix(read.csv(input$data_file$datapath, row.names = 1))
        result <- auto_fit_glasso(x = mat)

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

      } else {
        stop("Unsupported file type.")
      }
    }

    # ----- clamp top_k safely -----
    k <- input$top_k
    gene_count <- length(result$features)

    if (k > gene_count) {
      updateNumericInput(session, "top_k", value = gene_count)
      return(result)   # stop observer!
    }
    if (k < 1) {
      updateNumericInput(session, "top_k", value = 1)
      return(result)   # stop observer!
    }

    result
  })


  # =====================================================
  # 2. Graph reactive
  # =====================================================
  G_reactive <- reactiveVal(NULL)

  observeEvent(data_reactive(), {
    fit <- data_reactive()
    k   <- input$top_k

    # validate k again (race-condition safety)
    k <- min(max(k, 1), length(fit$features))

    G_reactive(new_GLgraph(fit, k, target = 1))
  })

  observeEvent(input$top_k, {
    req(data_reactive())

    fit <- data_reactive()
    k   <- input$top_k

    # clamp again (critical)
    gene_count <- length(fit$features)
    if (k > gene_count) {
      updateNumericInput(session, "top_k", value = gene_count)
      return()
    }
    if (k < 1) {
      updateNumericInput(session, "top_k", value = 1)
      return()
    }

    G_reactive(new_GLgraph(fit, k, target = 1))
  })


  # =====================================================
  # 3. Render network
  # =====================================================
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


  # =====================================================
  # 4. Toggle knockout when node clicked
  # =====================================================
  observeEvent(input$node_click, {
    req(G_reactive())

    clicked <- input$node_click

    # update G
    G_reactive(visnetwork_toggle_knock(G_reactive(), clicked))

    # redraw graph
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
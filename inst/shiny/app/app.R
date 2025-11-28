library(shiny)
library(visNetwork)
library(getNC)

ui <- fluidPage(

  tags$head(tags$style(HTML("
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

      #main_content {
        margin-left: 280px;
        padding: 10px;
      }
  "))),

  titlePanel("Interactive Gene Network (GLgraph)"),

  # ---------------- Sidebar ----------------
  div(
    id = "sidebar",

    h4("Settings"),

    numericInput(
      "glasso_features",
      "Maximum genes for Glasso computation:",
      value = 200,
      min = 20,
      max = 5000,
      step = 10
    ),

    numericInput(
      "top_k",
      "Number of genes to draw in network:",
      value = 10,
      min = 1,
      max = 200,
      step = 1
    ),

    textInput(
      "target_gene",
      "Target gene:",
      value = "",
      placeholder = "Enter gene name"
    ),

    actionButton("apply_target", "Apply Target Gene"),

    fileInput(
      "data_file",
      "Upload raw matrix / Seurat:",
      accept = c(".csv", ".rds")
    ),

    helpText("Leave empty to use the default dataset."),

    h4("Conditional Prediction"),
    verbatimTextOutput("cond_mean"),
    verbatimTextOutput("cond_var")
  ),

  # ---------------- Main content ----------------
  div(
    id = "main_content",
    fluidRow(
      visNetworkOutput("network", height = "700px")
    )
  )
)



server <- function(input, output, session) {

  # =====================================================
  # 1. DATA & GLASSO FIT
  # =====================================================
  data_reactive <- reactive({

    req(input$glasso_features)

    if (is.null(input$data_file)) {

      result <- auto_fit_glasso(
        x = NULL,
        nfeatures = input$glasso_features
      )

    } else {

      ext <- tools::file_ext(input$data_file$name)

      if (ext == "csv") {

        df <- read.csv(input$data_file$datapath, 
                        row.names = 1, 
                        check.names = FALSE)

        mat <- as.matrix(df)

        # FORCE numeric storage
        mode(mat) <- "numeric"   # or storage.mode(mat) <- "double"
        result <- auto_fit_glasso(x = mat, nfeatures = input$glasso_features)

      } else if (ext == "rds") {

        obj <- readRDS(input$data_file$datapath)

        if (inherits(obj, "Seurat")) {
          mat <- preprocess_seurat_data(obj)
        } else if (is.matrix(obj)) {
          mat <- obj
        } else {
          stop("RDS must contain a matrix or Seurat object.")
        }

        result <- auto_fit_glasso(x = mat, nfeatures = input$glasso_features)

      } else {
        stop("Unsupported file type.")
      }
    }

    # ensure glasso_features <= available available genes
    nf <- input$glasso_features
    gene_count <- length(result$features)

    if (nf > gene_count) {
      updateNumericInput(session, "glasso_features", value = gene_count)
    }

    result
  })


  # =====================================================
  # 2. VALID TARGET GENE — with APPLY BUTTON
  # =====================================================
  target_gene_reactive <- reactiveVal(NULL)

  observeEvent(input$apply_target, {
    fit <- data_reactive()
    genes <- fit$features
    tg <- input$target_gene

    if (!(tg %in% genes)) {
      showNotification("Invalid gene name — using first gene.", type="warning")
      tg <- genes[1]
      updateTextInput(session, "target_gene", value = tg)
    }

    target_gene_reactive(tg)

    # immediately rebuild graph with new target
    k <- min(max(input$top_k, 1), length(fit$features))
    G_reactive(new_GLgraph(fit, k, target = tg))
  })


  # =====================================================
  # 3. GLgraph builder
  # =====================================================
  G_reactive <- reactiveVal(NULL)

  # new dataset loaded → build initial graph
  observeEvent(data_reactive(), {
    fit <- data_reactive()

    # default target is first gene
    genes <- fit$features
    target_gene_reactive(genes[1])
    updateTextInput(session, "target_gene", value = genes[1])

    k <- min(max(input$top_k, 1), length(fit$features))

    G_reactive(new_GLgraph(fit, k, target = genes[1]))
  })

  # only top_k changed → update graph (same target)
  observeEvent(input$top_k, {
    req(data_reactive())
    fit <- data_reactive()

    k <- min(max(input$top_k, 1), length(fit$features))
    tg <- target_gene_reactive()

    G_reactive(new_GLgraph(fit, k, target = tg))
  })


  # =====================================================
  # 4. Render network
  # =====================================================
  output$network <- renderVisNetwork({
    req(G_reactive())

    visnetwork_from_GLgraph(G_reactive()) %>%
      visNetwork::visLayout(randomSeed = 1) %>%
      visNetwork::visEvents(
        click = "function(params) {
          if (params.nodes.length > 0) {
            Shiny.setInputValue('node_click', params.nodes[0], {priority:'event'});
          }
        }"
      )
  })


  # =====================================================
  # 5. Toggle knockout state
  # =====================================================
  observeEvent(input$node_click, {
    req(G_reactive())
    G_reactive(visnetwork_toggle_knock(G_reactive(), input$node_click))
  })


  # =====================================================
  # 6. CONDITIONAL PREDICTION
  # =====================================================
  compute_conditional <- reactive({

    req(G_reactive())
    fit <- data_reactive()
    G   <- G_reactive()

    target <- target_gene_reactive()
    genes  <- fit$features

    knocked_nodes <- vapply(G$nodes, function(n)
      if (n$knocked) genes[n$index] else NA_character_, character(1))

    knocked_set <- knocked_nodes[!is.na(knocked_nodes)]

    out <- predict_knockout_from_fit(
      fit,
      target = target,
      knocked = knocked_set,
      use_original_units = TRUE
    )
    return(out)
  })

  output$cond_mean <- renderPrint({
    if (is.null(compute_conditional())) return("No result yet.")
    cat("Conditional Mean:\n")
    print(compute_conditional()$mean_cond)
  })

  output$cond_var <- renderPrint({
    if (is.null(compute_conditional())) return("No result yet.")
    cat("Conditional Variance:\n")
    print(compute_conditional()$cov_cond)
  })

}


shinyApp(ui, server)

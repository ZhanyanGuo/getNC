library(shiny)
library(visNetwork)
library(getNC)
library(plotly)

# ---------- helper: info icon with bootstrap popover ----------
label_with_info <- function(label, info) {
  tags$span(
    label,
    tags$span(
      icon("info-circle", class = "text-primary"),
      style = "margin-left:6px; cursor:pointer;",
      `data-toggle` = "popover",
      `data-trigger` = "hover",
      `data-container` = "body",
      title = label,
      `data-content` = info
    )
  )
}

# ================================================================
#                             UI
# ================================================================
ui <- fluidPage(

  # enable popovers
  tags$script("$(function(){$('[data-toggle=\"popover\"]').popover();});"),

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

  titlePanel("Click to Select Knockouts"),

  # ---------------- Sidebar ----------------
  div(
    id = "sidebar",

    h4("Settings"),

    numericInput(
      "glasso_features",
      label_with_info(
        "Maximum genes for Glasso computation",
        "Top variable genes chosen for sparse Gaussian graphical model estimation. Higher values give richer networks but cost more time."
      ),
      value = 200, min = 20, max = 5000, step = 10
    ),

    numericInput(
      "top_k",
      label_with_info(
        "Number of genes to draw",
        "The number of strongest absolute corrolated partner genes to display around the target in the network diagram."
      ),
      value = 10, min = 1, max = 200, step = 1
    ),

    textInput(
      "target_gene",
      label_with_info(
        "Target gene",
        "The gene whose conditional distribution and knockouts will be analyzed."
      ),
      value = "", placeholder = "Enter gene name"
    ),

    actionButton("apply_target", "Apply Target Gene"),

    fileInput(
      "data_file",
      label_with_info(
        "Upload raw matrix / Seurat",
        "Upload a CSV matrix (cells × genes) or a Seurat RDS object. Leave empty to use PBMC_small."
      ),
      accept = c(".csv", ".rds")
    ),

    helpText("Leave empty to use the default dataset."),

    h4("Knockout Prediction"),
    verbatimTextOutput("mean"),
    verbatimTextOutput("var")
  ),

  # ---------------- Floating Density Panel ----------------
  absolutePanel(
    id = "density_panel",
    top = 80, right = 20, width = 550, height = 500,
    draggable = TRUE,
    style = "background-color:white; border:1px solid #ccc; padding:10px; z-index:2000; overflow:auto;",
    h4("Density Plots (One More Knockout)"),
    plotlyOutput("density_mean", height = "220px"),
    plotlyOutput("density_var", height = "220px")
  ),

  # ---------------- Main content ----------------
  div(
    id = "main_content",
    fluidRow(
      visNetworkOutput("network", height = "700px")
    )
  ),

  tags$script(HTML("
  $(document).ready(function(){
    $('[data-toggle=\"popover\"]').popover(); 
  });
"))

)

# ================================================================
#                           SERVER
# ================================================================
server <- function(input, output, session) {

  # --------------------------------------------------------------
  # 1. DATA & GLASSO
  # --------------------------------------------------------------
  data_reactive <- reactive({
    req(input$glasso_features)

    if (is.null(input$data_file)) {
      # default dataset
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
        mode(mat) <- "numeric"

        result <- auto_fit_glasso(
          x = mat,
          nfeatures = input$glasso_features,
          min_genes = 0   # easier for test data
        )
        showNotification("CSV loaded successfully.", type="message")

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
      }
    }

    # ensure top K fits #genes
    nf <- input$glasso_features
    gene_count <- length(result$features)
    if (nf > gene_count) {
      updateNumericInput(session, "glasso_features", value = gene_count)
    }

    result
  })


  # --------------------------------------------------------------
  # 2. TARGET GENE LOGIC
  # --------------------------------------------------------------
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

    k <- min(max(input$top_k, 1), length(genes))
    G_reactive(new_GLgraph(fit, k, target = tg))
  })


  # --------------------------------------------------------------
  # 3. GLgraph
  # --------------------------------------------------------------
  G_reactive <- reactiveVal(NULL)

  observeEvent(data_reactive(), {
    fit <- data_reactive()
    genes <- fit$features
    target_gene_reactive(genes[1])
    updateTextInput(session, "target_gene", value = genes[1])

    k <- min(max(input$top_k, 1), length(genes))
    G_reactive(new_GLgraph(fit, k, target = genes[1]))
  })

  observeEvent(input$top_k, {
    req(data_reactive())
    fit <- data_reactive()
    tg <- target_gene_reactive()
    k <- min(max(input$top_k, 1), length(fit$features))

    G_reactive(new_GLgraph(fit, k, target = tg))
  })


  # --------------------------------------------------------------
  # 4. Render Network
  # --------------------------------------------------------------
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


  # --------------------------------------------------------------
  # 5. Toggle Knockout
  # --------------------------------------------------------------
  observeEvent(input$node_click, {
    req(G_reactive())
    G_reactive(visnetwork_toggle_knock(G_reactive(), input$node_click))
  })


  # --------------------------------------------------------------
  # 6. Conditional Prediction
  # --------------------------------------------------------------
  compute_conditional <- reactive({
    req(G_reactive())
    fit <- data_reactive()
    G   <- G_reactive()
    genes <- fit$features
    target <- target_gene_reactive()

    knocked_nodes <- vapply(G$nodes, function(n)
      if (n$knocked) genes[n$index] else NA_character_,
      character(1))

    if (all(is.na(knocked_nodes))) return(NULL)

    knocked_set <- knocked_nodes[!is.na(knocked_nodes)]

    predict_knockout_from_fit(
      fit,
      target = target,
      knocked = knocked_set,
      use_original_units = TRUE
    )
  })

  output$mean <- renderPrint({
    if (is.null(compute_conditional())) return("No result yet.")
    cat("Conditional Mean:\n")
    print(compute_conditional()$mean_cond)
  })

  output$var <- renderPrint({
    if (is.null(compute_conditional())) return("No result yet.")
    cat("Conditional Variance:\n")
    print(compute_conditional()$cov_cond)
  })


  # --------------------------------------------------------------
  # 7. Density Plots
  # --------------------------------------------------------------
  density_reactive <- reactive({
    req(G_reactive(), data_reactive())

    fit <- data_reactive()
    G   <- G_reactive()
    target <- target_gene_reactive()
    genes  <- fit$features

    knocked_nodes <- vapply(G$nodes, function(n)
      if (n$knocked) genes[n$index] else NA_character_,
      character(1))

    knocked_set <- knocked_nodes[!is.na(knocked_nodes)]

    plot_partner_knockout_densities_dual(
      fit,
      target = target,
      knocked = knocked_set,
      k = input$top_k,
      use_original_units = TRUE
    )
  })

  output$density_mean <- renderPlotly({
    d <- density_reactive()
    if (is.null(d)) return(NULL)
    d$mean_plot
  })

  output$density_var <- renderPlotly({
    d <- density_reactive()
    if (is.null(d)) return(NULL)
    d$var_plot
  })
}

# run app
shinyApp(ui, server)

# UI template generated by AI. With reactions and Server costumized.
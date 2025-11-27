library(shiny)
library(visNetwork)
library(getNC)

ui <- fluidPage(
  titlePanel("Interactive Gene Network (GLgraph)"),
  
  fluidRow(
    visNetworkOutput("network", height = "600px")
  )
)

server <- function(input, output, session) {
  
  fit <- auto_fit_glasso(x = NULL)
  G0  <- new_GLgraph(fit, k = 10, target = 1)

  # keep G reactive
  G_reactive <- reactiveVal(G0)
  
  # initial graph
  output$network <- renderVisNetwork({
    visnetwork_from_GLgraph(G_reactive()) %>%
      visNetwork::visEvents(click = "function(params) { 
      if (params.nodes.length > 0) { 
        Shiny.setInputValue('node_click', params.nodes[0], {priority: 'event'}); 
        } 
      }"
      )
  })

  observeEvent(input$node_click, {
    clicked <- input$node_click
    if (is.null(clicked)) return()

    G_reactive(visnetwork_toggle_knock(G_reactive(), clicked))

    output$network <- renderVisNetwork({
      visnetwork_from_GLgraph(G_reactive()) %>%
        visNetwork::visEvents(
          click = "function(params) { 
            if (params.nodes.length > 0) { 
              Shiny.setInputValue('node_click', params.nodes[0], {priority: 'event'}); 
            } 
          }"
        )
    })
  })
}


shinyApp(ui, server)
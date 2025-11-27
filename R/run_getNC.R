#' Launch the getNC Shiny Application
#' @export
run_getNC <- function() {
  app_dir <- system.file("shiny/app", package = "getNC")
  if (app_dir == "" || !dir.exists(app_dir)) {
    stop("Shiny app not found; please reinstall the getNC package.")
  }
  shiny::runApp(app_dir, display.mode = "normal")
}

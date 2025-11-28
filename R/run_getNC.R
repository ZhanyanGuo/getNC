#' Launch the getNC Interactive Shiny Application
#'
#' This function starts the \strong{getNC} Shiny application, an interactive
#' graphical tool for exploring Gaussian graphical models (glasso), selecting
#' gene knockouts directly from a network visualization, and viewing
#' conditional density plots for predicted perturbation effects.
#'
#' The app provides:
#' \itemize{
#'   \item Interactive \strong{GLgraph} visualization with selectable gene knockouts.
#'   \item Automatic processing of raw matrices or Seurat objects.
#'   \item Conditional mean/variance predictions under single or multiple knockouts.
#'   \item Interactive 3D density plots of predicted conditional distributions.
#' }
#'
#' By default, the Shiny UI loads using the PBMC example dataset bundled with
#' the package. Users can optionally upload their own data in CSV or RDS format.
#'
#' @details
#' This function locates the Shiny application directory installed with the
#' \pkg{getNC} package and runs it using \code{shiny::runApp()}. If the app
#' directory cannot be found (e.g., incomplete installation), an error is thrown.
#'
#' Typical usage:
#' \preformatted{
#'     library(getNC)
#'     run_getNC()
#' }
#'
#' The application runs in the user's default web browser unless otherwise
#' configured via Shiny options.
#'
#' @return
#' Launches the Shiny application and returns the Shiny object created by
#' \code{shiny::runApp()}. The function is normally called for its side effect
#' of opening the interactive UI.
#'
#' @seealso
#' \code{\link{auto_fit_glasso}}, \code{\link{new_GLgraph}},
#' \code{\link{predict_knockout_from_fit}},
#' \code{\link{plot_partner_knockout_densities_dual}}
#'
#' @export
run_getNC <- function() {
  app_dir <- system.file("shiny/app", package = "getNC")
  if (app_dir == "" || !dir.exists(app_dir)) {
    stop("Shiny app not found; please reinstall the getNC package.")
  }
  shiny::runApp(app_dir, display.mode = "normal")
}

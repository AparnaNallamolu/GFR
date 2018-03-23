#' BFR INTERFACE
#'
#' Interface
#'
#' This is a shiny app
#'
#' @export
runInterface <- function() {
  appDir <- system.file("shiny", "GFRInterface", package = "GFR")
  if (appDir == "") {
    Error("Could not find the app directory. Try re-installing `GFR`.")
  }

  shiny::runApp(appDir, display.mode = "normal")
}

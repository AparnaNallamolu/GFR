#' BFR INTERFACE
#'
#' Interface
#'
#' This is a shiny app
#'
#' @export
runInterface <- function() {
  appDir <- system.file("shiny-examples", "BFRInterface", package = "BFR")
  if (appDir == "") {
    Error("Could not find the app directory. Try re-installing `BFR`.")
  }

  shiny::runApp(appDir, display.mode = "normal")
}

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
    stop("Could not find the app directory. Try re-installing `GFR` package.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

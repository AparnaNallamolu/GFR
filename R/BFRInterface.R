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
    stop("Could not find example directory. Try re-installing `BFR`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

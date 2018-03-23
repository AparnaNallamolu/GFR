#' BFR INTERFACE
#'
#' Interface
#'
#' This is a shiny app
#'
#' @export
runBFRInterface <- function() {
  appDir <- system.file("shiny", "BFRInterface", package = "GFR")
  if (appDir == "") {
    Error("Could not find the app directory. Try re-installing `GFR`.")
  }

  shiny::runApp(appDir, display.mode = "normal")
}


#' FFR INTERFACE
#'
#' Interface
#'
#' This is a shiny app
#'
#' @export
runFFRInterface <- function() {
  appDir <- system.file("shiny", "FFRInterface", package = "GFR")
  if (appDir == "") {
    Error("Could not find the app directory. Try re-installing `GFR`.")
  }

  shiny::runApp(appDir, display.mode = "normal")
}

#' Launch App function
#'
#' Function launching the crispr app..;
#'
#'
#' @return None
#'
#' @examples
#' if (interactive()) crispr_app()
#' @export

crispr_app <- function() {
  library(shinyBS)
  shinyApp(ui = ui_crispr_app, server = server_crispr_app)
}

#' @export
crispr_app_profvis <- function() {
  library(shinyBS)
  profvis::profvis(runApp(shinyApp(ui = ui_crispr_app, server = server_crispr_app)))
  #profvis::profvis(shinyAppDir("/home/cbenoit/Documents/Gitlab/crispr-app/R"))
  #profvis::profvis(runApp("/home/cbenoit/Documents/Gitlab/crispr-app/R"))
}

# crispr_app <- function() {  appDir <- system.file("app", package = "CRISPRApp");shiny::runApp(appDir, display.mode = "normal")}
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
  if(interactive()) shinyApp(ui = ui_crispr_app, server = server_crispr_app)
}

# crispr_app <- function() {  appDir <- system.file("app", package = "CRISPRApp");shiny::runApp(appDir, display.mode = "normal")}
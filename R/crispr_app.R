#' Launch App function
#'
#' Function launching the crispr app..;
#'
#'
#' @return None
#'
#' @examples
#' crispr_app()
#' @export
crispr_app <- function() {  appDir <- system.file("app", package = "CRISPRApp");shiny::runApp(appDir, display.mode = "normal")}
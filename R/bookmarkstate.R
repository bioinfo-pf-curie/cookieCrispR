ui <- function(req) {
  fluidPage(
    textInput("txt", "Text"),
    checkboxInput("chk", "Checkbox"),
    bookmarkButton()
  )
}
server <- function(input, output, session) {
  # observe({
  #   # Trigger this observer every time an input changes
  #   reactiveValuesToList(input)
  #   session$doBookmark()
  # })
}

shinyApp(ui, server, enableBookmarking = "server")
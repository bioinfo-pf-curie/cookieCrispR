###save App state to binary data
saveState <- function(filename,reactives,separators,input) {
#saveState <- function(session,filename) {
  isolate({
    
    print("save inputs")
    #r_inputs <<- lapply(reactiveValuesToList(input), unclass)
    #r_inputs <- reactiveValuesToList(input)
    r_inputs <<- input
    print("save reactives")
    # #r_reactives <- reactiveValuesToList(reactives)
    r_reactives <<- reactives
    print("save separators")
    # #r_separators <- reactiveValuesToList(separators)
    r_separators <<- separators
    save(list = c("r_inputs", "r_reactives", "r_separators"),  file = filename)

  })
}

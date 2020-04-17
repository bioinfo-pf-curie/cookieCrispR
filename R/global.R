###save App state to binary data
saveState <- function(filename,reactives,separators,input,output) {
#saveState <- function(filename,reactives,separators,input) {
#saveState <- function(session,filename) {
  isolate({
    #session <- session
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
    #save(list = c("LiveInputs", "r_data" , "separators"),  file = filename)
    #environment(environment)
    print("save outputs")
    r_outputs <<- output
    #save.image(file = filename)
    save(list = c("r_inputs", "r_reactives", "r_separators","r_outputs"),  file = filename)
    #save(list = c("r_inputs"),  file = filename)
    
  })
}


# load previous state if available
# ip_inputs <- paste0("r_state")
# ip_data <- paste0("r_data")
# ip_module <- paste0("r_module")
# if (exists("r_state") && exists("r_data")) {
#   r_data <- do.call(reactiveValues, r_data)
#   lapply(names(r_state),
#          function(x) session$sendInputMessage(x, list(value = r_state[[x]]))
#   )
#   rm(r_data, r_state, envir = .GlobalEnv)
# } else {
#   r_state <- list()
#   r_data <- init_state(reactiveValues())
# }



# ssign(ip_inputs, lapply(reactiveValuesToList(input), unclass), envir = .GlobalEnv)
# assign(ip_data, reactiveValuesToList(r_data), envir = .GlobalEnv)
# assign(ip_module, r_module, envir = .GlobalEnv)



#######################################
# Load previous state
#######################################
# observe({
#   inFile <- input$uploadState
#   if(!is.null(inFile)) {
#     isolate({
#       tmpEnv <- new.env()
#       load(inFile$datapath, envir=tmpEnv)
#       if (exists("r_data", envir=tmpEnv, inherits=FALSE)){
#         assign(ip_data, tmpEnv$r_data, envir=.GlobalEnv)
#       }
#       if (exists("r_state", envir=tmpEnv, inherits=FALSE)) {
#         assign(ip_inputs, tmpEnv$r_state, envir=.GlobalEnv)
#         lapply(names(r_state),
#                function(x) session$sendInputMessage(x, list(value = r_state[[x]]))
#         )
#       }
#       if (exists("r_module", envir=tmpEnv, inherits=FALSE)){
#         assign(ip_module, tmpEnv$r_module, envir=.GlobalEnv)
#       }
#       #assign(ip_dump, lubridate::now(), envir = .GlobalEnv)
#       rm(tmpEnv)
#     })
#   }
# })


#enableBookmarking(store = "url")

#### example Datasets PATH ##########

metadata_path <- system.file("extdata", "SampleDescriptiondatatest.txt", package = "CRISPRApp")
counts_path <- system.file("extdata", "global_counts_table_datatest.csv", package = "CRISPRApp")
essential_path <- system.file("extdata", "essentials.csv", package = "CRISPRApp")
non_essential_path <- system.file("extdata", "non_essentials_datatest.csv", package = "CRISPRApp")




#######


# saveBook < function(state = state) {
#   state$values$sampleplan <<- reactives$sampleplan
#   # state$reactives$sampleplanGood <<- reactives$sampleplanGood
#   # state$reactives$sampleplanRaw <<- reactives$sampleplanRaw
#   # state$separators$counts <<- separators$counts
#   # state$separators$sampleplan <<- separators$sampleplan
# }

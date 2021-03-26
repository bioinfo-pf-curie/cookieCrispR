#' @title Perform DE analysis ui side
#'
#' @description
#'
#' @param id Module's id.
#' @param label Button's label.
#' @param icon Button's icon.
#' @param ... Arguments passed to \code{\link{actionButton}}
#'
#' @return a \code{\link[shiny]{reactiveValues}} containing the data selected under slot \code{data}
#' and the name of the selected \code{data.frame} under slot \code{name}.
#' @export
#'
#'
#' @importFrom htmltools tagList tags singleton
#' @importFrom shinydashboard infoBoxOutput valueBoxOutput renderValueBox valueBox infoBox
#' @importFrom shiny NS actionButton icon uiOutput
#' @importFrom shinyWidgets updatePickerInput pickerInput


CRISPRDeaModUI <- function(id)  {
  ns <- NS(id)
  fluidPage(
    tags$head(
      tags$style(type='text/css', ".span161 { width: 850px; }"),
    ),
    tags$head(tags$style(type = "#boxPopUp1 .modal-body{ min-height:550px}")),
    br(),
    fluidPage(
      fluidRow(
        #div(class = "span161",
        tabsetPanel(type = "pills",id = ns("DEGtabs"),
                    tabPanel("RUN DEA",
                             #fluidPage(
                            br(),
                            fluidRow(
                            column(width = 6,
                                   pickerInput(ns("celline"),width = '100%',
                                       label = "Select a cell line to perform analysis ",
                                       #choices =  sampleplanmodel$table$Cell_line,
                                       choices = NULL,
                                       options = pickerOptions(
                                         actionsBox = TRUE,
                                         title = "Select cell line",
                                         liveSearch = TRUE,
                                         liveSearchStyle = "contains",
                                       ),
                                       selected = NULL,
                                       multiple = FALSE,
                                       choicesOpt = NULL,
                                       inline = FALSE)),
                            column(width = 6,
                                   pickerInput(ns("comptype"),width = '100%',
                                         label = "Select a type of comparison",
                                         choices = c("Intra-Treatment","Inter-Treatment"),
                                         options = pickerOptions(
                                           actionsBox = TRUE,
                                           title = "Select comparison",
                                           liveSearch = TRUE,
                                           liveSearchStyle = "contains",
                                         ),
                                  selected = "Intra-Treatment",
                                  multiple = FALSE,
                                  choicesOpt = NULL,
                                  inline = FALSE)),
                                  column(width = 6,
                                  br(),
                                  actionButton(ns("Build"),"Build Model",width = '100%')),
                                  column(width = 6,
                                  conditionalPanel(condition = 'input.comptype == "Intra-Treatment"' ,ns = NS(id),
                                                   pickerInput(ns("Treatlevel"),width = '100%',
                                                   "Select a treatment to perform intra comparison",choices = NULL,
                                                    options = pickerOptions(
                                                     actionsBox = TRUE,
                                                     title = "Select treament value",
                                                     liveSearch = TRUE,
                                                     liveSearchStyle = "contains"))
                                                   ),
                                  conditionalPanel(condition = 'input.comptype == "Inter-Treatment"' ,ns = NS(id),
                                  br(),
                                  actionButton(ns("newcomparison"),"Set up a new comparison",width='100%'))),
                                  bsModal(id = ns("compmodal"),"Set up groups two compare",
                                          trigger = "comptype",
                                          column(width = 6,
                                              conditionalPanel(condition = 'output.conditionnalTreatment != "1"' ,ns = NS(id),
                                                 pickerInput(ns("TreatlevelInter1"),width = '100%',
                                                      "Group1 Treatment type",choices = NULL,
                                                      options = pickerOptions(
                                                        actionsBox = TRUE,
                                                        title = "Select Treatment",
                                                        liveSearch = TRUE,
                                                        liveSearchStyle = "contains"))),
                                              conditionalPanel(condition = 'output.conditionnalMut != "1"' ,ns = NS(id),
                                                 pickerInput(ns("Mut1"),
                                                             "Group1 sample mutation type",choices = NULL,
                                                             options = pickerOptions(
                                                               actionsBox = TRUE,
                                                               title = "Select mutation type",
                                                               liveSearch = TRUE,
                                                               liveSearchStyle = "contains")))),
                                          column(width = 6, 
                                                 conditionalPanel(condition = 'output.conditionnalTreatment != "1"' ,ns = NS(id),
                                                 pickerInput(ns("TreatlevelInter2"),
                                                             "Group2 Treatment type",choices = NULL,
                                                             options = pickerOptions(
                                                               actionsBox = TRUE,
                                                               title = "Select Treament",
                                                               liveSearch = TRUE,
                                                               liveSearchStyle = "contains"))),
                                                 conditionalPanel(condition = 'output.conditionnalMut != "1"' ,ns = NS(id),
                                                 pickerInput(ns("Mut2"),
                                                             "Group2 sample mutation type",choices = NULL,
                                                             options = pickerOptions(
                                                               actionsBox = TRUE,
                                                               title = "Select mutation type",
                                                               liveSearch = TRUE,
                                                               liveSearchStyle = "contains")))
                                                 ),
                                          column(width = 12,
                                          textOutput(ns("Intercompresume1"))),
                                          br(),
                                          column(width = 4,
                                          span(textOutput(ns("Intercompresume2")), style="color:red")),
                                          column(width = 4,
                                          textOutput(ns("Intercompresume3"))),
                                          column(width = 4,
                                          span(textOutput(ns("Intercompresume4")), style="color:red"))
                                          ),
                                  tags$head(tags$style("#compmodal .modal-footer{ display:none}"))
                                  ),
                             br(),
                             br()
                    ), # end of first tabs
                    tabPanel("Explore results",
                             br(),
                             tagList(
                               #fluidRow(column(width = 12,infoBoxOutput(ns("selected_line")))),
                               fluidRow(column(width = 12,uiOutput(ns("selected_line")))),
                               fluidRow(box(title = span(icon("cogs"), "Parameters"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                   status = "success",width= 12,
                                   #fluidRow(
                                     column(width = 12,pickerInput(ns("ExploreIntra"),label= "Select a comparison to explore", selected = NULL,
                                                                  multiple = FALSE,choicesOpt = NULL,inline = FALSE,choices = NULL,
                                                                  options = pickerOptions(
                                                                    title = "Select comparison",
                                                                    liveSearch = TRUE,
                                                                    liveSearchStyle = "contains"))),
                                     column(width = 6,numericInput(ns("FCT"),"LogFC threshold", min = 0, max = 10 , value = 1, step = 0.05)),
                                     column(width = 6,numericInput(ns("PvalsT"),"adjusted P values threshold", min = 0, max = 1 , value = 0.05, step = 0.01)))
                               ),
                               fluidRow(uiOutput(ns('features_value_box'))),
                               fluidRow(box(title = span(icon("chart-bar"),"DEA Figures"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                   status = "success",width = 12,
                                   #fluidRow(
                                     column(width =12,girafeOutput(ns("Pvals_distrib"))),
                                   #),
                                   br(),
                                   br(),
                                   #fluidRow(
                                     column(width = 12,
                                                   pickerInput(ns("GeneVolcano"),"Select genes to annotate on volcano",
                                                               selected = NULL,
                                                               multiple = TRUE,
                                                               choicesOpt = NULL,
                                                               inline = FALSE,
                                                               choices = NULL,
                                                               options = pickerOptions(
                                                                 title = "Select genes to annotate",
                                                                 liveSearch = TRUE,
                                                                 liveSearchStyle = "contains",
                                                                 actionsBox = TRUE
                                                               ))
                                   ),#),
                                   #fluidRow(
                                    column(width = 12,
                                                   plotOutput(ns("Volcano"))),
                                    column(width = 12,
                                                   conditionalPanel("input.GeneVolcano != undifined && input.GeneVolcano.length <= 1", ns = ns,
                                                                    textOutput(ns('boxplots_error'))),
                                                   conditionalPanel("input.GeneVolcano != undifined && input.GeneVolcano.length >=2", ns = ns,
                                                                    plotOutput(ns("boxplots"))))
                                   #)
                               )
                               ),# end of fluidRow
                               fluidRow(box(title = span(icon("arrow-circle-down"),"DEA Tables"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                   status = "success",width= 12,
                                   #fluidRow(
                                   tags$head(tags$style(".butt{background-color:#2E8B57;}")),
                                   # column(width =12,
                                   #        h4("All genes :"),
                                   #        DT::dataTableOutput(ns('results_table')),
                                   #        downloadButton(ns("resdl"),"All genes",class = "butt")),
                                   column(width = 6,
                                          br(),
                                          h4("Upp regulated genes :"),
                                          DT::dataTableOutput(ns('up_table')),
                                          downloadButton(ns("updl"),"Up-regulated")),
                                   column(width = 6,
                                          br(),
                                          h4("Down regulated genes :"),
                                          #br(),
                                          DT::dataTableOutput(ns('down_table')),
                                          downloadButton(ns("downdl"),"Down-regulated"))
                               )#,
                             ) # end of box
                             ) # end of Taglist
                    ), # end of second tab
                    tabPanel("RRAscores",
                             fluidRow(
                             column(width = 6,
                             br(),
                             pickerInput(ns("ExploreIntra2"),label= "Select a comparison to compute", selected = NULL,
                                         multiple = FALSE,choicesOpt = NULL,inline = FALSE,choices = NULL,
                                         options = pickerOptions(
                                           title = "Select comparison",
                                           liveSearch = TRUE,
                                           liveSearchStyle = "contains"),width = '100%')),
                             column(width = 6,
                                    numericInput(ns("n_perm"),"Please select a number of permutations for rra scores computing",value = 200, min = 100, max = 4000),
                                    "Increasing this value will also increase both ",
                                    "results precision and computing time",
                                    br()
                             )),
                             br(),
                             fluidRow(column(width = 12,
                                    actionButton(ns("computeRRA"),"Launch permutations",width = '100%'))),
                             br(),
                             #numericInput("n_perm","Select a number of permutation"),
                             fluidRow(column(width = 12,DT::dataTableOutput(ns("RRAscores")))),
                             br(),
                             fluidRow(column(width = 12,downloadButton(ns("scoresdl"),"Download RRA scores")))
                    ) # end of 3 tab
        )
        #) # end of div
      )
    ) # end of fluidRow
  ) # end of FLuidPage
}

#' @param input,output,session standards \code{shiny} server arguments.
#' @param header Does the file have a Header
#' @param sep What is the file separator
#'
#' @export
#'
#' @title Perform DE analysis server side
#'
#' @importFrom shiny showModal modalDialog observeEvent reactiveValues callModule observe icon
#' @importFrom htmltools tags HTML
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors sizeFactors
#' @importFrom edgeR DGEList
#' @import ggiraph
#' @import ggplot2
#' @importFrom shinydashboard renderInfoBox infoBox
#' @importFrom shiny renderUI modalDialog observeEvent reactiveValues callModule observe icon
#' @import limma
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr gather
#' @importFrom shinyWidgets updatePickerInput pickerInput pickerOptions
#' @import dplyr
#' @importFrom tictoc tic toc

#CRISPRDeaModServer <- function(input, output, session,sampleplan = NULL, var = NULL,
CRISPRDeaModServer <- function(input, output, session,sampleplan = NULL,
                               norm_data = NULL) {
  
  ### Define reactives #############
  req(sampleplan)
  req(norm_data)
  ns <- session$ns
  
  reactives <- reactiveValues(design = NULL, formula = NULL, contrast = NULL, selectedcomp = NULL,
                              selectedFC = NULL,selectedPvalT = NULL,GeneVolcano = NULL)
  sampleplanmodel <- reactiveValues(table = NULL)
  results <- reactiveValues(res = NULL, up = NULL, down = NULL,nsignfc = NULL,v = NULL,boxplots = NULL,
                            scores = NULL,old_res = "NULL")
  
  observeEvent(c(input$celline,sampleplan$table),priority = -1,{
    req(sampleplan$table)
    sampleplanmodel$table <- sampleplan$table %>%
      filter(Cell_line == input$celline)
  })
  
  observeEvent(c(input$comptype,input$newcomparison), {
    if(input$comptype == "Inter-Treatment"){
    toggleModal(session=session,"compmodal", toggle = "open")
    }
  })
  
  observeEvent(c(sampleplanmodel$table,input$comptype),{
    updatePickerInput(session=session,"Treatlevel",
                      choices = unique(as.character(sampleplanmodel$table[,"Treatment"])),
                      selected = unique(as.character(sampleplanmodel$table[,"Treatment"]))[1])
  })
  observeEvent(c(sampleplan$table),{
    updatePickerInput(session=session,"celline",
                      choices = unique(sampleplan$table$Cell_line),selected = unique(sampleplan$table$Cell_line)[1])
  })
  observeEvent(c(sampleplanmodel$table,input$comptype),{
    if(length(unique(sampleplanmodel$table$Treatment)) >= 2){
    updatePickerInput(session=session,"TreatlevelInter1",
                      choices = unique(as.character(sampleplanmodel$table[,"Treatment"])),
                      selected = unique(as.character(sampleplanmodel$table[,"Treatment"]))[1])
    updatePickerInput(session=session,"TreatlevelInter2",
                      choices = unique(as.character(sampleplanmodel$table[,"Treatment"])),
                      selected = unique(as.character(sampleplanmodel$table[,"Treatment"]))[2])
    }
  })
  observeEvent(c(sampleplanmodel$table,input$comptype),{
    if(length(unique(sampleplanmodel$table$SupplementaryInfo)) >= 2){
    updatePickerInput(session=session,"Mut1",
                      choices = unique(as.character(sampleplanmodel$table[,"SupplementaryInfo"])),
                      selected = unique(as.character(sampleplanmodel$table[,"SupplementaryInfo"]))[1])
    updatePickerInput(session=session,"Mut2",
                      choices = unique(as.character(sampleplanmodel$table[,"SupplementaryInfo"])),
                      selected = unique(as.character(sampleplanmodel$table[,"SupplementaryInfo"]))[2])
    }
  })
  #observeEvent(c(input$TreatlevelInter1,input$TreatlevelInter2,input$Mut1,input$Mut2),{
    output$Intercompresume1 <- renderText({
      paste0("You are about to compare these samples : ")
    })
    output$Intercompresume2 <- renderText({
      if(length(unique(sampleplanmodel$table$Treatment)) >= 2 & length(unique(sampleplanmodel$table$SupplementaryInfo)) >= 2){
      paste0(input$TreatlevelInter1,"_",input$Mut1)
      } else if (length(unique(sampleplanmodel$table$Treatment)) == 1  & length(unique(sampleplanmodel$table$SupplementaryInfo)) >= 2){
        paste0(input$Mut1)
      } else if (length(unique(sampleplanmodel$table$Treatment)) >=2  & length(unique(sampleplanmodel$table$SupplementaryInfo)) == 1){
        paste0(input$TreatlevelInter1)
      }
    })
    output$Intercompresume3 <- renderText({
      "vs"
    })
    output$Intercompresume4 <- renderText({
      if(length(unique(sampleplanmodel$table$Treatment)) >= 2 & length(unique(sampleplanmodel$table$SupplementaryInfo)) >= 2){
      paste0(input$TreatlevelInter2,"_",input$Mut2)
      } else if (length(unique(sampleplanmodel$table$Treatment)) == 1  & length(unique(sampleplanmodel$table$SupplementaryInfo)) >= 2){
      paste0(input$Mut2)
      } else if (length(unique(sampleplanmodel$table$Treatment)) >=2  & length(unique(sampleplanmodel$table$SupplementaryInfo)) == 1){
      paste0(input$TreatlevelInter2)
      }
    })
    output$Intercompresume5 <- renderText({
      "samples"
    })
    
    output$conditionnalTreatment <- reactive({
      as.character(length(unique(sampleplanmodel$table$Treatment)))
    })
    output$conditionnalMut <- reactive({
      as.character(length(unique(sampleplanmodel$table$SupplementaryInfo)))
    })
    outputOptions(output, "conditionnalTreatment", suspendWhenHidden = FALSE) 
    outputOptions(output, "conditionnalMut", suspendWhenHidden = FALSE)  
  

  observeEvent(input$Build,{
    if(!is.null(sampleplanmodel$table)){
            if(!is.null(norm_data$data)){
              req(input$celline)
                data <- norm_data$data$samples %>%
                  filter(Cell_line == input$celline)
                Tnums <- as.numeric(unique(gsub("T","",data$Timepoint)))
                T0 <- min(Tnums[order(Tnums)])
                minT <- min(Tnums[order(Tnums)][-1])
                maxT <- max(Tnums)
                if(input$comptype == "Intra-Treatment"){
                  data$Intra <- paste0(data$Treatment,"_",data$Timepoint)
                  design <- model.matrix(~ 0 + Intra, data = data)
                  colnames(design) <- make.names(gsub("Intra","",colnames(design)))
                } else if(input$comptype == "Inter-Treatment"){
                  data$Inter <- paste0(data$Treatment,"_",data$SupplementaryInfo,"_",data$Timepoint)
                  design <- model.matrix(~ 0 + Inter, data = data)
                  colnames(design) <- make.names(gsub("Inter","",colnames(design)))
                }
                if(input$comptype == "Intra-Treatment"){
                  contrast <- purrr::map(glue::glue("{input$Treatlevel}_T{minT:maxT}-{input$Treatlevel}_T{T0}"),~makeContrasts(contrasts = .x,levels=design))
                  reactives$contrast <- contrast
                  reactives$design <- design
                } else if (input$comptype == "Inter-Treatment"){
                  if (length(unique(sampleplanmodel$table$SupplementaryInfo)) > 1 & length(unique(sampleplanmodel$table$Treatment)) > 1){
                    gluelist <- data.frame(gluelistcol = unique(glue::glue("{input$TreatlevelInter1}_{input$Mut1}_{data$Timepoint}-{input$TreatlevelInter2}_{input$Mut2}_{data$Timepoint}")))
                    missing1 <- unique(glue::glue("{input$TreatlevelInter1}_{input$Mut1}_{data$Timepoint}"))
                    missing2 <- unique(glue::glue("{input$TreatlevelInter2}_{input$Mut2}_{data$Timepoint}"))
                    missing2 <- missing2[!(missing2 %in% data$Inter)]
                    missing1 <- missing1[!(missing1 %in% data$Inter)]
                    for (missing in missing2){
                      gluelist <- gluelist %>% filter(!grepl(missing,gluelistcol))
                    }
                    for (missing in missing1){
                      gluelist <- gluelist %>% filter(!grepl(missing,gluelistcol))
                    }
                    gluelist <- as.list(gluelist$gluelistcol)
                  } else if (length(unique(sampleplanmodel$table$SupplementaryInfo)) > 1 & length(unique(sampleplanmodel$table$Treatment)) == 1){
                  gluelist <- data.frame(gluelistcol = unique(glue::glue("{unique(sampleplanmodel$table$Treatment)}_{input$Mut1}_{data$Timepoint}-{unique(sampleplanmodel$table$Treatment)}_{input$Mut2}_{data$Timepoint}")))
                  missing1 <- unique(glue::glue("{unique(sampleplanmodel$table$Treatment)}_{input$Mut1}_{data$Timepoint}"))
                  missing2 <- unique(glue::glue("{unique(sampleplanmodel$table$Treatment)}_{input$Mut2}_{data$Timepoint}"))
                  missing2 <- missing2[!(missing2 %in% data$Inter)]
                  missing1 <- missing1[!(missing1 %in% data$Inter)]
                  for (missing in missing2){
                    gluelist <- gluelist %>% filter(!grepl(missing,gluelistcol))
                  }
                  for (missing in missing1){
                    gluelist <- gluelist %>% filter(!grepl(missing,gluelistcol))
                  } 
                  gluelist <- as.list(gluelist$gluelistcol)

                  } else if (length(unique(sampleplanmodel$table$SupplementaryInfo)) == 1 & length(unique(sampleplanmodel$table$Treatment)) > 1){
                  gluelist <- data.frame(gluelistcol = unique(glue::glue("{input$TreatlevelInter1}_{unique(sampleplanmodel$table$SupplementaryInfo)}_{data$Timepoint}-{input$TreatlevelInter2}_{unique(sampleplanmodel$table$SupplementaryInfo)}_{data$Timepoint}")))
                  missing1 <- unique(glue::glue("{input$TreatlevelInter1}_{unique(sampleplanmodel$table$SupplementaryInfo)}_{data$Timepoint}"))
                  missing2 <- unique(glue::glue("{input$TreatlevelInter2}_{unique(sampleplanmodel$table$SupplementaryInfo)}_{data$Timepoint}"))
                  missing2 <- missing2[!(missing2 %in% data$Inter)]
                  missing1 <- missing1[!(missing1 %in% data$Inter)]
                  for (missing in missing2){
                    gluelist <- gluelist %>% filter(!grepl(missing,gluelistcol))
                  }
                  for (missing in missing1){
                    gluelist <- gluelist %>% filter(!grepl(missing,gluelistcol))
                  } 
                  gluelist <- as.list(gluelist$gluelistcol)
                  }
                  if (length(gluelist) != 0){
                  contrast <- purrr::map(gluelist,~makeContrasts(contrasts = .x,levels=design))
                  reactives$contrast <- contrast
                  reactives$design <- design
                  } else {
                    showModal(modalDialog(HTML(
                      paste0("<b> These samples would be suitable : </b></br>",paste0(c(missing1,missing2),collapse = "</br>"))),
                      title = "Missing required samples for this comparisons",
                      footer = tagList(
                        modalButton("Got it"))
                    ))
                  }
                }
              }
    }
  })
  
  ################ Compute Model new block ##########################
  observeEvent(c(norm_data$data$counts,
                 reactives$contrast),priority = 10,{
                   complist <- lapply(1:length(reactives$contrast), function(x){colnames(reactives$contrast[[x]])}) %>% unlist()
                   if(!(TRUE %in% c(complist %in% names(concatenated$results))) | is.null(concatenated$results) ){
                   if (!is.null(norm_data$data$counts) && !is.null(reactives$design)){
                     withProgress(message = 'Computing differential analysis', value = 0.5,{
                     incProgress(0.3,detail = "Filtering low expressed guides")
                     counts <- norm_data$data$counts[,colnames(norm_data$data$counts)%in%rownames(reactives$design)]
                     
                     # Remove control guides RNA for differential analysis
                     #counts <- filter(as.data.frame(counts), str_detect(Gene,"Non-Targeting"))
                     
                     ## Supprime car difficile à calibrer ############ 
                    ### Filtre sur nombre minimal de counts
                     kept <- which(rowSums(edgeR::cpm(counts) >= 1) >= 2)
                     counts <- counts[kept,]
                     
                     #incProgress(0.3,detail = "voomWithQualityWeights")
                     #results$v <- limma::voomWithQualityWeights(counts, design = reactives$design, normalize.method = "none", span = 0.5, plot = FALSE)
                     incProgress(0.3,detail = "voom")
                     results$v <- limma::voom(counts, design = reactives$design, normalize.method = "none", span = 0.5, plot = FALSE)
                     
                     res_fit <- limma::lmFit(results$v, method = "ls")
                     incProgress(0.3,detail = "fitting model")
                     res_eb <- eBayes(res_fit, robust = FALSE)
                     fit <- purrr::map(reactives$contrast, ~contrasts.fit(res_fit, contrasts = .x))
                     res_eb <- purrr::map(fit, ~eBayes(.x, robust = FALSE))
                     incProgress(0.3,detail = "Formating results")
                     tab <- purrr::map(res_eb,~process_res(.x,sgRNA_annot = norm_data$data$genes))
                     names(tab) <- lapply(tab, function(x){print(paste0(input$celline," || ",unique(as.character(x$term))))})
                     results$res <- tab
                     setProgress(1)
                     })
                   }
                   } else { 
                     showModal(modalDialog(HTML(
                     "<b>Your results are available on the other outlets. </b></br>"),
                     title = "This comparison was already computed !",
                     footer = tagList(
                       modalButton("Got it"))))
                    }
                 }) # end of observer
  ###########################################################################
 
  concatenated <- reactiveValues(resultsIntraNames = NULL, resultsInterNames = NULL,results = NULL)
  observeEvent(results$res,{
    req(results$res)
    if(input$comptype == "Intra-Treatment"){
    concatenated$resultsIntraNames <- c(concatenated$resultsIntraNames,names(results$res))
    } else if(input$comptype == "Inter-Treatment"){
    concatenated$resultsInterNames <- c(concatenated$resultsInterNames,names(results$res))
    }
    concatenated$results <- c(concatenated$results,results$res)
  })
  
  createLink <- function(val) {
    #sprintf('<a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s" target="_blank" class="btn btn-primary">Info</a>',val)
    sprintf('<a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s" target="_blank" class="btn btn-primary">Info</a>',gsub("_[0-9]","",gsub("sg","",val)))
    
  }
  
  observeEvent(c(concatenated$results),priority = 10,{
    if(length(concatenated$results) > 1){
    updatePickerInput(session = session,"ExploreIntra",
                      selected = c(concatenated$resultsInterNames,concatenated$resultsIntraNames)[1],
                      choices = list(
                        Inter_comparisons = concatenated$resultsInterNames,
                        Intra_comparisons = concatenated$resultsIntraNames
                      )
                      )
    updatePickerInput(session = session,"ExploreIntra2",
                      selected = c(concatenated$resultsInterNames,concatenated$resultsIntraNames)[1],
                      choices = list(
                             Inter_comparisons = concatenated$resultsInterNames,
                             Intra_comparisons = concatenated$resultsIntraNames
                             )
                      )
    } else {
      updatePickerInput(session = session,"ExploreIntra",
                        selected = c(concatenated$resultsInterNames,concatenated$resultsIntraNames)[1],
                        choices = names(concatenated$results)
      )
      updatePickerInput(session = session,"ExploreIntra2",
                        selected = c(concatenated$resultsInterNames,concatenated$resultsIntraNames)[1],
                        choices = names(concatenated$results)
      )
    }
  })
  ################ Compute Model old block ##########################
  observeEvent(c(concatenated$results,
                 input$FCT,
                 input$PvalsT,input$DEGtabs,input$ExploreIntra),ignoreInit = TRUE,{
                   if(input$DEGtabs == "Explore results"){
                   if(!is.null(concatenated$results)){
                   req(input$ExploreIntra)
                   res <- concatenated$results[[input$ExploreIntra]]
                   nsign <- length(which(res$adj_p.value < input$PvalsT))
                   results$nsignfc <- length(which(res$adj_p.value < input$PvalsT & abs(res$estimate) > input$FCT))
                   up <- which(res$adj_p.value_enrich < input$PvalsT)
                   down <- which(res$adj_p.value_dep < input$PvalsT)
                   #print(head(res$Gene))
                   res$ENSEMBL <- createLink(res$Gene)
                   print('end of DEG')
                   results$up <- res[up,]
                   results$down <- res[down,]
                   results$restable <- res
                   } else {
                   showModal(modalDialog(
                       HTML(
                         "<b> Please perform differential analysis first</b></br>
         Select variables in the RUN DEA menu then click on the build button.
                       "),
                       title = "Missing previous step !",
                       footer = tagList(
                         modalButton("Got it"))
                     ))
                   }
                   }
                 }) # end of observer
  
  #n_perm <- 2000
  # alpha_thr <- 0.3
  # observeEvent(c(input$DEGtabs,norm_data$data$genes),{
  #   req(norm_data$data$genes)
  #   if(input$DEGtabs == "RRAscores") {
  #     if(!is.null(concatenated$results)){
  #       if((length(concatenated$results) != length(results$old_res)) | length(concatenated$results) == 1){
  #         showModal(modalDialog(
  #           title = "Please select a number of permutations for rra scores computing",
  #           numericInput(ns("n_perm"),"Select a number of permutation",value = 200, min = 100, max = 4000),
  #           br(),
  #           "Increasing this value will also increase both results precision and computing time",
  #           footer = tagList(
  #             actionButton(ns("gotit"),"Launch permutations")
  #           )))
  #       observeEvent(input$gotit,{
  #       req(input$n_perm)
  #       removeModal()
  #       results$old_res <- concatenated$results
  #       n <- length(concatenated$results)
  #       names <- names(concatenated$results)
  #       sgRNAannot <- norm_data$data$genes
  #       print(sgRNAannot$Gene)
  # 
  #       sgRNAannot <- filter(sgRNAannot, str_detect(Gene,paste(c("Non-Targeting","negative_control"),collapse = "|")))
  #       
  #       withProgress(message = paste0('Computing per gene RRA score from DEA results (nperms = ',input$n_perm,')'), value = 0.5,{
  #         results$scores <- NULL
  #         
  #         #### Compute en simultane toutes les comps
  #         scores <- lapply(seq_along(concatenated$results),function(x){
  #           incProgress(as.numeric(1/n),detail = paste0("processing ",names[x]))
  #           computed_scores <- compute_score_RRA(concatenated$results[[x]],alpha_thr = alpha_thr)
  #           print(paste0("Launching RRA adj pval with ",input$n_perm," permutations..."))
  #           compute_RRA_pval(guide_res = concatenated$results[[x]],
  #                            gene_res = computed_scores,
  #                            n_guides = 6,
  #                            non_target =  sgRNAannot,
  #                            alpha_thr = alpha_thr,
  #                            n_perm = as.integer(input$n_perm))
  #         })
  # 
  #         ### Compute seulement la comp active
  #         names(scores) <- names
  #         print("resutsscores")
  #         results$scores <- scores
  #         print(head(results$scores))
  #         setProgress(1)
  #       })
  #       })# end of progress
  #     } else {
  #       showModal(modalDialog(HTML(
  #        "<b> Please perform differential analysis first or rerun with </b></br>
  #        To do so select at least two variables for the comparison then click on the build button.
  #       "),
  #         title = "Missing previous step !",
  #         footer = tagList(
  #           modalButton("Got it"))
  #       ))
  #     #})
  #       } # end of if
  #     }
  #   }
  # }) # end of observer
  
  #n_perm <- 2000
  alpha_thr <- 0.3
  #observeEvent(c(input$DEGtabs,norm_data$data$genes),{
  observeEvent(c(input$computeRRA),{
    req(input$computeRRA)
    req(norm_data$data$genes)
    #req(input$ExploreIntra2)
    req(input$n_perm)
    if(!is.null(concatenated$results)){
    if(!(input$ExploreIntra2 %in% names(results$scores))){
            results$old_res <- concatenated$results
            n <- length(concatenated$results)
            names <- names(concatenated$results)
            sgRNAannot <- norm_data$data$genes

            sgRNAannot <- filter(sgRNAannot, str_detect(Gene,paste(c("Non-Targeting","negative_control"),collapse = "|")))

            withProgress(message = paste0('Computing per gene RRA score from DEA results (nperms = ',input$n_perm,')'), value = 0.5,{
              results$scores <- NULL
              incProgress(as.numeric(1/n),detail = paste0("processing ",input$ExploreIntra2))
              computed_scores <- compute_score_RRA(concatenated$results[[input$ExploreIntra2]],alpha_thr = alpha_thr)
              print(paste0("Launching RRA adj pval with ",input$n_perm," permutations..."))
              scores <- compute_RRA_pval(guide_res = concatenated$results[[input$ExploreIntra2]],
                                         gene_res = computed_scores,
                                         n_guides = 6,
                                         non_target =  sgRNAannot, 
                                         alpha_thr = alpha_thr,
                                         n_perm = as.integer(input$n_perm))

              results$scores[[input$ExploreIntra2]] <- scores
              print(head(results$scores))
              setProgress(1)
          })# end of progress
    } else {
      showModal(modalDialog(
        "Do you want to remove the previous computation ?",
        "",
        "You'll have to re-click on the Launch Permutations button with your desired number of perm",
        title = "This scores were already computed",
        footer = tagList(
          actionButton(ns("Yes"),label = "Yes"),
          actionButton(ns('No'),label = "No"))
      ))
      #removeModal()
      }
      } else {
          showModal(modalDialog(HTML(
            "<b> Please perform differential analysis first</b></br>
         Select variables in the RUN DEA menu then click on the build button.
        "),
            title = "Missing previous step !",
            footer = tagList(
              modalButton("Got it"))
          ))
        } # end of if
  }) # end of observer
  
  observeEvent(input$Yes,{
    req(input$Yes)
    results$scores[[input$ExploreIntra2]] <- NULL
    removeModal()
  })
  observeEvent(input$No,{
    req(input$No)
    removeModal()
  })
  
  output$RRAscores <- DT::renderDataTable({
    #req(input$ExploreIntra2)
    #req(results$scores)
    req(results$scores[[input$ExploreIntra2]])
    scores <- results$scores[[input$ExploreIntra2]] %>% select(c("Gene","RRA_adjp","RRA_enrich_adjp","RRA_dep_adjp"))
    
    datatable(scores,rownames=FALSE,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE,initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
        "}")))
    
  })
  output$scoresdl <- downloadHandler(
    filename = function() {
      req(results$scores)
      req(input$ExploreIntra2)
      paste(input$ExploreIntra2,"_RRAscores", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      scores <- results$scores[[input$ExploreIntra2]] %>% select(c("Gene","RRA_adjp","RRA_enrich_adjp","RRA_dep_adjp"))
      #scores <- scores[order(scores$RRA_adjp),]
      write.csv(scores,file)
    }
  )

  output$Pvals_distrib <- renderGirafe({
    req(concatenated$results)
    req(input$ExploreIntra)
    res <- concatenated$results[[input$ExploreIntra]]
    #plot <- ggplot(data = res) + aes(x = `adj_p.value`) +
    plot <- ggplot(data = res) + aes(x = `p.value`) +
      geom_histogram_interactive(fill = "steelblue",breaks = seq(0, 1, length.out = 30))
    build  <- ggplot_build(plot)
    plot <- plot +  labs(title = paste0(input$ExploreIntra," : P values distribution"), x = "P values", y = "Occurences")# +
    ggiraph::girafe(code = {print(plot)})
  })
  # 
  observe({
    updatePickerInput("GeneVolcano", session = session, choices = rownames(norm_data$data$counts))
  })
  # 
  Volcano <- reactiveValues(plot = NULL)
  #observeEvent(concatenated$results,{
  observe({
    req(concatenated$results)
    req(input$ExploreIntra)
    res <- concatenated$results[[input$ExploreIntra]]
    tic("Ploting Volcano")
    ggplot <- ggplot(res, aes(x = estimate, y = -log10(p.value))) +
      ggtitle(colnames(reactives$contrast)) +
      scale_fill_gradient(low = "lightgray", high = "navy") +
      scale_color_gradient(low = "lightgray", high = "navy") +
      expand_limits(y = c(min(-log10(res$p.value)), 1)) +
      geom_point(data = res,
                 color = "grey", alpha = 0.5) +
      geom_point(data = subset(res, estimate > input$FCT),
                 color = "red", alpha = 0.5) +
      geom_point(data = subset(res, estimate < -input$FCT),
                 color = "blue", alpha = 0.5) +
      geom_point(data = subset(res, adj_p.value < input$PvalsT),
                 color = "green", alpha = 0.5) +
      geom_hline(yintercept = -log10(max(subset(res, adj_p.value < input$PvalsT)$p.value)), linetype = "dashed") +
      geom_vline(xintercept = c(-input$FCT, input$FCT), linetype = "dashed") +
      theme_linedraw() +
      theme(panel.grid = element_blank()) +
      xlab("Fold change (log2)") +
      ylab("-log10(P-Value)")
    Volcano$plot <- ggplot
    toc(log = TRUE)
  })
  # 
observeEvent(Volcano$plot,{
  output$Volcano <- renderPlot({
    req(Volcano$plot)
    req(input$ExploreIntra)
    req(concatenated$results)
    res <- concatenated$results[[input$ExploreIntra]]
    tic("Rendering Volcano...")
    ggplot <- Volcano$plot +
      geom_point(data = subset(res,sgRNA %in% input$GeneVolcano),
                 color = "purple", alpha = 0.6) +
      ggrepel::geom_text_repel(
        data = subset(res,sgRNA %in% input$GeneVolcano),
        aes(label = sgRNA),
        size = 5,
        force = 2,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )
    return(ggplot)
    toc(log = TRUE)
  })
})
  # 
  observeEvent(input$GeneVolcano,{
    if(length(input$GeneVolcano) >1){
      req(norm_data$data)
      req(sampleplanmodel$table)
      req(input$ExploreIntra)
      withProgress(message = 'Drawing boxplots with selected genes', value = 0.5,{
        incProgress(0.3)
        groups_table <- sampleplanmodel$table
        groups_table$Samples <- rownames(groups_table)
        if(input$comptype == "Intra-Treatment"){
          group1 <- stringr::str_split_fixed(gsub(paste0(input$Treatlevel,"_"),"",input$ExploreIntra),"-",n=2)[,1]
          group2 <- stringr::str_split_fixed(gsub(paste0(input$Treatlevel,"_"),"",input$ExploreIntra),"-",n=2)[,2]
          groups_table <- groups_table[,c("Treatment","Timepoint","Samples")] %>%
            filter(Treatment == input$Treatlevel)
         } else if(input$comptype == "Inter-Treatment"){
           if (length(unique(sampleplanmodel$table$SupplementaryInfo)) > 1 & length(unique(sampleplanmodel$table$Treatment)) > 1){
             group1 <- gsub(paste0(input$Mut1,"_"),"",gsub(paste0(input$TreatlevelInter1,"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,1]))
             group2 <- gsub(paste0(input$Mut2,"_"),"",gsub(paste0(input$TreatlevelInter2,"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,2]))
             groups_table <- groups_table[,c("Treatment","Timepoint","Samples","SupplementaryInfo")] %>%
               filter(Treatment == c(input$TreatlevelInter1,input$TreatlevelInter2)) %>%
               filter(SupplementaryInfo == c(input$Mut1,input$Mut2))  %>%
               mutate(COMP = paste0(Treatment,"_",SupplementaryInfo))
           } else if (length(unique(sampleplanmodel$table$SupplementaryInfo)) > 1 & length(unique(sampleplanmodel$table$Treatment)) == 1){
             group1 <- gsub(paste0(input$Mut1,"_"),"",gsub(paste0(unique(sampleplanmodel$table$Treatment),"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,1]))
             group2 <- gsub(paste0(input$Mut2,"_"),"",gsub(paste0(unique(sampleplanmodel$table$Treatment),"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,2]))
             groups_table <- groups_table[,c("Treatment","Timepoint","Samples","SupplementaryInfo")] %>%
               filter(Treatment %in% unique(sampleplanmodel$table$Treatment)) %>%
               filter(SupplementaryInfo == c(input$Mut1,input$Mut2))  %>%
               mutate(COMP = paste0(Treatment,"_",SupplementaryInfo))
           } else if (length(unique(sampleplanmodel$table$SupplementaryInfo)) == 1 & length(unique(sampleplanmodel$table$Treatment)) > 1){
             group1 <- gsub(paste0(unique(sampleplanmodel$table$SupplementaryInfo),"_"),"",gsub(paste0(input$TreatlevelInter1,"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,1]))
             group2 <- gsub(paste0(unique(sampleplanmodel$table$SupplementaryInfo),"_"),"",gsub(paste0(input$TreatlevelInter2,"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,2]))
             groups_table <- groups_table[,c("Treatment","Timepoint","Samples","SupplementaryInfo")] %>%
               filter(Treatment %in% c(input$TreatlevelInter1,input$TreatlevelInter2)) %>%
               filter(SupplementaryInfo == unique(sampleplanmodel$table$SupplementaryInfo))  %>%
               mutate(COMP = paste0(Treatment,"_",SupplementaryInfo))
          }
         }

        boxplotdata <- results$v$E[which(rownames(results$v$E) %in% input$GeneVolcano),]
        boxplotdata <- rbind(boxplotdata,colnames(boxplotdata))
        rownames(boxplotdata)[nrow(boxplotdata)] <- "Samples"
        incProgress(0.3)
        boxplotdata <- as.data.frame(t(boxplotdata)) %>%  gather(key = "GENE",value = "COUNTS", -Samples)
        boxplotdata$Samples <- as.character(boxplotdata$Samples)
        boxplotdata <- inner_join(boxplotdata,groups_table, by = "Samples")
        boxplotdata$COUNTS <- as.numeric(boxplotdata$COUNTS)
        if(input$comptype == "Intra-Treatment"){
            boxplotdata[,c(group1,group2)] <- as.character(boxplotdata[,"Timepoint"])
            results$boxplots <- ggplot(boxplotdata, aes(x=Timepoint, y=COUNTS,fill = Timepoint)) +
              geom_boxplot() +
              facet_grid(. ~ GENE) +
              geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.2,
                                                       seed = 1234),
                         pch=21,
                         # size = 2,
                         aes_string(fill="Timepoint"), show.legend = T)
            setProgress(1)
        } else if(input$comptype == "Inter-Treatment"){
          boxplotdata[,c(group1,group2)] <- as.character(boxplotdata[,"Timepoint"])
          results$boxplots <- ggplot(boxplotdata, aes(x=COMP, y=COUNTS,fill = COMP)) +
            geom_boxplot() +
            facet_grid(Timepoint ~ GENE) +
            geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.2,
                                                     seed = 1234),
                       pch=21,
                       # size = 2,
                       aes_string(fill="COMP"), show.legend = T)
          setProgress(1)
        }
      }) # end of progress
    }
  })
   
  output$boxplots <- renderPlot(results$boxplots)
  output$boxplots_error <- renderText({
    validate(
      need(length(input$GeneVolcano) >= 2, "Select at least two genes to draw boxplots...")
    )
  })
   

  output$up_table <- DT::renderDataTable({
    ups <- results$up %>%
      column_to_rownames("sgRNA") %>%
    #   select(c("estimate","p.value","adj_p.value_enrich","ENSEMBL"))
    # colnames(ups) <- c("logFC","p.value","adj_p.value_enrich","ENSEMBL")
    select(c("estimate","adj_p.value_enrich","ENSEMBL"))
    colnames(ups) <- c("logFC","adj_p.value_enrich","ENSEMBL")
    ups <- ups[order(ups$adj_p.value_enrich),]
    datatable(
      ups,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE,initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
        "}")))
  })

  output$updl <- downloadHandler(
    filename = function() {
      paste("DEA-UPS-PDX-Results", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      ups <- results$up %>%
        column_to_rownames("sgRNA") %>%
      #   select(c("estimate","p.value","adj_p.value_enrich"))
      # colnames(ups) <- c("logFC","p.value","adj_p.value_enrich")
      select(c("estimate","adj_p.value_enrich"))
      colnames(ups) <- c("logFC","adj_p.value_enrich")
      ups <- ups[order(ups$adj_p.value_enrich),]
      write.csv(ups, file)
    }
  )

  output$down_table <- DT::renderDataTable({
    down <- results$down %>%
      column_to_rownames("sgRNA") %>%
      select(c("estimate","adj_p.value_dep","ENSEMBL"))
    colnames(down) <- c("logFC","adj_p.value_dep","ENSEMBL")
    down <- down[order(down$adj_p.value_dep),]
    datatable(
      down,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE,initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
        "}")))
  })

  output$downdl <- downloadHandler(
    filename = function() {
      paste("DEA-DOWN-PDX-Results", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      down <- results$down %>%
        column_to_rownames("sgRNA")%>%
        select(c("estimate","adj_p.value_dep"))
      colnames(down) <- c("logFC","adj_p.value_dep")
      down <- down[order(down$adj_p.value_dep),]
      write.csv(down, file)
    }
  )
  # 
  output$featuress <-
    renderInfoBox({
      req(results$nsignfc)
      infoBox(
        "Number of features passing FC and Pval Filters",
        paste(results$nsignfc,"Among the ",nrow(results$res),"features pass the filters"),
        icon = icon("dna")
      )
    })
  
  output$selected_line <-
    renderUI({
      req(input$celline)
      fluidRow(
      infoBox(        
        "Selected Cell line : ",
        as.character(input$celline),
        icon = icon("dot-circle"),color = "purple",width = 12
      ))
    })

  output$upp_numbers <-
    renderValueBox({
      req(results$up)
      valueBox(
        as.character(nrow(results$up)),
        "Up-regulated features",
        icon = icon("dna"),color = "red"
      )
    })

  output$features_value_box <- renderUI({
    fluidRow(
      column(width = 6,
             valueBoxOutput(ns('down_numbers'),width =  12)),
      column(width = 6,
             valueBoxOutput(ns("upp_numbers"),width =  12)
      ))
  })

  output$down_numbers <-
    renderValueBox({
      req(results$down)
      valueBox(
        as.character(nrow(results$down)),
        "Down-regulated features",
        icon = icon("dna"),color = "blue"
      )
    })
   
  observeEvent(c(input$ExploreIntra,input$FCT,input$PvalsT,input$GeneVolcano),{
    reactives$selectedcomp <- input$ExploreIntra
    reactives$selectedFC <- as.numeric(input$FCT)
    reactives$selectedPvalT <- as.numeric(input$PvalsT) 
    reactives$GeneVolcano <- input$GeneVolcano
  })
  
  selected_comp_rra <- reactiveValues(list = NULL)
  observe({
    selected_comp_rra$list <- as.character(input$ExploreIntra2)
  })
  
  #observeEvent(c(selected_comp_rra$list,results$res, concatenated$results), ignoreInit = TRUE,{
  #observe({
  #print("selected comp rra")
  #print(selected_comp_rra)
  return(list(results=results,reactives=reactives,concatenated=concatenated,selected_comp_rra = selected_comp_rra ))
  #})
  
}

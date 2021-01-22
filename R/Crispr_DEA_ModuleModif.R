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


CRISPRDeaModUIMod <- function(id)  {
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
                                   pickerInput(ns("celline"),
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
                                   pickerInput(ns("comptype"),
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
                                                   pickerInput(ns("Treatlevel"),
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
                                  bsModal(id = ns("compmodal"),"Set up groups two compare",trigger = "comptype",
                                          column(width = 6,
                                              conditionalPanel(condition = 'output.conditionnalTreatment != "1"' ,ns = NS(id),
                                                 pickerInput(ns("TreatlevelInter1"),
                                                      "Group1 Treatment type",choices = NULL,
                                                      options = pickerOptions(
                                                        actionsBox = TRUE,
                                                        title = "Select Treatment",
                                                        liveSearch = TRUE,
                                                        liveSearchStyle = "contains"))),
                                              conditionalPanel(condition = 'output.conditionnalMutTreatment != "1"' ,ns = NS(id),
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
                    tabPanel("Figures",
                             br(),
                             tagList(
                               #fluidRow(column(width = 12,infoBoxOutput(ns("selected_line")))),
                               fluidRow(column(width = 12,uiOutput(ns("selected_line")))),
                               fluidRow(box(title = span(icon("cogs"), "Parameters"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                   status = "success",width= 12,
                                   #fluidRow(
                                     column(width = 6,pickerInput(ns("ExploreIntra"),label= "Select a comparison to explore", selected = NULL,
                                                                  multiple = FALSE,choicesOpt = NULL,inline = FALSE,choices = NULL,
                                                                  options = pickerOptions(
                                                                    title = "Select genes to annotate",
                                                                    liveSearch = TRUE,
                                                                    liveSearchStyle = "contains"))),
                                     column(width = 6,numericInput(ns("FCT"),"LogFC threshold", min = 0, max = 10 , value = 1, step = 0.05)),
                                     column(width = 6,pickerInput(ns("AdjMeth"),"Select an adjustment method", choices = c("BH","BY","holm","none"),
                                                                  multiple = FALSE,choicesOpt = NULL,inline = FALSE,selected = "BH",
                                                                  options = pickerOptions(
                                                                    title = "Select genes to annotate",
                                                                    liveSearch = TRUE,
                                                                    liveSearchStyle = "contains"))),
                                     column(width = 6,numericInput(ns("PvalsT"),"adjusted P values threshold", min = 0, max = 1 , value = 0.05, step = 0.01)))
                               ),
                               fluidRow(uiOutput(ns('features_value_box'))),
                               fluidRow(box(title = span(icon("chart-bar"),"DEA figures"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                   status = "success",width = 12,
                                   #fluidRow(
                                     column(width =12,girafeOutput(ns("Pvals_distrib"))),
                                   #),
                                   br(),
                                   br(),
                                   #fluidRow(
                                     column(width = 12,
                                                   pickerInput(ns("GeneVolcano"),"Select Genes to annotate on volcano",
                                                               selected = NULL,
                                                               multiple = TRUE,
                                                               choicesOpt = NULL,
                                                               inline = FALSE,
                                                               choices = NULL,
                                                               options = pickerOptions(
                                                                 title = "Select genes to annotate",
                                                                 liveSearch = TRUE,
                                                                 liveSearchStyle = "contains",
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
                               )# end of fluidRow
                             ) # end of Taglist
                    ), # end of second tab
                    tabPanel("Tables",
                             br(),
                             fluidRow(
                             box(title = span(icon("arrow-circle-down"),"Tables"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                 status = "success",width= 12,
                                 #fluidRow(
                                   tags$head(tags$style(".butt{background-color:#2E8B57;}")),
                                   column(width =12,
                                          h4("All genes :"),
                                          DT::dataTableOutput(ns('results_table')),
                                          downloadButton(ns("resdl"),"All genes",class = "butt")),
                                   column(width = 12,
                                          br(),
                                          h4("Upp regulated genes :"),
                                          DT::dataTableOutput(ns('up_table')),
                                          downloadButton(ns("uppdl"),"Up-regulated",class = "butt")),
                                   column(width = 12,
                                          br(),
                                          h4("Down regulated genes :"),
                                          #br(),
                                          DT::dataTableOutput(ns('down_table')),
                                          downloadButton(ns("downdl"),"Down-regulated",class = "butt"))
                                 )#,
                             ) # end of box
                    ),
                    tabPanel("RRAscores",
                             fluidRow(column(width =12,
                             pickerInput(ns("ExploreIntra2"),label= "Select a comparison to explore", selected = NULL,
                                         multiple = FALSE,choicesOpt = NULL,inline = FALSE,choices = NULL,
                                         options = pickerOptions(
                                           title = "Select genes to annotate",
                                           liveSearch = TRUE,
                                           liveSearchStyle = "contains")))),
                             fluidRow(DT::dataTableOutput(ns("RRAscores"))),
                             fluidRow(downloadButton(ns("scoresdl"),"Download RRA scores"))
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

CRISPRDeaModServerMod <- function(input, output, session,sampleplan = NULL, var = NULL,
                               norm_data = NULL) {
  
  ### Define reactives #############
  req(sampleplan)
  ns <- session$ns
  
  reactives <- reactiveValues(design = NULL, formula = NULL, contrast = NULL)
  sampleplanmodel <- reactiveValues(table = NULL)
  results <- reactiveValues(res = NULL, up = NULL, down = NULL,nsignfc = NULL,v = NULL,boxplots = NULL,
                            scores = NULL,old_res = "NULL")
  
  observeEvent(c(input$celline,sampleplan$table),{
    sampleplanmodel$table <- sampleplan$table %>%
      filter(Cell_line == input$celline)
  })
  
  observeEvent(c(input$comptype,input$newcomparison), {
    if(input$comptype == "Inter-Treatment"){
    toggleModal(session, "compmodal", toggle = "open")
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
    if(length(unique(sampleplanmodel$table$Sample_description)) >= 2){
    updatePickerInput(session=session,"Mut1",
                      choices = unique(as.character(sampleplanmodel$table[,"Sample_description"])),
                      selected = unique(as.character(sampleplanmodel$table[,"Sample_description"]))[1])
    updatePickerInput(session=session,"Mut2",
                      choices = unique(as.character(sampleplanmodel$table[,"Sample_description"])),
                      selected = unique(as.character(sampleplanmodel$table[,"Sample_description"]))[2])
    }
  })
  #observeEvent(c(input$TreatlevelInter1,input$TreatlevelInter2,input$Mut1,input$Mut2),{
    output$Intercompresume1 <- renderText({
      paste0("You are about to compare these samples : ")
    })
    output$Intercompresume2 <- renderText({
      if(length(unique(sampleplanmodel$table$Treatment)) >= 2 & length(unique(sampleplanmodel$table$Sample_description)) >= 2){
      paste0(input$TreatlevelInter1,"_",input$Mut1)
      } else if (length(unique(sampleplanmodel$table$Treatment)) == 1  & length(unique(sampleplanmodel$table$Sample_description)) >= 2){
        paste0(input$Mut1)
      } else if (length(unique(sampleplanmodel$table$Treatment)) >=2  & length(unique(sampleplanmodel$table$Sample_description)) == 1){
        paste0(input$TreatlevelInter1)
      }
    })
    output$Intercompresume3 <- renderText({
      "vs"
    })
    output$Intercompresume4 <- renderText({
      if(length(unique(sampleplanmodel$table$Treatment)) >= 2 & length(unique(sampleplanmodel$table$Sample_description)) >= 2){
      paste0(input$TreatlevelInter2,"_",input$Mut2)
      } else if (length(unique(sampleplanmodel$table$Treatment)) == 1  & length(unique(sampleplanmodel$table$Sample_description)) >= 2){
      paste0(input$Mut2)
      } else if (length(unique(sampleplanmodel$table$Treatment)) >=2  & length(unique(sampleplanmodel$table$Sample_description)) == 1){
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
      as.character(length(unique(sampleplanmodel$table$Sample_description)))
    })
    outputOptions(output, "conditionnalTreatment", suspendWhenHidden = FALSE)  
      

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
                  data$Inter <- paste0(data$Treatment,"_",data$Sample_description,"_",data$Timepoint)
                  design <- model.matrix(~ 0 + Inter, data = data)
                  colnames(design) <- make.names(gsub("Inter","",colnames(design)))
                } 
                save(list = c('design'),file = "~/coockiecrisprtestRDA/CRISPR_contrast.rda")
                if(input$comptype == "Intra-Treatment"){
                  contrast <- purrr::map(glue::glue("{input$Treatlevel}_T{minT:maxT}-{input$Treatlevel}_T{T0}"),~makeContrasts(contrasts = .x,levels=design))
                } else if (input$comptype == "Inter-Treatment"){
                  if (length(unique(sampleplanmodel$table$Sample_description)) > 1 & length(unique(sampleplanmodel$table$Treatment)) > 1){
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
                  } else if (length(unique(sampleplanmodel$table$Sample_description)) > 1 & length(unique(sampleplanmodel$table$Treatment)) == 1){
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

                  } else if (length(unique(sampleplanmodel$table$Sample_description)) == 1 & length(unique(sampleplanmodel$table$Treatment)) > 1){
                  gluelist <- data.frame(gluelistcol = unique(glue::glue("{input$TreatlevelInter1}_{unique(sampleplanmodel$table$Sample_description)}_{data$Timepoint}-{input$TreatlevelInter2}_{unique(sampleplanmodel$table$Sample_description)}_{data$Timepoint}")))
                  missing1 <- unique(glue::glue("{input$TreatlevelInter1}_{unique(sampleplanmodel$table$Sample_description)}_{data$Timepoint}"))
                  missing2 <- unique(glue::glue("{input$TreatlevelInter2}_{unique(sampleplanmodel$table$Sample_description)}_{data$Timepoint}"))
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
                  contrast <- purrr::map(gluelist,~makeContrasts(contrasts = .x,levels=design))
                }
                save(list = c('contrast'),file = "~/coockiecrisprtestRDA/CRISPR_contrast.rda")
                reactives$contrast <- contrast
                reactives$design <- design
              }
            }
  })
  

  ################ Compute Model new block ##########################
  # observeEvent(c(matrix$table,
  #                reactives$contrast),priority = 10,{
  observeEvent(c(norm_data$data$counts,
                 reactives$contrast),priority = 10,{
                   if (!is.null(norm_data$data$counts) && !is.null(reactives$design)){

                     
                     tictoc::tic("Voom transormation and fitting linear model..")
                     counts <- norm_data$data$counts[,colnames(norm_data$data$counts)%in%rownames(reactives$design)]
                     kept <- which(rowSums(edgeR::cpm(counts) >= 1) >= 2)
                     counts <- counts[kept,]
                     results$v <- limma::voomWithQualityWeights(counts, design = reactives$design, normalize.method = "none", span = 0.5, plot = FALSE)
                     res_fit <- limma::lmFit(results$v, method = "ls")
                     res_eb <- eBayes(res_fit, robust = FALSE)
                     fit <- purrr::map(reactives$contrast, ~contrasts.fit(res_fit, contrasts = .x))
                     res_eb <- purrr::map(fit, ~eBayes(.x, robust = FALSE))
                     
                     tab <- purrr::map(res_eb,~process_res(.x))
                     names(tab) <- lapply(tab, function(x){print(unique(as.character(x$term)))})
                     save(list = c('tab'),file = "~/coockiecrisprtestRDA/TAB_DEA.RData")
                     results$res <- tab
                     # save(list = c("tab"), file = "~/coockiecrisprtestRDA/results_res_new.rda")
                   } # end of if NULL
                   toc(log = TRUE)
                 }) # end of observer
  ###########################################################################
  
  createLink <- function(val) {
    #sprintf('<a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s" target="_blank" class="btn btn-primary">Info</a>',val)
    sprintf('<a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s" target="_blank" class="btn btn-primary">Info</a>',gsub("_[0-9]","",gsub("sg","",val)))
    
  }
  
  observeEvent(c(results$res,input$DEGtabs),{
    updatePickerInput(session = session,"ExploreIntra",
                      choices = names(results$res), selected = names(results$res)[1])
    updatePickerInput(session = session,"ExploreIntra2",
                      choices = names(results$res), selected = names(results$res)[1])
    # choices = list(
    #   Group1 = c(opt1 = "g11",
    #              opt2 = "g12",
    #              opt3 = "g13"),
    #   Group2 = c(opt1 = "g21"),
    #   Group3 = c(opt1 = "g31"),
    #   Group4 = c(opt1 = "g41", 
    #              opt2 = "g42",
    #              opt3 = "g43")
    # )
  })
  
  ################ Compute Model old block ##########################
  observeEvent(c(results$res,
                 input$FCT,
                 input$PvalsT,input$DEGtabs,input$ExploreIntra),ignoreInit = TRUE,{
                   if(input$DEGtabs == "Figures"){
                   req(input$ExploreIntra)
                   req(results$res)
                   res <- results$res[[input$ExploreIntra]]
                   #print(res)
                   nsign <- length(which(res$adj_p.value < input$PvalsT))
                   results$nsignfc <- length(which(res$adj_p.value < input$PvalsT & abs(res$estimate) > input$FCT))
                   up <- which(res$adj_p.value_enrich < input$PvalsT)
                   down <- which(res$adj_p.value_dep < input$PvalsT)
                   res$ENSEMBL <- createLink(res$gene)
                   print('end of DEG')
                   results$up <- res[up,]
                   results$down <- res[down,]
                   results$restable <- res
                   }
                 }) # end of observer
  
  n_perm <- 20000
  alpha_thr <- 0.3
  observeEvent(c(input$DEGtabs),{
    if(input$DEGtabs == "RRAscores") {
      if(!is.null(results$res)){
        if(results$res != results$old_res){
        results$old_res <- results$res
        print(results$res)
        withProgress(message = 'Computing per gene RRA score from DEA results', value = 0.5,{
          incProgress(0.3)
          results$scores <- lapply(results$res,function(x){compute_score_RRA(x,alpha_thr = alpha_thr)})
          setProgress(1)
        })
        }# end of progress
      } else {
        showModal(modalDialog(HTML(
          "<b>Please perform differential analysis first : </b></br>
         To do so select at least two variables for the comparison then click on the build button.
        "),
          title = "Variables infos",
          footer = tagList(
            modalButton("Got it"))
        ))
    } # end of if
    }
    
  }) # end of observer
  
  output$RRAscores <- DT::renderDataTable({
    req(input$ExploreIntra2)
    req(results$scores)
    datatable(
      results$scores[[input$ExploreIntra2]],rownames=FALSE,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE,initComplete = JS(
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
      write.csv(results$scores[[input$ExploreIntra2]],file)
    }
  )

  output$Pvals_distrib <- renderGirafe({
    req(results$res)
    res <- results$res[[input$ExploreIntra]]
    plot <- ggplot(data = res) + aes(x = `adj_p.value`) +
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
  observe({
    req(results$res)
    req(input$ExploreIntra)
    res <- results$res[[input$ExploreIntra]]
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
  output$Volcano <- renderPlot({
    req(Volcano$plot)
    req(input$ExploreIntra)
    res <- results$res[[input$ExploreIntra]]
    tic("Rendering Volcano...")
    ggplot <- Volcano$plot +
      geom_point(data = subset(res,gene %in% input$GeneVolcano),
                 color = "purple", alpha = 0.6) +
      ggrepel::geom_text_repel(
        data = subset(res,gene %in% input$GeneVolcano),
        aes(label = gene),
        size = 5,
        force = 2,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )
    return(ggplot)
    toc(log = TRUE)
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
        #print(head(groups_table))
        if(input$comptype == "Intra-Treatment"){
          group1 <- stringr::str_split_fixed(gsub(paste0(input$Treatlevel,"_"),"",input$ExploreIntra),"-",n=2)[,1]
          group2 <- stringr::str_split_fixed(gsub(paste0(input$Treatlevel,"_"),"",input$ExploreIntra),"-",n=2)[,2]
          groups_table <- groups_table[,c("Treatment","Timepoint","Samples")] %>%
            filter(Treatment == input$Treatlevel)
         } else if(input$comptype == "Inter-Treatment"){
           if (length(unique(sampleplanmodel$table$Sample_description)) > 1 & length(unique(sampleplanmodel$table$Treatment)) > 1){
             group1 <- gsub(paste0(input$Mut1,"_"),"",gsub(paste0(input$TreatlevelInter1,"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,1]))
             group2 <- gsub(paste0(input$Mut2,"_"),"",gsub(paste0(input$TreatlevelInter2,"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,2]))
             groups_table <- groups_table[,c("Treatment","Timepoint","Samples","Sample_description")] %>%
               filter(Treatment == c(input$TreatlevelInter1,input$TreatlevelInter2)) %>%
               filter(Sample_description == c(input$Mut1,input$Mut2))  %>%
               mutate(COMP = paste0(Treatment,"_",Sample_description))
           } else if (length(unique(sampleplanmodel$table$Sample_description)) > 1 & length(unique(sampleplanmodel$table$Treatment)) == 1){
             group1 <- gsub(paste0(input$Mut1,"_"),"",gsub(paste0(unique(sampleplanmodel$table$Treatment),"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,1]))
             group2 <- gsub(paste0(input$Mut2,"_"),"",gsub(paste0(unique(sampleplanmodel$table$Treatment),"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,2]))
             groups_table <- groups_table[,c("Treatment","Timepoint","Samples","Sample_description")] %>%
               filter(Treatment %in% unique(sampleplanmodel$table$Treatment)) %>%
               filter(Sample_description == c(input$Mut1,input$Mut2))  %>%
               mutate(COMP = paste0(Treatment,"_",Sample_description))
           } else if (length(unique(sampleplanmodel$table$Sample_description)) == 1 & length(unique(sampleplanmodel$table$Treatment)) > 1){
             group1 <- gsub(paste0(unique(sampleplanmodel$table$Sample_description),"_"),"",gsub(paste0(input$TreatlevelInter1,"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,1]))
             group2 <- gsub(paste0(unique(sampleplanmodel$table$Sample_description),"_"),"",gsub(paste0(input$TreatlevelInter2,"_"),"",stringr::str_split_fixed(input$ExploreIntra,"-",n=2)[,2]))
             groups_table <- groups_table[,c("Treatment","Timepoint","Samples","Sample_description")] %>%
               filter(Treatment %in% c(input$TreatlevelInter1,input$TreatlevelInter2)) %>%
               filter(Sample_description == unique(sampleplanmodel$table$Sample_description))  %>%
               mutate(COMP = paste0(Treatment,"_",Sample_description))
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
   
  output$results_table <- DT::renderDataTable({
    res <- bind_rows(results$up,results$down) %>%
      column_to_rownames("gene")%>%
      select(c("estimate","adj_p.value_dep","adj_p.value_enrich","ENSEMBL")) %>%
      mutate(adj_p.value = min(adj_p.value_dep,adj_p.value_enrich)) %>%
      select(c("estimate","adj_p.value","ENSEMBL"))
    colnames(res) <- c("logFC","adj_p.value","ENSEMBL")
    datatable(
      res,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE,initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
        "}")))
  })

  output$resdl <- downloadHandler(
    filename = function() {
      paste("DEA-PDX-Results", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      #write.csv(results$res, file)
      write.csv(results$res[[input$ExploreIntra]],file)
    }
  )

  output$up_table <- DT::renderDataTable({
    ups <- results$up %>%
      column_to_rownames("gene") %>%
      select(c("estimate","adj_p.value_enrich","ENSEMBL"))
    colnames(ups) <- c("logFC","adj_p.value_enrich","ENSEMBL")
    datatable(
      ups,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE,initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
        "}")))
  })

  output$updl <- downloadHandler(
    filename = function() {
      paste("DEA-UPPS-PDX-Results", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      ups <- results$up %>%
        column_to_rownames("gene") %>%
        select(c("estimate","adj_p.value_enrich"))
      colnames(ups) <- c("logFC","adj_p.value_enrich")
      write.csv(ups, file)
    }
  )

  output$down_table <- DT::renderDataTable({
    downs <- results$down %>%
      column_to_rownames("gene") %>%
      select(c("estimate","adj_p.value_dep","ENSEMBL"))
    colnames(downs) <- c("logFC","adj_p.value_dep","ENSEMBL")
    datatable(
      downs,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE,initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
        "}")))
  })

  output$downdl <- downloadHandler(
    filename = function() {
      paste("DEA-DOWN-PDX-Results", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      downs <- results$down %>%
        column_to_rownames("gene")%>%
        select(c("estimate","adj_p.value_dep"))
      colnames(downs) <- c("logFC","adj_p.value_dep")
      write.csv(downs, file)
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
        "Upp regulated features",
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
        "Down regulated features",
        icon = icon("dna"),color = "blue"
      )
    })
  # 
  return(results)
  
}

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
      tags$style(type='text/css', ".span161 { width: 850px; }")#,
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
                             #fluidPage(
                             box(title = "Creates DEG model",collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                 status = "primary",width= 12,
                                 #fluidPage(
                                 uiOutput(ns("vars_selUI")
                                 )
                                )
                                 #)
                             ), # end of box
                             #), # end of fluidRow
                             fluidRow(
                             box(title = "Comparison",collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                 status = "primary",width= 12,
                                 column(width=6,
                                        uiOutput(ns("Group1")),
                                        textOutput(ns("group1table")),
                                        pickerInput(
                                          ns("move1"),
                                          label = "Move samples to Group 2",
                                          choices = NULL,
                                          selected = NULL,
                                          multiple = TRUE,
                                          choicesOpt = NULL,
                                          inline = FALSE,
                                          options = pickerOptions(
                                            actionsBox = TRUE,
                                            title = "Select samples to move",
                                            liveSearch = TRUE,
                                            liveSearchStyle = "contains",
                                          )
                                        ),
                                        pickerInput(
                                          ns("remove1"),
                                          label = "remove samples from Group1",
                                          choices = NULL,
                                          options = pickerOptions(
                                            actionsBox = TRUE,
                                            title = "Select samples to remove",
                                            liveSearch = TRUE,
                                            liveSearchStyle = "contains",
                                          ),
                                          selected = NULL,
                                          multiple = TRUE,
                                          choicesOpt = NULL,
                                          inline = FALSE
                                        )),
                                 column(width = 6,uiOutput(ns("Group2")),
                                        textOutput(ns("group2table")),
                                        # selectInput(ns("remove2"),"remove samples from Group 1",multiple = TRUE,selected = NULL,choices = NULL),
                                        # selectInput(ns("move2"),"Move samples to Group 2",multiple = TRUE,selected = NULL,choices = NULL))
                                        pickerInput(
                                          ns("move2"),
                                          label = "Move samples to Group 1",
                                          choices = NULL,
                                          selected = NULL,
                                          multiple = TRUE,
                                          choicesOpt = NULL,
                                          inline = FALSE,
                                          options = pickerOptions(
                                            actionsBox = TRUE,
                                            title = "Select samples to move",
                                            liveSearch = TRUE,
                                            liveSearchStyle = "contains"
                                          )
                                        ),
                                        pickerInput(
                                          ns("remove2"),
                                          label = "remove samples from Group 2",
                                          choices = NULL,
                                          options = pickerOptions(
                                            actionsBox = TRUE,
                                            title = "Select samples to remove",
                                            liveSearch = TRUE,
                                            liveSearchStyle = "contains",
                                          ),
                                          selected = NULL,
                                          multiple = TRUE,
                                          choicesOpt = NULL,
                                          inline = FALSE
                                        )))),
                             br(),
                             br(),
                             actionButton(ns("Build"),"Build Model")
                    ), # end of first tabs
                    
                    tabPanel("Figures",
                             br(),
                             tagList(
                               box(title = span(icon("cogs"), "Parameters"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                   status = "success",width= 12,
                                   fluidRow(
                                   ),
                                   fluidRow(
                                     column(width = 6,selectInput(ns("AdjMeth"),"Select an adjustment method", choices = c("BH","none","BY","holm"), selected = "BH")),
                                     column(width = 6,numericInput(ns("PvalsT"),"adjusted P values threshold", min = 0, max = 1 , value = 0.05, step = 0.01))),
                                   fluidRow(
                                     column(width = 12,numericInput(ns("FCT"),"LogFC threshold", min = 0, max = 10 , value = 1, step = 0.05))#,
                                     #column(width = 6, infoBoxOutput(ns("featuress"),tags$style("#featuress {width:230px;}")))
                                   ),
                                   fluidRow(uiOutput(ns('features_value_box')))

                               ), # end of box
                               #fluidRow(
                               #column(width = 12 ,
                               box(title = span(icon("chart-bar"),"DEA figures"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                   status = "success",width = 12,
                                   # fluidRow(column(width =6,girafeOutput(ns("Pvals_distrib"))),
                                   # column(width =6,plotOutput(ns("scatter")))),
                                   fluidRow(column(width =12,girafeOutput(ns("Pvals_distrib")))),
                                   br(),
                                   br(),
                                   fluidRow(column(width = 12,
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
                                   )),
                                   fluidRow(column(width = 12,
                                                   plotOutput(ns("Volcano"))),
                                            column(width = 12,
                                                   conditionalPanel("input.GeneVolcano != undifined && input.GeneVolcano.length <= 1", ns = ns,
                                                                    #textOutput(ns('boxplots_error'))))
                                                                    textOutput(ns('boxplots_error'))),
                                                   conditionalPanel("input.GeneVolcano != undifined && input.GeneVolcano.length >=2", ns = ns,
                                                                    #plotOutput(ns("boxplots"))),
                                                                    plotOutput(ns("boxplots"))))
                                            
                                            
                                            
                                   )
                                   #girafeOutput(ns("Volcano"))))
                                   
                                   
                                   #)
                                   #)
                               )
                             ) # end of Taglist
                    ), # end of second tab
                    tabPanel("Tables",
                             br(),
                             box(title = span(icon("arrow-circle-down"),"Tables"),collapsible = TRUE, collapsed = FALSE,solidHeader = TRUE,
                                 status = "success",width= 12,
                                 #fluidRow(plotOutput(ns("Volcano"))),
                                 fluidRow(
                                   tags$head(tags$style(".butt{background-color:#2E8B57;}")),
                                   column(width =12,
                                          h4("All genes :",style="padding-left:20px"),
                                          br(),
                                          DT::dataTableOutput(ns('results_table')),
                                          downloadButton(ns("resdl"),"All genes",class = "butt")),
                                   #tags$br(),
                                   
                                   column(width = 12,
                                          h4("Upp regulated genes :",style="padding-left:20px"),
                                          br(),
                                          DT::dataTableOutput(ns('up_table')),
                                          downloadButton(ns("uppdl"),"Up-regulated",class = "butt")),
                                   #tags$br(),
                                   column(width = 12,
                                          h4("Down regulated genes :",style="padding-left:20px"),
                                          br(),
                                          DT::dataTableOutput(ns('down_table')),
                                          downloadButton(ns("downdl"),"Down-regulated",class = "butt"))
                                 )#,
                             ) # end of box
                    ),
                    tabPanel("RRAscores",
                    DT::dataTableOutput("RRAscores")
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

CRISPRDeaModServer <- function(input, output, session, matrix = NULL,sampleplan = NULL, var = NULL,
                               norm_data = NULL) {
  
  ### Define reactives #############
  req(sampleplan)
  #req(matrix)
  ns <- session$ns
  
  reactives <- reactiveValues(design = NULL, formula = NULL, contrast = NULL)
  groups <- reactiveValues(Group1 = NULL, Group2 = NULL)
  #sampleplanmodel <- reactiveValues(table = sampleplan$table)
  sampleplanmodel <- reactiveValues(table = NULL)
  results <- reactiveValues(res = NULL, up = NULL, down = NULL,nsignfc = NULL,v = NULL,boxplots = NULL,
                            scores = NULL)
  
  observe({
  sampleplanmodel$table <- sampleplan$table
  })
  
  observeEvent(input$remove1,{
    sampleplanmodel$table[input$remove1,input$var] <- "removed"
  })
  observeEvent(input$move1,{
    sampleplanmodel$table[input$move1,input$var] <- input$Group2sel
  })
  observeEvent(input$remove2,{
    sampleplanmodel$table[input$remove2,input$var] <- "removed"
  })
  observeEvent(input$move2,{
    sampleplanmodel$table[input$move2,input$var] <- input$Group1sel
  })
  
  observeEvent(c(input$var,
                 input$covar,
                 input$ok
  ),{
    if(!is.null(sampleplanmodel$table)){
      if(!is.null(input$var)){
        #if(!is.null(matrix$table)){
        if(!is.null(norm_data$data)){
            
          
          sampleplan <- sampleplanmodel$table
          print("sampleplanformula")
          #print(nrow(sampleplan))
          design.idx <- colnames(sampleplan)
          if(input$covar != ""){
            #print("with covar")
            vector <- c(input$var,input$covar)
            if(input$var != "Create your own groups"){
              formula <- as.formula(
                paste0('~0+',input$var,"+",paste0(input$covar,collapse = "+"),"+",
                       paste0(combn(vector,2,FUN = paste,collapse =":"),collapse = "+"))
              )} else if (input$var == "Create your own groups"){
                formula <- as.formula(
                  paste0('~0+',"personalisedGroup","+",paste0(input$covar,collapse = "+"),"+",
                         paste0(combn(vector,2,FUN = paste,collapse =":"),collapse = "+"))
                )} } else {
                  
                  if(input$var != "Create your own groups"){
                    #print("without covar")
                    formula <- as.formula(
                      paste0('~0+',as.character(input$var))
                    )} else if (input$var == "Create your own groups"){
                      print("without covar")
                      formula <- as.formula(
                        paste0('~0+',"personalisedGroup"))
                    }
                }
          reactives$formula <- formula
          print(reactives$formula)
        }}}
    
  })
  
  
  observeEvent(input$Build,{
    
    if(!is.null(sampleplanmodel$table)){
      if(!is.null(input$var)){
        if(!is.null(input$Group2sel)){
          if(!is.null(input$Group1sel)){
            #if(!is.null(matrix$table)){
            if(!is.null(norm_data$data)){
              if(!is.null(reactives$formula)){
                data <- sampleplanmodel$table
                if (input$var == "Create your own groups"){
                  data[groups$Group2,"personalisedGroup"] <- "Group2"
                  data[groups$Group1,"personalisedGroup"] <- "Group1"
                  completeVec <- complete.cases(data[,"personalisedGroup"])
                  data <- data[completeVec,]
                }
                ########### OLD ##########
                #mat <- matrix$table[,rownames(data)]
                #design <- model.matrix(reactives$formula, data=data)
                #design <- design[which(rownames(design) %in% colnames(mat)), ]
                
                ##### NEW ########
                design <- model.matrix(reactives$formula, data=norm_data$data$samples)
                #design <- model.matrix(~ 0 + Sample_description * Timepoint * Cell_line, data = norm_data$data$samples)
                design <- design[which(rownames(design) %in% colnames(norm_data$data$counts)), ]
                colnames(design) <- make.names(colnames(design))
                
                if (input$var == "Create your own groups"){
                  contrast <-makeContrasts(contrasts = paste0(paste0("personalisedGroup","Group1"),"-",(paste0("personalisedGroup","Group2"))) ,
                                           levels=design)
                } else {
                  contrast <-makeContrasts(contrasts = paste0(paste0(input$var,input$Group1sel),"-",(paste0(input$var,input$Group2sel))),
                                           levels=design)
                }
                save(list = c('contrast'),file = "~/coockiecrisprtestRDA/CRISPR_contrast.rda")
                reactives$contrast <- contrast
                reactives$design <- design
              }
            }
          }}}
    }
    
  })
  
  output$formula <- renderUI({
    if(!is.null(sampleplan$table)){
      print(Reduce(paste,deparse(reactives$formula)))
    } else {
      print("Provides a sampleplan first ")
    }
  })
  
  
  observeEvent(input$help1,{
    
    showModal(modalDialog(HTML(
      "<b>Variable of interest :</b></br>
    The variable that you want to use to create sample groups for the DE analysis </br></br></br>

    <b>Co-variables :</b></br>

    Select co-variables if you want that app to take into account their respective effects on genes' expression.

    "),
      title = "Variables infos",
      footer = tagList(
        modalButton("Got it"),
      )))
    
    
  })
  
observe({
  output$vars_selUI <- renderUI({
    tagList(
      column(width = 5,
               selectInput(ns("var"),"Variable of interest :",choices = c(var,"Create your own groups"),
                           multiple = FALSE, selected = var[1],width = "100%")),
      column(width = 5,
               selectInput(ns("covar"),"Covariables :",choices = c("None" = "",var),
                           multiple = FALSE, selectize = TRUE, selected = "",width = "100%")),
      column(width = 2,
               br(),
               actionButton(ns("help1"),"",icon = icon("info"),width = "100%"))
      #) # end of fluidRrow
    )
  })
})

observe({ 
  output$Group1 <- renderUI({
    tagList(
      if(!is.null(input$var)) {
        #if(length(input$var) != 0) {
        if(input$var != "Create your own groups"){
          selectInput(ns("Group1sel"),"Group 1", choices = na.omit(levels(sampleplan$table[,input$var])),
                      selected= na.omit(levels(sampleplanmodel$table[,input$var])[1]))
          
        }
      }
    )
  })
})
  observe({
    
    if(!is.null(input$var)) {
      if(input$var != "Create your own groups"){
        groups$Group1 <- rownames(sampleplanmodel$table[which(sampleplanmodel$table[,input$var] == input$Group1sel),])
      } else if (input$var == "Create your own groups"){
        
        #div(class = "#boxPopUp1",
        showModal(
          fluidPage(
            modalDialog(
              title = "Create two groups for differential analysis",
              fluidRow(
                column(width = 6,
                       # checkboxGroupInput(ns("createGroup1"), "select samples to add in group 2",
                       #                    choices = rownames(sampleplanmodel$table),
                       #                    selected = NULL)),
                       pickerInput(ns("createGroup1"), "select samples to add in group 1",
                                   choices = rownames(sampleplanmodel$table),
                                   selected = NULL,
                                   multiple = TRUE,
                                   choicesOpt = NULL,
                                   inline = FALSE,
                                   options = pickerOptions(
                                     actionsBox = TRUE,
                                     title = "Select samples to add",
                                     liveSearch = TRUE,
                                     liveSearchStyle = "contains"
                                   ))
                ),
                column(width = 6,
                       # checkboxGroupInput(ns("createGroup2"), "select samples to add in group 2",
                       #                                      choices = rownames(sampleplanmodel$table),
                       #                                      selected = NULL)
                       pickerInput(ns("createGroup2"), "select samples to add in group 2",
                                   choices = rownames(sampleplanmodel$table),
                                   selected = NULL,
                                   multiple = TRUE,
                                   choicesOpt = NULL,
                                   inline = FALSE,
                                   options = pickerOptions(
                                     actionsBox = TRUE,
                                     title = "Select samples to add",
                                     liveSearch = TRUE,
                                     liveSearchStyle = "contains"
                                   ))
                )),
              easyClose = TRUE,
              footer = tagList(
                modalButton(ns("Cancel")),
                actionButton(ns("ok"),"OK")
              )
            ) # end of fluidRow
          )
        )
        #) # end of div
      }
    }
  })
  
  observeEvent(input$ok,{
    groups$Group1 <- rownames(sampleplanmodel$table)[which(rownames(sampleplanmodel$table) %in% input$createGroup1)]
    groups$Group2 <- rownames(sampleplanmodel$table)[which(rownames(sampleplanmodel$table) %in% input$createGroup2)]
    removeModal()
  })
  
  
  observe({
    if(!is.null(groups$Group1)){
      updatePickerInput(session = session,
                        "remove1","remove samples from Group 1",selected = NULL,choices = groups$Group1,
                        pickerOptions(
                          actionsBox = TRUE,
                          title = "Select samples to remove",
                          header = "This is a title"
                        ))
      updatePickerInput(session = session,
                        "move1","Move samples to Group 2",selected = NULL,choices = groups$Group1)
      
    }
  })
  
  observe({
    if(!is.null(groups$Group2)){
      updatePickerInput(session = session,
                        "remove2","remove samples from Group 2",selected = NULL,choices = groups$Group2,
                        pickerOptions(
                          actionsBox = TRUE,
                          title = "Select samples to remove"
                        ))
      updatePickerInput(session = session,
                        "move2","Move samples to Group 1",selected = NULL,choices = groups$Group2)
      
    }
  })
  
  output$group1table <- renderText({
    groups$Group1
  })

  observe({  
  output$Group2 <- renderUI({
    tagList(
      if(input$var != "Create your own groups"){
        selectInput(ns("Group2sel"),"Group 2", choices  = na.omit(levels(sampleplan$table[,input$var])),
                    selected = na.omit(levels(sampleplanmodel$table[,input$var])[2]) )
      }
    )
    
  })
})

  observe({
    if(!is.null(input$var)) {
      if(input$var != "Create your own groups"){
        groups$Group2 <- rownames(sampleplanmodel$table[which(sampleplanmodel$table[,input$var] == input$Group2sel),])
      }
    }
  })
  
  output$group2table <- renderText({
    groups$Group2
  })
  
  observeEvent(c(input$Group1sel,input$Group2sel),ignoreInit = TRUE,{
    if(!is.null(input$Group1sel)){
      #if(length(input$Group1sel) != 0){
      if(input$Group1sel != "" & input$Group2sel != ""){
        if(input$Group1sel == input$Group2sel){
          showModal(modalDialog(p("You must select two different groups to make a comparison"),
                                title = "Identical groups",
                                footer = tagList(
                                  modalButton("Got it"),
                                )))
        }
      }
    }
  })
  
  observeEvent({input$var
    input$covar},{
      if(input$var %in% input$covar){
        validationModalModel(
          msg = "You can not use the same sample Plan column for variable and covariable ",
        )
        return(-1)
      }
    })
  
  ## Function def
  validationModalModel <- function(msg = "", title = "Model Error") {
    showModal(modalDialog(p(msg),
                          title = title,
                          footer = tagList(
                            modalButton("Dismiss"),
                          )))
    
  }
  
################ Compute Model new block ##########################
  # observeEvent(c(matrix$table,
  #                reactives$contrast),priority = 10,{
  observeEvent(c(norm_data$data$counts,
                 reactives$contrast),priority = 10,{

                   #if (!is.null(matrix$table) && !is.null(reactives$design)){
                   if (!is.null(norm_data$data$counts) && !is.null(reactives$design)){

                     tictoc::tic("Voom transormation and fitting linear model..")
                     #counts <- matrix$table[,colnames(matrix$table)%in%rownames(reactives$design)]
                     counts <- norm_data$data$counts[,colnames(norm_data$data$counts)%in%rownames(reactives$design)]
                     kept <- which(rowSums(edgeR::cpm(counts) >= 1) >= 2)
                     counts <- counts[kept,]
                     results$v <- limma::voomWithQualityWeights(counts, design = reactives$design, normalize.method = "none", span = 0.5, plot = FALSE)
                     res_fit <- limma::lmFit(results$v, method = "ls")
                     res_eb <- eBayes(res_fit, robust = FALSE)
                     
                     fit <- contrasts.fit(res_fit, contrasts = reactives$contrast)
                     res_eb <- eBayes(fit, robust = FALSE)
                     
                     save(list = c('fit','res_eb'),file = "~/coockiecrisprtestRDA/CRISPR_DEA.RData")
                     tab <- biobroom::tidy.MArrayLM(res_eb)
                     ## add unilateral pvalues
                     tab$p.value_dep <- pt(tab$statistic, df = res_eb$df.total[1], lower.tail = TRUE)
                     tab$p.value_enrich <- pt(tab$statistic, df = res_eb$df.total[1], lower.tail = FALSE)
                     ## compute FDR
                     tab$adj_p.value <- p.adjust(tab$p.value, method = input$AdjMeth)
                     tab$adj_p.value_dep <- p.adjust(tab$p.value_dep, method = input$AdjMeth)
                     tab$adj_p.value_enrich <- p.adjust(tab$p.value_enrich, method = input$AdjMeth)
                     tab <- tab[order(tab$adj_p.value),]
                     results$res <- tab
                     save(list = c("tab"), file = "~/coockiecrisprtestRDA/results_res_new.rda")
                     
                     # res <- topTable(fit, number=nrow(counts), adjust.method=input$AdjMeth)
                     # res <- res[order(res$adj.P.Val),]
                     # res$genes <- rownames(res)
                     #results$res <- res
                   } # end of if NULL
                   toc(log = TRUE)
                 }) # end of observer
###########################################################################
  
  createLink <- function(val) {
    sprintf('<a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s" target="_blank" class="btn btn-primary">Info</a>',val)
  }
  
  ################ Compute Model old block ##########################
  observeEvent(c(results$res,
                 input$FCT,
                 input$PvalsT),ignoreInit = TRUE,{

                   res <- results$res
                   nsign <- length(which(res$adj_p.value < input$PvalsT))
                   results$nsignfc <- length(which(res$adj_p.value < input$PvalsT & abs(res$estimate) > input$FCT))
                   #up <- which(res$adj_p.value_enrich < input$PvalsT & res$estimate > input$FCT)
                   #down <- which(res$adj_p.value_dep < input$PvalsT & res$estimate < -input$FCT)
                   up <- which(res$adj_p.value_enrich < input$PvalsT)
                   down <- which(res$adj_p.value_dep < input$PvalsT)
                   #res$t <- NULL
                   #res$P.Value <- NULL
                   #res$B <- NULL
                   #res$label <- NULL
                   res$ENSEMBL <- createLink(rownames(res))
                   print('end of DEG')
                   results$up <- res[up,]
                   results$down <- res[down,]
                   results$restable <- res
                 }) # end of observer

n_perm <- 20000
alpha_thr <- 0.3
observeEvent(c(input$DEGtabs),{
      if(input$DEGtabs == "RRAscores") {
        if(!is.null(results$res)){
        print(results$res)
        withProgress(message = 'Computing per gene RRA score from DEA results', value = 0.5,{
        incProgress(0.3)
          results$scores <- compute_score_RRA(results$res,alpha_thr = alpha_thr)
          scores <- results$scores
          save(list = c("scores"), file = "~/coockiecrisprtestRDA/RRA_scores_table.rda")
        setProgress(1)
        print(score_RRA)
       }) # end of progress
      } else {
      showModal(modalDialog(HTML(
        "<b>Please perform differential analysis first : </b></br>
         To do so select at least two variables for the comparison then click on the build button.
        "),
      title = "Variables infos",
      footer = tagList(
      modalButton("Got it"),
      )))
    }
  } # end of if
}) # end of observer

output$RRAscores <- DT::renderDataTable({
  results$scores
})
  
  
  output$Pvals_distrib <- renderGirafe({
    req(results$res)
    # plot <- ggplot(data = results$res) + aes(x = `P.Value`) +
    #   geom_histogram_interactive(fill = "steelblue",breaks = seq(0, 1, length.out = 20))
    # build  <- ggplot_build(plot)
    # plot <- plot +  labs(title = "P values distribution", x = "P values", y = "Occurences")# +
    # ggiraph::girafe(code = {print(plot)})
    #plot <- ggplot(data = results$res) + aes(x = `p.value`) +
    plot <- ggplot(data = results$res) + aes(x = `adj_p.value`) +
      geom_histogram_interactive(fill = "steelblue",breaks = seq(0, 1, length.out = 30))
    build  <- ggplot_build(plot)
    plot <- plot +  labs(title = "P values distribution", x = "P values", y = "Occurences")# +
    ggiraph::girafe(code = {print(plot)})

  })
  
  #input$GeneVolcano
  observe({
    #updatePickerInput("GeneVolcano", session = session, choices = rownames(matrix$table))
    updatePickerInput("GeneVolcano", session = session, choices = rownames(norm_data$data$counts))
    
  })
  
  Volcano <- reactiveValues(plot = NULL)
  observe({
    req(results$res)
    tic("Ploting Volcano")
    ggplot <- ggplot(results$res, aes(x = estimate, y = -log10(p.value))) +
      ggtitle(colnames(reactives$contrast)) +
      scale_fill_gradient(low = "lightgray", high = "navy") +
      scale_color_gradient(low = "lightgray", high = "navy") +
      expand_limits(y = c(min(-log10(results$res$p.value)), 1)) +
      geom_point(data = results$res,
                 color = "grey", alpha = 0.5) +
      geom_point(data = subset(results$res, estimate > input$FCT),
                 color = "red", alpha = 0.5) +
      geom_point(data = subset(results$res, estimate < -input$FCT),
                 color = "blue", alpha = 0.5) +
      geom_point(data = subset(results$res, adj_p.value < input$PvalsT),
                 color = "green", alpha = 0.5) +
      geom_hline(yintercept = -log10(max(subset(results$res, adj_p.value < input$PvalsT)$p.value)), linetype = "dashed") +
      geom_vline(xintercept = c(-input$FCT, input$FCT), linetype = "dashed") +
      theme_linedraw() +
      theme(panel.grid = element_blank()) +
      xlab("Fold change (log2)") +
      ylab("-log10(P-Value)")

    Volcano$plot <- ggplot
    toc(log = TRUE)
  })
  
  output$Volcano <- renderPlot({

    req(Volcano$plot)
    tic("Rendering Volcano...")
    ggplot <- Volcano$plot +
      #geom_point(data = subset(results$res,genes %in% input$GeneVolcano),
      geom_point(data = subset(results$res,gene %in% input$GeneVolcano),
                 color = "purple", alpha = 0.6) +
      ggrepel::geom_text_repel(
        #data = subset(results$res, adj.P.Val < input$PvalsT),
        #data = results$res[which(rownames(results$res) %in% input$GeneVolcano),],
        #data = subset(results$res,genes %in% input$GeneVolcano),
        data = subset(results$res,gene %in% input$GeneVolcano),
        #aes(label = results$res$label),
        #aes(label = genes),
        aes(label = gene),
        size = 5,
        force = 2,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )
    return(ggplot)
    toc(log = TRUE)
  })
  
  observeEvent(input$GeneVolcano,{

    if(length(input$GeneVolcano) >1){
      req(matrix$table)
      req(sampleplanmodel$table)
      groups_table <- sampleplanmodel$table
      groups_table$Samples <- rownames(groups_table)
      if(input$var != "Create your own groups"){
        groups_table <- groups_table[,c(input$var,"Samples")]
      } else {
        groups_table[groups$Group2,"personalisedGroup"] <- "Group2"
        groups_table[groups$Group1,"personalisedGroup"] <- "Group1"
        completeVec <- complete.cases(groups_table[,"personalisedGroup"])
        groups_table <- groups_table[completeVec,]
        groups_table <- groups_table[,c("personalisedGroup","Samples")]
      }
      boxplotdata <- results$v$E[which(rownames(results$v$E) %in% input$GeneVolcano),]
      boxplotdata <- rbind(boxplotdata,colnames(boxplotdata))
      rownames(boxplotdata)[nrow(boxplotdata)] <- "Samples"
      boxplotdata <- as.data.frame(t(boxplotdata)) %>%  gather(key = "GENE",value = "COUNTS", -Samples)
      boxplotdata$Samples <- as.character(boxplotdata$Samples)
      boxplotdata <- inner_join(boxplotdata,groups_table, by = "Samples")
      boxplotdata$COUNTS <- as.numeric(boxplotdata$COUNTS)
      if(input$var != "Create your own groups"){
        boxplotdata[,input$var] <- as.character(boxplotdata[,input$var])
        results$boxplots <- ggplot(boxplotdata, aes(x=GENE, y=COUNTS, fill = GENE)) +
          geom_boxplot() +
          geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.2,
                                                   seed = 1234),
                     pch=21,
                     aes_string(fill=input$var), show.legend = T)
      } else {
        boxplotdata[,"personalisedGroup"] <- as.character(boxplotdata[,"personalisedGroup"])
        results$boxplots <- ggplot(boxplotdata, aes(x=GENE, y=COUNTS, fill = GENE)) +
          geom_boxplot() +
          geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.2,
                                                   seed = 1234),
                     pch=21,
                     # size = 2,
                     aes_string(fill="personalisedGroup"), show.legend = T)
      }
    }
  })
  
  # output$boxplots <- renderPlot(results$boxplots)
  # output$boxplots_error <- renderText({
  #   validate(
  #     need(length(input$GeneVolcano) >= 2, "Select at least two genes to draw boxplots...")
  #   )
  # })
  
  output$results_table <- DT::renderDataTable({
    a <- results$restable %>%
      column_to_rownames("gene")
    #a$genes <- NULL
    datatable(
      a,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
  
  
  output$resdl <- downloadHandler(
    filename = function() {
      paste("DEA-PDX-Results", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      write.csv(results$res, file)
    }
  )
  
  output$up_table <- DT::renderDataTable({
    a <- results$up %>%
      column_to_rownames("gene")
    #a$genes <- NULL
    datatable(
      a,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE))})
  
  output$updl <- downloadHandler(
    filename = function() {
      paste("DEA-UPPS-PDX-Results", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      write.csv(results$up, file)
    }
  )
  
  
  output$down_table <- DT::renderDataTable({
    a <- results$down %>%
      column_to_rownames("gene")
    #a$genes <- NULL
    datatable(
      a,escape = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE))})
  
  output$downdl <- downloadHandler(
    filename = function() {
      paste("DEA-DOWN-PDX-Results", Sys.Date(), ".csv", sep=",")
    },
    content = function(file) {
      write.csv(results$down, file)
    }
  )
  
  output$featuress <-
    renderInfoBox({
      req(results$nsignfc)
      infoBox(
        "Number of features passing FC and Pval Filters",
        paste(results$nsignfc,"Among the ",nrow(results$res),"features pass the filters"),
        icon = icon("dna")
      )
    })
  
  output$upp_numbers <-
    #renderInfoBox({
    renderValueBox({
      req(results$up)
      # column(width = 12,
      #infoBox(
      valueBox(
        as.character(nrow(results$up)),
        "Upp regulated features",
        icon = icon("dna"),color = "red"
        # "Number of features passing FC and Pval Filters",
        # paste(nrow(results$up),"Upp uppregulated features"),
      )
      # )
    })
  
  output$features_value_box <- renderUI({
    fluidRow(
      column(width = 6,
             #infoBoxOutput(ns("down_numbers")),
             valueBoxOutput(ns('down_numbers'),width =  12)),
      column(width = 6,
             valueBoxOutput(ns("upp_numbers"),width =  12)
      ))
  })
  
  output$down_numbers <-
    #renderInfoBox({
    renderValueBox({
      
      req(results$down)
      #column(width = 12,
      valueBox(
        as.character(nrow(results$down)),
        "Down regulated features",
        icon = icon("dna"),color = "blue"
      )
      #)
    })
  
  #observeEvent(c(results$res,results$nsignfc),{
  #observe({
  #})
  return(results)
  #})
  
  
}

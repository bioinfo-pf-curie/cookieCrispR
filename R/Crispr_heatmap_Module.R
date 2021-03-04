#' @title Draw clusters heatmap from counts matrix UI side
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
#'

ClusteringUIMod <- function(id){
  
  ns <- NS(id)
  #ui <- shinyUI(
  fluidPage(
    tagList(
      tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 13px;width: 105px}
                 .selectize-dropdown { font-size: 12px; line-height: 13px; }
                 .form-group, .selectize-control {margin-left:-10px;max-height: 100px !important;}
                 .box-body {
          padding-bottom: 0px;
      }"),
      sidebarLayout(
        sidebarPanel(width = 12,
                     fluidRow(
                       box(title = 'Clustering data',collapsible = TRUE,collapsed = FALSE,width = NULL, status = "primary",solidHeader = TRUE,
                           #htmltools::h4('Data'),
                           #uiOutput(ns('data')),
                           fluidPage(
                             #checkboxInput(ns('showSample'),'Subset Data'),
                             #conditionalPanel('input.showSample == 1',ns = ns,hr(),uiOutput(ns('sample'))),
                             #uiOutput(ns('sample')),
                             #htmltools::hr(),htmltools::h4('Data Preprocessing'),
                             #column(width=4,selectizeInput(ns('transpose'),'Transpose',choices = c('No'=FALSE,'Yes'=TRUE),selected = FALSE)),
                             #column(width=4,selectizeInput(ns("transform_fun"), "Transform", c(Identity=".",Sqrt='sqrt',log='log',Scale='scale',Normalize='normalize',Percentize='percentize',"Missing values"='is.na10', Correlation='cor'),selected = '.')),
                             #uiOutput(ns('annoVars')),
                             pickerInput(ns('annoVar'),'Annotation',
                                                              choices = NULL,
                                                              selected=NULL,
                                                              multiple=TRUE,
                                                              options = pickerOptions(
                                                                actionsBox = TRUE,
                                                                title = "Select variables for annotation",
                                                                liveSearch = TRUE,
                                                                liveSearchStyle = "contains",
                                                              )),
                             #htmltools::br(),htmltools::hr(),htmltools::h4('Row dendrogram'),
                             column(width=12,sliderInput(ns("r"), "Number of row Clusters", min = 1, max = 15, value = 2)),
                             #column(width=4,numericInput("r", "Number of Clusters", min = 1, max = 20, value = 2, step = 1)),
                             #htmltools::br(),htmltools::hr(),htmltools::h4('Column dendrogram'),
                             column(width=12,sliderInput(ns("c"), "Number of column Clusters", min = 1, max = 15, value = 2)),
                             #column(width=4,numericInput("c", "Number of Clusters", min = 1, max = 20, value = 2, step = 1)),
                             
                           )
                       ) # end of fluidPage
                     ) # end of box and fluidRow
        ),
        #mainPanel(
        mainPanel(width = 12,
                  tabsetPanel(
                    tabPanel("Heatmaply",
                             #htmltools::tags$a(id = 'downloadData', class = paste("btn btn-default shiny-download-link",'mybutton'), href = "", target = "_blank", download = NA, icon("clone"), 'Download Heatmap as HTML'),
                             #htmltools::tags$head(htmltools::tags$style(".mybutton{color:white;background-color:blue;} .skin-black .sidebar .mybutton{color: green;}") ),
                             #plotly::plotlyOutput(ns("heatout"),height=paste0(500,'px')),
                             #column(width = 12,plotly::plotlyOutput(ns("heatout"),height="100%",width = "700px")),
                             fluidRow(
                               column(width = 12,plotly::plotlyOutput(ns("heatout"),height="100%",width = "100%"))),
                             column(width = 12,
                                    br(),
                                    br(),
                                    fluidRow(
                                      box(title = "Additionnal Parameters", collapsible = TRUE,
                                          collapsed = TRUE, status = "primary", width = NULL, solidHeader = TRUE,
                                          fluidPage(
                                            #htmltools::br(),htmltools::hr(),
                                            #htmltools::h4('Additional Parameters'),
                                            column(3,checkboxInput(ns('showColor'),'Color')),
                                            column(3,checkboxInput(ns('showMargin'),'Layout')),
                                            column(3,checkboxInput(ns('showDendo'),'Dendrogram')),
                                            htmltools::hr(),
                                            conditionalPanel('input.showColor==1',ns = ns,
                                                             htmltools::hr(),
                                                             htmltools::h4('Color Manipulation'),
                                                             uiOutput(ns('colUI')),
                                                             sliderInput(ns("ncol"), "Set Number of Colors", min = 1, max = 256, value = 256),
                                                             checkboxInput(ns('colRngAuto'),'Auto Color Range',value = TRUE),
                                                             conditionalPanel('!input.colRngAuto',ns=ns,uiOutput(ns('colRng')))
                                            ),
                                            conditionalPanel('input.showDendo==1',ns = ns,
                                                             htmltools::hr(),
                                                             htmltools::h4('Dendrogram Manipulation'),
                                                             column(width=12,selectInput(ns('dendrogram'),'Dendrogram Type',choices = c("both", "row", "column", "none"),selected = 'both')),
                                                             #htmltools::br(),htmltools::hr(),
                                                             htmltools::h4('Row dendrogram'),
                                                             #column(width=6,selectizeInput(ns("distFun_row"), "Distance method", c(Euclidean="euclidean",Maximum='maximum',Manhattan='manhattan',Canberra='canberra',Binary='binary',Minkowski='minkowski'),selected = 'euclidean')),
                                                             column(width=6,selectizeInput(ns("distFun_row"), "Distance method", c(Euclidean="euclidean"),selected = 'euclidean')),
                                                             column(width=6,selectizeInput(ns("hclustFun_row"), "Clustering linkage", c(Complete= "complete",Single= "single",Average= "average",Mcquitty= "mcquitty",Median= "median",Centroid= "centroid",Ward.D= "ward.D",Ward.D2= "ward.D2"),selected = 'complete')),
                                                             htmltools::br(),htmltools::hr(),htmltools::h4('Column dendrogram'),
                                                             #column(width=6,selectizeInput(ns("distFun_col"), "Distance method", c(Euclidean="euclidean",Maximum='maximum',Manhattan='manhattan',Canberra='canberra',Binary='binary',Minkowski='minkowski'),selected = 'euclidean')),
                                                             column(width=6,selectizeInput(ns("distFun_col"), "Distance method", c(Euclidean="euclidean"),selected = 'euclidean')),
                                                             column(width=6,selectizeInput(ns("hclustFun_col"), "Clustering linkage", c(Complete= "complete",Single= "single",Average= "average",Mcquitty= "mcquitty",Median= "median",Centroid= "centroid",Ward.D= "ward.D",Ward.D2= "ward.D2"),selected = 'complete')),
                                                             br(),
                                                             column(width=12,
                                                                    selectizeInput(ns("seriation"), "Seriation", c(OLO="OLO",GW="GW",Mean="mean",None="none"),selected = 'OLO'),
                                                                    sliderInput(ns('branches_lwd'),'Dendrogram Branch Width',value = 0.6,min=0,max=5,step = 0.1))
                                            ),
                                            conditionalPanel('input.showMargin==1',ns = ns,
                                                             htmltools::hr(),
                                                             htmltools::h4('Widget Layout'),
                                                             column(4,textInput(ns('main'),'Title','')),
                                                             column(4,textInput(ns('xlab'),'X Title','')),
                                                             column(4,textInput(ns('ylab'),'Y Title','')),
                                                             sliderInput(ns('row_text_angle'),'Row Text Angle',value = 0,min=0,max=180),
                                                             br(),
                                                             sliderInput(ns('column_text_angle'),'Column Text Angle',value = 45,min=0,max=180),
                                                             sliderInput(ns("l"), "Set Margin Width", min = 0, max = 200, value = 5),
                                                             sliderInput(ns("b"), "Set Margin Height", min = 0, max = 200, value = 5)
                                            )
                                            
                                          )# end of FluidPage
                                      ) # end of Box
                                    ) # end of fluidRow
                             ) # end of columnbox
                    ),
                    tabPanel("Data",
                             fluidRow(dataTableOutput(ns('tables')))
                    )
                  )
        ) # end of box
        # )#  end of shinyUI(
      ) # end of sidebarlayout
    ) # end of tagList
  ) # end of fluidPage
}




#' @param input,output,session standards \code{shiny} server arguments.
#' @param name Does the file have a Header
#' @param data 
#'
#' @export
#'
#'
#' @title Draw clusters heatmap from counts matrix server side
#' @importFrom shiny observeEvent reactiveValues callModule observe icon
#' @importFrom htmltools tags HTML
#' @importFrom plotly plotlyOutput layout renderPlotly
#' @importFrom dplyr mutate_at vars
#' @import heatmaply
#' @importFrom stats cor dist hclust
#' @importFrom xtable xtable
#' @importFrom tools file_path_sans_ext
#' @importFrom rmarkdown pandoc_available pandoc_self_contained_html
#' @importFrom viridisLite viridis
#' @importFrom viridis magma plasma inferno
#' @importFrom SummarizedExperiment assay


# choices = c('Vidiris (Sequential)'="viridis",
#             'Magma (Sequential)'="magma",
#             'Plasma (Sequential)'="plasma",
#             'Inferno (Sequential)'="inferno",

#

ClusteringServerMod <- function(input, output, session, data = NULL, metadata = NULL,printRows = FALSE, vst = FALSE) {
  
  ns <- session$ns
  
  #req(data$table)
  #req(metadata$table)
  reactives <- reactiveValues(metadata = NULL)
  #,variableFeatures = genefilter::rowVars(vst$vars))
  reactives2 <- reactiveValues(selData = NULL)
  
  observeEvent(data$table,{
  reactives2$selData <- data$table
  })
  
  #observeEvent(metadata$table,{
  observe({
    req(metadata$table)
    print('metadata observer')
    reactives$metadata <-  metadata$table
    print(class(metadata$table))
    print(head(metadata$table))
    print(colnames(metadata$table))
    print(class(colnames(metadata$table)))
    updatePickerInput(session = session, "annoVar", choices = colnames(metadata$table))
  })
  
  print("entering heatmap module...")
  
  #observeEvent(metadata$table,{
  observeEvent(reactives$metadata,{
  #observe({
    print("updated annovar")
    print(colnames(reactives$metadata))
    print(colnames(metadata$table))
    # updatePickerInput(session = session, "annovar", choices = colnames(metadata$table))
     #updatePickerInput(session = session, "annovar", choices = colnames(reactives$metadata))
   })
  
   
    observeEvent(reactives2$selData,{
      req(reactives2$selData)
  
    subdata <- reactiveValues(rows = nrow(reactives2$selData),
                              #cols = names(data$table)
                              cols = colnames(reactives2$selData)
    )
    })
    

    output$colUI<-renderUI({
      colSel='Vidiris'
      
      selectizeInput(inputId =ns("pal"), label ="Select Color Palette",
                     choices = c('Vidiris (Sequential)'="viridis",
                                 'Magma (Sequential)'="magma",
                                 'Plasma (Sequential)'="plasma",
                                 'Inferno (Sequential)'="inferno",
                                 'Magma (Sequential)'="magma",
                                 'Magma (Sequential)'="magma",
                                 
                                 'RdBu (Diverging)'="RdBu",
                                 'RdYlBu (Diverging)'="RdYlBu",
                                 'RdYlGn (Diverging)'="RdYlGn",
                                 'BrBG (Diverging)'="BrBG",
                                 'Spectral (Diverging)'="Spectral",
                                 
                                 'BuGn (Sequential)'='BuGn',
                                 'PuBuGn (Sequential)'='PuBuGn',
                                 'YlOrRd (Sequential)'='YlOrRd',
                                 'Heat (Sequential)'='heat.colors',
                                 'Grey (Sequential)'='grey.colors'),
                     selected=colSel)
    })
    
    observeEvent({reactives2$selData},{
      output$colRng=renderUI({
        
        rng=range(reactives2$selData,na.rm = TRUE)
        
        n_data = nrow(reactives2$selData)

        transform_fun <- "."
        min_min_range = ifelse(transform_fun=='cor',-1,-Inf)
        min_max_range = ifelse(transform_fun=='cor',1,rng[1])
        min_value = ifelse(transform_fun=='cor',-1,rng[1])
        max_min_range = ifelse(transform_fun=='cor',-1,rng[2])
        max_max_range = ifelse(transform_fun=='cor',1,Inf)
        max_value = ifelse(transform_fun=='cor',1,rng[2])

        a_good_step = 0.1 # (max_range-min_range) / n_data
        
        list(
          numericInput(ns("colorRng_min"), "Set Color Range (min)", value = min_value, min = min_min_range, max = min_max_range, step = a_good_step),
          numericInput(ns("colorRng_max"), "Set Color Range (max)", value = max_value, min = max_min_range, max = max_max_range, step = a_good_step)
        )
        
      })
    })
    
    
    #interactiveHeatmap<- reactive({
    interactiveHeatmap <- reactiveValues(plot = NULL)
    observeEvent(c(input$annoVar,reactives$metadata,reactives2$selData),priority = -1,{
       
      withProgress(message = 'Runing heatmap', value = 0.5, {
      req(reactives$metadata)
      req(reactives2$selData)
      data.in <- reactives2$selData
 
      ss_num =  sapply(data.in, is.numeric) # in order to only transform the numeric values

      if(input$colRngAuto){
        ColLimits=NULL
      }else{
        ColLimits=c(input$colorRng_min, input$colorRng_max)
      }
      
      if (length(input$distFun_row) != 0){
        if(input$distFun_row != "pearson"){
          distfun_row = function(x) stats::dist(x, method = input$distFun_row)
        } else {
          distfun_row <- function(x) {
            #print(head(x))
            return(1- factoextra::get_dist(x, method = "pearson", stand = FALSE))}
        }
      }
      if (length(input$distFun_col) != 0){
        if(input$distFun_col != "pearson"){
          distfun_col = function(x) stats::dist(x, method = input$distFun_col)
        } else {
          distfun_col <- function(x) {
            #print(head(x))
            return(1- factoextra::get_dist(x, method = "pearson", stand = FALSE))}
        }
      }

      print("hclust on row and col")
      #hclustfun_row = function(x) stats::hclust(x, method = input$hclustFun_row)
      hclustfun_col = function(x) stats::hclust(x, method = input$hclustFun_col)
      print("done")
      
      if(length(input$annoVar)>=1){
        
        samplesAnnot <- reactives$metadata[,input$annoVar, drop = F]
        data.in <- data.in[,colnames(data.in) %in% rownames(samplesAnnot)]
        print(ncol(data.in))
        print(nrow(samplesAnnot))
        print("samplesannoy")
        #print(head(samplesAnnot))
        print("datain")
        #print(head(data.in))
        
        p <- heatmaply::heatmaply(data.in,
                                  main = input$main,xlab = input$xlab,ylab = input$ylab,
                                  row_text_angle = input$row_text_angle,
                                  column_text_angle = input$column_text_angle,
                                  dendrogram = input$dendrogram,
                                  branches_lwd = input$branches_lwd,
                                  seriate = input$seriation,
                                  colors=eval(parse(text=paste0(input$pal,'(',input$ncol,')'))),
                                  distfun_row =  distfun_row,
                                  #hclustfun_row = hclustfun_row,
                                  distfun_col = distfun_col,
                                  hclustfun_col = hclustfun_col,
                                  k_col = input$c,
                                  k_row = input$r,
                                  col_side_colors = samplesAnnot,
                                  showticklabels = c(TRUE, printRows),
                                  limits = ColLimits) %>%
          plotly::layout(margin = list(l = input$l, b = input$b))
        
        #p$elementId <- NULL
        
        #return(p)
        #interactiveHeatmap$plot <- p
        
      } else {
        
        print("datain2")
        print(head(data.in))
        print(ncol(data.in))
        print(nrow(data.in))
        
        p <- heatmaply::heatmaply(data.in,
                                  main = input$main,xlab = input$xlab,ylab = input$ylab,
                                  row_text_angle = input$row_text_angle,
                                  column_text_angle = input$column_text_angle,
                                  dendrogram = input$dendrogram,
                                  branches_lwd = input$branches_lwd,
                                  seriate = input$seriation,
                                  colors=eval(parse(text=paste0(input$pal,'(',input$ncol,')'))),
                                  distfun_row =  distfun_row,
                                  #hclustfun_row = hclustfun_row,
                                  distfun_col = distfun_col,
                                  hclustfun_col = hclustfun_col,
                                  k_col = input$c,
                                  k_row = input$r,
                                  showticklabels = c(TRUE, printRows),
                                  limits = ColLimits) %>%
          plotly::layout(margin = list(l = input$l, b = input$b))
        
      }  
        
      print("hclust on row ")
      
        if (length(input$distFun_row) != 0){
          print("zded")
          if(input$distFun_row != "pearson"){
            # p$dendRow <- cutree(
            #   hclust(dist(data.in,method = input$distFun_row),method = input$hclustFun_row),
            #   k = input$r)
          } else {
            print("zeded")
            
            # p$dendRow <- cutree(
            #   hclust(1- factoextra::get_dist(data.in, method = "pearson", stand = FALSE),
            #          method = input$hclustFun_row),
            #   k = input$r)
          }
        }
      print("done")
        
        if (length(input$distFun_col) != 0){
          print("zdedd")
          
          if(input$distFun_col != "pearson"){
            # p$dendCol <- cutree(
            #   hclust(dist(t(data.in),method = input$distFun_col),method = input$hclustFun_col),
            #   k = input$c)
          } else {
            print("zddedd")
            
            # p$dendCol <- cutree(
            #   hclust(1- factoextra::get_dist(t(data.in), method = "pearson", stand = FALSE),
            #          method = input$hclustFun_col),
            #   k = input$c)
          }

        p$elementId <- NULL
        
        #return(p)
        interactiveHeatmap$plot <- p
      }
    })
    })
    
    observeEvent(interactiveHeatmap$plot,{
      output$heatout <- plotly::renderPlotly({
        interactiveHeatmap$plot
      })
    })
    
    output$tables=renderDataTable(reactives2$selData#,server = TRUE,filter='top',
    )
    
    #Clone Heatmap ----
    observeEvent({interactiveHeatmap$plot},{
      h<-interactiveHeatmap$plot
      
      l<-list(main = input$main,xlab = input$xlab,ylab = input$ylab,
              row_text_angle = input$row_text_angle,
              column_text_angle = input$column_text_angle,
              dendrogram = input$dendrogram,
              branches_lwd = input$branches_lwd,
              seriate = input$seriation,
              colors=paste0(input$pal,'(',input$ncol,')'),
              distfun_row =  input$distFun_row,
              hclustfun_row = input$hclustFun_row,
              distfun_col = input$distFun_col,
              hclustfun_col = input$hclustFun_col,
              k_col = input$c,
              k_row = input$r,
              limits = paste(c(input$colorRng_min, input$colorRng_max),collapse=',')
      )
      
      
      print('hey')
      l=data.frame(Parameter=names(l),Value=do.call('rbind',l),row.names = NULL,stringsAsFactors = FALSE)

      l[which(l$Value==''),2]='NULL'
      paramTbl=print(xtable::xtable(l),type = 'html',include.rownames=FALSE,print.results = FALSE,html.table.attributes = c('border=0'))

      
      h$width='100%'
      h$height='800px'
      s<-htmltools::tags$div(style="position: relative; bottom: 5px;",
                             htmltools::HTML(paramTbl),
                             htmltools::tags$em('This heatmap visualization was created using',
                                                htmltools::tags$a(href="https://github.com/yonicd/shinyHeatmaply/",target="_blank",'shinyHeatmaply'),
                                                Sys.time()
                             )
      )
      print("end of htmltools")
      
      output$downloadData <- downloadHandler(
        filename = function() {
          paste0("heatmaply-", strftime(Sys.time(),'%Y%m%d_%H%M%S'), ".html")
        },
        content = function(file) {
          libdir <- paste0(tools::file_path_sans_ext(basename(file)),"_files")
          
          htmltools::save_html(htmltools::browsable(htmltools::tagList(h,s)),file=file,libdir = libdir)
          
          if (!rmarkdown::pandoc_available()) {
            stop("Saving a widget with selfcontained = TRUE requires pandoc. For details see:\n",
                 "https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md")
          }
          
          rmarkdown::pandoc_self_contained_html(file, file)
          unlink(libdir, recursive = TRUE)
        }
      )
      print("end ov obssrver")
    })
    
  #}) # end of if !is.null(data)
  observeEvent(interactiveHeatmap$plot,ignoreInit = TRUE,{
    return(interactiveHeatmap$plot)
  })
  #} # end of if is nul seldata
}

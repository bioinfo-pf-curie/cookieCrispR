ui_crispr_app <- function(request){
  
  # get modes and themes for the ace editor
  modes <- shinyAce::getAceModes()
  themes <- shinyAce::getAceThemes()
  
  tagList(dashboardPage(skin = "green",
    dashboardHeader(title = "Cookie CRISPR"#,
    ),  # End of Header
    dashboardSidebar(
      sidebarUserPanel("Institut Curie",
                       subtitle = "Bioinformatic Platform",
                       #image = "https://upload.wikimedia.org/wikipedia/commons/7/77/Logo_Institut_Curie.jpg"
                       image = "data:image/jpeg;base64,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"
      ),
      sidebarMenu(id ="sidebarmenu",
        menuItem("Data Input", tabName = "DataInput"),
        menuItem("Descriptive analysis", tabName = "Descriptive_data_analysis",
                 menuSubItem("Raw distributions","Rawdist"),
                 menuSubItem("Temporal evolution","Tev"),
                 menuSubItem("ROC Curves","Roc")#,
                 #menuSubItem("Clustering","Clustering"),
                 #menuSubItem('Compare conditions','CompCond')
                 ),
        menuItem("Statistical analysis",tabName = "Statistical_analysis"),
        menuItem("Clustering", tabName = "Clustering"),
        menuItem("Report", tabName = "Report")
      )
    ),
    dashboardBody(
      tags$head(
        tags$style(
          HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              opacity:0.91;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
          )
        )),
      fluidPage(
        use_cicerone(),
      tabItems(
        tabItem(tabName = "DataInput",
                useShinyjs(),
                #useShinyalert(),
                infoBoxOutput("Totguidenumber",width = 6),
                valueBoxOutput("Depth",width = 6),
                br(),
                fluidRow(
                    box(title = p('Inputs',actionButton("startCicerone0",label=NULL,icon = icon("info-circle")),
                                  actionButton("exampledatasetbut",label=NULL,icon = icon("teeth"))),
                    width = 12,
                    status = "success",solidHeader = TRUE,
                    #title="Inputs",
                        tabsetPanel(type = "tabs",
                                tabPanel("Uploads",
                                         #fluidPage(
                                         #fluidRow(
                                           column(width=6,div(id = "sample_plandiv",fileInput("sample_plan","Sample infos"))),
                                           column(width=6,
                                                  div(id = "orderUIdiv",uiOutput("orderUI"))),
                                           fluidRow(),
                                           #div(id = 'genelistdiv',
                                               column(width=6,
                                                      div(id = "essentialdiv",fileInput("essential","Essential genes"))),
                                               column(width=6,
                                                      div(id = "nonessentialdiv",fileInput("nonessential","Other genes"))),
                                               #),
                                           #fluidRow(
                                          column(width = 6,
                                                 div(id = "countsdiv",fileInput("counts","Global counts"))),
                                          column(width = 6,
                                                 radioButtons("screentype","Screening type :", choices = c("negative","positive"),
                                                              selected = "negative")),
                                           # column(width = 12,fileInput(inputId = "restore", accept = ".rda", label = "Restore Previous analysis",
                                           #                            buttonLabel=list(icon("angle-double-up"))))
                                           ),
                                tabPanel("Help",
                                         uiOutput("Datahelptext"),
                                         fluidRow(
                                           column(width = 4,
                                                  downloadButton("DlTestSplan","DL sample plan example",style = "width:100%;"),#, class = "butt"),
                                                  tags$head(tags$style(".butt{background-color:#add8e6;} .butt{color: #337ab7;}"))
                                                  ),
                                           column(width = 4,downloadButton("DlTestCounts","DL counts matrix example",style = "width:100%;")),#, class = "butt")),
                                           column(width = 4,downloadButton("DlTesGuideList","DL Genes list example",style = "width:100%;"))#, class = "butt"))
                                           )#,
                                         #downloadButton("state_save_sc","Save State as .rda",style = "visibility: hidden;"),
                                         #downloadButton("exit_and_save","Save State as .rda",style = "visibility: hidden;")
                                )
                    )# end of tabset
                    
                    )# end of box
                    #) # end of div
                    ), 
                fluidRow(
                  box(
                    width = 12, status = "success", solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE,
                    title = "Counts table", id ='countstablebox',
                    column(width = 12,
                           div(id = "removegenesdiv",pickerInput("removegenes", "Remove genes for further analysis",
                                       choices = NULL,
                                       selected = NULL,
                                       multiple = TRUE,
                                       choicesOpt = NULL,
                                       inline = FALSE,
                                       options = pickerOptions(
                                         actionsBox = TRUE,
                                         title = "Select genes you want to remove",
                                         liveSearch = TRUE,
                                         liveSearchStyle = "contains"
                                       )))),
                    tabsetPanel(id="countstabset",
                      tabPanel("Rawcounts",br(),DT::dataTableOutput("counts_table")),
                      tabPanel("log10(cpm)",br(),DT::dataTableOutput("normalized_counts_table"))
                    )
                  )),
                fluidRow(
                  box(
                    width = 12, status = "success", solideHeader = TRUE, collapsed = FALSE,collapsible = TRUE,
                    title = "Sample plan", id = "sampleplanbox",
                    column(width = 12,
                           div(id = "removesamplesdiv",pickerInput("removesamples", "Remove samples for further analysis",
                                choices = NULL,
                                selected = NULL,
                                multiple = TRUE,
                                choicesOpt = NULL,
                                inline = FALSE,
                                options = pickerOptions(
                                  actionsBox = TRUE,
                                  title = "Select samples you want to remove",
                                  liveSearch = TRUE,
                                  liveSearchStyle = "contains"
                                )))),
                    column(width = 12,div(style = 'overflow-x: scroll',DT::dataTableOutput("sample_plan_table"))),
                  ),
                  column(width = 12,
                         div(id = 'correlationsAnnotdiv',pickerInput(inputId = "correlationsAnnot","Choose variables to construct annotations",width = '100%',
                                     choices = NULL,
                                     selected = NULL,
                                     multiple = TRUE,
                                     choicesOpt = NULL,
                                     inline = FALSE,
                                     options = pickerOptions(
                                       actionsBox = TRUE,
                                       title = "Select multiple conditions here",
                                       liveSearch = TRUE,
                                       liveSearchStyle = "contains",
                                     ))),
                         plotOutput("correlation_heatmap")),
                  br(),br(),
                  column(width = 6,
                         br(),
                         downloadButton("dlcorrelation_heatmap","Download correlations heatmaps",
                                        style = "width:100%;")),
                  column(width = 6,
                         br(),
                         downloadButton("dlcorrelation_coefficients","Download correlations coefficients",style = "width:100%;")),
                )),
        tabItem(tabName = "Rawdist",
                fluidRow(
                  box(collapsible = TRUE, collapsed = FALSE,
                      width = 12,status = "success",solidHeader = TRUE,
                      title="Read counts",
                      column(width=12,
                             plotOutput("read_number"),
                             downloadButton("dlreadnumber","Download read numbers plot",style = "width:100%;"),
                             br(),br()
                      ))),
                fluidRow(
                   box(collapsible = TRUE, collapsed = FALSE,
                       width = 12, status = "success", solidHeader = TRUE,
                       title = "Normalized log_cpm distributions",
                       column(width = 12,
                       plotOutput("boxplot_all", width = "100%", height = 600),
                       ),
                       column(width = 6,plotOutput("boxplot_noness"),
                              ),
                       column(width = 6,plotOutput("boxplot_ess"),
                              ),
                       br(),
                       column(width = 12,downloadButton("dlbox_all","Download Boxplots ",style = "width:100%;"),br(),br()),
                       # br()
                )),
                fluidRow(
                  box(collapsible = TRUE, collapsed = FALSE,
                      width = 12, status = "success", solidHeader = TRUE,
                      title = "Counts distributions for essential and non essential gene",
                      column(width=6,plotOutput("essential_distribs", width = "100%", height = 600)),
                      column(width=6,plotOutput("nonessential_distribs", width = "100%", height = 600)),
                      downloadButton("splited_distribs","Download distributions per gene categories",style = "width:100%;"),br(),br()
                  ))
        ),
        tabItem("Tev",
                fluidRow(
                      column(width=6,plotOutput("diff_box_all", width = "100%", height = 600),br()),
                      column(width=6,plotOutput("diff_box_ess", width = "100%", height = 600),br()),
                      downloadButton("dldiffboxes","Download difference to zero boxes",style = "width:100%;"),
                      br()
                )
            ),
            tabItem("Roc",
                    fluidRow(
                      column(width = 12,
                         checkboxInput("labels","Print AUC labels on plots",value = FALSE),
                         div(style = 'overflow-x: scroll',plotOutput("roc")),
                         br(),
                         downloadButton("dlROC","Download ROC plots",style = "width:100%;"),#class = "butt"
                         br(),br(),
                         DT::dataTableOutput('auc'),
                         br(),
                         downloadButton("dlauc","Download AUCs table",style = "width:100%;")#, class = "butt"
                  )
                )
            ),
            tabItem("Clustering",
                  fluidRow(
                    column(width = 12,infoBoxOutput("InfoCompHeatmap", width = 12)),
                    br(),
                      column(width = 12,ClusteringUIMod(id = "heatmapID"))
                  )
        ),
        tabItem("CompCond",
                fluidRow(
                column(width=8,pickerInput(inputId = "conditionreference1","Choose conditions to compare",width = '100%',
                                                     choices = NULL,
                                                     selected = NULL,
                                                     multiple = TRUE,
                                                     choicesOpt = NULL,
                                                     inline = FALSE,
                                                     options = pickerOptions(
                                                       actionsBox = TRUE,
                                                       title = "Select multiple conditions here",
                                                       liveSearch = TRUE,
                                                       liveSearchStyle = "contains",
                                                     ))),
                column(width=4,checkboxInput("splitcelline","Split boxplots by cell line ?"))),
                fluidRow(
                  column(width=8,
                         pickerInput(inputId = "selectguidescomp","Annotate guide on boxplots",width = '100%',
                                                    choices = NULL,
                                                    selected = NULL,
                                                    multiple = FALSE,
                                                    choicesOpt = NULL,
                                                    inline = FALSE,
                                                    options = pickerOptions(
                                                      actionsBox = TRUE,
                                                      title = "Select guides here",
                                                      liveSearch = TRUE,
                                                      liveSearchStyle = "contains",
                                                    ))),
                                column(width=4,pickerInput(inputId = "selecttimepointscomp","Remove timepoints from boxplots data",width = '100%',
                                                            choices = NULL,
                                                            selected = NULL,
                                                            multiple = TRUE,
                                                            choicesOpt = NULL,
                                                            inline = FALSE,
                                                            options = pickerOptions(
                                                              actionsBox = TRUE,
                                                              title = "Select timepoints conditions here",
                                                              liveSearch = TRUE,
                                                              liveSearchStyle = "contains",
                                                            )))),# end of fluidRow
                br(),
                fluidRow(column(width =12,girafeOutput("positive_boxplots", width = "90%", height = "90%")))
        ),
        tabItem(tabName = "Statistical_analysis",
                CRISPRDeaModUI(id = "DEA")),
        tabItem(tabName = "Report",
                          h2("About this report"),
                          h4("This content has been loaded from the template report `.Rmd` file. Please edit it at your best convenience!"),
                          h4("Some <!–html_preserve–> text may appear on this preview. Do not worry it will disappear when generating and saving the report."),            icon = icon("pencil"),
                          h1("Report Editor"),
                          fluidRow(
                             column(width = 12,
                                box(width = '100%',
                                title = "markdown options", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                                column(width = 6,radioButtons("rmd_dl_format", label = "Choose Format:", c("HTML" = "html", "R Markdown" = "rmd"), inline = TRUE),
                                textInput("report_title", "Title: "),
                                textInput("report_author", "Author: "),
                                radioButtons("report_toc", "Table of Contents", choices = list("Yes" = "true", "No" = "false"))
                                ),
                                column(width = 6,
                                       fluidRow(
                                              checkboxGroupInput("include",label = "Add/remove sections in the report",
                                                                 choices = c("Data Overview" = "DataOverview","read_numbers",
                                                                             "density_ridges","temporal_evolution",
                                                                             "Volcano plots","SessionInfo")),
                                       br(),
                                       br(),
                                       br(),
                                       radioButtons("report_ns", "Number sections", choices = list("Yes" = "true", "No" = "false")),
                                       pickerInput("report_theme", "Theme",
                                                   choices = list("Default" = "default", "Cerulean" = "cerulean",
                                                                  "Journal" = "journal", "Flatly" = "flatly",
                                                                  "Readable" = "readable", "Spacelab" = "spacelab",
                                                                  "United" = "united", "Cosmo" = "cosmo"),multiple = FALSE,
                                                   choicesOpt = NULL,
                                                   inline = FALSE,
                                                   options = pickerOptions(
                                                     actionsBox = TRUE,
                                                     liveSearch = TRUE,
                                                     liveSearchStyle = "contains"
                                                   )))
                                       )
                                    )
                                )
                          ),
                          fluidRow(
                            column(3,
                                   actionButton("updatepreview_button", "Update report",class = "btn btn-primary"),p()
                            ),
                            column(3, downloadButton("saveRmd", "Generate & Save",class = "btn btn-success"))
                          ),
                          tabBox(
                            width = NULL,
                            id="report_tabbox",
                            tabPanel("Report preview",
                                     icon = icon("file-text"),
                                     htmlOutput("knitDoc"),
                                     aceEditor("acereport_rmd", mode="markdown",theme = "solarized_light",autoComplete = "live",
                                               value="_Initialization of the_ `CRISPRApp` _report generation..._",
                                               placeholder = "You can enter some code and text in R Markdown format",
                                               height="800px"),
                                     #box(
                                     #   title = "editor_options", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                                     #   column(width = 12,
                                              checkboxInput("enableAutocomplete", "Enable AutoComplete", TRUE),
                                              #conditionalPanel(
                                               # "input.enableAutocomplete",
                                                #wellPanel(
                                                  checkboxInput("enableLiveCompletion", "Live auto completion", TRUE),
                                                  checkboxInput("enableRCompletion", "R code completion", TRUE),
                                               # )
                                    #          ),
                                              pickerInput("mode", "Mode: ", choices=modes, selected="markdown",
                                                          multiple = FALSE,
                                                          choicesOpt = NULL,
                                                          inline = FALSE,
                                                          options = pickerOptions(
                                                            actionsBox = TRUE,
                                                            liveSearch = TRUE,
                                                            liveSearchStyle = "contains"
                                                          )),
                                              pickerInput("theme", "Theme: ", choices=themes, selected="solarized_light",
                                                          multiple = FALSE,
                                                          choicesOpt = NULL,
                                                          inline = FALSE,
                                                          options = pickerOptions(
                                                            actionsBox = TRUE,
                                                            liveSearch = TRUE,
                                                            liveSearchStyle = "contains"
                                                          ))
                                              
                                     #  ))
                                     ),
                            tabPanel(
                              "About the app", icon = icon("institution"),
                              includeMarkdown(system.file("extdata", "about.md",package = "CRISPRApp")),
                              hr(),
                            )
                          )
                ) # end of tabBox
                ))
      )
  ),
  tags$footer(
    wellPanel(
      HTML(
        '
      <style>
      .footer {
        background: #0ea743;
        font-family: "serif", cursive;
        font-weight: 500;
        line-height: 1.1;
        color: #ffffff;
        height : 150px;
      }
      </style>

      <div class="footer">
      </br>
      <p align="center" width="4">Documentation available on the <a href="https://gitlab.curie.fr/r-shiny/bioshiny/blob/Master/Applications_tutorials/Per_App_Documentation/COOKIE_CRISPR/COOKIE_CRISP-R_DOC.html" style="color:#041518;"> bioshiny git page </a></p>
	<br/>
      <p align="center" width="4">Developped by: <a href="mailto: benoitclement.sand@gmail.com" style="color:#041518;">Clément BENOIT</a> & <a href="mailto: pierre.gestraud@curie.fr" style="color:#041518;">Pierre Gestraud</a></p>
      <p align="center" width="4"> For the <a href ="https://science.curie.fr/plateformes/criblage-genetique-crispr-crisprit/" style="color:#041518;">CRISPR IT Platform of the Institut Curie</a></p>
      </div>  
      '
      )
    )),
    tags$script(src = "imgModal.js")
  )
}
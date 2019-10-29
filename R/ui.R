dashboardPage(
  dashboardHeader(title = "CRISPR Screening"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Input", tabName = "DataInput"),
      menuItem("Descriptive data analysis", tabName = "Descriptive_data_analysis"),
      menuItem("Results", tabName = "Results")


    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "DataInput",
              fluidRow(
                valueBoxOutput("Totguidenumber"),
                valueBoxOutput("Depth")
              ),
              fluidRow(
                box(
                  width = 12,status = "info",solidHeader = TRUE,
                  title="Inputs",
                  #column(width=4,fileInput("list","Enter your essential and non essential genes list")),
                  fluidRow(
                    column(width=3,fileInput("sample_plan","Sample infos")),
                    column(width=3,fileInput("counts","Global counts")),
                    column(width=3,fileInput("essential","Essential genes")),
                    column(width=3,fileInput("nonessential","Non Essential genes"))
                  ))),
              
              fluidRow(
                box(
                  width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE,
                  title = "Counts table",
                  div(style = 'overflow-x: scroll', DT::dataTableOutput("counts_table"))
                )),
              fluidRow(
                box(
                  width = 12, status = "info", solideHeader = TRUE, collapsed = FALSE,collapsible = TRUE,
                  title = "Sample plan",
                  div(style = 'overflow-x: scroll',DT::dataTableOutput("sample_plan_table")),
                ))),
      tabItem(tabName = "Descriptive_data_analysis",

              fluidRow(
                box(collapsible = TRUE, collapsed = TRUE,
                    width = 12,status = "info",solidHeader = TRUE,
                    title="Read counts",
                    column(width=12,
                           #div(style = 'overflow-x: scroll',DT::dataTableOutput("read_number"))
                           plotOutput("read_number"),
                           downloadButton("dlreadnumber","Download read numbers plot")
                    ))),
              fluidRow(
                box(collapsible = TRUE, collapsed = TRUE,
                    width = 12, status = "info", solidHeader = TRUE,
                    title = "Normalized log_cpm distributions",
                    plotOutput("boxplot_all", width = "100%", height = 600),
                    downloadButton("dlbox_all","Download Boxplots")
                )),
              fluidRow(
                box(collapsible = TRUE, collapsed = TRUE,
                    width = 12,status = "info",solidHeader = TRUE,
                    title="Density ridges",
                    plotOutput("density_ridge", width = "100%", height = 600),
                    downloadButton("dldensity_ridge","Download Densitity ridges plots")
                )),
              fluidRow(
                box(collapsible = TRUE, collapsed = TRUE,
                    width = 12, status = "info", solidHeader = TRUE,
                    title = "Counts distributions for essential and non essential gene",
                    column(width=6,plotOutput("essential_distribs", width = "100%", height = 600)),
                    column(width=6,plotOutput("nonessential_distribs", width = "100%", height = 600)),
                    downloadButton("splited_distribs","Download distributions per gene categories")
                ))
              ),
      tabItem(tabName = "Results",
              box(title = "Comparison settings : ", solidHeader = TRUE, collapsible = TRUE, status = "info", width = 12,
              column(width=12,uiOutput("orderUI")),
              fluidRow()),
              #column(width=6,selectInput(inputId = "conditionreference2", "Choose a condition :", choices = c())),
              #column(width=6,selectInput(inputId = "conditionreference1", "Choose a second condition compare with :", choices = c()))),
              fluidRow(
                box(collapsible = TRUE, collapsed = TRUE,
                    width = 12, status = "info", solidHeader = TRUE,
                    title = "Difference to initial timepoint",
                    column(width=6,plotOutput("diff_box_all", width = "100%", height = 600)),
                    column(width=6,plotOutput("diff_box_ess", width = "100%", height = 600)),
                    downloadButton("dldiffboxes","Download difference to zero boxes")
                ),
                box(collapsible = TRUE, collapsed = TRUE,
                    width = 12, status = "info", solidHeader = TRUE,
                    title = "ROC",
                div(style = 'overflow-x: scroll',plotOutput("roc")),
                downloadButton("dlROC","Download ROC plots")
                ),
                )
      )


               
      )

    )
  )


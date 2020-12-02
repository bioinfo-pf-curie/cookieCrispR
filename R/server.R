#' Server function
#'
#' Function containing all the server side of the Shiny App, must be called inside the run_app() function 
#'
#'
#' @param input input
#' @param output output
#' @param session session
#'
#' @return None
#' @importFrom shinyjs runjs
#' @importFrom shinyAce updateAceEditor aceEditor aceAutocomplete
#' @import knitr
#' @importFrom xfun embed_files
#' @import kableExtra
#' @importFrom dplyr first
#' @importFrom gridExtra grid.arrange

server_crispr_app <- function(input, output, session) {
    
  ######### Global options ##############
  session$onSessionEnded(stopApp)
  options(shiny.sanitize.errors = TRUE,shiny.maxRequestSize = 3000*1024^2)
  session$allowReconnect(TRUE)

  ##### Usefull variables #############
  reactives <- reactiveValues(sampleplan = NULL,sampleplanGood = FALSE, sampleplanRaw = NULL,
                              joined = NULL,countsRaw = NULL, counts = NULL)
  separators <- reactiveValues(counts = NULL,sampleplan = NULL)
  
  ##### Upload files and datatable construction ####################
  ## counts table    
  observeEvent(c(input$counts,
                 separators$counts),{
      req(input$counts)
      inFile <- input$counts

      if (is.null(separators$counts)){
        counts <- read.table(inFile$datapath, sep = ",", header = TRUE)
      } else {
        counts <- read.table(inFile$datapath, sep = separators$counts, header = TRUE)
      }
      if (!("X" %in% colnames(counts))){
        print("a")
      showModal(modalDialog(p(""),
                              title = "Incorrect count matrix format",
                              tagList(h6('Check the help section of the app')),
                              footer = tagList(
                                modalButton("Got it"))
      ))
      } else{
        print("b")
      counts <- dplyr::rename(counts, sgRNA = .data$X)
      counts <- dplyr::select(counts, -.data$sequence)
      reactives$countsRaw <- counts
      }
})
  
observeEvent(reactives$countsRaw,{
  
      counts  <- reactives$countsRaw
      
      annot_sgRNA <- dplyr::select(counts, .data$sgRNA, Gene = .data$gene)
      # 
      counts <- gather(counts, value = "count", key = "Sample_ID", -.data$sgRNA, -.data$gene)
      counts <- dplyr::mutate(counts, Sample_ID = gsub(".R[1-9].fastq","",.data$Sample_ID))
      counts <- counts %>%
      dplyr::group_by(.data$Sample_ID) %>%
         dplyr::mutate(cpm = 1e6 * .data$count / sum(.data$count), log_cpm = log10(1 + .data$cpm)) %>%
         dplyr::ungroup()
     
     reactives$counts <- counts
      
  })  # end of observer
    
    ## samples annotations
observeEvent(c(input$sample_plan,
              separators$sampleplan),{#input$restore

      # req(input$sample_plan)
      # inFile <- input$sample_plan
      # samples <- read_delim(inFile$datapath, delim = "|",  col_names = c("Sample_ID", "sample","condition"))
      # samples <- samples %>%
      #   separate(sample, into = c("date","rep","Cell_line","day"), remove = FALSE) %>%
      #   mutate(day = as.factor(.data$Timepoint)) %>%
      #   mutate(condition = as.factor(.data$Treatment)) %>%
      #   mutate(Timepoint_num = as.numeric(gsub("DAY","",.data$Timepoint))) %>% 
      #   mutate(day = fct_reorder(.f = factor(.data$Timepoint), .x = .data$Timepoint_num))
      # return(samples)
      
      req(input$sample_plan)
      inFile <- input$sample_plan
      # samples <- read.table(inFile$datapath, sep = input$FSs,header = TRUE)
      if (is.null(separators$sampleplan)){
      samples <- read.table(inFile$datapath, sep = ";", header = TRUE)
      } else {
      samples <- read.table(inFile$datapath, sep = separators$sampleplan, header = TRUE)
      }
      reactives$sampleplanRaw <- samples

}) # end of observer
    
observeEvent(reactives$sampleplanRaw,{    
      
      if(reactives$sampleplanGood == TRUE){
        
      print(head(reactives$sampleplanRaw))
      print("Arranging sample plan")
      samples <- reactives$sampleplanRaw %>%
      #   #separate(sample, into = c("date","rep","Cell_line","day"), remove = FALSE) %>%
          mutate(Timepoint = as.factor(.data$Timepoint)) %>%
          mutate(Treatment = as.factor(.data$Treatment)) %>%
          mutate(Timepoint_num = as.numeric(gsub("[^0-9.-]", "", .data$Timepoint))) %>% 
          mutate(Timepoint = fct_reorder(.f = factor(.data$Timepoint), .x = .data$Timepoint_num))
      print("DONE")
      
      reactives$sampleplan <- samples
      }
    
}) # end of observer
    
    ## Join counts and annotations
    #joined <- reactive({
observeEvent(c(reactives$counts,reactives$sampleplan,
               input$timepoints_order),{
  
      if(!is.null(reactives$counts) & !is.null(reactives$sampleplan)){
      samples <- reactives$sampleplan
      counts <- reactives$counts
      
      if (TRUE %in% unique(!(unique(counts$Sample_ID) %in% samples$Sample_ID))){
      if(input$sidebarmenu == "DataInput"){
       showModal(modalDialog(p(""),
                    title = "Missing samples in provided sampleplan",
                    tagList(h6('The folowing samples :'),
                            h6(paste0(unique(counts[!(counts$Sample_ID %in% samples$Sample_ID),]$Sample_ID),collapse = " | "),style="color:red"),
                            h6("are missing in the sampleplan."),
                            p(),
                            h6("Please check that Samples IDs are correctly matching between ssplan and count matrix files, unless they will be removed for the following analysis")),
                    footer = tagList(
                      modalButton("Got it"),
                    ))
       )
      }
      counts <- counts  %>%
          filter(Sample_ID %in% samples$Sample_ID)
      }
      if(!(is.null(input$timepoints_order))){
      counts <- full_join(counts, samples) %>%
        #mutate(day = factor(.data$Timepoint, levels = input$timepoints_order))
        mutate(Timepoint = factor(.data$Timepoint, levels = input$timepoints_order))
      } else {
      counts <- full_join(counts, samples) %>%
          #mutate(day = factor(.data$Timepoint, levels = input$timepoints_order))
          mutate(Timepoint = factor(.data$Timepoint, level = unique(.data$Timepoint)))
      }
      
      print("joined_runed")
      reactives$joined <- counts
      }# end of if
}) # End of observer    
    #timepoints <- reactiveValues(a = reactives$joined$day)
    
    ess_genes <- reactive({
      req(input$essential)
      inFile <- input$essential
      ess_genes <- read.table(inFile$datapath, header=FALSE)
      return(ess_genes)
    })
    
    non_ess_genes <- reactive({
      req(input$nonessential)
      inFile <- input$nonessential
      non_ess <- read.table(inFile$datapath, header = FALSE)
      return(non_ess)
    })
    
    #Compute difference to 0
    diff_t0 <- reactive({
      req(input$timepoints_order)
      withProgress(message = 'Difference to initial timepoint calculation', value = 0.5, {
      counts <- reactives$joined
      firstpoint <- input$timepoints_order[[1]]
      counts$Timepoint <- relevel(as.factor(counts$Timepoint), ref = firstpoint)
      fin <- counts %>%
        group_by(sgRNA, Cell_line, Replicate) %>%
        arrange(Timepoint) %>%
        mutate(diff = log_cpm - first(log_cpm)) %>% 
        ungroup()
     
      print("DONE")
      return(fin)
      })
    })

######### Plots and tables outputs ####################################
observe({
  print(input$restore)
  
  if (!is.null(reactives$counts)){
    output$counts_table <- DT::renderDataTable({
      DT::datatable(reactives$counts, rownames = FALSE)
      
    })
}
}) # end of observer

observe({
    print(input$restore)
    output$sample_plan_table <- DT::renderDataTable({
      if (!is.null(reactives$sampleplan)){
      print(class(reactives$sampleplan))
      sample_plan <- reactives$sampleplan
      sample_plan <- dplyr::select(sample_plan, -.data$Timepoint_num)
      return(DT::datatable(sample_plan,rownames = FALSE))
      }
    })
}) # end of observer
       
    output$joined <- DT::renderDataTable({
      diff_t0()
    })
    
    read_number <- reactive({

      counts <- as.data.frame(reactives$joined) %>%
        group_by(.data$Sample_ID, .data$Replicate, .data$Cell_line, .data$Timepoint) %>%
        summarise(total = sum(.data$count)) %>% 
        as.data.frame() %>%
        ggplot(aes(x = .data$Timepoint, y = .data$total)) +
        geom_col(position = position_dodge()) + facet_wrap(vars(.data$Replicate, .data$Cell_line), nrow = 1) +
        labs(title = "Number of reads per sample", xlab ="Timepoint", ylab ="Total counts")
      #scale_x_continuous(breaks = seq(min(counts$day),max(counts$day)))
      return(plot(counts))
      #return(counts)
    })
    
    output$read_number <- renderPlot({
      plot(read_number())
    })
    
    output$dlreadnumber <- downloadHandler(
      filename = function(){
        paste("Read_number_",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(read_number())
        dev.off()
      }
    )
    
    boxplot_all <- reactive({
      counts <- reactives$joined
      counts %>% 
        ggplot(aes(x = .data$Timepoint, y = .data$log_cpm, fill = .data$Replicate)) + geom_boxplot() + facet_grid(. ~ .data$Cell_line) + 
        labs(title = "Distribution of normalized log-cpm by sample", subtitle = "All guides")
    })
    
    output$boxplot_all <- renderPlot({
      plot(boxplot_all())
    })
    
    output$dlbox_all <- downloadHandler(
      filename = function(){
        paste("Boxplots_all_samples",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(boxplot_all())
        dev.off()
      }
    )
    
    density_ridge <- reactive({
      counts <- reactives$joined
      withProgress(message = 'Calculating density ridges', value = 0.5, {
      incProgress(0.3)
      counts <-  counts %>%
        ggplot(aes(x = .data$log_cpm, y = .data$Timepoint)) +
        geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
        facet_wrap(vars(.data$Cell_line, .data$Replicate), ncol = 1, strip.position = "right") +
        labs(title = "Distribution of normalzed log-cpm by sample", subtitle = "All guides")
      return(plot(counts))
      }) # end of progress
    })
    
    output$density_ridge <- renderPlot({
      plot(density_ridge())
    })
    
    output$dldensity_ridge <- downloadHandler(
      filename = function(){
        paste("Density_ridges_",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(density_ridge())
        dev.off()
      }
    )
    
    essential_distribs <- reactive({
      counts <- reactives$joined
      ess_genes <- ess_genes()
      withProgress(message = 'Calculating density ridges', value = 0.5, {
        incProgress(0.3)
      counts <- filter(counts, .data$gene %in% ess_genes$V1)
      counts_plot <- ggplot(counts, aes(x = .data$log_cpm, y = .data$Timepoint)) + 
        geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
        facet_wrap(vars(.data$Cell_line, .data$Replicate), ncol = 1, strip.position = "right") +
        scale_fill_viridis_c() +
        labs(title = "Distribution of normalised log-cpms", subtitle = "Essential genes")
      return(plot(counts_plot))
      })
    })
    
    output$essential_distribs <- renderPlot({
      plot(essential_distribs())
    })
    
    nonessential_distribs <- reactive({
      counts <- reactives$joined
      ess_genes <- non_ess_genes()
      counts <- counts %>%
        filter(.data$gene %in% ess_genes$V1)
      
      counts_plot <- ggplot(counts, aes(x = .data$log_cpm, y = .data$Timepoint)) +
        geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
        facet_wrap(vars(.data$Cell_line, .data$Replicate), ncol = 1, strip.position = "right") +
        scale_fill_viridis_c() +
        labs(title = "Distribution of normalised log-cpms", subtitle = "Non Essential genes")
      
      return(plot(counts_plot))
    })
    
    output$nonessential_distribs <- renderPlot({
      plot(nonessential_distribs())
    })
    
    output$splited_distribs <- downloadHandler(
      filename = function(){
        paste("Splitted_distributions_",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(essential_distribs())
        plot(nonessential_distribs())
        dev.off()
      }
    )
    
    output$orderUI <- renderUI({
      counts <- levels(reactives$sampleplan$Timepoint)
      print(counts)
      orderInput(inputId = 'timepoints', label = 'Re-order your timepoints here :', items = counts, as_source = F)
    })
    
    ############## NEGATIV SCREENING ################
    diff_boxes <- reactiveValues(diff_box_ess = NULL,diff_box_all = NULL)
    #diff_box_all <- reactive({
    observe({
    firstpoint <- input$timepoints_order[[1]]
      
      print("diff box all")
      diff_boxes$diff_box_all <- diff_t0() %>%
        filter(.data$Timepoint != firstpoint) %>%
        ggplot(aes(x = .data$Timepoint, y = .data$diff, fill = .data$Replicate)) + geom_boxplot() + facet_grid(.data$Treatment ~ .data$Cell_line) +
        #ggplot(aes(x = .data$Timepoint, y = .data$diff, fill = .data$Replicate)) + geom_boxplot() + facet_grid(paste0(".data$",input$FacetChoice,collapse = " ~ ")) +
        ylab(paste0("diff_",firstpoint)) + 
        labs(title = paste0("Boxplots of log fold change from ", firstpoint ," - all guides"))
      
      print('DONE')
      #plot(diff_boxes$diff_box_all)
    })
    
    output$diff_box_all <- renderPlot({
      #plot(diff_box_all())
      plot(diff_boxes$diff_box_all)
    })
    
    # diff_box_ess <- reactive({
    observe({
      firstpoint <- input$timepoints_order[[1]]
      ess_genes <- ess_genes()
      
      diff_boxes$diff_box_ess <-  diff_t0() %>% 
        #select(-condition, -log_cpmt0, -diff) %>%
        filter(.data$Timepoint != !!firstpoint) %>%
        filter(.data$gene %in% ess_genes[,1]) %>%
        ggplot(aes(x = .data$Timepoint, y= .data$diff, fill = .data$Replicate)) + geom_boxplot() + facet_grid(.~ .data$Cell_line) +
        ylab(paste0("diff_",firstpoint)) + 
        labs(title = paste0("Boxplots of log fold change from ", firstpoint ," - essential genes's guides"))
      
      #plot(diff_boxes$diff_box_ess)
    })
    
    output$diff_box_ess <- renderPlot({
      #plot(diff_box_ess())
      plot(diff_boxes$diff_box_ess)
    })
    
    output$dldiffboxes <- downloadHandler(
      filename = function(){
        paste("Diff_boxes_",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(diff_box_ess())
        plot(diff_box_all())
        dev.off()
      }
    )
    
#################### ROC ####################
ROC <- reactiveValues(plot = NULL,AUC = NULL)
    
#observeEvent(input$sidebarmenu,{
observe({
  if (input$sidebarmenu == "Pscreen"){
    if(is.null(input$essential) | is.null(input$nonessential)){
      showModal(
        modalDialog(tagList(h3("You must provide essentials and non essentials genes list to perform positive screening")),
               footer = tagList(
               modalButton("Got it"))
               )
        )
    } else {
    ess_genes <- ess_genes()
      non_ess_genes <- non_ess_genes()
      #print(head(diff_t0()))
      withProgress(message = 'Calculating ROC curves', value = 0.5, {
      d <- diff_t0() %>% select(.data$sgRNA, .data$Cell_line, .data$Replicate, .data$Timepoint, .data$gene, .data$Treatment, .data$diff) %>%
      #d <- as.data.frame(diff_t0()) %>% select(.data$sgRNA, .data$Cell_line, .data$Replicate, .data$Timepoint, .data$gene, .data$Treatment, .data$diff) %>%
        group_by(.data$Timepoint, .data$Treatment, .data$Cell_line, .data$Replicate) %>%
        arrange(.data$diff) %>% 
        mutate(type = case_when(
          .data$gene %in% ess_genes[,1] ~ "+",
          .data$gene %in% non_ess_genes[,1] ~ "-", 
          TRUE ~ NA_character_)) %>%
        filter(!is.na(.data$type)) %>%
        mutate(TP = cumsum(.data$type == "+") / sum(.data$type == "+"), FP = cumsum(.data$type == "-") / sum(.data$type == "-")) %>%
        ungroup()
      print("DONE")
      
      print("Calulating AU ROC Curves")
      by_rep <- split(d, f= d$Replicate)
      by_rep_cl <- unlist(lapply(X = by_rep, FUN = function(x){split(x, f = x$Cell_line)}),recursive = FALSE)
      by_rep_cl_treat <- unlist(lapply(X = by_rep_cl, FUN = function(x){split(x, f = x$Treatment)}),recursive = FALSE)
      by_rep_cl_treat_tpt <- unlist(lapply(X = by_rep_cl_treat, FUN = function(x){split(x, f = x$Timepoint)}),recursive = FALSE)
      
      AUC_values <- lapply(X= by_rep_cl_treat_tpt, FUN = function(x){
        print(x)
        orders <- order(x$FP)
        print(sum(diff(x$FP[orders])*rollmean(x$TP[orders],2)))
      }
      )
      print(AUC_values)
      
      AUC_table <- data.frame(Replicate = NA, Cell_line = NA, Treatment = NA, Timepoint = NA, AUC = NA)
      #for (table in by_rep_cl_treat_tpt){
      for (tablenum in 1:length(by_rep_cl_treat_tpt)){
        vars <- names(by_rep_cl_treat_tpt)[tablenum]
        print(vars)
        row <- unlist(strsplit(vars, split =".",fixed = TRUE))
        row <- c(row,as.character(AUC_values[tablenum]))
        print(row)
        AUC_table <- rbind(AUC_table,row)
      }
      ROC$AUC <- AUC_table
      print("DONE")

      ROC$plot <- d %>% ggplot(aes(x = FP, y = TP, color = Timepoint)) + geom_abline(slope = 1, lty = 3) + geom_line() + facet_grid(Treatment + Cell_line ~ Replicate) + coord_equal()
      })
    }}
})
    
output$roc <- renderPlot({
   req(ROC$plot)
   plot(ROC$plot)
})
    
    output$dlROC <- downloadHandler(
      filename = function(){
        paste("ROC_plots",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        #plot(ROC())
        plot(ROC$plot)
        dev.off()
      }
    )

output$auc <- renderDataTable({
   req(ROC$AUC)
   ROC$AUC
})    

output$dlauc <- downloadHandler(
  filename = function(){
    paste("AUC_table",Sys.Date(),".pdf",sep="")
  },
  content = function(file){
    write.table(x = ROC$plot,file = file,sep = ",", quote = FALSE, row.names = FALSE)
  }
)
    
    ########################## Observers ############################
  observe({
      
      updateSelectInput(session,"conditionreference2",
                        choices = reactives$sampleplan$condition)
      
      updateSelectInput(session,"conditionreference1",
                        choices = reactives$sampleplan$condition)
    })
    
    ########################################################
    ################## Help section ############
    
    output$Totguidenumber <- renderInfoBox(
      infoBox(
        "Total guides number :",
        as.numeric(nrow(reactives$countsRaw))
      )
    )
    
    mean_depth <- reactive({
      req(reactives$counts)
      depth <- reactives$counts %>% 
        group_by(Sample_ID) %>% 
        summarise(m = mean(count)) %>% 
        summarise(m_all = mean(m))
      return(round(depth$m_all))
    })
    
    output$Depth <- renderInfoBox(
      infoBox(
        "Average sequencing depth across samples :",
        mean_depth()
      )
    )
    
    output$Datahelptext <- renderUI({HTML(
      "

<ol> 
<li>Sample infos file :
<br/><br/>
<B>Format description :</B>
<br/>
<ul>
<li>A csv file using ';' or ',' as field separator, each line of the file must respects the following format specifications :<br/>
Sample_ID;Sample_description;Cell_line;Timepoint;Treatment;Replicate</li>
<li>Column names must respect Upper and lower case. </li>
<li>Values in the table must not contain spaces, use '_' instead. </li>
</ul>
<br/>
<B>For example :</B>
<br/>
D308R1;Ctrl_R1;HEK;T0;Ref_R1;Replica_1
<br/>
<br/>
<br/>
</li><li>Counts table file :</li>
<br/>
<B>Format description :</B>
<br/>
<ul>
<li>A csv formated text file using ; or , as field separator. </li>
<li>The first line of the file is a header, it contains samples Sample_IDs as colnames and two supplementary columns called 'gene
' and 'sequence'. </li>
<li>AGuides names' are rownames. </li>
<li>AValues are read counts.</li>
</ul>
<br/>
<B>For example :</B><br/>
\"\";\"D115R01\";\"D115R02\";\"gene\";\"sequence\"<br/>
\"sgITGB8_1\";339;379;\"ITGB8\";\"AAAACACCCAGGCTGCCCAG\"
<br/>
<br/>
<br/>
</li><li>Genes' lists</li><br/>
<B>Format description :</B>
<br/>
Essentials and non essentials genes' list file just contains one column with all the genes of the list.<br/>
<br/>
<B>For example:</B><br/>

Gene1<br/>
Gene2<br/>
<br/>

...


</ol>

  
  
  
  "
    )})
    
    metadata_path <- system.file("extdata", "SampleDescriptiondatatest.txt", package = "CRISPRApp")
    counts_path <- system.file("extdata", "global_counts_table_datatest.csv", package = "CRISPRApp")
    essential_path <- system.file("extdata", "essentials.csv", package = "CRISPRApp")
    non_essential_path <- system.file("extdata", "non_essentials_datatest.csv", package = "CRISPRApp")
    
    output$DlTestCounts <- downloadHandler(
      filename = function() {
        paste("COOKIE_CRISPR_COUNTS_EXAMPLE-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        file.copy(from = counts_path, to = file)
      }
    )
    
    output$DlTestSplan <- downloadHandler(
      filename = function() {
        paste("COOKIE_CRISPR_SAMPLEPLAN_EXAMPLE-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        file.copy(from = metadata_path, to = file)
      }
    )
    
    output$DlTesGuideList <- downloadHandler(
      filename = function() {
        paste("COOKIE_CRISPR_GENESLIST_EXAMPLE-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        file.copy(from = essential_path, to = file)
        #file.copy(from = non_essential_path , to = file)
      }
    )
    
    ## Modal Dialogue 
    
    dataModal <- function(failed = FALSE) {
        modalDialog(p(""),
                    title = "Set up inputs files formats specifications",
                    selectInput("FSc","Field separator (counts)",choices=c(";",","), selected = separators$counts),
                    selectInput("FSs","Field separator (sampleplan)",choices=c(";",","), selected = separators$sampleplan),
                    footer = tagList(
                      actionButton("apply",label = "Apply this parameters",icon = icon("check")),
                      modalButton("Dismiss"),
                    )
            )
    
    }
    
    observeEvent(input$settings,{
    showModal(dataModal())
    })
    
    observeEvent(input$apply,{
      separators$sampleplan <- as.character(input$FSs)
      separators$counts <- as.character(input$FSc)
      removeModal()
    })
    
    observeEvent(c(reactives$sampleplanRaw),priority = 10,{
    
    colnames <- colnames(reactives$sampleplanRaw)
    if( (!"Sample_ID" %in% colnames|
        !"Sample_description" %in% colnames|
        !"Cell_line" %in% colnames|
        !"Timepoint" %in% colnames|
        !"Treatment" %in% colnames|
        !"Replicate" %in% colnames)  ){
      showModal(modalDialog(tagList(h3('Colnames must contain these values :',style="color:red"),
                                    h4("| Sample_ID | Sample_description | Cell_line | Timepoint | Treatment | Replicate |"),
                                    #p(),
                                    h4("colnames must respect majuscules",style="color:red"),
                                    p(),
                                    h3("Your current file colnames are : ",strong ="bold"),
                                    #h4(paste(colnames,collapse ="  "))
                                    h4(paste("| ",paste(colnames,collapse =" | ")," |"))
                                    
                                    ),
                            title = "Anormal sample plan columns naming",
                            footer = tagList(
                              modalButton("Return to input tab"),
                          )))
      
    reactives$sampleplanGood <- FALSE
    } else {
    reactives$sampleplanGood <- TRUE
    }
       
    })
    
########### Save State ############
    observeEvent(c(input$init2,
                   input$init)
                 ,priority =10,ignoreInit = TRUE,{
      
      cat("save data \n")
      saveState(filename = "/tmp/WorkingEnvironment.rda",
                 reactives= reactives,
                 separators = separators,
                 input = input)
    })
    
    output$exit_and_save <- downloadHandler(
      filename = function() {
        paste0("COOKIE_CRISPR_rState_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".rda")
      },
      content = function(file) {
        #saveState(filename)
        file.copy(from = "/tmp/WorkingEnvironment.rda", to = file)
        file.remove("/tmp/WorkingEnvironment.rda")
        stopApp("COOCKIE CRISPR closed")
      })
    
    observeEvent(c(input$init2),
                 ignoreInit = TRUE, {
                   shinyjs::runjs("$('#exit_and_save')[0].click();")
                 })
    
    observeEvent(c(input$init),
                 ignoreInit = TRUE, {
      shinyjs::runjs("$('#state_save_sc')[0].click();")
    })
    
    output$state_save_sc <- downloadHandler(
      filename = function() {
        paste0("COOKIE_CRISPR_rState_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".rda")
      },
      content = function(file) {
        file.copy(from = "/tmp/WorkingEnvironment.rda", to = file)
        file.remove("/tmp/WorkingEnvironment.rda")
      }
    )
    
##########"## Restore state ############
observeEvent(input$restore,priority = 10,{
      withProgress(message = 'Loading analysis state', value = 0.5, {
       inFile <- input$restore
       if(!is.null(inFile)) {
         isolate({
           tmpEnv <- new.env()
           load(inFile$datapath, envir=tmpEnv)
            if (exists("r_reactives", envir=tmpEnv, inherits=FALSE)) {#
              print("load reactives")
              reactives$sampleplan <- tmpEnv$r_reactives$sampleplan
              reactives$counts <- tmpEnv$r_reactives$counts
            }
           incProgress(0.3)
            if (exists("r_separators", envir=tmpEnv, inherits=FALSE)){
              print("load separators")
              separators$counts <- tmpEnv$r_separators$counts
              separators$sampleplan <- tmpEnv$r_separators$sampleplan
            }
            incProgress(0.3)
            if (exists("r_inputs", envir=tmpEnv, inherits=FALSE)){
            print("load inputs")
              input <- tmpEnv$r_inputs
              lapply(names(input),
                    function(x) session$sendInputMessage(x, list(value = input[[x]]))
              )
             print(input)
         }
           rm(tmpEnv)
         }) #end of isolate
       }
       setProgress(1)
})
}) # end of observer
  
################### Report Section ################

# server report editor ---------------------------------------------------------
### yaml generation
rmd_yaml <- reactive({
  paste0("---",
         "\ntitle: '", input$report_title,
         "'\nauthor: '", input$report_author,
         "'\ndate: '", Sys.Date(),
         "'\noutput:\n  html_document:\n    toc: ", input$report_toc, "\n    number_sections: ", input$report_ns, "\n    theme: ", input$report_theme, "\n---\n\n",collapse = "\n")
})

### loading report template
# update aceEditor module
observe({
  # loading rmd report from disk
  inFile <- system.file("extdata", "reportTemplate.Rmd",package = "CRISPRApp")
  isolate({
    if(!is.null(inFile) && !is.na(inFile)) {
      
      rmdfilecontent <- paste0(readLines(inFile),collapse="\n")
      
      shinyAce::updateAceEditor(session, "acereport_rmd", value = rmdfilecontent)
    }
  })
})

### ace editor options
observe({
  autoComplete <- if(input$enableAutocomplete) {
    if(input$enableLiveCompletion) "live" else "enabled"
  } else {
    "disabled"
  }
  updateAceEditor(session, "acereport_rmd", autoComplete = autoComplete,theme=input$theme, mode=input$mode)
  # updateAceEditor(session, "plot", autoComplete = autoComplete)
})

#Enable/Disable R code completion
rmdOb <- aceAutocomplete("acereport_rmd")
observe({
  if(input$enableRCompletion) {
    rmdOb$resume()
  } else {
    rmdOb$suspend()
  }
})

## currently not working as I want with rmarkdown::render, but can leave it like this - the yaml will be taken in the final version only
output$knitDoc <- renderUI({
  input$updatepreview_button
  return(
    withProgress(
      message = "Updating the report in the app body",
      detail = "This can take some time",
      {
        # temporarily switch to the temp dir, in case you do not have write
        # permission to the current working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        tmp_content <- paste0(rmd_yaml(),input$acereport_rmd,collapse = "\n")
        incProgress(0.5, detail = "Rendering report...")
        htmlpreview <- knit2html(text = tmp_content, fragment.only = TRUE, quiet = TRUE)
        incProgress(1)
        isolate(HTML(htmlpreview))
      })
  )
})

## Download Report

output$saveRmd <- downloadHandler(
  filename = function() {
    if(input$rmd_dl_format == "rmd") {
      "report.Rmd"
    } else {
      "report.html"
    }
  },
  content = function(file) {
    # knit2html(text = input$rmd, fragment.only = TRUE, quiet = TRUE))
    tmp_content <- paste0(rmd_yaml(),input$acereport_rmd,collapse = "\n")
    # input$acereport_rmd
    if(input$rmd_dl_format == "rmd") {
      cat(tmp_content,file=file,sep="\n")
    } else {
      # write it somewhere too keeping the source
      # tmpfile <- tempfile()
      # file.create(tmpfile)
      # fileConn<- file(tempfile())
      # writeLines(tmp_content, fileConn)
      # close(fileConn)
      if(input$rmd_dl_format == "html") {
        
        # temporarily switch to the temp dir, in case you do not have write
        # permission to the current working directory
        # owd <- setwd(tempdir())
        # on.exit(setwd(owd))
        
        cat(tmp_content,file="/tmp/tempreport.Rmd",sep="\n")
        rmarkdown::render(input = "/tmp/tempreport.Rmd",
                          output_file = file,
                          "html_document",
                          # fragment.only = TRUE,
                          quiet = TRUE)
      }
    }
  })

} # end of Server

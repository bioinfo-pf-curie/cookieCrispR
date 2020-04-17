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
server_crispr_app <- function(input, output, session) {
    
  ######### Global options ##############
  session$onSessionEnded(stopApp)
  options(shiny.sanitize.errors = TRUE)
  session$allowReconnect(TRUE)

  ##### Usefull variables #############
  reactives <- reactiveValues(sampleplan = NULL,sampleplanGood = FALSE, sampleplanRaw = NULL, joined = NULL)
  separators <- reactiveValues(counts = NULL,sampleplan = NULL)
  
  ##### Upload files and datatable construction ####################
  ## counts table    
  observeEvent(c(input$counts,
                 separators$counts),{  
  #counts <- reactive({
      
      # req(input$counts)
      # inFile <- input$counts
      # counts <- read_delim(inFile$datapath, delim = input$FSc)
      # counts <- dplyr::rename(counts, sgRNA = .data$X1)
      # counts <- dplyr::select(counts, -.data$sequence)
      # 
      # annot_sgRNA <- dplyr::select(counts, .data$sgRNA, Gene = .data$gene)
      # 
      # counts <- gather(counts, value = "count", key = "barcode", -.data$sgRNA, -.data$gene)
      # counts <- dplyr::mutate(counts, barcode = str_remove(.data$barcode, ".R1.fastq"))
      # 
      # counts <- counts %>%
      #   dplyr::group_by(.data$barcode) %>%
      #   dplyr::mutate(cpm = 1e6 * .data$count / sum(.data$count), log_cpm = log10(1 + .data$cpm)) %>%
      #   dplyr::ungroup()
      # 
      # return(list(counts, annot_sgRNA))
      
      req(input$counts)
      inFile <- input$counts

      if (is.null(separators$counts)){
        counts <- read.table(inFile$datapath, sep = ";", header = TRUE)
      } else {
        counts <- read.table(inFile$datapath, sep = separators$counts, header = TRUE)
      }
      print(colnames(counts))
      counts <- dplyr::rename(counts, sgRNA = .data$X)
      counts <- dplyr::select(counts, -.data$sequence)
      
      
      reactives$counts <- counts
      
  })  # end of observer
    
    ## samples annotations
observeEvent(c(input$sample_plan,
              separators$sampleplan),{#input$restore

      # req(input$sample_plan)
      # inFile <- input$sample_plan
      # samples <- read_delim(inFile$datapath, delim = "|",  col_names = c("barcode", "sample","condition"))
      # samples <- samples %>%
      #   separate(sample, into = c("date","rep","clone","day"), remove = FALSE) %>%
      #   mutate(day = as.factor(.data$day)) %>%
      #   mutate(condition = as.factor(.data$condition)) %>%
      #   mutate(Timepoint_num = as.numeric(gsub("DAY","",.data$day))) %>% 
      #   mutate(day = fct_reorder(.f = factor(.data$day), .x = .data$Timepoint_num))
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
      samples <- reactives$sampleplanRaw %>%
      #   #separate(sample, into = c("date","rep","clone","day"), remove = FALSE) %>%
          mutate(Timepoint = as.factor(.data$Timepoint)) %>%
          mutate(Treatment = as.factor(.data$Treatment)) %>%
          mutate(Timepoint_num = as.numeric(gsub("[^0-9.-]", "", .data$Timepoint))) %>% 
          mutate(Timepoint = fct_reorder(.f = factor(.data$Timepoint), .x = .data$Timepoint_num))
      
      reactives$sampleplan <- samples
      }
    
}) # end of observer
    
    ## counts and annotations
    joined <- reactive({
      counts <- reactives$counts
      samples <- reactives$sampleplan
      counts <- full_join(counts, samples) %>%
        mutate(day = factor(.data$day, levels = input$timepoints_order))
      return(counts)
    })
    
    #timepoints <- reactiveValues(a = joined()$day)
    
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
      counts <- reactives$counts
      firstpoint <- input$timepoints_order[[1]]
      
      all <- joined() %>%
        select(.data$sgRNA, .data$clone, .data$rep, .data$day, .data$log_cpm, .data$gene, .data$condition)
      
      t0 <- all %>%  
        filter(.data$day == firstpoint)  %>%
        mutate(log_cpmt0 = .data$log_cpm) %>%
        mutate(log_cpm = NULL) %>%
        mutate(day = NULL) %>%
        mutate(gene = NULL) %>%
        mutate(condition = NULL)
      
      all <- inner_join(all,t0,c("sgRNA","clone","rep")) %>%
        mutate(diff = .data$log_cpm - .data$log_cpmt0) %>%
        mutate(log_cpm = NULL) %>%
        mutate(condtion = NULL)
      
      fin <- inner_join(joined() %>% 
                          mutate(gene=NULL) %>%
                          mutate(day=NULL) %>%
                          mutate(condition = NULL), all, by=c("sgRNA","clone","rep"))
      return(fin)
    })
    #
    
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
      counts <- joined()
      counts <- counts %>%
        group_by(.data$barcode, .data$rep, .data$clone, .data$day) %>%
        summarise(total = sum(.data$count)) %>% 
        as.data.frame() %>%
        ggplot(aes(x = .data$day, y = .data$total)) +
        geom_col(position = position_dodge()) + facet_wrap(vars(.data$rep, .data$clone), nrow = 1) +
        labs(title = "Number of reads by sample") #+
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
      counts <- joined()
      counts %>% 
        ggplot(aes(x = .data$day, y = .data$log_cpm, fill = .data$rep)) + geom_boxplot() + facet_grid(. ~ .data$clone) + 
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
      counts <- joined()
      
      counts <-  counts %>%
        ggplot(aes(x = .data$log_cpm, y = .data$day)) +
        geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
        facet_wrap(vars(.data$clone, .data$rep), ncol = 1, strip.position = "right") +
        labs(title = "Distribution of normalzed log-cpm by sample", subtitle = "All guides")
      return(plot(counts))
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
      counts <- joined()
      ess_genes <- ess_genes()
      counts <- filter(counts, .data$gene %in% ess_genes$V1)
      counts_plot <- ggplot(counts, aes(x = .data$log_cpm, y = .data$day)) + 
        geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
        facet_wrap(vars(.data$clone, .data$rep), ncol = 1, strip.position = "right") +
        scale_fill_viridis_c() +
        labs(title = "Distribution of normalised log-cpms", subtitle = "Essential genes")
      return(plot(counts_plot))
    })
    
    output$essential_distribs <- renderPlot({
      plot(essential_distribs())
    })
    
    nonessential_distribs <- reactive({
      counts <- joined()
      ess_genes <- non_ess_genes()
      counts <- counts %>%
        filter(.data$gene %in% ess_genes$V1)
      
      counts_plot <- ggplot(counts, aes(x = .data$log_cpm, y = .data$day)) +
        geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
        facet_wrap(vars(.data$clone, .data$rep), ncol = 1, strip.position = "right") +
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
    
    diff_box_all <- reactive({
      firstpoint <- input$timepoints_order[[1]]
      
      diff_box_all <- diff_t0() %>%
        filter(.data$day != firstpoint) %>%
        ggplot(aes(x = .data$day, y = .data$diff, fill = .data$rep)) + geom_boxplot() + facet_grid(.data$condition ~ .data$clone) +
        #ggplot(aes(x = day, y = diff_DAY1, fill = rep)) + geom_boxplot() + facet_grid(.~ clone) +
        ylab(paste0("diff_",firstpoint)) + 
        labs(title = paste0("Boxplots of log fold change from ", firstpoint ," - all guides"))
      
      plot(diff_box_all)
    })
    
    output$diff_box_all <- renderPlot({
      plot(diff_box_all())
    })
    
    
    diff_box_ess <- reactive({
      
      firstpoint <- input$timepoints_order[[1]]
      ess_genes <- ess_genes()
      
      diff_box_ess <-  diff_t0() %>% 
        #select(-condition, -log_cpmt0, -diff) %>%
        filter(.data$day != !!firstpoint) %>%
        filter(.data$gene %in% ess_genes[,1]) %>%
        ggplot(aes(x = .data$day, y= .data$diff, fill = .data$rep)) + geom_boxplot() + facet_grid(.~ .data$clone) +
        ylab(paste0("diff_",firstpoint)) + 
        labs(title = paste0("Boxplots of log fold change from ", firstpoint ," - essential genes's guides"))
      
      plot(diff_box_ess)
    })
    
    
    
    output$diff_box_ess <- renderPlot({
      plot(diff_box_ess())
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
    
    
    ROC <- reactive({
      
      ess_genes <- ess_genes()
      non_ess_genes <- non_ess_genes()
      
      d <- diff_t0() %>% select(.data$sgRNA, .data$clone, .data$rep, .data$day, .data$log_cpm, .data$gene, .data$condition, .data$diff) %>%
        group_by(.data$day, .data$condition, .data$clone, .data$rep) %>%
        arrange(.data$diff) %>% 
        mutate(type = case_when(
          .data$gene %in% ess_genes[,1] ~ "+",
          .data$gene %in% non_ess_genes[,1] ~ "-", 
          TRUE ~ NA_character_)) %>%
        filter(!is.na(.data$type)) %>%
        mutate(TP = cumsum(.data$type == "+") / sum(.data$type == "+"), FP = cumsum(.data$type == "-") / sum(.data$type == "-")) %>%
        ungroup()
      
      write.csv(diff_t0(),"~/to_shiny.csv")
      
      d <- d %>% ggplot(aes(x = .data$FP, y = .data$TP, color = .data$day)) + geom_abline(slope = 1, lty = 3) + geom_line() + facet_grid(.data$condition + .data$clone ~ .data$rep) + coord_equal()
      #   # 
      return(d)
    })
    
    
    
    output$roc <- renderPlot({
      ROC()
    })
    
    output$dlROC <- downloadHandler(
      filename = function(){
        paste("ROC_plots",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(ROC())
        dev.off()
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
        as.numeric(nrow(reactives$counts))
      )
    )
    
    output$Depth <- renderInfoBox(
      infoBox(
        "Average sequencing depth across samples :",
        round(sum(reactives$counts$count)/(nrow(reactives$sampleplan) + nrow(reactives$counts)))
        
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
<a href='https://gitlab.curie.fr/r-shiny/bioshiny/blob/devel/example_datasets/Crispr/SampleDescription' target='_blank'>Click here to download the Sample Description example file</a>
<br/>
<br/>
</li><li>Counts table file :</li>
<br/>
<B>Format description :</B>
<br/>
<ul>
<li>A csv formated text file using ; or , as field separator. </li>
<li>The first line of the file is a header, it contains samples barcodes as colnames and two supplementary columns called 'gene
' and 'sequence'. </li>
<li>AGuides names' are rownames. </li>
<li>AValues are read counts.</li>
</ul>
<br/>
<B>For example :</B><br/>
\"\";\"D115R01\";\"D115R02\";\"gene\";\"sequence\"<br/>
\"sgITGB8_1\";339;379;\"ITGB8\";\"AAAACACCCAGGCTGCCCAG\"<br/>
<br/>
<a href='https://gitlab.curie.fr/r-shiny/bioshiny/blob/devel/example_datasets/Crispr/global_counts_table.csv' target='_blank'>Click here to download the Counts table example file</a>
<br/>
<br/><br/>
</li><li>Genes' lists</li><br/>
<B>Format description :</B>
<br/>
Essentials and non essentials genes' list file just contains one column with all the genes of the list.<br/>
<br/>
<B>For example:</B><br/>

Gene1<br/>
Gene2<br/>
<br/>
<a href='https://gitlab.curie.fr/r-shiny/bioshiny/blob/devel/example_datasets/Crispr/essentials.csv' target='_blank'>Click here to download the essentials genes'list example file</a>

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
    
    ###### Save State ############
    
    observeEvent(c(input$exit_and_save,
                   input$init)
                 ,priority =10,ignoreInit = TRUE,{
      
      cat("save data \n")
      saveState(filename = "WorkingEnvironment.rda",
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
        file.copy(from = "WorkingEnvironment.rda", to = filename)
        file.remove("WorkingEnvironment.rda")
        stopApp("COOCKIE CRISPR closed")
      })
    
    observeEvent(input$init, {
      runjs("$('#state_save_sc')[0].click();")
    })
    
    output$state_save_sc <- downloadHandler(
      filename = function() {
        paste0("COOKIE_CRISPR_rState_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".rda")
      },
      content = function(file) {
        file.copy(from = "WorkingEnvironment.rda", to = file)
        file.remove("WorkingEnvironment.rda")
      }
    )
    
    ## Restore state 
    

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
  

}

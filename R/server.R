#' Server function
#'
#' Function containing all the server side of the Shiny App, must be called inside the run_app() function 
#'
#' @param input input
#' @param output output
#' @param session session
#'
#' @return None

server_crispr_app <- function(input, output, session) {
    
######### Global options ##############
session$onSessionEnded(stopApp)
options(shiny.sanitize.errors = TRUE,shiny.maxRequestSize = 3000*1024^2)
session$allowReconnect(TRUE)

##### Usefull variables #############
reactives <- reactiveValues(sampleplan = NULL,sampleplanGood = FALSE, sampleplanRaw = NULL,
                              joined = NULL,countsRaw = NULL, counts = NULL,
                            annot_sgRNA = NULL, norm_data =  NULL,norm_cpm =  NULL, guidelist = NULL,
                            genelist =  NULL,sample = NULL,checkcoherence = TRUE, normalize = TRUE,
                            checkcountscols = TRUE,interactive_boxplots = NULL,diff_t0 = NULL,
                            control_sgRNA = NULL)
geneslists <- reactiveValues(essential = NULL, non_essential = NULL)

############ Cicerone ############

guide0 <- Cicerone$
  new(id = "ciceroneGuide0")$
  step(el = "sample_plandiv",
       title ="Import your sample plan here",
       position = "top"#,
  )$
  step(
    el = "orderUIdiv",
    title = "Re-order timepoints",
    HTML("The app try to automaticaly figure out the correct order for your timepoints. If you see errors you can use this button to re-order them")
  )$
  step(el = 'essentialdiv',title = "Import your essentials' genes list here")$
  step(el = 'nonessentialdiv',title = "Import your non essentials' genes list here")$
  step(
    el = "screentype",
    title = 'Select the screen type corresponding to your data',
    HTML("The main difference is the normalization method used")
  )$
  step(el = 'countsdiv',title = "Import your counts table here")$
  step(el = "[data-value='Help']",
       title = "More information about input files formats is available here",
       is_id = FALSE)$
  #step(el = "Counts table",
  step(el = "countstablebox",
       title = "Here is an overview of your counts data")$
  step("[data-value='log10(cpm)']",title = "Here is an overview of your log2cpm(counts)",is_id = FALSE)$
  step("removegenesdiv",title = "You have got the possibility to remove genes for further analysis here")$
  step("sampleplanbox",title = "Here is an overview of your sampleplan")$  
  step("removesamplesdiv",title = "You have got the possibility to remove samples for further analysis here")$
  step("correlation_heatmap", title = "Sample to sample correlation heatmap")$
  step("correlationsAnnotdiv", title = "Use this button to annotate the correlation heatmap")$
  step("dlcorrelation_heatmap",title = "Click here to download the heatmap")$
  step("dlcorrelation_coefficients",title = "Click here to download correlation coeffcients")



  


 
    # HTML("Two types of comparison are available in the app </br><li><b>Intra Treatment comparison :</b> 
    #   For the samples with the same Treatment value, will compute differential analysis for each timepoint against the initial timepoint </li>
    #   <li><b>Inter Treatment comparison :</b> Set up a comparison using <i>Treatment</i> and <i>Supplementary info</i> values. Then run the differential analysis. You must have samples with equivalent timeponts available in your data</li>")
  #)#$

observeEvent(input$startCicerone0, {
  withProgress(message = 'Loading example dataset...', value = 0.5, {
  metadata_path <- system.file("extdata", "SampleDescriptiondatatest.txt", package = "CRISPRApp")
  reactives$sampleplan <- read.table(metadata_path, sep = ";",header = TRUE) %>%
    mutate(Timepoint = as.factor(.data$Timepoint)) %>%
    mutate(Treatment = as.factor(.data$Treatment)) %>%
    mutate(Timepoint_num = as.numeric(gsub("[^0-9.-]", "", .data$Timepoint))) %>% 
    mutate(Timepoint = fct_reorder(.f = factor(.data$Timepoint), .x = .data$Timepoint_num))
  reactives$samples <- as.character(reactives$sampleplan$Sample_ID)
  
  counts_path <- system.file("extdata", "global_counts_table_datatest.csv", package = "CRISPRApp")
  counts <- as.data.frame(fread(counts_path, sep = ",", header = TRUE,fill = TRUE))
  counts <- counts %>% rename(X = V1)
  counts <- dplyr::rename(counts, sgRNA = .data$X)
  counts <- dplyr::select(counts, -.data$sequence)
  reactives$guidelist <- as.character(unique(counts$sgRNA))
  reactives$genelist <- as.character(unique(counts$gene))
  reactives$countsRaw <- counts %>% select(sgRNA,gene, everything())
  
  essential_path <- system.file("extdata", "essentials.csv", package = "CRISPRApp")
  ess_genes <- read.table(essential_path, header=FALSE)
  ess_genes <- ess_genes %>% rename(X = V1)
  geneslists$essential <- ess_genes
  non_essential_path <- system.file("extdata", "non_essentials_datatest.csv", package = "CRISPRApp")
  non_ess_genes <- read.table(non_essential_path, header=FALSE)
  non_ess_genes <- non_ess_genes %>% rename(X = V1)
  geneslists$non_essential <- non_ess_genes
  
  print("guideinit")
  guide0$init()$start()
  })
})

##################### END OF CICERONE ##########################

##### Upload files and datatable construction ####################
## counts table  
precheck <- reactiveValues(counts = FALSE,sampleplan = FALSE,Fs = NULL)
observeEvent(c(input$counts,
                 input$Fsc),{
      print("checking file separator in counts matrix ...")
      req(input$counts)
      inFile <- input$counts
      if(grepl(".xls",inFile$name) == FALSE){
      #req(input$Fsc)
      semicolon <- FALSE
      comma <- FALSE
      tabssplan <- FALSE
      #counts <- read.table(inFile$datapath, sep = "d", header = TRUE,fill = TRUE)
      counts <- as.data.frame(fread(inFile$datapath, sep = "d", header = TRUE,fill = TRUE))
      for(col in 1:ncol(counts)){
        if (TRUE %in% grepl(";",counts[,col])){
          semicolon <- TRUE
          precheck$Fs <- ";" 
        }
        if (TRUE %in% grepl(",",counts[,col])){
          comma <- TRUE
          precheck$Fs <- "," 
        }
        if (TRUE %in% grepl("\t",counts[,col])){
          tabssplan <- TRUE
          precheck$Fs <- "\t"
        }
      }
    if(semicolon ==  TRUE & comma == TRUE){
        showModal(modalDialog(p(""),
                              title = h4(HTML("<b>Both semicolon and comma</b> detected in your counts matrix file"),style="color:red"),
                              h6('Use a unique field separator'),
                              footer = tagList(
                                modalButton("Got it"))
        ))
      } else if(semicolon ==  TRUE & tabssplan == TRUE){
        showModal(modalDialog(p(""),
                              title = h4(HTML("<b>Both semicolon and tabulation</b> detected in your counts matrix file"),style="color:red"),
                              h6('Use a unique field separator'),
                              footer = tagList(
                                modalButton("Got it"))
        ))
      } else if(comma ==  TRUE & tabssplan == TRUE){
          showModal(modalDialog(p(""),
                                title = h4(HTML("<b>Both comma and tabulation</b> detected in your counts matrix file"),style="color:red"),
                                h6('Use a unique field separator'),
                                footer = tagList(
                                  modalButton("Got it"))
          ))
      } else {
        precheck$counts <- TRUE
      }
      }else{precheck$counts <- TRUE}
})
  
observeEvent(precheck$counts,{
  req(input$counts)
  inFile <- input$counts
  if(precheck$counts == TRUE){
    if(grepl(".xls",inFile$name) == FALSE){
        print("reading file...")
        #counts <- read.table(inFile$datapath, sep = precheck$Fs, header = TRUE)
        counts <- as.data.frame(fread(inFile$datapath, sep = precheck$Fs, header = TRUE,fill = TRUE))
        print(precheck$Fs)
      if (!("V1" %in% colnames(counts))){
        if (paste0("V",ncol(counts)) %in% colnames(counts)) {
          colnames(counts) <- c("X",colnames(counts)[1:(ncol(counts)-1)])
        } else {
          counts <- tibble::rownames_to_column(counts,"X")
        }
      } else {
        #counts <- tibble::column_to_rownames(counts,"V1")
        counts <- counts %>% rename(X = V1)
     }
    } else {
      counts <- readxl::read_excel(inFile$datapath)
      
      if ("...1" %in% colnames(counts)){
        colnames(counts)[1] <- "X"
        print("...1")
      } else{colnames(counts) <- c("X",colnames(counts)[-length(colnames(counts))])}
    }
    counts <- dplyr::rename(counts, sgRNA = .data$X)
    counts <- dplyr::select(counts, -.data$sequence)
    reactives$guidelist <- as.character(unique(counts$sgRNA))
    reactives$genelist <- as.character(unique(counts$gene))
    reactives$countsRaw <- counts %>% select(sgRNA,gene, everything())
    precheck$counts <- FALSE
  }
})

observeEvent(reactives$genelist,priority = -1,{
   req(reactives$genelist)
   withProgress(message = 'Updating genes filter...', value = 0.5, {
   updatePickerInput(session, "removegenes",choices = as.character(reactives$genelist))
   setProgress(1)
   })
})


observeEvent(c(input$removeguides,input$removegenes,input$removesamples,reactives$countsRaw),ignoreInit = TRUE,{
  withProgress(message = 'Applying data filters ...', value = 0.5, {
  reactives$normalize <- TRUE
  if (!is.null(input$removeguides) | !is.null(input$removegenes) | !is.null(input$removesamples)){
  tic("data filters observer")
  reactives$checkcountscols <- FALSE
  if (!is.null(input$removeguides) | !is.null(input$removegenes)){
  incProgress(0.3)
  selectedcountsRaw <-   reactives$countsRaw %>%
    #filter(!(sgRNA %in% input$removeguides)) %>%
    filter(!(gene %in% input$removegenes))
  } else {
    selectedcountsRaw  <- reactives$countsRaw
  }
  if(!is.null(input$removesamples)){
    incProgress(0.3)
    selectedcountsRaw <- selectedcountsRaw[,!(colnames(selectedcountsRaw) %in% input$removesamples)]
  }
  toc(log = TRUE)
  } else {
    tic("data filters observer")
    selectedcountsRaw <- reactives$countsRaw
    toc(log = TRUE)
  }
    reactives$selectedcountsRaw <- selectedcountsRaw
    read_number$toplot <- TRUE
    dists$toplot <- TRUE
    differencetoT0$toplot <- TRUE
    ROC$toplot <- TRUE
  setProgress(1)
  })
})

observeEvent(c(input$sample_plan),priority = 10,{
                 print("checking file separator in sample plan file ...")
                 req(input$sample_plan)
                 inFile <- input$sample_plan
                 if(grepl(".xls",inFile$name) == TRUE){
                  samples <- readxl::read_excel(inFile$datapath)
                  reactives$sampleplanRaw <- samples
                  reactives$samples <- samples$Sample_ID
                 } else {
                 semicolonssplan <- FALSE
                 commassplan <- FALSE
                 tabssplan <- FALSE
                 #ssplan <- read.table(inFile$datapath, sep = "d", header = TRUE, fill = TRUE)
                 ssplan <- as.data.frame(fread(file = inFile$datapath, sep = "d", header = TRUE, fill = TRUE,nThread=2))
                 for(col in colnames(ssplan)){
                   
                   if (TRUE %in% grepl(";",ssplan[,col])){
                     semicolonssplan <- TRUE
                   }
                   if (TRUE %in% grepl(",",ssplan[,col])){
                     commassplan <- TRUE
                   }
                   if (TRUE %in% grepl("\t",ssplan[,col])){
                     tabssplan <- TRUE
                   }
                 }
                 if(semicolonssplan ==  TRUE && commassplan == FALSE){
                  #samples <- read.table(inFile$datapath, sep = ";", header = TRUE, fill = TRUE)
                  samples <- as.data.frame(fread(inFile$datapath, sep = ";", header = TRUE, fill = TRUE))
                  samples <- samples[which(!(samples[,1] %in% c(""," ",NA))),]
                  reactives$sampleplanRaw <- samples
                  reactives$samples <- as.character(samples$Sample_ID)
                 } else if (commassplan == TRUE && semicolonssplan ==  FALSE ){
                   samples <- as.data.frame(fread(inFile$datapath, sep = ",", header = TRUE, fill = TRUE))
                   samples <- samples[which(!(samples[,1] %in% c(""," ",NA))),]
                   #samples <- read.table(inFile$datapath, sep = ",", header = TRUE, fill = TRUE)
                   reactives$sampleplanRaw <- samples
                   reactives$samples <- as.character(samples$Sample_ID)
                 } else if (commassplan == FALSE && semicolonssplan ==  FALSE && tabssplan == TRUE ){
                   #samples <- read.table(inFile$datapath, sep = "\t", header = TRUE, fill = TRUE)
                   samples <- as.data.frame(fread(inFile$datapath, sep = "\t", header = TRUE, fill = TRUE))
                   samples <- samples[which(!(samples[,1] %in% c(""," ",NA))),]
                   reactives$sampleplanRaw <- samples
                   reactives$samples <- as.character(samples$Sample_ID)
                 }  else if(semicolonssplan ==  TRUE && commassplan == TRUE){
                     showModal(modalDialog(p(""),
                                           title = h4(HTML("Both <b>Comma</b>  ans <b>semicolon </b> detected in your sample plan file"),style="color:red"),
                                           h6("Check for file separators homogeneity"),
                                           footer = tagList(
                                             modalButton("Got it"))
                     ))
                 }
                 }
      print("end of sampleplan check sep")
})

observeEvent(c(reactives$sampleplanRaw,input$removesamples),{   
  
      if(reactives$sampleplanGood == TRUE){
      print("Arranging sample plan")
      samples <- reactives$sampleplanRaw
      samples$Sample_ID <- as.character(samples$Sample_ID)
      samples$Treatment <- as.character(samples$Treatment)
      samples$Cell_line <- as.character(samples$Cell_line)
      samples$Replicate <- as.character(samples$Replicate)
      samples$SupplementaryInfo <- as.character(samples$SupplementaryInfo)
      samples$Timepoint <- as.character(samples$Timepoint)
      
      samples <- samples  %>% replace_with_na_all(~.x == "<NA>")
      
      samples <- samples %>% filter(!(Sample_ID %in% input$removesamples)) %>% 
        replace_na(list(Sample_ID = "unknown", Cell_line = "unknown", Timepoint= "unknown", Treatment= "unknown", Replicate= "unknown", SupplementaryInfo = "unknown"))
      
      samples$Timepoint <- as.factor(samples$Timepoint)
      
      samples <- samples %>% 
          mutate(Timepoint = as.factor(.data$Timepoint)) %>%
          mutate(Treatment = as.factor(.data$Treatment)) %>%
          mutate(Timepoint_num = as.numeric(gsub("[^0-9.-]", "", .data$Timepoint))) %>% 
          mutate(Timepoint = fct_reorder(.f = factor(.data$Timepoint), .x = .data$Timepoint_num))
      print("DONE")
      read_number$toplot <- TRUE
      differencetoT0$toplot <- TRUE
      dists$toplot <- TRUE
      ROC$toplot <- TRUE
      reactives$sampleplan <- samples
      }
    
}) # end of observer
    
## Join counts and annotations
observeEvent(c(reactives$selectedcountsRaw,reactives$sampleplan,
               input$timepoints_order),{
      if(!is.null(reactives$selectedcountsRaw) & !is.null(reactives$sampleplan)){
        samples <- reactives$sampleplan
        counts <- reactives$selectedcountsRaw
        
        tic("tic calculating cpm")
        withProgress(message = 'Calculating cpm', value = 0.5,{
          incProgress(0.3,detail = "Formating data...")
        reactives$annot_sgRNA <- dplyr::select(counts, .data$sgRNA, Gene = .data$gene)
        counts <- gather(counts, value = "count", key = "Sample_ID", -.data$sgRNA, -.data$gene)
        counts <- dplyr::mutate(counts, Sample_ID = gsub(".R[1-9].fastq","",.data$Sample_ID))
        sumcols <- sum(counts$count)
        incProgress(0.3,detail = "Computing cpm...")
        counts <- counts %>%
          dplyr::group_by(.data$Sample_ID) %>%
          dplyr::mutate(cpm = 1e6 * .data$count / as.numeric(sumcols), log_cpm = log10(1 + cpm),log2_cpm = log2(1+cpm))  %>%
          dplyr::ungroup()
        
        print(head(counts))
        #toc(log=TRUE)
        })
      
      if(reactives$checkcoherence ==  TRUE){
      if (TRUE %in% unique(!(unique(counts$Sample_ID) %in% samples$Sample_ID))){
      if(input$sidebarmenu == "DataInput" && reactives$checkcountscols == TRUE){
       showModal(modalDialog(p(""),
                    title = "Missing samples in provided sampleplan",
                    tagList(h6('The folowing samples :'),
                            h6(paste0(unique(counts[!(counts$Sample_ID %in% samples$Sample_ID),]$Sample_ID),collapse = " | "),style="color:red"),
                            h6("are missing in the sampleplan."),
                            p(),
                            h6("Check that Samples IDs are correctly matching between ssplan and count matrix files, unless they will be removed for the following analysis")),
                    footer = tagList(
                      modalButton("Got it")
                    ))
       )
      }
      counts <- counts  %>%
          filter(Sample_ID %in% samples$Sample_ID)
      } # end of check unique
      } # end of check coherence
      reactives$counts <- counts
      } # end of if is null

})

# #########################################################################################################
# ########################################## Sample to sample correlations ################################"
# #########################################################################################################
# 

observeEvent(reactives$sampleplan,{
  req(reactives$sampleplan)
  cn <- colnames(reactives$sampleplan %>% select(-Sample_ID,-Timepoint_num))
  updatePickerInput(session=session,"correlationsAnnot",choices = cn)
})

correlation_reactives <- reactiveValues(coefficients = NULL, plot = NULL)
observeEvent(c(reactives$counts,input$correlationsAnnot),{

  req(reactives$counts)
  req(reactives$samples)
  counts <- reactives$counts

    counts <- reactives$counts %>%
      dplyr::select(sgRNA, Sample_ID,log2_cpm) %>%
      spread(key = Sample_ID, value = log2_cpm) %>%
      as.data.frame()

    counts <- counts %>%
      column_to_rownames("sgRNA")
    
    samples <- reactives$sampleplan
    #save(counts,samples,file = "~/sampletosamplecorrelations.rda")

    spearman_corr <- cor(x=counts, method = 'spearman')
    dd <- as.dist((1-spearman_corr)/2)
    hc <- hclust(dd)
    spearman_corr <- spearman_corr[hc$order,hc$order]
    
    spearman_corr <- signif(spearman_corr, digits = 3)
    correlation_reactives$coefficients <- spearman_corr

    spearman_corr <- as.matrix(spearman_corr)
    cn <- colnames(spearman_corr)
    rn <- rownames(spearman_corr)
    samples <- samples %>% column_to_rownames("Sample_ID")
    
    top_annotation <- NULL
    if(length(input$correlationsAnnot) > 0){
      samples <- samples[cn,input$correlationsAnnot]
      if (length(input$correlationsAnnot) >= 2 ){
        print("b")
        top_annotation <- HeatmapAnnotation(
          df = samples,
          which = "col",
          show_legend = TRUE)
      } else {
        print("ddede")
        top_annotation <- HeatmapAnnotation(
          Annot_name = samples,
          which = "col",
          show_legend = TRUE)
        names(top_annotation) <- input$correlationsAnnot
        
      }
    }
    
    print('plot')
    print(top_annotation)
    save(spearman_corr,samples,file = "~/sampletosamplecorrelations.rda")
    
    correlation_reactives$plot  <- Heatmap(spearman_corr,
            name = "coefficient",
            column_title = "spearman correlation on log2(cpm +1) counts",
            cluster_rows = FALSE,cluster_columns = FALSE,
            top_annotation = top_annotation
            )
})

output$correlation_heatmap <- renderPlot({
  req(correlation_reactives$plot)
  correlation_reactives$plot
})

output$dlcorrelation_heatmap <- downloadHandler(
  filename = function(){
    paste("Correlations_heatmap_",Sys.Date(),".pdf",sep="")
  },
  content = function(file){
    pdf(file = file)
    #plot(correlation_reactives$plot)
    draw(correlation_reactives$plot)
    dev.off()
  }
)

output$dlcorrelation_coefficients <- downloadHandler(
  filename = function(){
    paste("Correlation_coefficients_",Sys.Date(),".csv",sep="")
  },
  content = function(file){
    req(correlation_reactives$coefficients)
    write.table(x = correlation_reactives$coefficients,file = file,sep = ",", quote = FALSE, row.names = FALSE)
  }
)

ctrlterm <- reactiveValues(term = NULL)
observeEvent(reactives$annot_sgRNA,{
  if(input$screentype == "negative"){
  req(reactives$annot_sgRNA)
  annot_sgRNA <- reactives$annot_sgRNA

  control_sgRNA <- filter(annot_sgRNA, str_detect(Gene,paste(c("Non-Targeting","negative_control"),collapse = "|")))
  
  if(nrow(control_sgRNA) == 0){
    showModal(modalDialog(
      title = "No control guides found in data",
      HTML("By default the app search guides annoted as Non-Targeting or negative_control </br></br>
            Maybe controls in your data are annoted in a different manner </br></br>"),
      textInput("customctrl",label = "Enter here the character string defining the control guides in your data (case sensitive)",width = "100%"),br(),
      HTML("YOU CAN ENTER A RANDOM CHARACTER STRING IF YOU PLAN TO PROCESS A POSITIVE SRCEENING </br></br>"),
      footer = tagList(
        actionButton("useit","Use this term to search for control guides",width = "100%")
      ))
    )
    observeEvent(input$useit,{
      req(input$customctrl)
      req(input$useit)
      control_sgRNA <- filter(annot_sgRNA,str_detect(Gene,
                                                     paste(c("Non-Targeting","negative_control",input$customctrl),collapse = "|")))
      ctrlterm$term <- input$customctrl
      if(nrow(control_sgRNA) == 0){
        control_sgRNA <- NULL
      }
      reactives$control_sgRNA <- control_sgRNA
      removeModal()
    })
  }
  if(nrow(control_sgRNA) == 0){
    control_sgRNA <- NULL
  }
  reactives$control_sgRNA <- control_sgRNA
  } # end of screentype check
})

observeEvent(c(reactives$counts,reactives$sampleplan,input$sidebarmenu,input$countstabset,reactives$control_sgRNA),priority = 10,{
      samples <- reactives$sampleplan
      if(!is.null(reactives$sampleplan) & !is.null(reactives$counts) & input$sidebarmenu != "DataInput" & reactives$normalize == TRUE | input$countstabset=="log10(cpm)"  & reactives$normalize == TRUE){

      reactives$normalize <- "DONE"
      withProgress(message = 'Data normalization edgeR', value = 0.5, {
      annot_sgRNA <- reactives$annot_sgRNA
      
      counts <- reactives$counts %>% 
           dplyr::select(sgRNA, Sample_ID,count) %>%
           spread(key = Sample_ID, value = count) %>%
           as.data.frame() 

      counts <- counts %>%
           column_to_rownames("sgRNA")

      if(input$screentype == "negative"){
      norm_data <- sg_norm(counts,
                             sample_annot = column_to_rownames(samples, "Sample_ID")[colnames(counts), ],
                             sgRNA_annot = annot_sgRNA, control_sgRNA = reactives$control_sgRNA$sgRNA)
      } else {
        norm_data <- sg_norm(counts,
                             sample_annot = column_to_rownames(samples, "Sample_ID")[colnames(counts), ],
                             sgRNA_annot = annot_sgRNA, control_sgRNA = annot_sgRNA$sgRNA)
      }
      reactives$norm_data <- norm_data
      
      incProgress(0.3)

      norm_cpm <- cpm(norm_data$counts, log = TRUE)
      reactives$norm_cpm <- norm_cpm
      setProgress(1)
    })
      
    withProgress(message = 'Data transformation for distributions', value = 0.5, {
        req(reactives$sampleplan)
        #req(reactives$selectedcountsRaw)
        req(reactives$counts)
        samples <- reactives$sampleplan
        counts <- reactives$counts
        
        if(!(is.null(input$timepoints_order))){
        counts <- full_join(counts, samples) %>%
          mutate(Timepoint = factor(.data$Timepoint, levels = input$timepoints_order)) %>%  mutate(Gene = gene)
        } else {
        counts <- full_join(counts, samples) %>%
            mutate(Timepoint = factor(.data$Timepoint, level = unique(.data$Timepoint))) %>%  mutate(Gene = gene)
        }
      reactives$joined <- counts %>%  mutate(Gene = gene)
      joined <- reactives$joined 
      save(joined, file = "~/crispr_joined.rda")
      print("joining done")
    setProgress(1)
    })
      #setProgress(1)
      #}) # end of progress
      }# end of if
    #} # end of else
}) # End of observer    

     observe({ 
      req(input$essential)
      inFile <- input$essential
      ess_genes <- read.table(inFile$datapath, header=FALSE)
      if ("V1" %in% colnames(ess_genes)){
        ess_genes <- ess_genes %>% rename(X = V1)
      }
      geneslists$essential <- ess_genes
    })
    
    observe({ 
      req(input$nonessential)
      inFile <- input$nonessential
      non_ess <- read.table(inFile$datapath, header = FALSE)
      if ("V1" %in% colnames(non_ess)){
        non_ess <- non_ess %>% rename(X = V1)
      }
      return(non_ess)
      geneslists$non_essential <- non_ess
    })

#########################################################################################################
#################################### Negative screening #################################################
#########################################################################################################
differencetoT0 <- reactiveValues(toplot = TRUE)
observeEvent(input$sidebarmenu,{    
      if (input$sidebarmenu ==  "Tev"){
        if(differencetoT0$toplot == TRUE){
      req(input$timepoints_order)
      req(reactives$joined)
      withProgress(message = 'Difference to initial timepoint calculation', value = 0.5, {
      counts <- reactives$joined
      firstpoint <- input$timepoints_order[[1]]
      counts$Timepoint <- relevel(as.factor(counts$Timepoint), ref = firstpoint)
      
      fin <- counts %>%
        group_by(sgRNA, Cell_line, Replicate) %>%
        arrange(Timepoint) %>%
        mutate(diff = log_cpm - first(log_cpm)) %>% 
        ungroup()
     
      reactives$diff_t0 <- fin
      })
      differencetoT0$toplot <- FALSE
      }
      }
})

######### Plots and tables outputs ####################################
    output$counts_table <- DT::renderDataTable({
      print("selectedCountsRaw table observer")
      DT::datatable(reactives$selectedcountsRaw, rownames = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
    output$normalized_counts_table <- DT::renderDataTable({
      DT::datatable(reactives$joined %>% select(sgRNA,Gene,Sample_ID,log_cpm) %>% 
                      tidyr::spread(key=Sample_ID,value=log_cpm),
                    rownames = FALSE,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    }) # end of datatable

observeEvent(reactives$samples,{
req(reactives$samples)
withProgress(message = 'Updating samples filter...', value = 0.5, {
updatePickerInput(session=session,"removesamples",choices = reactives$samples)
setProgress(1)
})
})
    
observe({
    output$sample_plan_table <- DT::renderDataTable({
      if (!is.null(reactives$sampleplan)){
      print("sampleplan table observer")
      sample_plan <- reactives$sampleplan
      sample_plan <- dplyr::select(sample_plan, -.data$Timepoint_num)
      return(DT::datatable(sample_plan,rownames = FALSE))
      }
    })
}) # end of observer
       
    read_number <- reactiveValues(plot = NULL,toplot = TRUE)
    observeEvent(input$sidebarmenu,{    
        if (input$sidebarmenu ==  "Rawdist"){
          if(read_number$toplot == TRUE){
      req(reactives$counts)
      req(reactives$sampleplan)
      print("read number observer")
      withProgress(message = 'updating data', value = 0.5, {
        incProgress(0.3)
      if(!(is.null(input$timepoints_order))){
        counts <- full_join(reactives$counts, reactives$sampleplan) %>%
          mutate(Timepoint = factor(.data$Timepoint, levels = input$timepoints_order))
      } else {
        print(reactives$sampleplan$Timepoint)
        counts <- full_join(reactives$counts, reactives$sampleplan) %>%
          mutate(Timepoint = factor(.data$Timepoint, levels = unique(.data$Timepoint)))
      }
        incProgress(0.3)
      counts <- as.data.frame(counts) %>%
        group_by(.data$Sample_ID, .data$Replicate, .data$Cell_line, .data$Timepoint) %>%
        summarise(total = sum(.data$count)) %>% 
        as.data.frame() %>%
        ggplot(aes(x = .data$Timepoint, y = .data$total)) +
        geom_col(position = position_dodge()) + facet_wrap(vars(.data$Replicate, .data$Cell_line), nrow = 1) +
        labs(title = "Number of reads per sample", xlab ="Timepoint", ylab ="Total counts")
      read_number$plot <- counts
      setProgress(1)
      print("end of progress updating data")
      })
      read_number$toplot <- FALSE
    } # enf of toplot if
    } # end of sidebarif
  })
    
    output$read_number <- renderPlot({
      plot(read_number$plot)
    })
    
    output$dlreadnumber <- downloadHandler(
      filename = function(){
        paste("Read_number_",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(read_number$plot)
        dev.off()
      }
    )
    
    dists <- reactiveValues(boxall = NULL,boxess = NULL, boxnoness =  NULL,density = NULL,toplot = TRUE)
    observeEvent(c(input$sidebarmenu,reactives$joined),{
      req(reactives$joined)
      if (input$sidebarmenu ==  "Rawdist"){
        if(dists$toplot == TRUE){
      print("boxplots all")
      withProgress(message = 'Rendering boxplots all', value = 0.5,{

      if(length(unique(reactives$joined$Cell_line)) >= 2){
      plot <- reactives$joined %>%
        ggplot(aes(x = .data$Timepoint, y = .data$log_cpm, fill = .data$Replicate)) + geom_boxplot() + facet_grid(. ~ .data$Cell_line) +
        labs(title = "Distribution of normalized log-cpm by sample", subtitle = "All guides")
      } else{
      plot <- reactives$joined %>%
          ggplot(aes(x = .data$Timepoint, y = .data$log_cpm, fill = .data$Replicate)) + geom_boxplot() +
          labs(title = "Distribution of normalized log-cpm by sample", subtitle = "All guides")
      }
      dists$boxall <- plot
      #if(!is.null(ess_genes())){
      if(!is.null(geneslists$essential)){
      if(length(unique(reactives$joined$Cell_line)) >= 2){
        plot <- reactives$joined %>% 
          #filter(gene %in% as.character(ess_genes()$X)) %>%
          filter(gene %in% as.character(geneslists$essential$X)) %>%
          ggplot(aes(x = .data$Timepoint, y = .data$log_cpm, fill = .data$Replicate)) + geom_boxplot() + facet_grid(. ~ .data$Cell_line) +
          labs(title = "Distribution of normalized log-cpm by sample", subtitle = "essential guides")
      } else{
        plot <- reactives$joined %>% 
          filter(gene %in% as.character(geneslists$essential$X)) %>%
          #filter(gene %in% as.character(ess_genes()$X)) %>%
          ggplot(aes(x = .data$Timepoint, y = .data$log_cpm, fill = .data$Replicate)) + geom_boxplot() +
          labs(title = "Distribution of normalized log-cpm by sample", subtitle = "essential guides")
      }
      dists$boxess <- plot
      }
      #if(!is.null(non_ess_genes())){
      if(!is.null(geneslists$non_essential)){
      if(length(unique(reactives$joined$Cell_line)) >= 2){
        plot <- reactives$joined %>%
          #filter(gene %in% c("Non-Targeting","negative_control",as.character(non_ess_genes()$X))) %>%
          filter(gene %in% c("Non-Targeting","negative_control",as.character(geneslists$non_essential$X))) %>%
          ggplot(aes(x = .data$Timepoint, y = .data$log_cpm, fill = .data$Replicate)) + geom_boxplot() + facet_grid(. ~ .data$Cell_line) +
          labs(title = "Distribution of normalized log-cpm by sample", subtitle = "Non essential and control guides")
      } else{
        plot <- reactives$joined %>%
          filter(gene %in% c("Non-Targeting","negative_control",as.character(geneslists$non_essential$X))) %>%
          #filter(gene %in% c("Non-Targeting","negative_control",as.character(non_ess_genes()$X))) %>%
          ggplot(aes(x = .data$Timepoint, y = .data$log_cpm, fill = .data$Replicate)) + geom_boxplot() +
          labs(title = "Distribution of normalized log-cpm by sample", subtitle = "Non essential and control guides")
      }
      dists$boxnoness <- plot
      }
      setProgress(1)
      })
      dists$toplot <- FALSE
        }
      }
    })

    output$boxplot_all <- renderPlot({
      plot(dists$boxall)
    })
    
    output$dlbox_all <- downloadHandler(
      filename = function(){
        paste("Boxplots_all_guides",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(dists$boxall)
        plot(dists$boxess)
        plot(dists$boxnoness)
        dev.off()
      }
    )
    
    output$boxplot_ess <- renderPlot({
      plot(dists$boxess)
    })
    
    output$boxplot_noness <- renderPlot({
      plot(dists$boxnoness)
    })
    
    essential_distribs <- reactive({
      req(reactives$joined)
      req(geneslists$essential)
      counts <- reactives$joined
      ess_genes <- geneslists$essential
      #ess_genes <- ess_genes()
      withProgress(message = 'Calculating density ridges', value = 0.5, {
        incProgress(0.3)
  
      counts <- filter(counts, .data$gene %in% ess_genes$X)
      
      subtitle <- "Essential genes"
      
      counts_plot <- ggplot(counts, aes(x = .data$log_cpm, y = .data$Timepoint)) + 
        geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
        facet_wrap(vars(.data$Cell_line, .data$Replicate), ncol = 1, strip.position = "right") +
        scale_fill_viridis_c() +
        labs(title = "Distribution of normalised log-cpms", subtitle = subtitle)
      return(plot(counts_plot))
      })
    })
    
    output$essential_distribs <- renderPlot({
      plot(essential_distribs())
    })
    
    nonessential_distribs <- reactive({

      req(reactives$joined)
      #req(non_ess_genes())
      req(geneslists$non_essential)
      counts <- reactives$joined
      non_ess_genes <- geneslists$non_essential
      #non_ess_genes <- non_ess_genes()
      print("non ess and control selection")
      
      counts <- counts %>%
        filter(.data$gene %in% c(as.character(non_ess_genes$X),"Non-Targeting","negative_control"))
      
      subtitle <- 'Not enssential and control genes'
      
      if(length(counts$Cell_line) > 1 && length(counts$Replicate) > 1){
      counts_plot <- ggplot(counts, aes(x = .data$log_cpm, y = .data$Timepoint)) +
        geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
        facet_wrap(vars(.data$Cell_line, .data$Replicate), ncol = 1, strip.position = "right") +
        scale_fill_viridis_c() +
        labs(title = "Distribution of log-cpms", subtitle = subtitle)
      } else {
        counts_plot <- ggplot(counts, aes(x = .data$log_cpm, y = .data$Timepoint)) +
          geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
          facet_wrap(vars(.data$Cell_line, .data$Replicate), ncol = 1, strip.position = "right") +
          scale_fill_viridis_c() +
          labs(title = "Distribution of log-cpms", subtitle = subtitle)
      }
      
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
      timepoints <- levels(reactives$sampleplan$Timepoint)
      print(timepoints)
      orderInput(inputId = 'timepoints', label = 'Re-order your timepoints here :', items = timepoints, as_source = F)
    })
    
    ############## DIFF BOXES TEV CONTROL ################
    diff_boxes <- reactiveValues(diff_box_ess = NULL,diff_box_all = NULL)

    observeEvent(reactives$diff_t0,{
    req(input$timepoints_order[[1]])
    req(reactives$diff_t0)
    firstpoint <- input$timepoints_order[[1]]
    #if(is.null(input$essential) | is.null(input$nonessential)){
    if(is.null(geneslists$essential) | is.null(geneslists$non_essential)){
      showModal(
        modalDialog(tagList(h3("You must provide essentials and non essentials genes list to perform positive screening")),
                    footer = tagList(
                      modalButton("Got it"))
        )
      )
    } else {
    #non_ess_genes <- non_ess_genes()
    non_ess_genes <- geneslists$non_essential
    diff_t0 <- reactives$diff_t0 %>%
        filter(.data$Timepoint != firstpoint)
      
      diff_t0 <- diff_t0 %>%
        filter(.data$gene %in% c(as.character(non_ess_genes$X),"Non-Targeting","negative_control"))

      diff_boxes$diff_box_all <- diff_t0 %>%
        filter(.data$Timepoint != firstpoint) %>%
        filter(.data$gene %in% c(as.character(non_ess_genes$X),"Non-Targeting","negative_control")) %>%
        ggplot(aes(x = .data$Timepoint, y = .data$diff, fill = .data$Replicate)) + geom_boxplot() + facet_grid(.data$Treatment ~ .data$Cell_line) +
        ylab(paste0("diff_",firstpoint)) + 
        labs(title = paste0("Boxplots of log fold change from ", firstpoint ," - Non essential and control genes"))
      
      print('DONE')
    }
    })
    
    output$diff_box_all <- renderPlot({
      plot(diff_boxes$diff_box_all)
    })
    
    observe({
      subtitle <- ' - Essential genes'
      req(reactives$diff_t0)
      firstpoint <- input$timepoints_order[[1]]
      #if(is.null(input$essential) | is.null(input$nonessential)){
      if(is.null(geneslists$essential) | is.null(geneslists$non_essential)){
      } else {
      #ess_genes <- ess_genes()
      ess_genes <- geneslists$essential
      diff_boxes$diff_box_ess <-  reactives$diff_t0 %>% 
        filter(.data$Timepoint != !!firstpoint) %>%
        filter(.data$Gene %in% ess_genes[,1]) %>%
        ggplot(aes(x = .data$Timepoint, y= .data$diff, fill = .data$Replicate)) + geom_boxplot() +  facet_grid(.data$Treatment ~ .data$Cell_line) +
        ylab(paste0("diff_",firstpoint)) + 
        labs(title = paste0("Boxplots of log fold change from ", firstpoint ,paste0(subtitle,"'s guides")))
      }
    })
    
    output$diff_box_ess <- renderPlot({
      plot(diff_boxes$diff_box_ess)
    })
    
    output$dldiffboxes <- downloadHandler(
      filename = function(){
        paste("Diff_boxes_",Sys.Date(),".pdf",sep="")
      },
      content = function(file){
        pdf(file = file)
        plot(diff_boxes$diff_box_ess)
        plot(diff_boxes$diff_box_all)
        dev.off()
      }
    )
#################### ROC ####################
ROC <- reactiveValues(plot = NULL,AUC = NULL,toplot = TRUE, d = NULL)

observeEvent(c(input$sidebarmenu,reactives$joined),{
  if (input$sidebarmenu ==  "Roc"){
    if(ROC$toplot == TRUE){
    req(reactives$joined)
    if(is.null(geneslists$essential) | is.null(geneslists$non_essential)){
    #if(is.null(input$essential) | is.null(input$nonessential)){
      showModal(
        modalDialog(tagList(h3("You must provide essentials and non essentials genes list to perform positive screening")),
               footer = tagList(
               modalButton("Got it"))
               )
        )
    } else {
    # ess_genes <- ess_genes()
    # non_ess_genes <- non_ess_genes()
    ess_genes <- geneslists$essential
    non_ess_genes <- geneslists$non_essential
 
    req(input$timepoints_order)
    
    counts <- reactives$joined %>% filter(!(str_detect(gene,paste(c("Non-Targeting","negative_control"),collapse = "|"))))
         
        firstpoint <- input$timepoints_order[[1]]
        counts$Timepoint <- relevel(as.factor(counts$Timepoint), ref = firstpoint)

        fin <- counts %>%
          group_by(sgRNA, Cell_line, Replicate) %>%
          arrange(Timepoint) %>%
          mutate(diff = log_cpm - first(log_cpm)) %>% 
          ungroup()
      #########################################################
        
      withProgress(message = 'Calculating ROC curves', value = 0.5, {
      d <- fin %>% select(.data$sgRNA, .data$Cell_line, .data$Replicate, .data$Timepoint, .data$Gene, .data$Treatment, .data$diff) %>%
        group_by(.data$Timepoint, .data$Treatment, .data$Cell_line, .data$Replicate) %>%
        arrange(.data$diff) %>% 
        mutate(Gene = as.character(Gene)) %>%
        mutate(type = case_when(
          .data$Gene %in% ess_genes[,1] ~ "+",
          .data$Gene %in% non_ess_genes[,1] ~ "-", 
          TRUE ~ NA_character_)) %>%
        filter(!is.na(.data$type)) %>%
        mutate(TP = cumsum(.data$type == "+") / sum(.data$type == "+"), FP = cumsum(.data$type == "-") / sum(.data$type == "-")) %>%
        ungroup()
      
      print("Calulating AU ROC Curves")
      by_rep <- split(d, f= d$Replicate)
      by_rep_cl <- unlist(lapply(X = by_rep, FUN = function(x){split(x, f = x$Cell_line)}),recursive = FALSE)
      by_rep_cl_treat <- unlist(lapply(X = by_rep_cl, FUN = function(x){split(x, f = x$Treatment)}),recursive = FALSE)
      by_rep_cl_treat_tpt <- unlist(lapply(X = by_rep_cl_treat, FUN = function(x){split(x, f = x$Timepoint)}),recursive = FALSE)
      
      AUC_values <- lapply(X= by_rep_cl_treat_tpt, FUN = function(x){
        orders <- order(x$FP)
        print(sum(diff(x$FP[orders])*rollmean(x$TP[orders],2)))
      }
      )
      AUC_table <- data.frame(Replicate = NA, Cell_line = NA, Treatment = NA, Timepoint = NA, AUC = NA)
      for (tablenum in 1:length(by_rep_cl_treat_tpt)){
        vars <- names(by_rep_cl_treat_tpt)[tablenum]
        row <- unlist(strsplit(vars, split =".",fixed = TRUE))
        row <- c(row,as.character(AUC_values[tablenum]))
        AUC_table <- rbind(AUC_table,row)
      }
      AUC_table <- AUC_table[-1,]
      AUC_table$AUC <- signif(as.numeric(AUC_table$AUC),digits =2)
      AUC_table$AUC <- as.character(AUC_table$AUC)
      ROC$AUC <- AUC_table

      ROC$d <- left_join(d,AUC_table, by = c("Cell_line","Replicate","Treatment","Timepoint"))
      })
        ROC$toplot <- FALSE
    }
    }
    }
  })
      
    observeEvent(c(ROC$d,input$labels),{
      req(ROC$d)
      if(input$labels == TRUE){
        ROC$plot <- ROC$d %>% ggplot(aes(x = FP, y = TP, color = Timepoint)) + 
          geom_abline(slope = 1, lty = 3) + 
          xlim(0,1.3) +
          geom_dl(aes(label = AUC), method = "last.polygons")  +
          geom_line() + facet_grid(Treatment + Cell_line ~ Replicate) + 
          coord_equal()
      } else if(input$labels == FALSE){
      ROC$plot <- ROC$d %>% ggplot(aes(x = FP, y = TP, color = Timepoint)) + 
        geom_abline(slope = 1, lty = 3) + 
        xlim(0,1.3) +
        #geom_dl(aes(label = AUC), method = "last.polygons")  +
        geom_line() + facet_grid(Treatment + Cell_line ~ Replicate) + 
        coord_equal()
      }
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
  plot(ROC$plot)
  dev.off()
}
)

output$auc <- DT::renderDataTable({
   DT::datatable(ROC$AUC, filter = "top",
   rownames = FALSE,
   options = list(
    paging =TRUE,
    pageLength =  5 ,
    scrollX= TRUE
  )
  )
})

output$dlauc <- downloadHandler(
  filename = function(){
    paste("AUC_table",Sys.Date(),".csv",sep="")
  },
  content = function(file){
    write.table(x = ROC$AUC[input[["auc_rows_all"]], ],file = file,sep = ",", quote = FALSE, row.names = FALSE)
  }
)

########################## Observers ############################
observeEvent(input$sidebarmenu,{
  if(input$sidebarmenu == "CompCond"){
  updatePickerInput(session,"conditionreference1",
                    choices = as.character(unique(reactives$sampleplan$Treatment)))
  updatePickerInput(session = session,'selecttimepointscomp',
                    choices = as.character(unique(reactives$sampleplan$Timepoint)))
  }
})
observeEvent(input$sidebarmenu,{
  if(input$sidebarmenu == "CompCond"){
  withProgress(message = 'retrieving selected guides', value = 0.5, {
  updatePickerInput(session = session,"selectguidescomp",
                    choices = as.character(unique(reactives$selectedcountsRaw$sgRNA)))
  })
  }
})
## Boxplots 

observeEvent(c(input$splitcelline,input$conditionreference1,reactives$joined,input$selectguidescomp,input$selecttimepointscomp),{
  req(input$conditionreference1)
  if(length(input$conditionreference1) >= 2){
    counts <- reactives$joined %>% 
      filter(Treatment %in% input$conditionreference1) %>%
      filter(!(Timepoint %in% input$selecttimepointscomp))
    counts$Treatment <- as.character(counts$Treatment)
    counts$sgRNA <- as.character(counts$sgRNA)
    counts$Timepoint <- as.character(counts$Timepoint)

    if(input$splitcelline == TRUE){
      interactive_boxplots <- counts %>%
        #ggplot(aes(x = .data$Treatment, y = .data$log_cpm, fill = .data$Treatment)) +
        ggplot(aes(x = .data$Treatment, y = .data$log10_cpm, fill = .data$Treatment)) +
        geom_boxplot_interactive(outlier.shape = NA) + 
        geom_point_interactive(data = subset(counts, sgRNA %in% input$selectguidescomp),
                               #aes(x = .data$Treatment, y = .data$log_cpm, color = .data$Timepoint,tooltip = .data$sgRNA),alpha =0.8) +
                               aes(x = .data$Treatment, y = .data$log10_cpm, color = .data$Timepoint,tooltip = .data$sgRNA),alpha =0.8) +
        facet_grid(.~Cell_line)
    } else {
      interactive_boxplots <- counts %>%
        #ggplot(aes(x = .data$Treatment, y = .data$log_cpm, fill = .data$Treatment)) +
        ggplot(aes(x = .data$Treatment, y = .data$log10_cpm, fill = .data$Treatment)) +
        geom_boxplot_interactive(outlier.shape = NA) + theme(axis.title.x=element_blank(),axis.text.x=element_blank()) +
        geom_point_interactive(data = subset(counts, sgRNA %in% input$selectguidescomp),
                   aes(x = .data$Treatment, y = .data$log10_cpm, color = .data$Timepoint,tooltip = .data$sgRNA),alpha =0.8)
                   #aes(x = .data$Treatment, y = .data$log_cpm, color = .data$Timepoint,tooltip = .data$sgRNA),alpha =0.8)

    }
    reactives$interactive_boxplots <- girafe(ggobj = interactive_boxplots)
  }
})

output$positive_boxplots <- renderGirafe({
  reactives$interactive_boxplots
})

################################################################################
################################# HEATMAP ######################################"
################################################################################

observeEvent(DEA$selected_comp$list,{
  req(DEA$selected_comp$list)
})

output$InfoCompHeatmap <- renderInfoBox({
  infoBox(title = "Selected comparison from RRA scores outlet",
          value = as.character(DEA$selected_comp$list))

})

DeaToClustGenes <- reactiveValues(list = NULL)
#observeEvent(c(DEA$results$scores,DEA$selected_comp$list),priority = 1,{
#observe({
observeEvent(DEA$results$scores[[as.character(DEA$selected_comp$list)]],{
  print("quering rra score genes list")
  #req(DEA$results$scores)
  req(DEA$results$scores[[as.character(DEA$selected_comp$list)]])
  #req(input$ScoreThres)
  req(DEA$selected_comp$list)
   names <- DEA$results$scores[[as.character(DEA$selected_comp$list)]]# %>%
        #filter(RRA_dep_score < 10^-17 | RRA_enrich_score < 10^-17)
        #filter(RRA_dep_adjp < 0.002 | RRA_enrich_adjp < 0.002)
        #filter(RRA_adjp < input$ScoreThres)
   
   #remove control guides
   names <- filter(names, !(str_detect(Gene,paste(c("Non-Targeting","negative_control"),collapse = "|"))))
   if(nrow(names) > 100){
     names <- names[1:100,]
     names <- names$Gene
     DeaToClustGenes$list <- names
   } else {
     print("notenougthgenes for heatmap")
   }
})
# 
ClustData_ess <- reactiveValues(table = NULL)
ClustData_non_ess <- reactiveValues(table = NULL)
ClustData <- reactiveValues(table = NULL)
ClustMetadata <- reactiveValues(table = NULL)
observeEvent(c(input$sidebarmenu,
               DeaToClustGenes$list,
               reactives$sampleplan,
               #reactives$norm_data$counts,
               reactives$norm_cpm,
               DEA$selected_comp$list),{
  if (input$sidebarmenu == "Clustering"){
    print("runing clustering module")
    req(reactives$norm_data$genes)
    req(reactives$sampleplan)
    #req(reactives$norm_data$counts)
    req(reactives$norm_cpm)
    #if(!is.null(DeaToClustGenes$list)){
    if(!is.null(DEA$results$scores[[DEA$selected_comp$list]])){
    withProgress(message = 'Filtering data for clustering', value = 0.5, {
      
    ClustMetadataa <-  column_to_rownames(reactives$sampleplan,"Sample_ID")
    sgRNAannot <- reactives$norm_data$genes
    
    ################# NORM PIERRE ######
    ClustDataa <- reactives$norm_cpm
    #ClustDataa <- reactives$norm_data$counts
    ClustDataa <- tibble::rownames_to_column(as.data.frame(ClustDataa),"sgRNA")
    ClustDataa <- ClustDataa %>%
      left_join(sgRNAannot, by = "sgRNA")
    ClustDataa <- ClustDataa %>% filter(Gene %in% as.character(DeaToClustGenes$list))
    
    ClustDataa <- as.data.frame(ClustDataa)
    if(ncol(ClustDataa) > 0 && nrow(ClustDataa) > 0){
 
      ClustDataa <- column_to_rownames(ClustDataa,"sgRNA") %>% select(-c(Gene))
      ClustMetadataa <- ClustMetadataa[rownames(ClustMetadataa) %in% colnames(ClustDataa),]
      ClustMetadataa <- ClustMetadataa[rownames(ClustMetadataa) %in% colnames(ClustDataa),]
      ClustDataa <- ClustDataa[,colnames(ClustDataa) %in% rownames(ClustMetadataa)]
      ClustMetadata$table <- ClustMetadataa
      ClustData$table <- ClustDataa
    }
    })
    } else {
    # 
    showModal(modalDialog(HTML(
        "To do so go through all the Statistical analysis steps"),
        br(),
        title = "Compute RRA scores first ",
        footer = tagList(
        modalButton("Got it"))
    ))
    }
  }
})

heatmap <- callModule(ClusteringServerMod, id = "heatmapID", session = session,
                      data = ClustData , metadata =  ClustMetadata, printRows = FALSE)

observeEvent(ClustData$table,ignoreInit = TRUE,{
  if(nrow(ClustData$table) <= 3){
        showModal(modalDialog(HTML(
        "Rerun a comparison on the rra scores tab"),
        br(),
        title = "There is not enough sgRNA passing filters to perform heatmap",
        footer = tagList(
          modalButton("Got it")
        )))
  }
})

######################################################################################################
######################################## DEA #########################################################
######################################################################################################

DEAMetadata <- reactiveValues(table = NULL)
DEAnormdata <- reactiveValues(data = NULL)
observeEvent(c(reactives$sampleplan,input$sidebarmenu),{
  if(!is.null(reactives$sampleplan)){
    if(input$sidebarmenu=="Statistical_analysis"){

    DEAMetadata$table <- reactives$sampleplan %>%
       filter(!(Sample_ID %in% input$removesamples)) %>%
       column_to_rownames("Sample_ID") %>%
       mutate(Cell_line = as.character(Cell_line)) %>%
       select(c("Cell_line","Timepoint","Treatment","SupplementaryInfo")) 
    }
  }
})
observeEvent(c(reactives$norm_data,input$sidebarmenu),{
  if(!is.null(reactives$norm_data)){
    if(input$sidebarmenu=="Statistical_analysis"){
    DEAnormdata$data <- reactives$norm_data
    }
  }
})

observeEvent(input$sidebarmenu,{
  if(input$sidebarmenu=="Rawdist" | input$sidebarmenu=="Tev" | input$sidebarmenu=="Roc" | input$sidebarmenu == "Clustering" | input$sidebarmenu == "CompCond"){
    #if(is.null(input$counts) | is.null(input$sample_plan)){
    if(is.null(reactives$counts) | is.null(reactives$sampleplan)){
      showModal(modalDialog(
        title = "Upload both count matrix and sampleplan first",
        footer = tagList(
          modalButton("Got it")
        )))
  }}
})

DEA <- callModule(CRISPRDeaModServer, "DEA", session = session,
                  norm_data = DEAnormdata,
                  sampleplan = DEAMetadata, 
                  #ctrlterm = ctrlterm$term,
                  ctrlterm = ctrlterm,
                  # ess_genes=ess_genes,
                  # non_ess_genes = non_ess_genes)
                  # ess_genes=geneslists$essential,
                  # non_ess_genes = geneslists$non_essential)
                  ess_genes= geneslists,
                  non_ess_genes = geneslists)

observeEvent(DEA$concatenated$results,ignoreInit = TRUE,{
req(DEA$concatenated$results)
updatePickerInput(session = session, 'volcanoslist',choices = as.character(names(DEA$concatenated$results)))
})

observeEvent(input$sidebarmenu,priority = -1,{
  if(input$sidebarmenu=="Statistical_analysis"){
  if(is.null(DEAnormdata$data) | is.null(DEAMetadata$table)){
    if(input$screentype == "negative"){
    if(is.null(reactives$control_sgRNA)){
          showModal(modalDialog(
            title = "No control guides found in data",
            HTML("Have you used the appropriate term to find it ?"),
            footer = tagList(
              modalButton("Got it")
            ))
          )
  } else {
    print("hoooo")
  showModal(modalDialog(
    title = "Upload both count matrix and sampleplan first",
    footer = tagList(
      modalButton("Got it")
    )))
     }}
    }}
})
#########################################################################
############################ Help section ###############################
#########################################################################
    
output$Totguidenumber <- renderInfoBox(
    infoBox(
      "Total guides number :",
      as.numeric(nrow(reactives$countsRaw)),
      icon = shiny::icon("arrows-alt-h",lib = "font-awesome")
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
      mean_depth(),
    icon = shiny::icon("arrows-alt-v",lib = "font-awesome")
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
<li>An xlsx or a csv file using 'commas'  'semicolons' or 'tabulations' as field separator, each line of the file must respects the following format specifications :<br/>
Sample_ID;Cell_line;Timepoint;Treatment;Replicate;SupplementaryInfo</li>
<li>Column names must respect Upper and lower case. </li>
<li>Values in the table must not contain spaces, use '_' instead. </li>
<li>Values in the table must not contain dots </li>
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
<li>An xlsx or a csv formated text file using 'commas'  'semicolons' or 'tabulations' as field separator. </li>
<li>The first line of the file is a header, it contains samples Sample_IDs as colnames and two supplementary columns called 'gene
' and 'sequence'. </li>
<li>Guides names' are rownames. </li>
<li>Values are read counts.</li>
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
})
    
    observeEvent(c(reactives$sampleplanRaw),priority = 10,{
    colnames <- colnames(reactives$sampleplanRaw)
    if( (!"Sample_ID" %in% colnames|
        !"SupplementaryInfo" %in% colnames|
        !"Cell_line" %in% colnames|
        !"Timepoint" %in% colnames|
        !"Treatment" %in% colnames|
        !"Replicate" %in% colnames)  ){
      showModal(modalDialog(tagList(h3('Colnames must contain these values :',style="color:red"),
                                    h4("| Sample_ID | Cell_line | Timepoint | Treatment | Replicate | SupplementaryInfo"),
                                    h4("colnames must respect majuscules",style="color:red"),
                                    p(),
                                    h3("Your current file colnames are : ",strong ="bold"),
                                    h4(paste("| ",paste(colnames,collapse =" | ")))
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
                 separators = input$Fsc,
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
              input$Fsc <- tmpEnv$r_input$Fsc
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

shinyjs::hide(id = "acereport_rmd")
shinyjs::hide(id= "enableAutocomplete")
shinyjs::hide(id= "enableLiveCompletion")
shinyjs::hide("enableRCompletion")
shinyjs::hide(id= "mode")
shinyjs::hide(id= "theme")
#shinyjs::hide(id= "editor_options")
### ace editor options
observe({
  autoComplete <- if(input$enableAutocomplete) {
    if(input$enableLiveCompletion) "live" else "enabled"
  } else {
    "disabled"
  }
  updateAceEditor(session, "acereport_rmd", autoComplete = autoComplete,theme=input$theme, mode=input$mode)
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
        #tmp_content <- rmd_yaml()
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
    #tmp_content <- rmd_yaml()
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

#' Server function
#'
#' Function containing all the server side of the Shiny App, must be called inside the run_app() function 
#'
#'
#' @return None
#'
#' @examples
#'
#' @export
#' @import CRISPRApp
crispr_server <- function(input, output, session) {

##### Upload files and datatable construction ####################

counts <- reactive({

  req(input$counts)
  inFile <- input$counts
  counts <- read_delim(inFile$datapath, delim = input$FS)
  counts <- dplyr::rename(counts, sgRNA = X1)
  counts <- dplyr::select(counts, -sequence)

  annot_sgRNA <- dplyr::select(counts, sgRNA, Gene = gene)

  counts <- gather(counts, value = count, key = barcode, -sgRNA, -gene)
  counts <- dplyr::mutate(counts, barcode = str_remove(barcode, ".R1.fastq"))

  counts <- counts %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(cpm = 1e6 * count / sum(count), log_cpm = log10(1 + cpm)) %>%
    dplyr::ungroup()

  return(list(counts,annot_sgRNA))

})

sample_plan <- reactive({

  req(input$sample_plan)
  inFile <- input$sample_plan
  samples <- read_delim(inFile$datapath, delim = "|",  col_names = c("barcode", "sample","condition"))
  samples <- samples %>%
  separate(sample, into = c("date","rep","clone","day"), remove = FALSE) %>%
    mutate(day = as.factor(day)) %>%
    mutate(condition = as.factor(condition)) %>%
    mutate(day_num = as.numeric(gsub("DAY","",day))) %>% 
    mutate(day = fct_reorder(.f = factor(day), .x = day_num))
  return(samples)
})

joined <- reactive({
  counts <- counts()[[1]]
  samples <- sample_plan()
  counts <- full_join(counts, samples) %>%
  mutate(day = factor(day, levels = input$timepoints_order))
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
   counts <- counts()[[1]]
   firstpoint <- input$timepoints_order[[1]]
   
   all <- joined() %>%
     select(sgRNA,clone, rep, day, log_cpm, gene, condition)

   t0 <- all %>%  
     filter(day == firstpoint)  %>%
     mutate(log_cpmt0 = log_cpm) %>%
     mutate(log_cpm = NULL) %>%
     mutate(day = NULL) %>%
     mutate(gene = NULL) %>%
     mutate(condition = NULL)
   
   all <- inner_join(all,t0,c("sgRNA","clone","rep")) %>%
     mutate(diff = log_cpm - log_cpmt0) %>%
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
output$counts_table <- DT::renderDataTable({
  DT::datatable(counts()[[1]],rownames=FALSE)
  })

output$sample_plan_table <- DT::renderDataTable({
  sample_plan <- sample_plan()
  sample_plan <- subset(sample_plan, select = -c(day_num))
  return(DT::datatable(sample_plan,rownames = FALSE))
})

output$joined <- DT::renderDataTable({
  diff_t0()
})



read_number <- reactive({
  counts <- joined()
  counts <- counts %>%
    group_by(barcode, rep, clone, day) %>%
    summarise(total = sum(count)) %>% 
    as.data.frame() %>%
     ggplot(aes(x = day, y = total)) +
     geom_col(position = position_dodge()) + facet_wrap(vars(rep, clone), nrow = 1) +
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
     ggplot(aes(x = day, y = log_cpm, fill = rep)) + geom_boxplot() + facet_grid(. ~ clone) + 
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
   ggplot(aes(x = log_cpm, y = day, fill = ..x..)) +
    geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
    facet_wrap(vars(clone, rep), ncol = 1, strip.position = "right") +
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
  counts <- filter(counts, gene %in% ess_genes$V1)
  counts_plot <- ggplot(counts, aes(x = log_cpm, y = day, fill = ..x..)) +
    geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
    facet_wrap(vars(clone, rep), ncol = 1, strip.position = "right") +
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
    filter(gene %in% ess_genes$V1)
  
  counts_plot <- ggplot(counts, aes(x = log_cpm, y = day, fill = ..x..)) +
    geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
    facet_wrap(vars(clone, rep), ncol = 1, strip.position = "right") +
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
  counts <- levels(sample_plan()$day)
  print(counts)
  orderInput(inputId = 'timepoints', label = 'Re-order your timepoints here :', items = counts, as_source = F)
})



############## NEGATIV SCREENING ################

diff_box_all <- reactive({
  firstpoint <- input$timepoints_order[[1]]

    diff_box_all <- diff_t0() %>%
    filter(day != firstpoint) %>%
    ggplot(aes(x = day, y = diff, fill = rep)) + geom_boxplot() + facet_grid(condition ~ clone) +
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
    filter(day != !!firstpoint) %>%
    filter(gene %in% ess_genes[,1]) %>%
    ggplot(aes(x = day, y= diff, fill = rep)) + geom_boxplot() + facet_grid(.~ clone) +
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
   
  d <- diff_t0() %>% select(sgRNA,clone, rep, day,log_cpm, gene, condition, diff) %>%
    group_by(day,condition,clone,rep) %>%
      arrange(diff) %>% 
      mutate(type = case_when(gene %in% ess_genes[,1] ~ "+",gene %in% non_ess_genes[,1] ~ "-", TRUE ~ NA_character_)) %>%
     filter(!is.na(type)) %>%
     mutate(TP = cumsum(type == "+") / sum(type == "+"), FP = cumsum(type == "-") / sum(type == "-")) %>%
     ungroup()
  
  write.csv(diff_t0(),"~/to_shiny.csv")

    d <- d %>% ggplot(aes(x = FP, y = TP, color = day)) + geom_abline(slope = 1, lty = 3) + geom_line() + facet_grid(condition + clone ~ rep) + coord_equal()
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
                    choices = sample_plan()$condition)
  
  updateSelectInput(session,"conditionreference1",
                    choices = sample_plan()$condition)
})





########################################################
################## Help section ############


output$Totguidenumber <- renderInfoBox(
  infoBox(
    "Total guides number :",
    as.numeric(nrow(counts()[[1]]))
  )
)

output$Depth <- renderInfoBox(
  infoBox(
    "Average sequencing depth across samples :",
    sum(counts()[[1]]$count)/(nrow(sample_plan()) + nrow(counts()[[1]]))
    
  )
)

output$Datahelptext <- renderUI({HTML(
  "

<ul> 
<li>Sample infos file :
<br/><br/>
<B>Format description :</B>
<br/>
A text file, each line of the file must respects the following format specifications :<br/>
SampleBarcode|SampleCreationDate-Replicat-Souche-Day|Condition
<br/>
<B>For example :</B>
<br/>
D115R13|191015-Replica1-WT-DAY14|TREATED
<br/>
<a href='https://gitlab.curie.fr/r-shiny/bioshiny/blob/devel/example_datasets/Crispr/SampleDescription' target='_blank'>Click here to download the Sample Description example file</a>
<br/>
<br/>
</li><li>Counts table file :</li>
<br/>
<B>Format description :</B>
<br/>
A csv formated text file using ; or , as field separator. The first line of the file is a header, it contains samples barcodes as colnames and two supplementary columns called 'gene
' and 'sequence'. Guides names' are rownames. Values are read counts.
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
<a href='https://gitlab.curie.fr/r-shiny/bioshiny/blob/devel/example_datasets/Crispr/essentials.csv' target='_blank'>Click here to download the essentials genes'list example file</a>
<br/>

<B>For example:</B><br/>

Gene1<br/>
Gene2<br/>
...


</ul>

  
  
  
  "
)})



}


#shinyApp(server = crispr_server())
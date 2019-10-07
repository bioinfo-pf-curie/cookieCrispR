
function(input, output, session) {

require(shiny)
require(shinyjs)
require(shinydashboard)
library(tidyverse)
library(dplyr)
require(ggridges)

##### Upload files and datatable construction ####################"

counts <- reactive({

  req(input$counts)
  inFile <- input$counts
  counts <- read_delim(inFile$datapath, delim = ";")
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
  samples <- read_delim(inFile$datapath, delim = "|",  col_names = c("barcode", "sample"))
  samples <- samples %>%
  separate(sample, into = c("rep", "clone", "day"), remove = FALSE) %>%
    mutate(day = as.factor(day)) %>%
    mutate(day_num = as.numeric(gsub("J","",day))) %>% 
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



## compute diff to t0
diff_t0 <- reactive({

  
  req(input$timepoints_order)
  counts <- counts()[[1]]
  # t0 <- joined() %>%
  #   select(sgRNA, clone, gene ,rep, day, log_cpm) %>%
  #   spread(key = day, value = log_cpm) %>%
  #   mutate_at(vars(J4:J16), ~ . - J1) %>%
  #   mutate(J1 = 0) %>%
  #   gather(key = day, value = diff_J1, -sgRNA, -clone, -rep, -gene)
  # 
  # 
  # 
  firstpoint <- input$timepoints_order[[1]]
  joined <- joined()
  t0 <- joined %>%
    select(sgRNA,clone, rep, day, log_cpm, gene) %>%
    spread(key = day, value = log_cpm) %>%
    mutate_at(vars(input$timepoints_order[[2]]:input$timepoints_order[[length(input$timepoints_order)]]), ~ . - !!(as.name(firstpoint))) %>%
    mutate(!!firstpoint := 0) %>%
     gather(key = day, value = diff_J1, -sgRNA, -clone, -rep, -gene) %>%
     mutate(day = factor(day, levels = levels(joined$day)))

  
  t0 <- inner_join(counts, t0, by="sgRNA")
  return(t0)
})
  



######### Plots and tables outputs ####################################
output$counts_table <- DT::renderDataTable({

  DT::datatable(counts()[[1]],rownames=FALSE)

  })

output$sample_plan_table <- DT::renderDataTable({

  sample_plan <- sample_plan()
  sample_plan <- subset(sample_plan,select = -c(day_num))
  return(DT::datatable(sample_plan,rownames = FALSE))

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

  counts <- counts %>%
    filter(gene %in% ess_genes$V1)
  
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




diff_box_all <- reactive({
  
  
  firstpoint <- input$timepoints_order[[1]]
  
  diff_box_all <- diff_t0() %>%
    filter(day != firstpoint) %>%
    ggplot(aes(x = day, y = diff_J1, fill = rep)) + geom_boxplot() + facet_grid(.~ clone) +
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
    filter(day != !!firstpoint) %>%
    filter(gene.x %in% ess_genes$V1) %>%
    ggplot(aes(x = day, y= diff_J1, fill = rep)) + geom_boxplot() + facet_grid(.~ clone) +
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

########################## Observers ############################









########################################################




}





require(shiny)
require(shinyjs)
require(shinydashboard)
require(shinyjqui)
#
options(shiny.maxRequestSize = 1000*1024^2)
#

# 
# 
#   
# require(shiny)
# require(shinyjs)
# require(shinydashboard)
# library(tidyverse)
# library(dplyr)
# require(ggridges)
#   
# ##### Upload files and datatable construction ####################"
#   
# 
#     inFile <- input$counts
#     counts <- read_delim(inFile$datapath, delim = ";")
#     counts <- dplyr::rename(counts, sgRNA = X1)
#     counts <- dplyr::select(counts, -sequence)
#     
#     annot_sgRNA <- dplyr::select(counts, sgRNA, Gene = gene)
#     
#     
#     counts <- gather(counts, value = count, key = barcode, -sgRNA, -gene)
#     counts <- dplyr::mutate(counts, barcode = str_remove(barcode, ".R1.fastq"))
#     
#     counts <- counts %>%
#       dplyr::group_by(barcode) %>%
#       dplyr::mutate(cpm = 1e6 * count / sum(count), log_cpm = log10(1 + cpm)) %>%
#       dplyr::ungroup()
#     
#     return(list(counts,annot_sgRNA))
#     
#   })
#   
#   sample_plan <- reactive({
#     
#     req(input$sample_plan)
#     inFile <- input$sample_plan
#     samples <- read_delim(inFile$datapath, delim = "|",  col_names = c("barcode", "sample"))
#     samples <- samples %>%
#       separate(sample, into = c("rep", "clone", "day"), remove = FALSE)
#     samples <- dplyr::mutate(samples, day = factor(day, levels = c("J1", "J4", "J7", "J9", "J11", "J14", "J16")))
#     return(samples)
#   })
#   
#   joined <- reactive({
#     
#     counts <- counts()[[1]]
#     samples <- sample_plan()
#     counts <- full_join(counts, samples) %>%
#       mutate(day_num = as.numeric(gsub("J","",day))) %>% 
#       mutate(day = fct_reorder(.f = factor(day), .x = day_num))
#     
#     return(counts)
#     
#   })
#   
#   
#   
#   #timepoints <- reactiveValues(a = joined()$day)
#   
#   
#   
#   
#   ess_genes <- reactive({
#     
#     req(input$essential)
#     inFile <- input$essential
#     ess_genes <- read.table(inFile$datapath, header=FALSE)
#     return(ess_genes)
#     
#     
#   })
#   
#   non_ess_genes <- reactive({
#     
#     req(input$nonessential)
#     inFile <- input$nonessential
#     non_ess <- read.table(inFile$datapath, header = FALSE)
#     return(non_ess)
#     
#     
#   })
#   
#   
#   # 
#   # output$order <- renderPrint({
#   #   print(input$timepoints_order)
#   #   })
#   # 
#   output$orderUI <- renderUI({
#     
#     counts <- levels(joined()$day)
#     #print(counts$day)
#     orderInput(inputId = 'timepoints', label = 'Re-order your timepoints here', items = c(1,2), as_source = F)
#     
#   })
#   #})
#   
#   
#   
#   
#   ######### Plots and tables outputs ####################################
#   output$counts_table <- DT::renderDataTable({
#     
#     counts()[[1]]
#     
#   })
#   
#   output$sample_plan_table <- DT::renderDataTable({
#     
#     sample_plan()
#     
#   })
#   
#   
#   
#   
#   output$read_number <- renderPlot({
#     #output$read_number <- DT::renderDataTable({
#     
#     
#     
#     
#     counts <- joined()
#     
#     counts <- counts %>%
#       group_by(barcode, rep, clone, day) %>%
#       summarise(total = sum(count)) %>% 
#       as.data.frame() %>%
#       ggplot(aes(x = day, y = total)) +
#       geom_col(position = position_dodge()) + facet_wrap(vars(rep, clone), nrow = 1) +
#       labs(title = "Number of reads by sample") #+
#     #scale_x_continuous(breaks = seq(min(counts$day),max(counts$day)))
#     
#     return(plot(counts))
#     
#     #return(counts)
#   })
#   
#   
#   
#   
#   output$boxplot_all <- renderPlot({
#     
#     
#     counts <- joined()
#     counts %>% 
#       ggplot(aes(x = day, y = log_cpm, fill = rep)) + geom_boxplot() + facet_grid(. ~ clone) + 
#       labs(title = "Distribution of normalized log-cpm by sample", subtitle = "All guides")
#     
#   })
#   
#   output$density_ridge <- renderPlot({
#     
#     
#     
#     counts <- joined()
#     
#     
#     counts <-  counts %>%
#       ggplot(aes(x = log_cpm, y = day, fill = ..x..)) +
#       geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
#       facet_wrap(vars(clone, rep), ncol = 1, strip.position = "right") +
#       labs(title = "Distribution of normalzed log-cpm by sample", subtitle = "All guides")
#     return(plot(counts))
#   })
#   
#   output$essential_distribs <- renderPlot({
#     
#     counts <- joined()
#     
#     ess_genes <- ess_genes()
#     
#     counts <- counts %>%
#       filter(gene %in% ess_genes$V1)
#     
#     counts_plot <- ggplot(counts, aes(x = log_cpm, y = day, fill = ..x..)) +
#       geom_density_ridges(alpha = 0.6, show.legend = FALSE, fill = "gray50") +
#       facet_wrap(vars(clone, rep), ncol = 1, strip.position = "right") +
#       scale_fill_viridis_c() +
#       labs(title = "Distribution of normalised log-cpms", subtitle = "Essential genes")
#     
#     return(plot(counts_plot))
#   })
#   
#   
#  
# 

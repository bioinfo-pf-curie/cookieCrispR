require(shiny)
require(shinyjs)
require(shinydashboard)
library(tidyverse)
library(dplyr)
require(ggridges)
require(shinyjqui)
# #
# options(shiny.maxRequestSize = 1000*1024^2)
# #
# 
# #
# #
# #
# folder <- "/bioinfo/users/cbenoit/Documents/Rshiny/ShinyApps/CRISPRApp_packaged_2/testdatasets"
# #
# # ##### Upload files and datatable construction ####################"
# 
# #counts
# counts <- read_delim(file.path(folder,"global_counts_table_truncated.csv"), delim = ";")
#      counts <- dplyr::rename(counts, sgRNA = X1)
#      counts <- dplyr::select(counts, -sequence)
#      annot_sgRNA <- dplyr::select(counts, sgRNA, Gene = gene)
#      counts <- gather(counts, value = count, key = barcode, -sgRNA, -gene)
#      counts <- dplyr::mutate(counts, barcode = str_remove(barcode, ".R1.fastq"))
#      counts <- counts %>%
#        dplyr::group_by(barcode) %>%
#        dplyr::mutate(cpm = 1e6 * count / sum(count), log_cpm = log10(1 + cpm)) %>%
#        dplyr::ungroup()
# #########################
# 
# #sample_plan##################
# 
# samples <- read_delim(file.path(folder,"D115.sampleDescription.txt"), delim = "|",  col_names = c("barcode", "sample"))
# samples <- samples %>%
# separate(sample, into = c("rep", "clone", "day"), remove = FALSE)
# samples <- dplyr::mutate(samples, day = factor(day, levels = c("J1", "J4", "J7", "J9", "J11", "J14", "J16")))
# 
# ################ Joined ###################
# joined <- full_join(counts, samples) %>%
#        mutate(day_num = as.numeric(gsub("J","",day))) %>%
#        mutate(day = fct_reorder(.f = factor(day), .x = day_num))
# 
# ################# ess_genes #################
# ess_genes <- read.table(file.path(folder,"essential-genes.txt"), header=FALSE)
# 
# # timepoints_order <- c("J1", "J4", "J7", "J9", "J11", "J14", "J16")
# #
# #
# #
# #
# # t0 <- joined %>%
# #   select(sgRNA, clone, gene,rep, day, log_cpm) %>%
# #   spread(key = day, value = log_cpm) %>%
# #   mutate_at(vars(J4:J16), ~ . - J1) %>%
# #   mutate(J1 = 0) %>%
# #   gather(key = day, value = diff_J1, -sgRNA, -clone, -rep,-gene)
# 
# 
# 
# timepoints_order <- c("J1", "J4", "J7", "J9", "J11", "J14", "J16")
# 
# firstpoint <- timepoints_order[[1]]
# t0 <- joined %>%
#   select(sgRNA,clone, rep, day, log_cpm) %>%
#   spread(key = day, value = log_cpm) %>%
#   mutate_at(vars(timepoints_order[[2]]:timepoints_order[[length(timepoints_order)]]), ~ . - !!(as.name(firstpoint))) %>%
#   mutate(!!firstpoint := 0) %>%
#    gather(key = day, value = diff_J1, -sgRNA, -clone, -rep) #%>%
#   # mutate(day = factor(day, levels = levels(counts$day)))
# 
# 
# 
# 
# t0 <- inner_join(counts, t0, by="sgRNA")
# 
# 
# diff_box_ess <-  t0 %>%
#   filter(day != !!firstpoint) %>%
#   filter(gene %in% ess_genes$V1) %>%
#   ggplot(aes(x = day, y = diff_J1, fill = rep)) + geom_boxplot() + facet_grid(.~ clone) +
#   labs(title = "Boxplots of log fold change from J1 - essential genes")
# 
# plot(diff_box_ess)
# 


guide_plot <- function(counts, gene_id, normd0 = TRUE) {
  d <- filter(counts, gene == gene_id)
  if (normd0) {
    d0 <- group_by(d, sgRNA, rep, clone) %>% filter(day == "J1") %>% transmute(cpm_d0 = norm_cpm)
    d <- full_join(d, d0) %>% mutate(norm_cpm = norm_cpm - cpm_d0)
  }
  g <- ggplot(d, aes(x = day, y = norm_cpm, color = sgRNA, group = sgRNA)) + geom_point() + geom_line(linetype = 2) + ggtitle(gene_id) + facet_grid(rep ~ clone)
  return(g)
}



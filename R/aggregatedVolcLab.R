# load("~/aggregated_volcano.rda")
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# #library(tidyverse)
# res <- res %>% select(Gene, estimate)
# res_aggregated <- aggregate(res$estimate,list(res$Gene),mean)
# colnames(res_aggregated) <- c("Gene","meanlogFC")
# selected_scores <- selected_scores %>% select(Gene,RRA_dep_score)
# joined <- left_join(res_aggregated,selected_scores, by = "Gene")
# 
# ess_genes <- read.table("/home/cbenoit/Bureau/DATATEST COOCKIE/CRISPR/Mine/Essentials.txt")
# non_ess_genes <- read.table("/home/cbenoit/Bureau/DATATEST COOCKIE/CRISPR/Mine/Non essential.txt")
# ess_genes$V2 <- 'essential genes'
# non_ess_genes$V2 <- 'non essential genes'
# colnames(ess_genes) <- c("Gene","Type")
# colnames(non_ess_genes) <- c("Gene","Type")
# genetypes <- rbind(ess_genes,non_ess_genes)
# joined <- left_join(joined,genetypes, by = c("Gene")) %>%
#   replace_na(list(Type ="Other"))
# 
# ggplot <- ggplot(joined, aes(x = meanlogFC, y = -log10(RRA_dep_score), colour = factor(Type),alpha = factor(Type),
#                              size = factor(Type))) +
#   expand_limits(y = c(min(-log10(joined$RRA_dep_score)), 1)) +
#   scale_colour_manual(name = NULL, 
#                       values =c('Other'='black','essential genes'='red', "non essential genes" = "blue"), 
#                       labels = c('Other','essential genes',"non essential genes")) +
#   theme(legend.position = "top",
#         legend.direction = "horizontal") +
#   scale_size_manual(values =c('Other'=0.5,'essential genes'=1.3, "non essential genes" = 1.3)) +
#   scale_alpha_manual(values =c('Other'=0.4,'essential genes'=0.6, "non essential genes" = 0.6)) + 
#   guides(size = FALSE, alpha =FALSE)
# 
# plot(ggplot)
# 

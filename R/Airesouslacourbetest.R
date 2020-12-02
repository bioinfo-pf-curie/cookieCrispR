rm(list = ls())
load("~/coockie_crispr_d.rda")
library(zoo)
library(ggplot2)
library(dplyr)
library(tidyverse)

plot <- d %>% ggplot(aes(x = FP, y = TP, color = Timepoint)) + geom_abline(slope = 1, lty = 3) + geom_line() + facet_grid(Treatment + Cell_line ~ Replicate) + coord_equal()
plot(plot)

### Splitting data #####
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
AUC_table <- AUC_table[-1,]
head(AUC_table)

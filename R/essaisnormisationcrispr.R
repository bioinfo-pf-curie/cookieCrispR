# 
# load("~/c.rda")
# library(stringr)
# library(DESeq2)
# library(edgeR)
# library(dplyr)
# 
# control_sgRNA <- filter(sgRNAannot, str_detect(Gene,"Non-Targeting"))
# control_sgRNA <- sgRNAannot[str_detect(sgRNAannot$Gene,"Non-Targeting"),]
# control_sgRNA <- sgRNAannot[str_detect(sgRNAannot$Gene,"BCAS2"),]
# control_sgRNA
# 
# norm_data <- sg_norm(ClustData %>% select(-c(gene,sgRNA)),sgRNAannot= sgRNAannot)

###save App state to binary data
saveState <- function(filename,reactives,separators,input) {
#saveState <- function(session,filename) {
  isolate({
    
    print("save inputs")
    #r_inputs <<- lapply(reactiveValuesToList(input), unclass)
    #r_inputs <- reactiveValuesToList(input)
    r_inputs <<- input
    print("save reactives")
    # #r_reactives <- reactiveValuesToList(reactives)
    r_reactives <<- reactives
    print("save separators")
    # #r_separators <- reactiveValuesToList(separators)
    r_separators <<- separators
    save(list = c("r_inputs", "r_reactives", "r_separators"),  file = filename)

  })
}

###### Normalisation function #################
sg_norm <- function(counts, control_sgRNA = NULL, sample_annot = NULL, sgRNA_annot = NULL, ...){
  
  # contruct DGEList
  dge_all <- DGEList(counts = counts, samples = sample_annot, genes = sgRNA_annot)
  
  # compute normFactors on control if present
  if (!is.null(control_sgRNA)) {
    norm_factors <- calcNormFactors(counts[control_sgRNA, , drop = FALSE], ...)
  } else norm_factors <- calcNormFactors(counts, ...)
  
  # set normFactors to all data
  dge_all$samples$norm.factors <- norm_factors
  # set lib.size to number of control reads
  dge_all$samples$lib.size <- colSums(counts[control_sgRNA, , drop = FALSE])
  
  return(dge_all)
}

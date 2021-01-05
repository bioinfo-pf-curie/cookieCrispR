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
  save(list = c("dge_all"), file = "~/DGE_list.rda")
  
  
  # compute normFactors on control if present
  if (!is.null(control_sgRNA)) {
    print("control guides found")
    norm_factors <- calcNormFactors(counts[control_sgRNA, , drop = FALSE], ...)
  } else {
    norm_factors <- calcNormFactors(counts, ...)
    print("no control guides found")
  }
  
  # set normFactors to all data
  dge_all$samples$norm.factors <- norm_factors
  # set lib.size to number of control reads
  dge_all$samples$lib.size <- colSums(counts[control_sgRNA, , drop = FALSE])
  
  return(dge_all)
}




compute_model <- function(norm_data, design, cpm_thr = 1, sample_thr = 2,
                          normalize.method = "none", span = 0.5, ## voom options
                          method = "ls",  ## lmFit options
                          robust = FALSE,
                          voom_qw = FALSE,
                          plot_voom = FALSE
){
  kept <- which(rowSums(edgeR::cpm(norm_data) >= cpm_thr) >= sample_thr)
  if (!voom_qw) {
    v <- limma::voom(norm_data[kept, ], design = design, normalize.method = normalize.method, span = span, plot = plot_voom)
  } else 
    v <- limma::voomWithQualityWeights(norm_data[kept, ], design = design, normalize.method = normalize.method, span = span, plot = plot_voom)
  res_fit <- limma::lmFit(v, method = method)
  res_eb <- limma::eBayes(res_fit, robust = robust)
  return(list(v = v, res_eb = res_eb))
}

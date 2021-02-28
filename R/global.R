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


compute_RRA <- function(x, alpha_thr = 1){
  x <- x[x <= alpha_thr]
  l <- length(x)
  if (length(x)){
    pb <- pbeta(sort(x), 1:l, l + 1 - 1:l) ## ~Beta(k, n + 1 -k)
    score <- min(pb)
    w <- which.min(pb)
  } else {
    score <- 1
    w <- 1
  }
  res <- data.frame(score = score,  which = w)
  return(res)
}

compute_score_RRA <- function(object, alpha_thr = 1){
  RRA_pvalue <- tapply(object$p.value, object$Gene, compute_RRA, alpha_thr = alpha_thr) %>%
    # bind_rows(.id = "Gene") %>%
    imap_dfr(~mutate(.x, Gene = .y)) %>%
    dplyr::rename(RRA_score = score, RRA_which = which)
  RRA_pvalue_dep <- tapply(object$p.value_dep, object$Gene, compute_RRA, alpha_thr = alpha_thr) %>%
    # bind_rows(.id = "Gene") %>%
    imap_dfr( ~mutate(.x, Gene = .y)) %>%
    dplyr::rename(RRA_dep_score = score, RRA_dep_which = which)
  RRA_pvalue_enrich <- tapply(object$p.value_enrich, object$Gene, compute_RRA, alpha_thr = alpha_thr) %>%
    # bind_rows(.id = "Gene") %>%
    imap_dfr( ~mutate(.x, Gene = .y)) %>%
    dplyr::rename(RRA_enrich_score = score, RRA_enrich_which = which)
  object <- group_by(object, Gene) %>% nest
  res <- suppressMessages(reduce(list(object, RRA_pvalue, RRA_pvalue_dep, RRA_pvalue_enrich), full_join, by = "Gene")) %>%
    select(c("RRA_score","RRA_dep_score","RRA_enrich_score"))
  return(res)
}


compute_RRA_pval <- function(guide_res, gene_res, non_target, n_perm = 10000, n_guides = 6, alpha_thr = 1){
  ## create random table
  object_non_target <- filter_(guide_res, ~sgRNA %in% non_target$sgRNA)
  #object_non_target <- filter_(guide_res, ~Gene %in% non_target$Gene)
  random_tab <- data.frame(sgRNA = sample(non_target$sgRNA, n_guides * n_perm, replace = TRUE),
                           Gene = rep(paste0("random_non_target_", seq_len(n_perm)), each = n_guides),
                           p.value = sample(object_non_target$p.value, size = n_guides * n_perm, replace = TRUE),
                           p.value_dep = sample(object_non_target$p.value_dep, size = n_guides * n_perm, replace = TRUE),
                           p.value_enrich = sample(object_non_target$p.value_enrich, size = n_guides * n_perm, replace = TRUE)
  )
  res_random <- compute_score_RRA(random_tab, alpha_thr = alpha_thr)
  gene_res <- mutate_(gene_res,
                      RRA_pvalue = ~map_dbl(RRA_score, ~((1 + sum(. >= res_random$RRA_score)) / (1 + n_perm))),
                      RRA_dep_pvalue = ~map_dbl(RRA_dep_score, ~((1 +sum(. >= res_random$RRA_dep_score)) / (1 + n_perm))),
                      RRA_enrich_pvalue = ~map_dbl(RRA_enrich_score, ~((1 +sum(. >= res_random$RRA_enrich_score)) / (1 + n_perm)))
  )
  gene_res$RRA_adjp <- p.adjust(gene_res$RRA_pvalue, method = "BH")
  gene_res$RRA_dep_adjp <- p.adjust(gene_res$RRA_dep_pvalue, method = "BH")
  gene_res$RRA_enrich_adjp <- p.adjust(gene_res$RRA_enrich_pvalue, method = "BH")
  return(gene_res)
}


process_res <- function(object, sgRNA_annot, fdr_method = "BH"){
  ## tidy results
  tab <- biobroom::tidy.MArrayLM(object)
  
  ## add unilateral pvalues
  tab$p.value_dep <- pt(tab$statistic, df = object$df.total[1], lower.tail = TRUE)
  tab$p.value_enrich <- pt(tab$statistic, df = object$df.total[1], lower.tail = FALSE)
  
  ## compute FDR
  tab$adj_p.value <- p.adjust(tab$p.value, method = fdr_method)
  tab$adj_p.value_dep <- p.adjust(tab$p.value_dep, method = fdr_method)
  tab$adj_p.value_enrich <- p.adjust(tab$p.value_enrich, method = fdr_method)
  tab <- tab[order(tab$adj_p.value),]
  
  tab <- left_join(tab, sgRNA_annot, by = c("gene" = "sgRNA"))
  tab <- rename_(tab, sgRNA = "gene")
  
  return(tab)
}


pickerSelectOptions <- function(choices, selected = NULL, choicesOpt = NULL, maxOptGroup = NULL) {
  print("This is the modified pickerSelectOptions")
  if (is.null(choicesOpt)) choicesOpt <- list()
  l <- sapply(choices, length)
  if (!is.null(maxOptGroup)) maxOptGroup <- rep_len(x = maxOptGroup, length.out = sum(l))
  cs <- cumsum(l)
  m <- matrix(data = c(1, cs[-length(l)] + 1, cs), ncol = 2)
  namechoice <- names(choices)
  tagList(lapply(1:length(choices), function(i) {
    label <- namechoice[i]
    choice <- choices[[i]]
    if (is.list(choice)) {
      optionTag <- list(
        label = htmlEscape(label, TRUE),
        pickerSelectOptions4(
          choice, selected,
          choicesOpt = lapply(
            X = choicesOpt,
            FUN = function(j) {
              j[m[i, 1]:m[i, 2]]
            }
          )
        )
      )
      if (!is.null(maxOptGroup))
        optionTag[["data-max-options"]] <- maxOptGroup[i]
      optionTag <- dropNulls1(optionTag)
      do.call(tags$optgroup, optionTag)
    } else {
      if (length(choicesOpt) == 0) {
        optionTag <- list(
          value = choice, if (is.null(label)) HTML(NULL) else HTML(htmlEscape(label)),
          selected = if (any(choice == selected)) "selected" else NULL
        )
      } else {
        optionTag <- list(
          value = choice, if (is.null(label)) HTML(NULL) else HTML(htmlEscape(label)),
          style = choicesOpt$style[i],
          `data-icon` = choicesOpt$icon[i],
          `data-subtext` = choicesOpt$subtext[i],
          `data-content` = choicesOpt$content[i],
          disabled = if (!is.null(choicesOpt$disabled[i]) && choicesOpt$disabled[i]) "disabled",
          selected = if (any(choice == selected)) "selected" else NULL
        )
      }
      # optionTag$attribs <- c(optionTag$attribs, list(if (choice %in% selected) " selected" else ""))
      optionTag <- dropNulls1(optionTag)
      do.call(tags$option, optionTag)
    }
  }))
}
run_enrich <- function(
  contrast,                # Comparison as specified in 'contrast' column in DE results
  DE_direction = "either", # 'either', 'up' (LFC>0), or 'down' (LFC<0)
  DE_res,                  # DE results df
  cat_map,                 # Functional category to gene mapping
  p_DE = 0.05,
  lfc_DE = 0,
  p_enrich = 0.05,
  q_enrich = 0.2
) {
  fcontrast <- contrast
  
  ## Filter the DE results, if needed: only take up- or downregulated
  if (DE_direction == "up") DE_res <- DE_res %>% filter(log2FoldChange > 0)
  if (DE_direction == "down") DE_res <- DE_res %>% filter(log2FoldChange < 0)
  
  ## Create a vector with DEGs
  DE_genes <- DE_res %>%
    filter(padj < p_DE,
           abs(log2FoldChange) > lfc_DE,
           contrast == fcontrast) %>%
    pull(gene_id)
  
  cat(fcontrast, " // DE Direction:", DE_direction, " // Nr DE genes: ", length(DE_genes))
  
  ## Run the enrichment analysis
  if (length(DE_genes) > 1) {
    enrich_res <- enricher(gene = DE_genes,
                           TERM2GENE = cat_map,
                           pAdjustMethod = "BH",
                           pvalueCutoff = 1,
                           qvalueCutoff = 1) %>%
      as.data.frame(.) %>%
      mutate(sig = ifelse(p.adjust < p_enrich & qvalue < q_enrich, 1, 0),
             contrast = fcontrast,
             DE_direction = DE_direction) %>%
      select(contrast,
             DE_direction,
             category = ID,
             numDEInCat = Count,
             GeneRatio,
             BgRatio,
             p = pvalue,
             padj = p.adjust,
             sig,
             description = Description,
             gene_ids = geneID)
    
    row.names(enrich_res) <- NULL
    
    cat(" // Nr enriched pathways:", sum(enrich_res$sig), "\n")
    return(enrich_res)
  } else {
    cat("\n")
  }
}

## Packages
if (!require("ggforce", quietly = TRUE)) install.packages("ggforce")
library(ggforce)

## GO analysis wrapper function
GO_wrap <- function(
  contrast_id,                       # A DE contrast as specified in the 'contrast' column in the 'DE_res' df
  DE_res,                            # Df with DE results from DESeq2
  GO_map,                            # Df with one GOterm-to-gene relation per row,
                                     # with columns 'gene_id' and 'go_term'
  gene_lens,                         # Df with gene lengths,
                                     # with columns 'gene_id' and 'length'
  DE_direction = "either",           # 'either' (= both together), 'up' (LFC>0), or 'down' (LFC>0)
  min_in_cat = 2, max_in_cat = Inf,  # Min. & max. nr of total terms in GO category 
  min_DE_in_cat = 2,                 # Min. nr. DE genes in GO term for a term to be significant
  p_DE = 0.05,                       # Adj. p-value threshold for DE significance
  lfc_DE = 0,                        # LFC threshold for DE significance
  ontologies = c("BP", "MF", "CC"),  # GO ontologies to consider
  filter_no_descrip = TRUE,          # Remove GO categories with no description
  rm_padj_na = TRUE,                 # Whether to remove genes with NA for `padj`
  ...
) {
  cat("\n-------------\nStarting analysis for contrast:", contrast_id, "\n")
  
  DE_vec <- get_DE_vec(contrast_id,
                       DE_res,
                       DE_direction = DE_direction,
                       rm_padj_na = rm_padj_na,
                       ...)
  
  GO_df <- GO_run(contrast_id,
                  DE_vec,
                  GO_map,
                  gene_lens,
                  min_in_cat = min_in_cat,
                  max_in_cat = max_in_cat,
                  ontologies = ontologies)
  
  return(GO_df)
}

## Function to run a GO analysis with goseq
## (Helper function, use GO_wrap to run the analysis)
GO_run <- function(contrast_id, DE_vec, GO_map, gene_lens,
                   DE_direction = "either",
                   min_in_cat = 2, max_in_cat = Inf,
                   min_DE_in_cat = 2,
                   ontologies = c("BP", "MF", "CC"),
                   filter_no_descrip = TRUE) {
  
  if (sum(DE_vec >= 2)) {
    ## Remove rows from gene length df not in the DE_vec
    fgene_lens <- gene_lens %>% filter(gene_id %in% names(DE_vec))
    if (nrow(fgene_lens) == 0) stop("Error: no rows left in gene length df, gene_id's likely don't match!")
    cat("- Nr genes removed from gene length df (no matching gene in DE vector):",
        nrow(gene_lens) - nrow(fgene_lens), "\n")
    
    ## Remove elements from DE_vec not among the gene lengths
    fDE_vec <- DE_vec[fgene_lens$gene_id %in% names(DE_vec)]
    if (length(fDE_vec) == 0) stop("Error: no entries left in DE vector, gene_id's likely don't match!")
    cat("- Nr genes removed from DE vector (no matching gene in gene length df):",
        length(DE_vec) - length(fDE_vec), "\n")
    cat("- Final length of DE vector:", length(fDE_vec), "\n")
    
    ## Check nr of genes with GO annotations
    ngenes_go <- length(intersect(GO_map$gene_id, names(fDE_vec)))
    cat("- Nr genes in DE vector with GO annotations:", ngenes_go, "\n")
    
    ## Check that gene lengths and contrast vector contain the same genes in the same order
    stopifnot(all(fgene_lens$gene_id == names(fDE_vec)))
    
    ## Probability weighting function based on gene lengths
    pwf <- nullp(DEgenes = fDE_vec,
                 bias.data = fgene_lens$length,
                 plot.fit = FALSE, genome = NULL, id = NULL)
    
    ## Run GO test
    GO_df <- goseq(pwf = pwf, gene2cat = GO_map, method = "Wallenius")
    
    ## Process GO results
    GO_df <- GO_df %>%
      filter(numDEInCat >= min_DE_in_cat,    # P-adjustment only for genes that were actually tested
             numInCat >= min_in_cat,         # Exclude very small categories
             numInCat <= max_in_cat,         # Exclude very large categories
             ontology %in% ontologies)       # Only select certain ontologies
    if (filter_no_descrip == TRUE) GO_df <- GO_df %>% filter(!is.na(term))
    
    GO_df <- GO_df %>%
      mutate(padj = p.adjust(over_represented_pvalue, method = "BH"),
             sig = ifelse(padj < 0.05 & numDEInCat >= min_DE_in_cat, 1, 0),
             contrast = contrast_id,
             DE_direction = DE_direction) %>%
      select(contrast, DE_direction,
             sig, p = over_represented_pvalue, padj,
             numDEInCat, numInCat,
             category, ontology, description = term)
    
    cat(contrast_id,
        "  DE dir.:", DE_direction,
        "  Nr DEG:", sum(DE_vec),
        "  Nr GO cat:", nrow(GO_df),
        "  Nr sig. (padj/p):", sum(GO_df$sig), "/", sum(GO_df$p < 0.05),
        "\n")
    
    return(GO_df)
    
  } else {
    cat("Contrast:", contrast_id, " -- 0 or 1 significant DE genes, skipping GO analysis\n")
  }
}

## Create named vector of DE genes (0s and 1s to indicate significance)
## for goseq analysis
## (Helper function, use GO_wrap to run the analysis)
get_DE_vec <- function(contrast_id,             # Focal comparison (contrast)
                       DE_res,                  # DE results df from DESeq2
                       DE_direction = "either", # either / both / up / down
                       rm_padj_na = TRUE,       # Whether to remove genes with NA for `padj`
                       p_DE = 0.05,             # padj threshold for DE
                       lfc_DE = 0) {            # LFC threshold for DE
  
  if (DE_direction == "up") DE_res <- DE_res %>% filter(log2FoldChange > 0)
  if (DE_direction == "down") DE_res <- DE_res %>% filter(log2FoldChange < 0)
  
  ## Create df for focal contrast, indicate which genes are significant
  fDE <- DE_res %>%
    filter(contrast == contrast_id) %>% # Select focal contrast
    mutate(sig = ifelse(padj < p_DE & abs(log2FoldChange) > lfc_DE, 1, 0)) %>%
    arrange(gene_id)
  cat("- Nr unique genes in DE results:", length(unique(fDE$gene_id)), "\n")
  
  # Exclude genes with NA adj-p-val - those were not tested
  if (rm_padj_na == TRUE) {
    fDE <- fDE %>% filter(!is.na(padj))
    cat("- Nr unique genes after removing those with NA for padj:",
        length(unique(fDE$gene_id)), "\n")
  } else {
    fDE <- fDE %>% mutate(sig = ifelse(is.na(sig), 0, 1))
  }
  
  DE_vec <- fDE$sig
  names(DE_vec) <- fDE$gene_id
  
  return(DE_vec)
}

## Function to plot the GO results
GO_plot <- function(GO_res, contrasts,
                    DE_directions = c("up", "down", "either"),
                    x_var = "contrast", facet_var = NULL,
                    label_count = TRUE, label_count_size = 2,
                    xlabs = NULL, ylabsize = 9,
                    title = NULL) {
  
  vars <- c("category", "ontology", "description",
            "contrast", "DE_direction", "padj")
  
  GO_sel <- GO_res %>%
    select(any_of(vars)) %>%
    filter(contrast %in% contrasts,
           DE_direction %in% DE_directions) %>%
    #select(any_of(vars)) %>%
    ## Pivot wider and then longer to include all terms in all contrasts
    pivot_wider(names_from = contrast, values_from = padj) %>%
    pivot_longer(cols = any_of(contrasts),
                 names_to = "contrast", values_to = "padj") %>%
    left_join(GO_res %>% select(contrast, category, DE_direction, numDEInCat, sig),
              by = c("contrast", "category", "DE_direction")) %>%
    ## No labels if not significant
    mutate(numDEInCat = ifelse(padj >= 0.05, NA, numDEInCat)) %>%
    mutate(contrast = sub("padj_", "", contrast),
           padj = ifelse(sig == FALSE, NA, padj),
           padj_log = -log10(padj)) %>% 
    ## Only take GO categories with at least one significant contrast
    filter(category %in% (filter(., sig == TRUE) %>% pull(category))) %>%
    ## Only take contrast with at least one significant category
    filter(contrast %in% (filter(., sig == TRUE) %>% pull(contrast))) %>%
    arrange(padj_log) %>% 
    mutate(description = paste0(category, " - ", description),
           description = str_trunc(description, width = 45),
           description = ifelse(is.na(description), category, description))
  #description = fct_inorder(description))
  
  p <- ggplot(GO_sel) +
    aes(x = .data[[x_var]],
        y = str_trunc(description, width = 40),
        fill = padj_log) +
    geom_tile(stat = "identity", size = 0.25, color = "grey80") +
    scale_fill_viridis_c(option = "D", na.value = "grey95") +
    scale_y_discrete(position = "right") +
    labs(fill = "-log10\n(adj. p)",
         title = title) +
    theme_minimal() +
    theme(legend.position = "left",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = ylabsize),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5))
  
  ## Add a count of the nr of DE genes in each GO category
  if (label_count == TRUE)
    p <- p + geom_label(aes(label = numDEInCat),
                        fill = "grey95", size = label_count_size)
  
  ## Faceting - ggforce::facet_col will keep tile heights constant
  if (is.null(facet_var)) {
    p <- p + facet_col(vars(ontology), scales = "free_y", space = "free")
  } else {
    p <- p + facet_grid(rows = vars(ontology),
                        cols = vars(!!sym(facet_var)),
                        scales = "free_y", space = "free_y",
                        switch = "y")
  }
  
  ## Add x-axis label
  if (!is.null(xlabs)) p <- p + scale_x_discrete(labels = xlabs)
  
  ## Print the final figure
  print(p)
}


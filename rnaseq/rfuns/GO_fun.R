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
  DE_direction = "either",           # 'either' (= both together), 'up', 'down', or 'both' (= both separately)
                                     # 'up' means LFC>0, 'down' means LFC<0
  min_in_cat = 2, max_in_cat = Inf,  # Min. & max. nr of total terms in GO category 
  min_DE_in_cat = 2,                 # Min. nr. DE genes in GO term for a term to be significant
  p_DE = 0.05,                       # Adj. p-value threshold for DE significance
  lfc_DE = 0,                        # LFC threshold for DE significance
  ontologies = c("BP", "MF", "CC"),  # GO ontologies to consider
  filter_no_descrip = TRUE,          # Remove GO categories with no description
  ...
  ) {
  
  if (DE_direction %in% c("either", "up", "down")) {
    DE_vec <- get_DE_vec(contrast_id, DE_res, DE_direction = DE_direction, ...)
    GO_df <- GO_run(contrast_id, DE_vec, GO_map, gene_lens,
                    min_in_cat = min_in_cat, max_in_cat = max_in_cat,
                    ontologies = ontologies)
    
  } else if (DE_direction == "both") {
    ## Up
    DE_vec_up <- get_DE_vec(contrast_id, DE_res, DE_direction = "up", ...)
    GO_df_up <- GO_run(contrast_id, DE_vec_up, GO_map, gene_lens,
                       DE_direction = "up",
                       min_in_cat = min_in_cat, max_in_cat = max_in_cat,
                       ontologies = ontologies)
    ## Down
    DE_vec_down <- get_DE_vec(contrast_id, DE_res, DE_direction = "down", ...)
    GO_df_down <- GO_run(contrast_id, DE_vec_down, GO_map, gene_lens,
                         DE_direction = "down",
                         min_in_cat = min_in_cat, max_in_cat = max_in_cat,
                         ontologies = ontologies)
    ## Combine
    GO_df <- bind_rows(GO_df_up, GO_df_down)
  } else {
    stop("Error: DE_direction should be one of 'either', 'both', 'up', or 'down'")
  }
}

## Create named vector of DE genes (0s and 1s to indicate significane)
get_DE_vec <- function(contrast_id,
                       DE_res,
                       DE_direction = "either",
                       p_DE = 0.05,
                       lfc_DE = 0) {
  
  if (DE_direction == "up") DE_res <- DE_res %>% filter(log2FoldChange > 0)
  if (DE_direction == "down") DE_res <- DE_res %>% filter(log2FoldChange < 0)
  
  fDE <- DE_res %>%
    filter(contrast == contrast_id,  # Select focal contrast
           !is.na(padj)) %>%         # Exclude genes with NA adj-p-val - those were not tested
    mutate(sig = ifelse(padj < p_DE & abs(log2FoldChange) > lfc_DE, 1, 0)) %>%
    arrange(gene_id)
  
  DE_vec <- fDE$sig
  names(DE_vec) <- fDE$gene_id
  
  return(DE_vec)
}

## Function to run a GO analysis with goseq
GO_run <- function(contrast_id, DE_vec, GO_map, gene_lens,
                   DE_direction = "either",
                   min_in_cat = 2, max_in_cat = Inf,
                   min_DE_in_cat = 2,
                   ontologies = c("BP", "MF", "CC"),
                   filter_no_descrip = TRUE) {
  
  if(sum(DE_vec > 0)) {
    cat("\n-------------\nStarting analysis for contrast:", contrast_id, "\n")
    
    ## Remove rows from gene length df not in the DE_vec
    fgene_lens <- gene_lens %>% filter(gene_id %in% names(DE_vec))
    if(nrow(fgene_lens) == 0) stop("Error: no rows left in gene length df, gene_id's likely don't match!")
    cat("Removed ", nrow(gene_lens) - nrow(fgene_lens), "rows in gene length df - no matching gene in DE vector\n")
    
    ## Remove elements from DE_vec not among the gene lengths
    fDE_vec <- DE_vec[names(DE_vec) %in% fgene_lens$gene_id]
    if(length(fDE_vec) == 0) stop("Error: no entries left in DE vector, gene_id's likely don't match!")
    cat("Removed ", length(DE_vec) - length(fDE_vec), "rows in DE vector - no matching gene in gene length df\n")
    cat("Final length of DE vector:", length(fDE_vec), "\n")
    
    ## Check that gene lengths and contrast vector contain the same genes in the same order
    stopifnot(all(fgene_lens$gene_id == names(fDE_vec)))
    
    ## Probability weighting function based on gene lengths
    pwf <- nullp(DEgenes = fDE_vec,
                 bias.data = fgene_lens$length,
                 plot.fit = FALSE)
    
    ## Run GO test
    GO_df <- goseq(pwf = pwf, gene2cat = GO_map, method = "Wallenius")
    
    ## Process GO results
    GO_df <- GO_df %>%
      filter(numDEInCat > 0,                 # P-adjustment only for genes that were actually tested
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
    
    cat("Contrast:", contrast_id,
        "    DE direction:", DE_direction,
        "    Nr GO cats:", nrow(GO_df),
        "    Nr DEGs:", sum(DE_vec),
        "    Nr sign. GO:", sum(GO_df$sig), "\n")
    
    return(GO_df)
  
  } else {
    cat("Contrast:", contrast_id, " -- No significant DE genes\n")
  }
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

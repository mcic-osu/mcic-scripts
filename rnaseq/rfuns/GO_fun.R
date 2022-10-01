## Packages
if (!require("ggforce", quietly = TRUE)) install.packages("ggforce")
library(ggforce)

## Create named vector of DE genes
get_DE_vec <- function(contrast_id, DE_res, DE_direction = "either",
                       p_DE = 0.05, lfc_DE = 0) {
  
  if (DE_direction == "up") DE_res <- DE_res %>% filter(log2FoldChange > 0)
  if (DE_direction == "down") DE_res <- DE_res %>% filter(log2FoldChange < 0)
  
  fDE <- DE_res %>%
    filter(contrast == contrast_id,
           !is.na(padj)) %>%   # Exclude genes with NA adj-p-val - not tested
    mutate(sig = ifelse(padj < p_DE & abs(log2FoldChange) > lfc_DE, 1, 0)) %>%
    arrange(gene_id)
  
  DE_vec <- fDE$sig
  names(DE_vec) <- fDE$gene_id
  
  return(DE_vec)
}

GO_run <- function(contrast_id, DE_vec, GO_map, gene_lens,
                   DE_direction = "either",
                   min_in_cat = 2, max_in_cat = Inf,
                   ontologies = c("BP", "MF", "CC")) {
  
  if(sum(DE_vec > 0)) {
    ## Remove rows from gene length df not in the DE_vec
    fgene_lens <- gene_lens %>% filter(gene_id %in% names(DE_vec))
    
    ## Remove elements from DE_vec not among the gene lengths
    fDE_vec <- DE_vec[names(DE_vec) %in% fgene_lens$gene_id]
    
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
             numInCat >= min_in_cat,
             numInCat <= max_in_cat,
             ontology %in% ontologies) %>%     # Exclude very large categories
      mutate(padj = p.adjust(over_represented_pvalue, method = "BH"),
             sig = ifelse(padj < 0.05 & numDEInCat > 1, 1, 0),
             contrast = contrast_id,
             DE_direction = DE_direction) %>%
      select(contrast, DE_direction,
             sig, p = over_represented_pvalue, padj,
             numDEInCat, numInCat,
             category, ontology, description = term)
    
    cat("Contrast:", contrast_id,
        "    DE direction:", DE_direction,
        "    Nr of rows:", nrow(GO_df),
        "    Nr DE:", sum(DE_vec),
        "    Nr sign. GO:", sum(GO_df$sig), "\n")
    
    GO_df <- GO_df %>% filter(!is.na(description))
    cat("Nr sign. after removing GO categories without a description: ",
        sum(GO_df$sig), "\n")
    
    return(GO_df)
  } else {
    cat("Contrast:", contrast_id, " -- No significant DE genes\n")
  }
}

GO_wrap <- function(contrast_id, DE_res, GO_map, gene_lens,
                    DE_direction = "either",
                    min_in_cat = 2, max_in_cat = Inf,
                    ontologies = c("BP", "MF", "CC"),
                    ...) {
  
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

GO_plot <- function(GO_res, contrasts,
                    DE_directions = c("up", "down", "either"),
                    x_var = "contrast", facet_var = NULL,
                    label_count = TRUE, label_count_size = 2,
                    xlabs = NULL, ylabsize = 9,
                    title = NULL) {
  
  # GO_sel <- GO_res %>%
  #   filter(contrast %in% contrasts,
  #          DE_direction %in% DE_directions) %>%
  #   group_by(category) %>%
  #   filter(sum(sig) > 0) %>%
  #   mutate(cat_descrip = paste0(category, " - ", description),
  #          padj_log = ifelse(sig == TRUE, -log10(padj), NA))
  
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

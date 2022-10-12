run_kegg <- function(fcontrast, direction = "both",
                     DE_res, kegg_map,
                     p_tres_DE = 0.05, lfc_tres_DE = 0,
                     p_tres_enrich = 0.05, q_tres = 0.2) {
  
  ## Filter the DE results, if needed: only take up- or downregulated
  if (direction == "up") DE_res <- DE_res %>% filter(log2FoldChange > 0)
  if (direction == "down") DE_res <- DE_res %>% filter(log2FoldChange < 0)
  
  ## Create a vector with DEGs
  DE_genes <- DE_res %>%
    filter(padj < p_tres_DE,
           abs(log2FoldChange) > lfc_tres_DE,
           contrast == fcontrast) %>%
    pull(gene_id)
  
  cat("## Contrast:", fcontrast,
      " // Direction:", direction,
      " // Nr DE genes: ", length(DE_genes))
  
  ## Run the enrichment analysis
  if (length(DE_genes) >  1) {
    kegg_res <- enricher(gene = DE_genes,
                         TERM2GENE = kegg_map,
                         pAdjustMethod = "BH",
                         pvalueCutoff = 1,
                         qvalueCutoff = 1) %>%
      as.data.frame(.) %>%
      mutate(sig = ifelse(p.adjust < p_tres_enrich & qvalue < q_tres, 1, 0),
             contrast = fcontrast,
             direction = direction)
    row.names(kegg_res) <- NULL
    cat(" // Nr enriched pathways:", sum(kegg_res$sig), "\n")
    return(kegg_res)
  } else cat("\n")
}

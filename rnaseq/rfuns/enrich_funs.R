# Packages
if (!require("ggforce", quietly = TRUE)) install.packages("ggforce")
library(ggforce)


# CLUSTERPROFILER FUNCTIONS ----------------------------------------------------
run_enrich <- function(
  contrast,                  # Comparison as specified in 'contrast' column in DE results
  DE_direction = "either",   # 'either', 'up' (LFC>0), or 'down' (LFC<0)
  DE_res,                    # DE results df
  cat_map,                   # Functional category to gene mapping with columns: 1:category, 2:gene_id, and optionally 3:description
  p_DE = 0.05,               # Adj. p-value threshold for DE
  lfc_DE = 0,                # LFC threshold for DE
  min_DE_in_cat = 2,         # Min. nr. DE genes in GO term for a term to be significant
  p_enrich = 0.05,           # Adj. p-value threshold for enrichment
  q_enrich = 0.2,            # Q value threshold for enrichment
  filter_no_descrip = TRUE,  # Remove pathways with no description
  return_df = FALSE          # Convert results object to a simple dataframe (tibble)
) {

  fcontrast <- contrast

  # Filter the DE results, if needed: only take over- or underexpressed
  if (DE_direction == "up") DE_res <- DE_res |> filter(log2FoldChange > 0)
  if (DE_direction == "down") DE_res <- DE_res |> filter(log2FoldChange < 0)

  # Create a vector with DEGs
  DE_res <- DE_res |> filter(padj < p_DE,
                             contrast == fcontrast)
  if (lfc_DE != 0) DE_res <- DE_res |> filter(abs(log2FoldChange) > lfc_DE)
  DE_genes <- DE_res$gene_id
  
  # Check if genes are present multiple times -- this would indicate there are multiple contrasts
  if (any(duplicated(DE_genes))) {
    stop("ERROR: Duplicated gene IDs detected -- there are probably multiple contrasts remaining in your input df")
  }
  
  # Report
  cat(fcontrast, " // DE Direction:", DE_direction,
      " // Nr DE genes: ", length(DE_genes))

  # Prep term mappings - if there's a third column, make a term2name df as well
  term2gene <- cat_map[, 1:2]
  if (ncol(cat_map) > 2) {
    term2name <- cat_map[, c(1, 3)]
    if (filter_no_descrip == TRUE) term2name <- term2name[!is.na(term2name[[2]]), ]
  } else {
    term2name <- NA
  }

  # Run the enrichment analysis
  if (length(DE_genes) <= 1) {
    cat("\n")
    return(NULL)
  }
    
  res <- enricher(gene = DE_genes,
                  TERM2GENE = term2gene,
                  TERM2NAME = term2name,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  if (return_df == FALSE) {
    res <- res |>
      filter(p.adjust < p_enrich,
             qvalue < q_enrich,
             Count >= min_DE_in_cat)
    cat(" // Nr enriched pathways:", nrow(res), "\n")
  }
  
  if (return_df == TRUE) {
    
    res <- as_tibble(res) |>
      mutate(sig = ifelse(p.adjust < p_enrich &
                            qvalue < q_enrich &
                            Count >= min_DE_in_cat,
                          TRUE, FALSE),
             contrast = fcontrast,
             DE_direction = DE_direction) |>
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
    
    cat(" // Nr enriched pathways:", sum(res$sig), "\n")
  }
  
  return(res)
}

# Function to run a Gene Set Enrichment Analysis (gsea)
run_gsea <- function(
    contrast,               # Contrast ID present in a column 'contrast'
    DE_res,                 # Differential expression level
    cat_map,                # Gene set category map (e.g., GO)
    p_enrich = 0.05,        # Adj. p-value threshold for enrichment
    return_df = FALSE       # Convert results object to a simple dataframe (tibble)
  ) {
  
  term2name <- NA
  fcontrast <- contrast
  
  # Prep the df to later create a gene vector
  gene_df <- DE_res |>
    filter(contrast == fcontrast,
           !is.na(log2FoldChange)) |>
    arrange(desc(log2FoldChange))
  
  # Check if genes are present multiple times -- this would indicate there are multiple contrasts
  if (any(duplicated(gene_df$gene_id))) {
    stop("ERROR: Duplicated gene IDs detected -- there are probably multiple contrasts remaining in your input df")
  }
  
  # Create a vector with lfc's and gene IDs
  gene_vec <- gene_df$log2FoldChange
  names(gene_vec) <- gene_df$gene_id
  
  # Report
  n_DE <- sum(gene_df$padj < 0.05, na.rm = TRUE)
  message("Contrast: ", fcontrast, " // Nr DE genes: ", n_DE)
  
  # Prep term mappings - if there's a third column, make a term2name df as well
  term2gene <- cat_map[, 1:2]
  if (ncol(cat_map) > 2) term2name <- cat_map[, c(1, 3)]
  
  # Run the enrichment analysis
  gsea_res <- GSEA(geneList = gene_vec,
                   TERM2GENE = term2gene,
                   TERM2NAME = term2name,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   verbose = FALSE)
  
  # Report
  n_sig <- sum(gsea_res$p.adjust < p_enrich)
  message("Nr enriched pathways: ", n_sig, "\n")
  
  # Return a df, if requested
  if (return_df == FALSE) {
    gsea_res <- gsea_res |> filter(p.adjust < p_enrich)
  
  } else {
    gsea_res <- as_tibble(gsea_res) |> 
      mutate(sig = ifelse(p.adjust < p_enrich, TRUE, FALSE),
             contrast = fcontrast) |>
      rename(category = ID)
  }

  return(gsea_res)
}


# GOSEQ PACKAGE FUNCTIONS ------------------------------------------------------
# This function will run the GO analysis 3 times for each contrast:
# for either DE direction, for LFC>0 DEGs only, and LFC<0 DEGs only.
# Then it will save the results in a TSV and print a summary.
run_GO_all <- function(contrasts, DE_res, GO_map, gene_lens, ...) {

  # Either DE direction
  GO_ei <- map_dfr(.x = contrasts, .f = run_GO,
                   DE_res = DE_res, gene_lens = gene_len_df, GO_map = GO_map,
                   DE_direction = "either", ...)

  # DE up (LFC>0)
  GO_up <- map_dfr(.x = contrasts, .f = run_GO,
                   DE_res = DE_res, gene_lens = gene_len_df, GO_map = GO_map,
                   DE_direction = "up", ...)

  # DE down (LFC<0)
  GO_dn <- map_dfr(.x = contrasts, .f = run_GO,
                   DE_res = DE_res, gene_lens = gene_len_df, GO_map = GO_map,
                   DE_direction = "down", ...)

  # Combine results for different DE directions
  GO_combined <- bind_rows(GO_ei, GO_up, GO_dn) |>
    mutate(DE_direction = factor(DE_direction, levels = c("either", "up", "down")))

  # Summarize the results
  smr <- GO_combined |>
    group_by(contrast, DE_direction) |>
    summarize(nsig = sum(sig), .groups = "drop") |>
    pivot_wider(id_cols = contrast, names_from = DE_direction, values_from = nsig)
  print(smr)

  return(GO_combined)
}

# GO analysis wrapper function
run_GO <- function(
  contrast,                          # A DE contrast as specified in the 'contrast' column in the 'DE_res' df
  DE_res,                            # Df with DE results from DESeq2
  GO_map,                            # Df with one GOterm-to-gene relation per row, w/ columns 'gene_id' and 'go_term'
  gene_lens,                         # Df with gene lengths w/ columns 'gene_id' and 'length'
  DE_direction = "either",           # 'either' (= both together), 'up' (LFC>0), or 'down' (LFC>0)
  min_in_cat = 2, max_in_cat = Inf,  # Min. & max. nr of total terms in GO category
  min_DE_in_cat = 2,                 # Min. nr. DE genes in GO term for a term to be significant
  p_DE = 0.05,                       # Adj. p-value threshold for DE significance
  lfc_DE = 0,                        # LFC threshold for DE significance
  use_sig_column = NULL,             # Specify column in DE results with TRUE/FALSE indicating DE significance
  ontologies = c("BP", "MF", "CC"),  # GO ontologies to consider
  filter_no_descrip = TRUE,          # Remove GO categories with no description
  rm_padj_na = TRUE,                 # Whether to remove genes with NA for `padj`
  verbose = FALSE,
  ...
) {

  if (verbose == TRUE) {
    cat("\n-------------\nStarting analysis for contrast:", contrast, "\n")
  }
  
  DE_vec <- get_DE_vec(
      contrast,
      DE_res,
      DE_direction = DE_direction,
      rm_padj_na = rm_padj_na,
      verbose = verbose,
      use_sig_column = use_sig_column,
      ...
  )

  GO_df <- run_GO_internal(
    contrast,
    DE_vec,
    GO_map,
    gene_lens,
    DE_direction,
    min_DE_in_cat = min_DE_in_cat,
    min_in_cat = min_in_cat,
    max_in_cat = max_in_cat,
    ontologies = ontologies,
    filter_no_descrip = filter_no_descrip,
    verbose = verbose
  )

  return(GO_df)
}

# Function to run a GO analysis with goseq
# (Helper function, use run_GO to run the analysis)
run_GO_internal <- function(
  contrast, DE_vec, GO_map, gene_lens,
  DE_direction = "either",
  min_in_cat = 2,
  max_in_cat = Inf,
  min_DE_in_cat = 2,
  ontologies = c("BP", "MF", "CC"),
  filter_no_descrip = TRUE,
  verbose = FALSE
  ) {

  if (verbose == TRUE) {
    message()
    message("run_GO_internal function...")
    message("GO analysis settings:")
    message("   - min_in_cat: ", min_in_cat, "  // max_in_cat: ", max_in_cat)
    message("   - min_DE_in_cat: ", min_DE_in_cat)
    message("   - Ontologies: ", ontologies)
    message("   - Filter GO cats with no description: ", filter_no_descrip)
  }

  if (sum(DE_vec) >= 2) {
    # Remove rows from gene length df not in the DE_vec
    fgene_lens <- gene_lens |> filter(gene_id %in% names(DE_vec))
    if (nrow(fgene_lens) == 0) stop("Error: no rows left in gene length df, gene_id's likely don't match!")
    
    n_removed <- nrow(gene_lens) - nrow(fgene_lens)
    if (verbose == TRUE) message("- Nr genes removed from gene length df (not in DE vector): ", n_removed)

    # Remove elements from DE_vec not among the gene lengths
    fDE_vec <- DE_vec[names(DE_vec) %in% fgene_lens$gene_id]
    if (length(fDE_vec) == 0) stop("Error: no entries left in DE vector, gene_id's likely don't match!")
    
    n_removed <- length(DE_vec) - length(fDE_vec)
    if (verbose == TRUE) message("- Nr genes removed from DE vector (not in gene length df): ", n_removed)
    if (verbose == TRUE) message("- Final length of DE vector: ", length(fDE_vec))

    # Check nr of genes with GO annotations
    ngenes_go <- length(intersect(GO_map$gene_id, names(fDE_vec)))
    if (verbose == TRUE) message("- Nr genes in DE vector with GO annotations: ", ngenes_go)

    # Check that gene lengths and contrast vector contain the same genes in the same order
    if (!all(fgene_lens$gene_id == names(fDE_vec))) {
      message("Gene IDs in gene length df do not match the gene IDs in the DE vector")
      message("IDs in gene length df:")
      print(head(fgene_lens$gene_id))
      message("IDs in DE vector:")
      print(head(names(fDE_vec)))
      stop()
    }

    # Probability weighting function based on gene lengths
    pwf <- nullp(DEgenes = fDE_vec,
                 bias.data = fgene_lens$length,
                 plot.fit = FALSE, genome = NULL, id = NULL)

    # Run GO test
    GO_df <- suppressMessages(
      goseq(pwf = pwf, gene2cat = GO_map, method = "Wallenius")
    )

    # Process GO results
    GO_df <- GO_df |> filter(numDEInCat >= min_DE_in_cat)
    GO_df <- GO_df |> filter(numInCat >= min_in_cat)
    GO_df <- GO_df |> filter(numInCat <= max_in_cat)
    GO_df <- GO_df |> filter(ontology %in% ontologies)
    if (filter_no_descrip == TRUE) GO_df <- GO_df |> filter(!is.na(term))

    GO_df <- GO_df |>
      mutate(padj = p.adjust(over_represented_pvalue, method = "BH"),
             sig = ifelse(padj < 0.05 & numDEInCat >= min_DE_in_cat, TRUE, FALSE),
             contrast = contrast,
             DE_direction = DE_direction) |>
      select(contrast, DE_direction,
             sig, p = over_represented_pvalue, padj,
             numDEInCat, numInCat,
             category, ontology, description = term)

    cat(contrast,
        "  DE dir.:", DE_direction,
        "  Nr DEG:", sum(DE_vec),
        "  Nr GO cat:", nrow(GO_df),
        "  Nr sig.:", sum(GO_df$sig),
        "\n")

    return(GO_df)

  } else {
    cat(contrast,
        "  DE dir.:", DE_direction,
        "  Nr DEG:", sum(DE_vec),
        " -- skipping GO analysis\n")
  }
}

# Create named vector of DE genes (0s and 1s to indicate significance) for goseq analysis
# (Helper function, use run_GO to run the analysis)
get_DE_vec <- function(
  contrast,                # Focal comparison (contrast)
  DE_res,                  # DE results df from DESeq2
  DE_direction = "either", # 'either' / 'up' / 'down'
  p_DE = 0.05,             # padj threshold for DE
  lfc_DE = 0,              # LFC threshold for DE
  use_sig_column = NULL,   # Use a column in the DE results with TRUE/FALSE indicating DE significance
  rm_padj_na = TRUE,       # Whether to remove genes with NA for `padj`
  verbose = FALSE
  ) {

  fcontrast <- contrast

  # If we use a column with precomputed DE significance, don't use thresholds
  if (!is.null(use_sig_column)) {
    if (verbose == TRUE) message("Using column ", use_sig_column, " to find DE genes")
    
    colnames(DE_res)[grep(use_sig_column, colnames(DE_res))] <- "isDE"
    p_DE <- NULL
    lfc_DE <- NULL
  }

  if (verbose == TRUE) {
    message("DE settings:")
    if (!is.null(p_DE)) message("    - DE p value: ", p_DE)
    if (!is.null(lfc_DE)) message("    - Min DE LFC: ", lfc_DE)
    message("    - Remove genes with NA as the adj. p-value: ", rm_padj_na)
  }

  # Create df for focal contrast
  fDE <- DE_res |>
    filter(contrast == fcontrast) |>
    arrange(gene_id)
  
  # Check if genes are present multiple times -- this would indicate there are multiple contrasts
  if (any(duplicated(fDE$gene_id))) {
    stop("ERROR: Duplicated gene IDs detected -- there are probably multiple contrasts remaining in your input df")
  }

  # Indicate which genes are significant
  if (is.null(use_sig_column)) {
    if (lfc_DE != 0) {
      fDE <- fDE |>
        mutate(isDE = ifelse(padj < p_DE & abs(log2FoldChange) > lfc_DE, TRUE, FALSE))
    } else {
      fDE <- fDE |>
        mutate(isDE = ifelse(padj < p_DE, TRUE, FALSE))
    }
  }

  # Subset to up/down DEGs if needed
  if (DE_direction == "up")
    fDE <- fDE |> mutate(isDE = ifelse(log2FoldChange > 0, isDE, FALSE))
  if (DE_direction == "down")
    fDE <- fDE |> mutate(isDE = ifelse(log2FoldChange < 0, isDE, FALSE))

  # Report nr of genes
  n_genes <- length(unique(fDE$gene_id))
  if (verbose == TRUE) message("- Nr unique genes in DE results: ", n_genes)

  if (rm_padj_na == TRUE) {
    # Exclude genes with NA adj-p-val - those were not tested
    fDE <- fDE |> filter(!is.na(padj))
    n_genes <- length(unique(fDE$gene_id))
    if (verbose == TRUE) message("- Nr unique genes after removing padj=NA: ", n_genes)
  } else {
    # Otherwise, turn NAs to FALSE (not DE)
    fDE <- fDE |> mutate(isDE = ifelse(is.na(isDE), FALSE, isDE))
  }

  DE_vec <- fDE$isDE
  names(DE_vec) <- fDE$gene_id
  
  if (verbose == TRUE) print(table(DE_vec))
    
  return(DE_vec)
}


# PLOTTING FUNCTIONS -----------------------------------------------------------
# Function to plot the GO results
enrich_plot <- function(
  enrich_res,                   # Enrichment results
  contrasts,                    # One or more contrasts
  DE_directions = c("up", "down", "either"),  # One or more DE directions
  plot_ontologies = TRUE,       # Show BP/CC/MF ontologies separately (should be FALSE for KEGG results)
  x_var = "contrast",           # What to plot along the x-axis (refer to column in `enrich_res`)
  x_var_levels = NULL,          # Factor levels for x_var (for ordering)
  facet_var = NULL,             # What to facet by (default: no faceting)
  label_count = TRUE,           # Print nrDEInCat genes in the box
  label_count_size = 2,         # Size of nrDEInCat label in the box
  padj_tres = 1,                # Further subset categories: only those with padj below this value
  n_tres = 0,                   # Further subset categories: only those with nrDEInCat at or above this value
  xlabs = NULL,
  ylab_size = 10,
  xlab_size = 14,
  xlab_angle = 0,
  plot_title = NULL,
  facet_labeller = "label_value",
  just_df = FALSE               # If true, don't make plot, just return modified df
  ) {

  # Columns to select
  vars <- c("category", "ontology", "description",
            "contrast", "DE_direction", "padj")

  # Prep the df
  enrich_res <- enrich_res |>
    #select(any_of(vars)) |>
    filter(contrast %in% contrasts,
           DE_direction %in% DE_directions) |>
    mutate(numDEInCat = ifelse(padj >= 0.05, NA, numDEInCat),
           contrast = sub("padj_", "", contrast),
           padj = ifelse(sig == FALSE, NA, padj),
           padj_log = -log10(padj)) |>
    # Only take GO categories with at least one significant contrast
    filter(category %in% (filter(., sig == TRUE) |> pull(category))) |>
    # Only take GO categories with at least one significant contrast
    filter(category %in% (filter(., padj < padj_tres & numDEInCat >= n_tres) |>
                            pull(category))) |>
    # Only take contrasts with at least one significant category
    filter(contrast %in% (filter(., sig == TRUE) |> pull(contrast))) |>
    mutate(description = paste0(category, " - ", description),
           description = str_trunc(description, width = 45),
           description = ifelse(is.na(description), category, description)) |>
    arrange(padj_log)

  # Make sure all combinations of contrast, DE_dir, and GO cat. are present
  # A) Make a lookup table with GO terms
  go_lookup <- enrich_res |>
    select(any_of(c("category", "ontology", "description"))) |>
    distinct()
  # B) Make a df with all possible combinations of contrast, DE_dir, and GO cat.
  enrich_rows <- enrich_res |>
    select(contrast, DE_direction, category) |>
    complete(contrast, DE_direction, category) |>
    left_join(go_lookup, by = c("category"))
  # C) Merge this with the enrich_res
  enrich_res <- left_join(enrich_rows,
                          enrich_res |> select(-(any_of(c("ontology", "description")))),
                          by = c("contrast", "DE_direction", "category"),
                          multiple = "all") |>
    mutate(sig = ifelse(is.na(sig), FALSE, sig))

  # If requested, don't make a plot and return the df
  if (just_df == TRUE) return(enrich_res)

  # Make sure the x-axis levels are ordered correctly
  if (!is.null(x_var_levels))
    enrich_res[[x_var]] <- factor(enrich_res[[x_var]], x_var_levels)

  # Legend title with subscript
  legend_title <- expression("-Log"[10]*" P")

  # Make the plot
  p <- ggplot(enrich_res) +
    aes(x = .data[[x_var]],
        y = str_trunc(description, width = 40),
        fill = padj_log) +
    geom_tile(stat = "identity", linewidth = 0.25, color = "grey80") +
    scale_fill_viridis_c(option = "D", na.value = "grey95") +
    scale_y_discrete(position = "right") +
    labs(fill = legend_title,
         title = plot_title) +
    theme_minimal() +
    theme(legend.position = "left",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = xlab_size, angle = xlab_angle),
          axis.text.y = element_text(size = ylab_size),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5))

  # Add a count of the nr of DE genes in each category
  if (label_count == TRUE)
    p <- p + geom_label(aes(label = numDEInCat),
                        fill = "grey95", size = label_count_size)

  # Faceting
  if (plot_ontologies == TRUE) {
    if (is.null(facet_var)) {
      # ggforce::facet_col will keep tile heights constant
      p <- p + facet_col(vars(ontology),
                         scales = "free_y", space = "free")
    } else {
      p <- p + facet_grid(rows = vars(ontology),
                          cols = vars(!!sym(facet_var)),
                          scales = "free_y",
                          space = "free_y",
                          switch = "y",
                          labeller = facet_labeller)
    }
  }

  if (plot_ontologies == FALSE) {
    if (!is.null(facet_var)) {
      p <- p + facet_grid(cols = vars(!!sym(facet_var)),
                          scales = "free_y",
                          space = "free_y",
                          switch = "y")
    }
  }

  # Add x-axis label
  if (!is.null(xlabs)) p <- p + scale_x_discrete(labels = xlabs)

  return(p)
}


# GO dotplot
GO_dotplot <- function(df, type = "GO") {

  if (type == "GO") group_by <- "ontology" else group_by <- NULL

  p <- ggdotchart(df,
                  x = "description", y = "padj_log",
                  color = "padj_log",
                  sorting = "descending",                       # Sort value in descending order
                  add = "segments",                             # Add segments from y = 0 to dots
                  rotate = TRUE,                                # Rotate vertically
                  #group = group_by,                             # Order by groups
                  dot.size = 5,                                 # Large dot size
                  label = "numDEInCat",                         # Add nr DE genes as dot labels
                  font.label = list(color = "white", size = 9, vjust = 0.5),
                  ggtheme = theme_bw()) +                       # ggplot2 theme
    labs(y = "-log10(adj. p-value)", x = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_color_viridis_c(option = "D", na.value = "grey95") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
          plot.title = element_text(size = 15, face = "bold"),
          strip.text.y = element_text(angle = 270, face = "bold"),
          strip.placement = "outside",
          axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          legend.position = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank())

  if (type == "GO") {
    p <- p + facet_grid(ontology~contrast, space = "free", scales = "free")
  } else if (type == "KEGG") {
    p <- p + facet_wrap(vars(contrast), scales = "free_x", nrow = 1)
  }

  print(p)
}


# KEGG DATABASE FUNCTIONS ------------------------------------------------------

# Function to get the description (technically: 'Name') of a KEGG pathway,
# given its pathway ID ('ko' or 'map' IDs).
# Needs tryCatch because come IDs fail (example of a failing pathway ID: "ko01130"),
# and the function needs to keep running.
# Use like: `pw_descrip <- map_dfr(.x = kegg_ids, .f = kegg_descrip)`
get_kegg_descrip <- function(kegg_pathway) {
  message(kegg_pathway)
  tryCatch( {
    description <- keggGet(kegg_pathway)[[1]]$NAME
    return(data.frame(kegg_pathway, description))
  }, error = function(cond) {
    message("keggGet failure")
    return(NULL)
  }
  )
}

# Get the genes belonging to a certain KEGG pathway
get_pw_genes <- function(pathway_id) {
  print(pathway_id)

  pw <- keggGet(pathway_id)
  if (is.null(pw[[1]]$GENE)) return(NA)
  pw2 <- pw[[1]]$GENE[c(TRUE, FALSE)]
  pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = TRUE),
                       function(x) x[1]))

  return(pw2)
}

# Function to get a KEGG pathway (ko-ID) associated with a KEGG K-term
get_pathway <- function(K_term, outdir) {
  # K_term <- "K13449"
  cat("K_term:", K_term, " ")

  tryCatch(
    {
      kegg_info <- keggGet(K_term)
      pathway_df <- data.frame(pathway_description = kegg_info[[1]]$PATHWAY) |>
        rownames_to_column("pathway_id") |>
        mutate("K_term" = K_term)

      cat(" Nr of pathways:", nrow(pathway_df), "\n")

      if(nrow(pathway_df) > 0) {
        pathway_df_file <- file.path(outdir, paste0(K_term, ".txt"))
        write_tsv(pathway_df, pathway_df_file)
        return(pathway_df)
      } else {
        return(NULL)
      }
    },
    error = function(cond) {
      message("keggGet failure")
      message(cond)
      return(NULL)
    }
  )
}

get_NCBI_id <- function(geneID) {
  geneID_NCBI <- entrez_search(db = "gene", term = geneID)$ids
  message(geneID, " - ", geneID_NCBI)
  if (is_empty(geneID_NCBI)) geneID_NCBI <- NA
  return(data.frame(geneID, geneID_NCBI))
}

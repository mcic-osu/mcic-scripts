# Packages
if (! "tidyverse" %in% installed.packages()) install.packages("tidyverse")
if (! "ggforce" %in% installed.packages()) install.packages("ggforce")
if (! "colorspace" %in% installed.packages()) install.packages("colorspace")
if (! "BiocManager" %in% installed.packages()) install.packages("BiocManager")
if (! "clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")

# CLUSTERPROFILER FUNCTIONS ----------------------------------------------------
# run_enrich: function to run a (GO or KEGG) standard overrepresentation analysis
run_enrich <- function(
  contrast,                  # DE comparison as specified in 'contrast' column in DE results
  DE_direction = "either",   # DE direction: 'either' (all DE genes), 'up' (lfc > 0), or 'down' (lfc < 0)
  DE_res,                    # DE results df from 'extract_DE()', should have columns 'gene', 'contrast', 'lfc', 'isDE', 'padj'
  cat_map,                   # Functional category/term to gene mapping with columns: 1:category, 2:gene ID, and optionally 3:description, 4: ontology
  p_enrich = 0.05,           # Adj. p-value threshold for enrichment
  q_enrich = 0.2,            # Q value threshold for enrichment
  min_DE_in_cat = 2,         # Nr. DE genes threshold for enrichment: at least this number of genes in the ontology category should be DE
                             # (Occasionally, 'small' categories with 1 DE gene can have p-values below 0.05 -- this excludes those)
  min_cat_size = 5,          # Min. nr. of genes in a category (= clusterProfiler 'minGSSize' argument - NOTE: clPr default is 10)
  max_cat_size = 500,        # Min. nr. of genes in a category (= clusterProfiler 'maxGSSize' argument)
  filter_no_descrip = TRUE,  # Remove categories/terms with no description (at least for GO terms, these tend to be old/deprecated ones)
  exclude_nontested = TRUE,  # Exclude genes that weren't tested for DE from the 'universe' of genes
  p_DE = NULL,               # Adj. p-value threshold for DE (default: use 'isDE' column to determine DE status)
  lfc_DE = NULL,             # LFC threshold for DE (default: use 'isDE' column to determine DE status)
  return_df = FALSE          # Convert results object to a simple dataframe (tibble), instead of keeping the ClusterProfiler object
                             # Should be FALSE if you want to use the enrichPlot functions directly 
) {

  # Filter DE results to only get those for the focal contrast
  fcontrast <- contrast
  DE_res <- DE_res |> filter(contrast == fcontrast)
  
  # Rename cat_map columns
  colnames(cat_map)[1:2] <- c("category", "gene")
  if(ncol(cat_map) > 2) colnames(cat_map)[3] <- "description"
  if(ncol(cat_map) > 3) colnames(cat_map)[4] <- "ontology"
  
  # Get the background 'universe' of genes:
  # genes that were tested for DE *and* occur in the ontology category map
  # (Excluding the latter is equivalent to goseq's 'use_genes_without_cat=FALSE',
  # and this is done by default by ClusterProfiler -- but non-tested genes *are* included)
  if (exclude_nontested == TRUE) {
    universe <- DE_res |>
      filter(!is.na(padj),
             gene %in% cat_map$gene) |>
      pull(gene)
  } else {
    universe <- NULL
  }
  
  # Prep term mappings - term-to-gene
  term2gene <- cat_map |> dplyr::select(category, gene)
  
  # Prep term mappings - term-to-name (description)
  if (ncol(cat_map) > 2) {
    term2name <- cat_map |> dplyr::select(category, description)
    
    if (filter_no_descrip == TRUE) {
      n_before <- length(unique(term2name$category))
      term2name <- term2name[!is.na(term2name$description), ]
      n_removed <- n_before - length(unique(term2name$category))
      if (n_removed > 0) message("Note: removed ", n_removed, " terms with no description")
    }
  } else {
    term2name <- NA
  }
  
  # Filter the DE results, if needed: only take over- or underexpressed
  if (DE_direction == "up") DE_res <- DE_res |> filter(lfc > 0)
  if (DE_direction == "down") DE_res <- DE_res |> filter(lfc < 0)

  # Create a vector with DEGs
  if (is.null(p_DE) & is.null(lfc_DE)) DE_res <- DE_res |> filter(isDE)
  if (!is.null(p_DE)) DE_res <- DE_res |> filter(padj < p_DE)
  if (!is.null(lfc_DE)) DE_res <- DE_res |> filter(abs(lfc) > lfc_DE)
  DE_genes <- DE_res$gene
  
  # Check if genes are present multiple times -- this would indicate there are multiple contrasts
  if (any(duplicated(DE_genes))) {
    stop("ERROR: Found duplicated gene IDs: you probably have multiple 'contrasts' in your input df")
  }
  
  # Check nr of DE genes in the category map
  genes_in_map <- DE_genes[DE_genes %in% term2gene$gene]
  cat("Contrast: ", fcontrast, " // DE direction: ", DE_direction,
      " // Nr DE genes (with category assigned): ", length(DE_genes),
      " (", length(genes_in_map), ")",
      sep = "")
  if (length(genes_in_map) == 0) {
    message("\nERROR: None of the DE genes are in the GO/KEGG category dataframe ('cat_map')")
    cat("First gene IDs from DE results: ", head(DE_genes), "\n")
    cat("First gene IDs from cat_map: ", head(term2gene$gene), "\n")
    stop()
  }
  
  # Skip the enrichment analysis if there are too few genes
  if (length(DE_genes) <= 1) {
    cat("Skipping enrichment: too few DE genes\n")
    return(NULL)
  }
  
  # Run the enrichment analysis
  res <- enricher(gene = DE_genes,
                  TERM2GENE = term2gene,
                  TERM2NAME = term2name,
                  universe = universe,
                  minGSSize = min_cat_size,
                  maxGSSize = max_cat_size,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  if (return_df == FALSE) {
    res <- res |>
      filter(p.adjust < p_enrich,
             qvalue < q_enrich,
             Count >= min_DE_in_cat)
    cat(" // Nr enriched:", nrow(res), "\n")
  } else {
    res <- as_tibble(res) |>
      mutate(sig = ifelse(p.adjust < p_enrich &
                            qvalue < q_enrich &
                            Count >= min_DE_in_cat,
                          TRUE,
                          FALSE),
             contrast = fcontrast,
             DE_direction = DE_direction) |>
      dplyr::select(contrast,
                    DE_direction,
                    category = ID,
                    n_DE_in_cat = Count,
                    GeneRatio,
                    BgRatio,
                    padj = p.adjust,
                    sig,
                    description = Description,
                    gene_ids = geneID)
    
    if ("ontology" %in% colnames(cat_map)) {
      res <- cat_map |>
        dplyr::select(category, ontology) |>
        distinct(category, .keep_all = TRUE) |> 
        right_join(res, by = "category") |>
        relocate(ontology, .before = "gene_ids")
    }
    
    # Add mean & median LFC value
    w_lfc <- res |>
      separate_longer_delim(cols = gene_ids, delim = "/") |>
      left_join(dplyr::select(DE_res, gene, lfc),
                by = join_by("gene_ids" == "gene"),
                relationship = "many-to-many") |>
      summarize(mean_lfc = mean(lfc),
                median_lfc = median(lfc),
                .by = c("category", "contrast", "DE_direction"))
    res <- left_join(res, w_lfc, by = c("category", "contrast", "DE_direction"))
    
    # Add gene numbers and enrichment ratio
    res <- res |>
    separate_wider_delim(cols = c("GeneRatio", "BgRatio"), delim = "/",
                         names_sep = "_") |>
      mutate(n_DE = as.integer(GeneRatio_2),
             n_cat = as.integer(BgRatio_1),
             n_total = as.integer(BgRatio_2)) |>
      select(-GeneRatio_1, -GeneRatio_2, -BgRatio_1, -BgRatio_2) |>
      mutate(fold_enrich = (n_DE_in_cat / n_DE) / (n_cat / n_total))
    
    # Report
    cat(" // Nr enriched:", sum(res$sig), "\n")
  }
  
  return(res)
}

# Function to run a Gene Set Enrichment Analysis (gsea)
run_gsea <- function(
    contrast,               # Contrast ID, should be one of the values in column 'contrast' in DE_res
    DE_res,                 # Differential expression results, should have columns 'contrast', 'gene', 'lfc'
    cat_map,                # Gene-to-category map (e.g., GO or KEGG)
    p_enrich = 0.05,        # Adj. p-value threshold for enrichment
    return_df = FALSE       # Convert results object to a simple dataframe (tibble)
  ) {
  
  term2name <- NA
  fcontrast <- contrast
  
  # Rename cat_map columns
  colnames(cat_map)[1:2] <- c("category", "gene")
  if(ncol(cat_map) > 2) colnames(cat_map)[3] <- "description"
  if(ncol(cat_map) > 3) colnames(cat_map)[4] <- "ontology"
  
  # Prep the df to later create a gene vector
  gene_df <- DE_res |>
    filter(contrast == fcontrast, !is.na(lfc)) |>
    arrange(desc(lfc))
  
  # Check if genes are present multiple times -- this would indicate there are multiple contrasts
  if (any(duplicated(gene_df$gene))) {
    stop("ERROR: Duplicated gene IDs detected -- there are probably multiple contrasts remaining in your input df")
  }
  
  # Create a vector with lfc's and gene IDs
  lfc_vec <- gene_df$lfc
  names(lfc_vec) <- gene_df$gene
  
  # Prep term mappings - if there's a third column, make a term2name df as well
  term2gene <- cat_map[, 1:2]
  if (ncol(cat_map) > 2) term2name <- cat_map[, c(1, 3)]
  
  # Report & check
  n_DE <- sum(gene_df$padj < 0.05, na.rm = TRUE)
  genes_in_map <- names(lfc_vec)[names(lfc_vec) %in% term2gene[[2]]]
  cat("Contrast: ", fcontrast, "// Nr DE genes (nr in cat_map):", n_DE, "(", length(genes_in_map), ")")
  if (length(genes_in_map) == 0) {
    message("\nERROR: None of the DE genes are in the cat_map dataframe")
    cat("First gene IDs from DE results: ", head(names(lfc_vec)), "\n")
    cat("First gene IDs from cat_map: ", head(term2gene[[2]]), "\n")
    stop()
  }
  
  # Run the enrichment analysis
  gsea_res <- GSEA(geneList = lfc_vec,
                   TERM2GENE = term2gene,
                   TERM2NAME = term2name,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   verbose = FALSE,
                   eps = 0,
                   seed = TRUE)
  
  # Report
  n_sig <- sum(gsea_res$p.adjust < p_enrich)
  cat(" // Nr enriched: ", n_sig, "\n")
  
  # Return a df, if requested
  if (return_df == FALSE) {
    gsea_res <- gsea_res |> filter(p.adjust < p_enrich)
  } else {
    gsea_res <- as_tibble(gsea_res) |> 
      mutate(sig = ifelse(p.adjust < p_enrich, TRUE, FALSE),
             contrast = fcontrast) |>
      dplyr::select(contrast,
                    category = ID,
                    padj = p.adjust,
                    sig,
                    description = Description,
                    gene_ids = core_enrichment)
    
    if ("ontology" %in% colnames(cat_map)) {
      gsea_res <- cat_map |>
        dplyr::select(category, ontology) |>
        distinct(category, .keep_all = TRUE) |> 
        right_join(gsea_res, by = "category") |>
        relocate(ontology, .before = "gene_ids")
    }
    
    # Add mean & median LFC value
    w_lfc <- gsea_res |>
      separate_longer_delim(cols = gene_ids, delim = "/") |>
      left_join(dplyr::select(DE_res, gene, lfc),
                by = join_by("gene_ids" == "gene"),
                relationship = "many-to-many") |>
      summarize(mean_lfc = mean(lfc),
                median_lfc = median(lfc),
                .by = c("category", "contrast"))
    
    gsea_res <- left_join(gsea_res, w_lfc, by = c("category", "contrast"))
  }

  return(gsea_res)
}


# PLOTTING FUNCTIONS -----------------------------------------------------------
# Function to plot the enrichment results in a heatmap format
enrichplot <- function(
  enrich_df,                    # Enrichment results
  contrasts = NULL,             # One or more contrasts (default: all)
  DE_directions = NULL,         # One or more DE directions (default: all)
  x_var = "contrast",           # Column in enrich_df to plot along the x-axis
  fill_var = "padj_log",        # Column in enrich_df to vary fill color by
                                # (the default, padj_log, will be computed from padj)
  facet_var1 = NULL,            # Column in enrich_df to facet by
  facet_var2 = NULL,            # Second column in enrich_df to facet by
                                # When specifiying both facet_var1 and var2: var1=>rows, var2=>columns
  countlab = TRUE,              # Print nrDEInCat genes in the box
  countlab_size = 2,            # Size of nrDEInCat label in the box
  xlabs = NULL,                 # Manually provide a vector with x-axis labels
  xlab_size = 13,               # Size of x-axis labels
  ylab_size = 10,               # Size of y-axis labels (= categories)
  xlab_angle = 0,               # Angle of x-axis labels
  facet_labeller = "label_value", # Facet labelling
  add_cat_id = TRUE,             # Add ontology category ID to its name
  merge_directions = FALSE
) {
  
  # Check
  if (is.null(facet_var1) && !is.null(facet_var2)) {
    stop("ERROR: Only use facet_var2 when also using facet_var1")
  }
    
  # Select contrasts & DE directions
  if (is.null(contrasts)) contrasts <- unique(enrich_df$contrast)
  if (is.null(DE_directions)) DE_directions <- unique(enrich_df$DE_direction)
  
  # Prep the df
  enrich_df <- enrich_df |>
    filter(contrast %in% contrasts,
           DE_direction %in% DE_directions) |>
    mutate(n_DE_in_cat = ifelse(padj >= 0.05, NA, n_DE_in_cat),
           contrast = sub("padj_", "", contrast),
           fold_enrich = ifelse(sig == FALSE, NA, fold_enrich),
           mean_lfc = ifelse(sig == FALSE, NA, mean_lfc),
           median_lfc = ifelse(sig == FALSE, NA, median_lfc),
           n_DE = ifelse(sig == FALSE, NA, n_DE),
           padj = ifelse(sig == FALSE, NA, padj),
           padj_log = -log10(padj)) %>%
    # Only take GO categories with at least one significant contrast (as pre-specified in 'sig' column)
    filter(category %in% (filter(., sig == TRUE) |> pull(category))) %>%
    # Only take contrasts with at least one significant category
    filter(contrast %in% (filter(., sig == TRUE) |> pull(contrast))) |>
    arrange(padj_log)
  
  # Modify the ontology category description
  trunc_width <- 40
  enrich_df <- enrich_df |>
    mutate(
      # Capitalize the first letter
      description = paste0(toupper(substr(description, 1, 1)),
                           substr(description, 2, nchar(description))),
      # If there is no description, use the category ID
      description = ifelse(is.na(description), category, description)
    )
  if (add_cat_id) {
    enrich_df <- enrich_df |>
      mutate(description = paste0(category, " - ", description))
    trunc_width <- 50
  }
  enrich_df <- enrich_df |>
    mutate(description = str_trunc(description, width = trunc_width))
  
  # Make sure all combinations of contrast, DE_dir, and GO cat. are present
  # A) Make a lookup table with GO terms
  go_lookup <- enrich_df |>
    dplyr::select(any_of(c("category", "ontology", "description"))) |>
    distinct()
  # B) Make a df with all possible combinations of contrast, DE_dir, and GO cat.
  enrich_rows <- enrich_df |>
    dplyr::select(contrast, DE_direction, category) |>
    complete(contrast, DE_direction, category) |>
    left_join(go_lookup, by = c("category"))
  # C) Merge this with the enrich_df
  enrich_df <- left_join(enrich_rows,
                         enrich_df |> dplyr::select(-(any_of(c("ontology", "description")))),
                         by = c("contrast", "DE_direction", "category"),
                         multiple = "all") |>
    mutate(sig = ifelse(is.na(sig), FALSE, sig))
  
  # Merge across DE directions
  if (! "DE_direction" %in% c(x_var, facet_var1, facet_var2)) {
    message("Merging DE directions (showing only most significant)
             because DE direction is not included in any plotting variable")
    merge_directions <- TRUE
  }
  if (merge_directions) {
    enrich_df <- enrich_df |>
      group_by(category, contrast) |>
      arrange(padj, .by_group = TRUE) |>
      slice_head(n = 1)
  }
  
  # Legend title with subscript
  if (fill_var == "padj_log") {
    fill_name <- expression("-Log"[10]*" P")
  } else if (fill_var == "median_lfc") {
    fill_name <- "Median\nLFC"
  } else if (fill_var == "mean_lfc") {
    fill_name <- "Mean\nLFC"
  } else {
    fill_name <- fill_var
  }
  
  # X-label position
  if (is.null(facet_var1)) xlab_pos <- "top" else xlab_pos <- "bottom" 
  
  # Color scale
  if (fill_var %in% c("mean_lfc", "median_lfc")) {
    col_scale <- colorspace::scale_fill_continuous_divergingx(
      palette = "Tropic", mid = 0.0, na.value = "grey98", rev = TRUE
      #https://stackoverflow.com/questions/58718527/setting-midpoint-for-continuous-diverging-color-scale-on-a-heatmap
      )
  } else {
    col_scale <- scale_fill_viridis_c(option = "D", na.value = "grey97")
  }
  
  # Make the plot
  p <- ggplot(enrich_df) +
    aes(x = .data[[x_var]],
        y = description,
        fill = .data[[fill_var]]) +
    geom_tile(stat = "identity", linewidth = 0.25, color = "grey80") +
    col_scale +
    scale_x_discrete(position = xlab_pos) +
    scale_y_discrete(position = "right") +
    labs(fill = fill_name) +
    theme_minimal() +
    theme(legend.position = "left",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = xlab_size, angle = xlab_angle),
          axis.text.y = element_text(size = ylab_size),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5))
  
  # Add x-axis label
  if (!is.null(xlabs)) p <- p + scale_x_discrete(labels = xlabs)
  
  # Add a count of the nr of DE genes in each category
  if (countlab == TRUE) {
    p <- p + suppressWarnings(
      geom_label(aes(label = n_DE_in_cat), fill = "grey95", size = countlab_size)
      )
  }
  
  # Faceting
  if (!is.null(facet_var1)) {
    if (!is.null(facet_var2)) {
      # Facet_grid when there are 2 facet vars
      p <- p + facet_grid(rows = vars(.data[[facet_var1]]),
                          cols = vars(.data[[facet_var2]]),
                          scales = "free_y",
                          space = "free_y",
                          switch = "y",
                          labeller = facet_labeller)
    } else {
      if (facet_var1 != "ontology") {
        p <- p +
          ggforce::facet_wrap(facets = vars(.data[[facet_var1]]),
                              scales = "fixed", nrow = 1)
      } else {
        p <- p +
          ggforce::facet_col(facets = vars(.data[[facet_var1]]),
                             scales = "free_y", space = "free")
      }
    }
  }
  
  return(p)
}


# Cleveland dotplot of enrichment results
cdotplot <- function(
    enrich_df,               # Dataframe with enrichment results
    contrasts = NULL,        # One or more contrasts (default: all)
    DE_directions = NULL,    # One or more DE directions (default: all)
    x_var = "padj_log",      # Column in enrich_df to plot along the x axis ('padj_log' will be computed from 'padj')
    fill_var = "median_lfc", # Column in enrich_df to vary fill color by ('padj_log' will be computed from 'padj')
    label_var = "n_DE_in_cat", # Column in enrich_df with a number to add as a label in the circles
    facet_var1 = NULL,       # Column in enrich_df to facet by
    facet_var2 = NULL,       # Second column in enrich_df to facet by (e.g. 'ontology' for GO)
                             # (will use facet_grid())
    facet_to_columns = TRUE, # When only using one facet_var1, facets are columns (or rows)
    x_title = NULL,          # X-axis title
    ylab_size = 11,          # Size of y-axis labels (= category labels)
    add_cat_id = TRUE,       # Add category ID (e.g., 'GO:009539') to its description
    point_size = 6
) {
  
  # Select contrasts & DE directions
  if (is.null(contrasts)) contrasts <- unique(enrich_df$contrast)
  
  # Prep the df
  enrich_df <- enrich_df |>
    filter(sig == TRUE,
           contrast %in% contrasts) |>
    mutate(padj_log = -log10(padj))
  
  if (! is.null(DE_directions))
    enrich_df <- enrich_df |> filter(DE_direction %in% DE_directions)
  
  # Modify the category description
  enrich_df <- enrich_df |>
    mutate(
      # Capitalize the first letter
      description = paste0(toupper(substr(description, 1, 1)),
                           substr(description, 2, nchar(description))),
      # If there is no description, use the category ID
      description = ifelse(is.na(description), category, description)
    )
  trunc_width <- 40
  if (add_cat_id) {
    enrich_df <- enrich_df |>
      mutate(description = paste0(category, " - ", description))
    trunc_width <- 50
  }
  enrich_df <- enrich_df |>
    mutate(description = str_trunc(description, width = trunc_width))
  
  # Legend position and title
  if (x_var == fill_var) legend_pos <- "none" else legend_pos <- "top"
  if (fill_var == "median_lfc") {
    color_name <- "Median\nLFC"
  } else if (fill_var == "mean_lfc") {
    color_name <- "Mean\nLFC"
  } else if (fill_var == "padj_log") {
    color_name <- expression("-Log"[10]*" P")
  } else {
    color_name <- fill_var
  }
  
  # X-axis title
  # Log-transform the p-value
  if (x_var == "padj_log") {
    x_title <- expression("-Log"[10]*" P")
  } else if (x_var == "padj") {
    x_title <- "Adjusted p-value"
  } else if (x_var == "fold_enrich") {
    x_title <- "Fold enrichment"
  } else if (x_var == "median_lfc") {
    x_title <- "Median LFC"
  } else if (fill_var == "mean_lfc") {
    x_title <- "Mean LFC"
  } 
  
  # Color scale
  if (fill_var %in% c("mean_lfc", "median_lfc")) {
    #https://carto.com/carto-colors/
    col_scale <- colorspace::scale_color_continuous_divergingx(
      palette = "Tropic", mid = 0.0, na.value = "grey97",
      name = color_name, rev = TRUE
    )
  } else {
    col_scale <- scale_color_viridis_c(
      option = "D", na.value = "grey95", name = color_name,
      )
  }
  
  # Create the base plot
  p <- ggpubr::ggdotchart(
    enrich_df,
    x = "description",
    y = x_var,
    label = label_var,
    color = fill_var,
    sorting = "descending",       # Sort value in descending order
    add = "segments",             # Add segments from y = 0 to dots
    rotate = TRUE,                # Rotate vertically
    dot.size = point_size,
    font.label = list(color = "white", size = 9, vjust = 0.5),
    ggtheme = theme_bw()
  )
  
  # Formatting
  p <- p +
    labs(x = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.075, 0.075))) +
    col_scale +
    theme(legend.position = legend_pos,
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
          strip.text.y = element_text(angle = 270, face = "bold"),
          strip.placement = "outside",
          axis.title.x = element_text(
            size = 12,
            margin = margin(t = 0.5, b = 0.5, unit = "cm")
          ),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = ylab_size),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank())
  
  # Other formatting
  if (! is.null(x_title)) p <- p + labs(y = x_title)
  if (x_var %in% c("mean_lfc", "median_lfc")) {
    p <- p + geom_hline(yintercept = 0, color = "grey70", linewidth = 1)
  }
  
  # Faceting
  if (!is.null(facet_var1)) {
    if (!is.null(facet_var2)) {
      p <- p +
        facet_grid(rows = vars(.data[[facet_var1]]),
                   cols = vars(.data[[facet_var2]]),
                   space = "free_y", scales = "free_y")
    } else if (facet_to_columns == FALSE ) {
      p <- p +
        ggforce::facet_col(facets = vars(.data[[facet_var1]]),
                           scales = "free_y",
                           space = "free")
    } else {
      # Default: facet into columns with facet_row()
      p <- p +
        ggforce::facet_row(facets = vars(.data[[facet_var1]]),
                           scales = "free_x",
                           space = "free")
    }
  }
  
  return(p)
}

# Heatmap for DE genes in significant GO terms
# NOTE: This function requires (calls) the pheat() function from mcic-scripts/rnaseq/DE_funs.R
GO_pheat <- function(
    GO_cat,                    # GO category to plot
    GO_res,                    # GO results from run_enrich() or run_gsea() with return_df=TRUE
    count_mat,                 # Normalized GO matrix
    meta_df,                   # Metadata df
    annot_df = NULL,           # Annotation df: gene IDs as rownames, gene names/descriptions as the single column
    contrast = NULL,           # DE contrast as specified in 'contrast' column in GO_res 
    DE_direction = "either",   # Direction of DE ('either', 'up', or 'down')
    ...                        # Other args to pass to pheat()
    ) {
  
  # Rename before filtering
  fcontrast <- contrast
  DE_dir <- DE_direction
  
  # Filter the GO results
  fgo <- GO_res |> filter(category == GO_cat)
  if (!is.null(DE_direction)) fgo <- fgo |> filter(DE_direction == DE_dir)
  if (!is.null(fcontrast)) fgo <- fgo |> filter(contrast == fcontrast)
  
  # Get a vector with gene IDs
  fgenes <- fgo |> separate_longer_delim(cols = gene_ids, delim = "/") |> pull(gene_ids)
  
  # Report
  message("GO category: ", GO_cat, " (", length(fgenes), " DE genes)")
  
  # Prepare the title
  descrip <- fgo$description[1]
  p <- formatC(fgo$padj, format = "e", digits = 1)
  title <- paste0(GO_cat, " (", descrip, ")\n(DEGs in ", fcontrast, ", enrich-p=", p, ")")
  
  # Make the heatmap
  p <- pheat(genes = fgenes,
             count_mat = count_mat,
             meta_df = meta_df,
             annot_df = annot_df,
             nchar_gene = 60,
             main = title,
             ...)
  print(p)
}


# GOSEQ PACKAGE FUNCTIONS ------------------------------------------------------
# This function will run the GO analysis 3 times for each contrast:
# for either DE direction, for LFC>0 DEGs only, and LFC<0 DEGs only.
# Then it will save the results in a TSV and print a summary.
rgoseq_all <- function(contrasts, DE_res, GO_map, gene_lens, ...) {

  # Either DE direction
  GO_ei <- map_dfr(.x = contrasts, .f = rgoseq,
                   DE_res = DE_res, gene_lens = gene_len_df, GO_map = GO_map,
                   DE_direction = "either", ...)

  # DE up (LFC>0)
  GO_up <- map_dfr(.x = contrasts, .f = rgoseq,
                   DE_res = DE_res, gene_lens = gene_len_df, GO_map = GO_map,
                   DE_direction = "up", ...)

  # DE down (LFC<0)
  GO_dn <- map_dfr(.x = contrasts, .f = rgoseq,
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
rgoseq <- function(
  contrast,                          # A DE contrast as specified in the 'contrast' column in the 'DE_res' df
  DE_res,                            # Df with DE results from DESeq2
  GO_map,                            # Df with one GOterm-to-gene relation per row, w/ columns 'gene' and 'go_term'
  gene_lens,                         # Df with gene lengths w/ columns 'gene' and 'length'
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

  GO_df <- rgoseq_internal(
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
# (Helper function, use rgoseq to run the analysis)
rgoseq_internal <- function(
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
    message("rgoseq_internal function...")
    message("GO analysis settings:")
    message("   - min_in_cat: ", min_in_cat, "  // max_in_cat: ", max_in_cat)
    message("   - min_DE_in_cat: ", min_DE_in_cat)
    message("   - Ontologies: ", ontologies)
    message("   - Filter GO cats with no description: ", filter_no_descrip)
  }

  if (sum(DE_vec) >= 2) {
    # Remove rows from gene length df not in the DE_vec
    fgene_lens <- gene_lens |> filter(gene %in% names(DE_vec))
    if (nrow(fgene_lens) == 0) stop("Error: no rows left in gene length df, gene IDs likely don't match!")
    
    n_removed <- nrow(gene_lens) - nrow(fgene_lens)
    if (verbose == TRUE) message("- Nr genes removed from gene length df (not in DE vector): ", n_removed)

    # Remove elements from DE_vec not among the gene lengths
    fDE_vec <- DE_vec[names(DE_vec) %in% fgene_lens$gene]
    if (length(fDE_vec) == 0) stop("Error: no entries left in DE vector, gene IDs likely don't match!")
    
    n_removed <- length(DE_vec) - length(fDE_vec)
    if (verbose == TRUE) message("- Nr genes removed from DE vector (not in gene length df): ", n_removed)
    if (verbose == TRUE) message("- Final length of DE vector: ", length(fDE_vec))

    # Check nr of genes with GO annotations
    ngenes_go <- length(intersect(GO_map$gene, names(fDE_vec)))
    if (verbose == TRUE) message("- Nr genes in DE vector with GO annotations: ", ngenes_go)

    # Check that gene lengths and contrast vector contain the same genes in the same order
    if (!all(fgene_lens$gene == names(fDE_vec))) {
      message("Gene IDs in gene length df do not match the gene IDs in the DE vector")
      message("IDs in gene length df:")
      print(head(fgene_lens$gene))
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
      goseq(pwf = pwf,
            gene2cat = GO_map,
            method = "Wallenius",
            use_genes_without_cat = FALSE)
    )

    # Process GO results
    GO_df <- GO_df |> filter(n_DE_in_cat >= min_DE_in_cat)
    GO_df <- GO_df |> filter(numInCat >= min_in_cat)
    GO_df <- GO_df |> filter(numInCat <= max_in_cat)
    GO_df <- GO_df |> filter(ontology %in% ontologies)
    if (filter_no_descrip == TRUE) GO_df <- GO_df |> filter(!is.na(term))

    GO_df <- GO_df |>
      mutate(padj = p.adjust(over_represented_pvalue, method = "BH"),
             sig = ifelse(padj < 0.05 & n_DE_in_cat >= min_DE_in_cat, TRUE, FALSE),
             contrast = contrast,
             DE_direction = DE_direction) |>
      dplyr::select(contrast, DE_direction,
             sig, p = over_represented_pvalue, padj,
             n_DE_in_cat, numInCat,
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
# (Helper function, use rgoseq to run the analysis)
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
  fDE <- DE_res |> filter(contrast == fcontrast) |> arrange(gene)
  
  # Check if genes are present multiple times -- this would indicate there are multiple contrasts
  if (any(duplicated(fDE$gene))) {
    stop("ERROR: Duplicated gene IDs detected -- there are probably multiple contrasts remaining in your input df")
  }

  # Indicate which genes are significant
  if (is.null(use_sig_column)) {
    if (lfc_DE != 0) {
      fDE <- fDE |>
        mutate(isDE = ifelse(padj < p_DE & abs(lfc) > lfc_DE, TRUE, FALSE))
    } else {
      fDE <- fDE |>
        mutate(isDE = ifelse(padj < p_DE, TRUE, FALSE))
    }
  }

  # Subset to up/down DEGs if needed
  if (DE_direction == "up")
    fDE <- fDE |> mutate(isDE = ifelse(lfc > 0, isDE, FALSE))
  if (DE_direction == "down")
    fDE <- fDE |> mutate(isDE = ifelse(lfc < 0, isDE, FALSE))

  # Report nr of genes
  n_genes <- length(unique(fDE$gene))
  if (verbose == TRUE) message("- Nr unique genes in DE results: ", n_genes)

  if (rm_padj_na == TRUE) {
    # Exclude genes with NA adj-p-val - those were not tested
    fDE <- fDE |> filter(!is.na(padj))
    n_genes <- length(unique(fDE$gene))
    if (verbose == TRUE) message("- Nr unique genes after removing padj=NA: ", n_genes)
  } else {
    # Otherwise, turn NAs to FALSE (not DE)
    fDE <- fDE |> mutate(isDE = ifelse(is.na(isDE), FALSE, isDE))
  }

  DE_vec <- fDE$isDE
  names(DE_vec) <- fDE$gene
  
  if (verbose == TRUE) print(table(DE_vec))
  return(DE_vec)
}


# KEGG DATABASE FUNCTIONS ------------------------------------------------------
# Function to get the description (technically: 'Name') of a KEGG pathway,
# given its pathway ID ('ko' or 'map' IDs).
# Needs tryCatch because some IDs fail (example of a failing pathway ID: "ko01130"),
# and the function needs to keep running.
# Use like so: `pw_descrip <- map_dfr(.x = kegg_ids, .f = kegg_descrip)`
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

# Get the genes belonging to a certain KEGG pathway
get_pw_kos <- function(pathway_id) {
  cat(pathway_id, "... ")
  
  pw <- keggGet(pathway_id)
  cat(pw[[1]]$NAME[[1]])
  
  if (is.null(pw[[1]]$GENE)) {
    cat(" - No KOs found!\n")
    ko_df <- tibble(pathway = pathway_id,
                    pathway_nr = sub("[[:alpha:]]+", "", pathway_id),
                    KO = NA)
  } else {
    kos <- pw[[1]]$GENE[c(FALSE, TRUE)]
    kos <- unique(sub(".*KO:(K\\d+).*", "\\1", x = kos))
    cat(" - Nr of KOs:", length(kos), "\n")
    
    ko_df <- tibble(pathway = pathway_id,
                    pathway_nr = sub("[[:alpha:]]+", "", pathway_id),
                    KO = kos)
  }
  
  return(ko_df)
}

# Function to get a KEGG pathway (ko-ID) associated with a KEGG K-term
# (Example K_term: "K13449")
get_pathway <- function(K_term, outdir) {
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

# Get the KO numbers associated with a module,
# straight from the KEGG database
get_kos_in_mod <- function(idx, modules) {
  module <- names(modules)[idx]
  mod_vec <- keggGet(module)[[1]]$ORTHOLOGY
  mod_df <- data.frame(KO_id = names(mod_vec), KO_descrip = mod_vec,
                       module = module, module_descrip = modules[idx],
                       row.names = NULL)
  return(mod_df)
}

# Get NCBI ID for a gene ID
get_NCBI_id <- function(geneID) {
  geneID_NCBI <- entrez_search(db = "gene", term = geneID)$ids
  message(geneID, " - ", geneID_NCBI)
  if (is_empty(geneID_NCBI)) geneID_NCBI <- NA
  return(data.frame(geneID, geneID_NCBI))
}

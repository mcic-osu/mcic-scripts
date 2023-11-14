# Packages
if (!require(janitor)) install.packages("janitor")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggiraph)) install.packages("ggiraph")
if (!require(ggrepel)) install.packages("ggrepel")
if (!require(pheatmap)) install.packages("pheatmap")
if (!require(ggforce)) install.packages("ggforce")
if (!require(ggiraph)) install.packages("ggiraph")
if (!require(patchwork)) install.packages("patchwork")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(PCAtools)) BiocManager::install("PCAtools")
if (!require(ggvenn)) BiocManager::install("ggvenn")
library(tidyverse)

# Function to run a DE analysis with DESeq2
run_DE <- function(
  dds,
  design = NULL,
  subset_factor = NULL,
  subset_levels = NULL,
  extract_factor = NULL,
  count_df = NULL,
  minReplicatesForReplace = 7,
  ...
) {
  nsample_org <- ncol(dds)
  
  if (!is.null(design)) design(dds) <- design
  
  if (!is.null(subset_factor)) {
    dds <- dds[, dds[[subset_factor]] %in% subset_levels]
    message("\nAfter subsetting ", subset_factor, " to keep ", subset_levels, " only, ",
            ncol(dds), " out of ", nsample_org, " samples are left")
  }
  
  dds <- suppressMessages(
    DESeq(dds,
          minReplicatesForReplace = minReplicatesForReplace,
          ...)
  )
  
  if (! is.null(extract_factor)) {
    DE_res <- extract_DE_all(dds, fac = extract_factor, count_df = count_df)
    return(DE_res)
  } else {
    return(dds)
  }
}


# Extract DE results from a DESeq2 object
# To specify the contrast, either use 'fct=' and 'comp=', or 'contrasts='
extract_DE <- function(
    dds,                    # DESeq2 object
    fct = NULL,             # Focal factor (column in metadata with the levels in 'comp')
    comp = NULL,            # Vector of 2 with focal levels of factor 'fac'
    contrasts = NULL,       # Vector of 1 or more contrasts as listed in resultsNames(dds)
    sig_only = FALSE,       # Only return significant results
    p_tres = 0.05,          # Adj. p-value threshold for DE significance
    lfc_tres = 0,           # Log2-fold change threshold for DE significance
    mean_tres = 0,          # Mean expr. level threshold for DE significance
    count_df = NULL,        # Optional: df with normalized counts
    annot = NULL            # Optional: df with gene annotations
  ) {

  # Final contrast
  if (is.null(contrasts)) {
    stopifnot(!is.null(fct))
    stopifnot(!is.null(comp))
    fcontrast <- c(fct, comp)
  } else {
    fcontrast <- list(contrasts)
  }
  
  # Get DEseq results
  res <- results(dds,
                 contrast = fcontrast,
                 lfcThreshold = lfc_tres,
                 alpha = p_tres,
                 tidy = TRUE) |>
    dplyr::rename(gene = row,
                  lfc = log2FoldChange,
                  mean = baseMean) |>
    dplyr::select(-lfcSE, -stat) |>
    arrange(padj) |>
    as_tibble()
  
  if (is.null(contrasts)) {
    res <- res |>
      mutate(group1 = comp[1],
             group2 = comp[2],
             contrast = paste0(comp, collapse = "_"))
  }

  # Include mean normalized counts
  if (!is.null(count_df) & is.null(contrasts)) {
    fcount_df <- count_df |> dplyr::filter(.data[[fct]] %in% comp)
    
    group_means <- fcount_df |>
      group_by(gene, .data[[fct]]) |>
      summarize(mean = mean(count), .groups = "drop") |>
      pivot_wider(id_cols = gene, values_from = mean, names_from = all_of(fct))
    colnames(group_means)[2:3] <- c("mean1", "mean2")
    
    overall_means <- fcount_df |>
      group_by(gene) |>
      summarize(mean = mean(count), .groups = "drop")
    
    fcount_df <- dplyr::left_join(group_means, overall_means, by = "gene")
    
    res <- dplyr::left_join(dplyr::select(res, -mean), fcount_df, by = "gene")

    # Determine whether a gene is DE
    res <- res |>
      mutate(
        isDE = ifelse(padj < p_tres & abs(lfc) > lfc_tres & (mean1 > mean_tres | mean2 > mean_tres),
                      TRUE, FALSE)
        )
  } else {
    # Determine whether a gene is DE
    res <- res |>
      mutate(
        isDE = ifelse(padj < p_tres & abs(lfc) > lfc_tres & mean > mean_tres,
                      TRUE, FALSE)
        )
  }

  # Only keep significant genes
  if (sig_only == TRUE) res <- res |> dplyr::filter(isDE == TRUE)

  # Add gene annotation
  if (!is.null(annot)) res <- dplyr::left_join(res, annot, by = "gene")

  # Arrange by adj p-value, remove pvalue column
  res <- res |>
    dplyr::select(-pvalue) |> 
    arrange(padj, -abs(lfc))
  
  # Report
  n_sig <- sum(res$isDE, na.rm = TRUE)
  n_up <- nrow(dplyr::filter(res, lfc > 0 & isDE == TRUE))
  n_down <- nrow(dplyr::filter(res, lfc < 0 & isDE == TRUE))
  
  if (is.null(contrasts)) comp_str <- paste0(comp[1], " vs ", comp[2])
  if (!is.null(contrasts)) comp_str <- paste0(contrasts, collapse = " - ")
  
  message(comp_str, " - Nr DEG: ", n_sig, " (", n_up, "/", n_down, " up/down in group1)")
  
  # Return the final results df
  return(res)
}

# Run the extract_DE function for all pairwise comparisons of a factor
extract_DE_all <- function(dds, fct, count_df = NULL, annot = NULL) {
  lvls <- levels(colData(dds)[[fct]])
  combs <- as.list(as.data.frame(combn(lvls, 2)))
  
  extract_DE2 <- function(lvl1, lvl2) {
    extract_DE(dds = dds, fct = fct, comp = c(lvl1, lvl2),
               count_df = count_df, annot = annot)
  } 
  
  map_dfr(.x = combs, .f = ~extract_DE2(lvl1 = .x[1], lvl2 = .x[2]))
}

# Function to provide shrunken LFC estimates
shrink_lfc <- function(
    dds,                # DESeq object
    fct,                # focal factor (column name)
    comp,               # vector of two with focal factor levels (column values)
    lfc_tres = 1,       # If set to 0, original p-values will be used; otherwise, LFC-based s-values
    p_tres = 0.05,      # P-value threshold
    s_tres = 0.005      # S-value threshold
  ) {

  coef <- paste0(fct, "_", paste0(comp, collapse = "_vs_"))
  message("\nCoefficient: ", coef)

  if (!coef %in% resultsNames(dds)) {
    message("This contrast is not in resultsNames, so we need to rerun DESeq()")
    dds[[fct]] <- relevel(dds[[fct]], ref = comp[2])
    design(dds) <- as.formula(paste("~", fct))
    dds <- DESeq(dds)
  }

  res <- lfcShrink(dds,
                   coef = coef,
                   type = "apeglm",
                   lfcThreshold = lfc_tres)

  res <- res |>
    as.data.frame() |>
    as_tibble() |>
    rownames_to_column("gene") |>
    dplyr::rename(mean = baseMean,
                  lfc = log2FoldChange) |>
    dplyr::select(-lfcSE) |>
    mutate(group1 = comp[1],
           group2 = comp[2],
           contrast = paste0(comp, collapse = "_"))

  if (lfc_tres == 0) {
    res <- res |> mutate(isDE = ifelse(padj < p_tres, TRUE, FALSE))
  } else {
    res <- res |> mutate(isDE = ifelse(svalue < tres, TRUE, FALSE))
  }

  # Report
  if (lfc_tres == 0) {
    n_sig <- sum(res$padj < p_tres, na.rm = TRUE)
    message(comp[1], " vs ", comp[2], " - Nr DEGs (from padj):    ", n_sig)
  } else {
    n_sig <- sum(res$svalue < s_tres, na.rm = TRUE)
    message(comp[1], " vs ", comp[2], " - Nr DEGs (from s-value):   ", n_sig)
  }

  return(res)
}


# Function to normalize counts from a DESeq2 object
norm_counts <- function(
    dds,                     # DESeq2 object
    transform = "rlog",      # Normalization/transformation method: either 'vst', 'rlog', or 'lib_size'
    annot = NULL,            # Gene ID column should be named 'gene'
    return_matrix = FALSE,   # Don't transform to tidy format & don't add metadata
    blind = TRUE,            # Whether to take statistical design of dds into account
    design = NULL            # Use this statistical design (will force 'blind' to FALSE)
    ) {
  
  # Set the design
  if (!is.null(design)) {
    design(dds) <- formula(design)
    blind <- FALSE
  }

  # Normalize the counts
  if (transform == "vst") {
    count_mat <- assay(vst(dds, blind = blind))
  } else if (transform == "rlog") {
    # Suppress messages to avoid "vst is much faster transformation - message"
    count_mat <- suppressMessages(assay(rlog(dds, blind = blind)))
  } else if (transform == "lib_size") {
    dds <- estimateSizeFactors(dds)
    count_mat <- sweep(assay(dds), 2, sizeFactors(dds), "/")
  } else {
    stop("Unknown transformation method '", transform, "'")
  }

  # Stop here if there's no metadata
  if (return_matrix == TRUE) return(count_mat)

  # Get metadata from dds if not provided
  meta_df <- as.data.frame(colData(dds))
  meta_df$sample <- rownames(meta_df)
  rownames(meta_df) <- NULL

  # Get a long-format count df with metadata
  count_df <- as.data.frame(count_mat) |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "sample", values_to = "count") |>
    dplyr::left_join(meta_df, by = "sample")

  if (!is.null(annot)) count_df <- count_df |> dplyr::left_join(annot, by = "gene")
  
  return(count_df)
}


# Function to make an MA plot
pMA <- function(
    deseq_results,               # should be DE results df from extract_DE()
    rm_padj_na = TRUE,
    ptsize_nonsig = 1,
    ptsize_sig = 1,
    ptcol_sig = "blue",
    alpha_nonsig = 0.3,
    alpha_sig = 0.3,
    x_min = NA,
    contrast = NULL,
    x_breaks = c(1, 10, 100, 1000, 10000, 100000) # Use 'waiver()' to get auto-breaks
    ) {

  # Subset to results for a specific contrast
  if (!is.null(contrast)) {
    fcontrast <- contrast
    deseq_results <- deseq_results |> dplyr::filter(contrast %in% fcontrast)
  }

  # Remove genes with NA as the adj. p-value
  if (rm_padj_na == TRUE) {
    message("Removing genes with NA as the adj. pvalue...")
    deseq_results <- deseq_results |> dplyr::filter(!is.na(padj))
  }

  # Make all NA isDE value FALSE, so they will be plotted
  deseq_results <- deseq_results |> mutate(isDE = ifelse(is.na(isDE), FALSE, isDE))

  # Labels
  xlab <- "Mean of normalized counts"
  ylab <- expression("Log"[2]*"-fold change")

  # Create the plot
  p <- ggplot() +
    geom_point(data = dplyr::filter(deseq_results, isDE == FALSE),
               aes(x = mean,
                   y = lfc),
               shape = 21,
               fill = "grey80",
               color = "grey40",
               size = ptsize_nonsig,
               alpha = alpha_nonsig) +
    geom_point(data = dplyr::filter(deseq_results, isDE == TRUE),
               aes(x = mean, y = lfc),
               shape = 21,
               fill = ptcol_sig,
               color = ptcol_sig,
               size = ptsize_sig,
               alpha = alpha_sig) +
    scale_x_log10(labels = scales::comma_format(accuracy = 1),
                  breaks = x_breaks) +
    coord_cartesian(xlim = c(x_min, NA)) +
    geom_hline(yintercept = 0) +
    guides(color = "none") +
    labs(x = xlab,
         y = ylab) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, hjust = 0.5))

  if (!is.null(contrast)) {
    if (length(fcontrast) > 1) p <- p + facet_wrap(vars(contrast))
  }

  return(p)
}


# Volcano plot
pvolc <- function(
    DE_df,                        # Output from extract_DE()
    contrasts = "all",            # Subset contrasts using 'contrast' column in DE_df. Default: no subsetting
    sig_only = TRUE,              # Whether to plot all (FALSE) or only significant (TRUE) genes
    interactive = FALSE,          # Use ggiraph to make the plot interactive
                                  # NOTE 1: Will use columns 'gene', 'gene_name', 'gene_description'
                                  # NOTE 2: Need to call 'girafe(ggobj = p)' on output!
    add_ids = FALSE,              # Whether to add gene IDs with ggrepel
    id_column = "gene_name",      # Gene ID column to use when add_ids == TRUE
    colors = NULL,                # Manual colors for each contrast
    plot_grid = FALSE,
    grid_rows = NULL,
    grid_cols = NULL,
    facet_scales = "fixed",
    return_plot = FALSE           # Default is to print but not return the plot object
                                  # (Except when interactive == TRUE)
    ) {

  # Rename contrasts vector
  fcontrasts <- contrasts

  # Whether to return the plot object
  if (interactive == TRUE) return_plot <- TRUE
  
  # Remove genes with NA as the adj. pvalue
  DE_df <- DE_df |> dplyr::filter(!is.na(padj))

  # Remove non-significant genes
  if (sig_only == TRUE) DE_df <- DE_df |> dplyr::filter(isDE == TRUE)

  # Get log10 of padj to remove Infinite values
  DE_df <- DE_df |>
    mutate(padj = -log10(padj)) |>
    dplyr::filter(!is.infinite(padj))

  # Only take results for the the focal contrast
  if (!is.null(fcontrasts) && fcontrasts[1] != "all") {
    message("Selecting only the focal contrasts...")
    DE_df <- DE_df |> dplyr::filter(contrast %in% fcontrasts)
  }
  
  # Arrange
  DE_df <- DE_df |>
    group_by(contrast) |>
    dplyr::arrange(-padj, -abs(lfc)) |>
    ungroup()
  
  # Make sure the gene info columns exist
  if (interactive == TRUE) {
    if(! "gene" %in% colnames(DE_df)) {
      warning("Column gene does not exist")
      DE_df$gene <- NA
    }
    if(! "gene_name" %in% colnames(DE_df)) {
      warning("Column gene_name does not exist")
      DE_df$gene_name <- NA
    } 
    if(! "gene_description" %in% colnames(DE_df)) {
      warning("Column gene_description does not exist")
      DE_df$gene_description <- NA
    }
  }
  
  # Axis titles
  xlab <- expression("Log"[2]*"-fold change")
  ylab <- expression("-Log"[10]*" of adj. P")

  # Make the base plot
  p <- ggplot(DE_df) +
    aes(x = lfc, y = padj) +
    geom_point(data = dplyr::filter(DE_df, isDE == FALSE),
               fill = "grey80", color = "grey40",
               size = 2, shape = 21, alpha = 0.3)
  
  # Add points
  if (interactive == TRUE) {
    p <- p + ggiraph::geom_point_interactive(
      data = dplyr::filter(DE_df, isDE == TRUE),
      mapping = aes(fill = contrast,
                    tooltip = paste(gene, "\n", gene_name, "\n", gene_description)),
      size = 2, shape = 21, color = "grey20", alpha = 0.5
    )
  } else {
    p <- p + geom_point(
      data = dplyr::filter(DE_df, isDE == TRUE),
      mapping = aes(fill = contrast),
      size = 2, shape = 21, color = "grey20", alpha = 0.5
    )
  }
    
  p <- p +
    geom_vline(xintercept = 0, color = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = xlab, y = ylab) +
    guides(fill = "none") +
    theme(panel.grid.minor = element_blank())

  # Add color scheme
  if (is.null(colors)) p <- p + scale_fill_brewer(palette = "Dark2")
  if (! is.null(colors)) p <- p + scale_fill_manual(values = colors)

  # Add gene IDs
  if (add_ids == TRUE) {
    DE_df[[id_column]] <- str_trunc(DE_df[[id_column]], width = 35)
    
    p <- p + ggrepel::geom_text_repel(
      data = slice_head(DE_df, n = 10, by = contrast),
      mapping = aes(label = .data[[id_column]]),
      max.overlaps = 20,
      min.segment.length = 0.1,
      color = "grey50",
      fontface = "italic"
    )
  }
  
  # When no focal contrast is specified, show all with a facet
  if ((length(fcontrasts) > 1 || fcontrasts[1] == "all") & plot_grid == FALSE) {
    p <- p + facet_wrap(vars(contrast), nrow = 1, scales = facet_scales)
  }
  
  if (plot_grid == TRUE) {
    p <- p + facet_grid(rows = vars(!!sym(grid_rows)),
                        cols = vars(!!sym(grid_cols)),
                        scales = facet_scales)
  }
  
  print(p)
  if (return_plot == TRUE) return(p)
}

# Heatmap with ComplexHeatmap
cheatmap <- function(
    genes,
    count_mat,
    meta_df,
    groups = NULL,
    show_rownames = TRUE,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    scale = FALSE,
    logtrans = FALSE,
    heat_color_scheme = "blue_black_yellow",
    heat_colors = NULL,
    ...
    ) {
  
  suppressPackageStartupMessages(library(ComplexHeatmap))
  
  # Set up 'groups' (factors) by which to group samples
  if (!is.null(groups)) {
    # Arrange metadata according to the columns with included factors
    meta_df <- meta_df |>
      dplyr::select(all_of(groups)) |>
      arrange(across(all_of(groups)))
    
    # Set colors
    col_list <- list()
    palettes <- c("Dark2", "Set1", "Set2", "Set3")
    for (i in seq_along(groups)) {
      group <- groups[i]
      palette <- palettes[i]
      
      if (is.factor(meta[[group]])) levs <- levels(meta[[group]])
      if (!is.factor(meta[[group]])) levs <- unique(meta[[group]])
      
      colors <- suppressWarnings(
        RColorBrewer::brewer.pal(n = length(levs), name = palette)[1:length(levs)]
      )
      col_list[[group]] <- colors
      names(col_list[[group]]) <- levs
    }
    
    # Set final annotation
    annot <- HeatmapAnnotation(
      df = meta_df,
      col = col_list
    )
    
  } else {
    # Don't include metadata if no groups are provided
    meta_df <- NA
  }
  
  # Select and arrange count matrix
  fmat <- count_mat[match(genes, rownames(count_mat)),
                    match(rownames(meta_df), colnames(count_mat))]
  fmat <- as.matrix(fmat)
  
  # Scale matrix
  if (scale == TRUE) {
    fmat <- t(scale(t(fmat)))
  } else if (logtrans == TRUE) {
    fmat <- log10(fmat)
    fmat[fmat == -Inf] <- 0
  }
  
  # Define expression-level colors
  if (is.null(heat_colors)) {
    if (heat_color_scheme == "blue_yellow_red") {
      heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
    } else if (heat_color_scheme == "blue_black_yellow") {
      heat_colors <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
    } else {
      stop("'heat_color_scheme' should be 'blue_yellow_red' or 'blue_black_yellow', or 'heat_colors' should be defined")
    }
  }
  
  if (scale == TRUE) {
    # Set min. and max. expression counts for colors (avoid effects of outliers)
    min_count <- -3
    max_count <- 3
  } else {
    min_count <- quantile(fmat, probs = 0.01)
    max_count <- quantile(fmat, probs = 0.99)
  }
  my_breaks <- seq(min_count, max_count, length.out = 100)
  color_fun <- circlize::colorRamp2(my_breaks, heat_colors)
  
  # Truncate long names
  row.names(fmat) <- str_trunc(row.names(fmat), width = 20, ellipsis = "")
  
  # Legend name
  if (scale == TRUE) {
    heat_legend_name <- "Normalized\n& scaled\ngene count"
  } else if (logtrans == TRUE) {
    heat_legend_name <- "Normalized &\nlog-transformed\ngene count"
  } else {
    heat_legend_name <- "Normalized\ngene count"
  }
  
  # Create the heatmap
  hm <- Heatmap(
    matrix = fmat,
    top_annotation = annot,
    col = color_fun,
    name = heat_legend_name,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 10)
    )
  
  draw(hm, merge_legend = TRUE)
}

# Heatmap plot showing abundances for multiple/many genes
#TODO Allow to add gene annotation
pheat <- function(
    genes,                      # Vector of genes to include
    count_mat,                  # Count matrix (with gene names as rownames)
    meta_df,                    # Metadata, should have sample IDs as rownames!
    samples = NULL,             # Vector of samples to include
    groups = NULL,              # Column names from metadata to show as groups at top of heatmap
    show_rownames = TRUE,       # Whether to show row (gene) names
    show_colnames = FALSE,      # Whether to show column (sample) names
    cluster_rows = TRUE,        # Whether to do hierarchical clustering on the rows (genes)
    cluster_cols = FALSE,       # Whether to do hierarchical clustering on the rows (samples)
    logtrans = FALSE,           # Whether to log-transform the counts
    scale = "none",
    annotation_colors = NULL,
    id_labsize = 10,
    ...                         # Arguments to be passed to the pheatmap function
    ) {

  # Arrange metadata according to the columns with included factors
  if (!is.null(groups)) {
    meta_df <- meta_df |>
      dplyr::select(all_of(groups)) |>
      arrange(across(all_of(groups)))
  }

  # Select and arrange count matrix
  fcount_mat <- count_mat[match(genes, rownames(count_mat)),
                          match(rownames(meta_df), colnames(count_mat))]
  if (!is.null(samples)) fcount_mat <- fcount_mat[, match(samples, colnames(fcount_mat))]
  fcount_mat <- as.matrix(fcount_mat)
  
  # Don't include metadata if no groups are provided
  if (is.null(groups)) meta_df <- NA
  
  # Log-transform
  if (logtrans == TRUE) {
    fcount_mat <- log10(fcount_mat)
    fcount_mat[fcount_mat == -Inf] <- 0
  }
  
  # If few features are included, reduce the cell (row) height
  cellheight <- ifelse(length(genes) > 20, NA, 20)
  id_labsize <- ifelse(length(genes) > 40, 6, id_labsize)

  # Truncate long names
  row.names(fcount_mat) <- str_trunc(row.names(fcount_mat),
                                     width = 25, ellipsis = "")

  # Create the plot
  p <- pheatmap::pheatmap(
    fcount_mat,
    annotation_col = meta_df,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    annotation_colors = annotation_colors,
    cellheight = cellheight,
    fontsize = 9,
    fontsize_row = id_labsize,
    cex = 1,
    scale = scale,
    ...
    )

  return(p)
}


# Boxplot showing abundances for a single gene
#TODO - Allow for annotation column in count_df, instead of separate annot df
pbox <- function(
    gene,                      # Gene ID
    count_df,                  # Count df with column 'count' with counts to plot
    x,                         # What to plot along the x-axis (column from count_df)
    color_by = NULL,           # What to color the boxes and points by (column from count_df)
    facet_by = NULL,           # What to facet the plot by (column from count_df)
    facet_scales = "fixed",
    annot = NULL,              # Annotation: gene info for plot title
    colors = NULL,             # Manually provide colors
    log_scale = FALSE,         # Whether counts (y-axis) should be on a y-scale 
    ymin = NA,                 # Force a min. value for the y-axis
    xlab = NULL,
    save_plot = FALSE,
    plotdir = "results/figures/geneplots",
    return_plot = FALSE        # Whether to return the plot object
    ) {

  # Rename gene variable
  fgene <- gene

  # Process annotation
  if (!is.null(annot)) {
    annot <- annot |> janitor::clean_names()  
    id_col_name <- colnames(annot)[1] 
    
    if ("gene_description" %in% colnames(annot)) {
      g_descrip <- annot |>
        dplyr::filter(.data[[id_col_name]] == fgene) |>
        pull(gene_description) |>
        str_trunc(width = 50)
    } else {
      g_descrip <- NA
    }
    
    if ("gene_name" %in% colnames(annot)) {
      g_name <- annot |>
        dplyr::filter(.data[[id_col_name]] == fgene) |>
        pull(gene_name)
    } else {
      g_name <- NA
    }
    
    if (!is.na(g_name)) {
      ptitle <- g_name
      if (!is.na(g_descrip)) psub <- paste0(g_descrip, "\n", fgene)
      if (is.na(g_descrip)) psub <- fgene 
    } else {
      ptitle <- fgene
      if (!is.na(g_descrip)) psub <- paste0(g_descrip)
      if (is.na(g_descrip)) psub <- NULL
    }
  } else {
    ptitle <- fgene
    psub <- NULL
  }

  # Filter to contain only focal ID
  count_df <- count_df |> dplyr::filter(gene %in% fgene)

  # Base plot
  p <- ggplot(count_df) +
    aes(y = count, x = .data[[x]])
  if (!is.null(color_by)) p <- p + aes(color = .data[[color_by]])
  
  # Add boxplot
  p <- p + geom_boxplot(outlier.shape = NA)
  
  # Add points
  if (!is.null(color_by)) {
    # With color-aesthetic, use jitter-dodge and smaller points:
    p <- p + geom_point(
      position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
      size = 2.5
      )
  } else {
    # Withour color aesthetic
    p <- p + geom_point(
      position = position_jitter(width = 0.2, height = 0),
      size = 2.5, color = "grey40"
      )
  }
  
  # Facet
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(
      vars(!!sym(facet_by)),
      nrow = 1,
      scales = facet_scales,
      labeller = "label_both")
  }
    
  # Plot formatting
  p <- p +
    labs(title = ptitle,
         subtitle = psub,
         x = xlab) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold.italic"),
          plot.subtitle = element_text(hjust = 0.5, face = "plain"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm"))

  # No legend needed if color-aes is same as x-aes
  if (!is.null(color_by)) {
    if (x == color_by | facet_by == color_by) p <- p + guides(color = "none")
  }

  # Y scale and label
  if (log_scale == TRUE) {
    p <- p +
      scale_y_log10(labels = scales::comma) +
      labs(y = "Gene count (log10-scale)")
  } else {
    p <- p +
      scale_y_continuous(labels = scales::comma,
                         limits = c(ymin, NA),
                         expand = expansion(mult = c(0.03, 0.03))) +
      labs(y = "Gene count")
  }
  
  # Color scheme
  if (is.null(colors)) p <- p + scale_color_brewer(palette = "Set1")
  if (!is.null(colors)) p <- p + scale_color_manual(values = colors)

  # Save plot
  if (save_plot == TRUE) {
    dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
    plotfile <- file.path(plotdir, paste0(gene, ".png"))
    ggsave(filename = plotfile, plot = p,
           width = 6, height = 6, dpi = "retina")
  }

  print(p)
  if (return_plot == TRUE) return(p)
}


# Wrapper function to combine 4 single-gene boxplots
p4box <- function(
    genes,
    count_df,
    annot = NULL,
    x,
    color_by,
    colors = NULL,
    ...
    ) {
  
  suppressPackageStartupMessages(library(patchwork))

  p1 <- pbox(genes[1], count_df, annot, cols,
             theme_base_size = 11, ...) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(5.5, 1, 5.5, 5.5, "pt"))

  p2 <- pbox(genes[2], count_df, annot, cols,
             theme_base_size = 11, ...) +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(5.5, 1, 5.5, 5.5, "pt"))

  p3 <- pbox(genes[3], count_df, annot, cols,
             theme_base_size = 11, ...) +
    theme(plot.margin = margin(5.5, 1, 5.5, 5.5, "pt"))

  p4 <- pbox(genes[4], count_df, annot, cols,
             theme_base_size = 11, ...) +
    theme(axis.title.y = element_blank(),
          plot.margin = margin(5.5, 1, 5.5, 5.5, "pt"))

  p <- (p1 + p2) / (p3 + p4)
  print(p)
}

# Function to create a PCA the same way as the DESeq2 function
# (The DEseq function uses `prcomp()` under the hood)
pca_prcomp <- function(
    dds,                      # DESeq2 object
    transform = "rlog",       # Normalization procedure, either 'rlog' or 'vst'
    ntop = 500                # Top n most variable genes to include
) {
  
  message("Using the top ", ntop, " most highly variable genes...")
  
  # Normalize the data
  if (transform == "vst") dds_norm <- vst(dds, blind = TRUE)
  if (transform == "rlog") dds_norm <- rlog(dds, blind = TRUE)
  
  # Run the PCA
  rv <- rowVars(assay(dds_norm))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(dds_norm)[select, ]))
  
  # Prep data for plotting
  percent_var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
  pca_df <- as.data.frame(cbind(colData(dds_norm), pca$x))
  
  # Put everything in a list
  pca_res <- list(df = pca_df, percent_var = percent_var)
  
  return(pca_res)
}

# Function to run the PCA using the PCAtools `pca()` function
pca_pcatools <- function(
    dds = NULL,          # Input option 1: a DESeq object
    count_mat = NULL,    # Input option 2: A (normalized) count matrix    
    remove_prop = 0.1,   # Remove this proportion of variables (genes) with the lowest variance
    transform = "rlog",  # Type of normalization
    meta = NULL          # If the input is count_mat, use this arg to add metadata
) {
  
  # If input is a DESeq object, normalize and get metadata
  if (!is.null(dds)) {
    if (transform == "vst") count_mat <- assay(vst(dds, blind = TRUE))
    if (transform == "rlog") count_mat <- assay(rlog(dds, blind = TRUE))
    meta <- colData(dds)
  }
  
  # Run the PCA
  PCAtools::pca(count_mat, metadata = meta, removeVar = remove_prop)
}

# Function to create a regular PCA plot from PCAtools PCA results
pca_plot <- function(
    pca_res,                      # PCA results after running `PCA_run()`
    x = "PC1",                    # Principal component to plot on the x-axis
    y = "PC2",                    # Principal component to plot on the y-axis
    color_by = NULL,              # Vary point color by this variable from the metadata
    fill_by = NULL,               # Vary point fill by this variable from the metadata
    shape_by = NULL,              # Vary point shape by this variable from the metadata
    ellipse_by = NULL,            # Draw ellipses around groups of points by this variable
    add_ids = FALSE,              # Add sample names to points TRUE/FALSE
    title = NULL,                 # Add a title as a string; "NULL" is no title
    aspect_ratio = FALSE,         # When TRUE, aspect ratio will be according to % variation for each PC
    pt_size = 5,                  # Point size
    pt_shape = NULL,              # Point shape (defaults: 16 or 21, the latter when using 'fill=')
    interactive = FALSE,          # Use ggiraph to make the plot interactive
                                  # NOTE: need to call 'girafe(ggobj = p)' on output!
    pc_var = TRUE                 # Indicate % of variation in axis titles
) {
  
  # Prep the dataframe
  if (class(pca_res) == "pca") {
    # Extract the data from the PCAtools object
    meta <- as.data.frame(pca_res$metadata)
    meta$sample <- rownames(meta)
    df <- as.data.frame(pca_res$rotated)
    df$sample <- rownames(df)
    df <- dplyr::left_join(df, meta, by = "sample")
    
    # Percent variation explained by each PC
    percent_var <- round(pca_res$variance, 2)
    
  } else {
    # Expecting a list from the pca_prcomp() function
    df <- pca_res$df
    df$sample <- rownames(df)
    
    # Percent variation explained by each PC
    percent_var <- pca_res$percent_var
  }
  
  # Sample names
  if (add_ids == TRUE) names <- pca_res$yvars
  if (add_ids == FALSE) names <- 0
  
  # Percent variation explained by each PC, and aspect ratio
  x_pct <- percent_var[as.integer(sub("PC", "", x))]
  y_pct <- percent_var[as.integer(sub("PC", "", y))]
  if (aspect_ratio == TRUE) asp <- y_pct / x_pct else asp <- 1
  
  # Axis labels
  if (pc_var == TRUE) {
    x_lab <- paste0(x, " (", x_pct, "% of variance)")
    y_lab <- paste0(y, " (", y_pct, "% of variance)")
  } else {
    x_lab <- x
    y_lab <- y
  }
  
  # Create the base plot
  p <- ggplot(data = df) +
    aes(x = .data[[x]], y = .data[[y]])
  
  # Color aesthethic
  if (!is.null(color_by)) {
    p <- p +
      aes(color = .data[[color_by]]) +
      scale_color_brewer(palette = "Dark2")
  }
  
  # Fill aesthetic
  if (!is.null(fill_by)) {
    p <- p +
      aes(fill = .data[[fill_by]]) +
      scale_fill_brewer(palette = "Dark2")
    if (is.null(shape_by)) pt_shape <- 21
  } else {
    if (is.null(shape_by)) pt_shape <- 16
  }
  
  # Shape aesthetic
  if (!is.null(shape_by)) {
    n_shapes <- length(unique(df[[shape_by]]))
    
    if (!is.null(fill_by)) {
      if (n_shapes > 5) warning("Too many unique values in ", shape_by, " (can only plot 5 hollow shapes)")
      shapes <- c(21:25)[1:n_shapes]
    } else {
      shapes <- c(15:25)[1:n_shapes]
    }
    
    p <- p +
      aes(shape = .data[[shape_by]]) +
      scale_shape_manual(values = shapes) +
      guides(fill = guide_legend(override.aes = list(shape = 21)))
  }
  
  # Add points
  if (interactive == FALSE) {
    
    # Non-interactive
    if (is.null(shape_by)) {
      p <- p + geom_point(size = pt_size, shape = pt_shape)
    } else {
      if (!is.null(shape_by)) p <- p + geom_point(size = pt_size)
    }
  } else {
    
    # Interactive
    if (is.null(shape_by)) {
      p <- p + ggiraph::geom_point_interactive(
        aes(tooltip = sample), size = pt_size, shape = pt_shape
        )
    } else {
      p <- p + ggiraph::geom_point_interactive(
        aes(tooltip = sample), size = pt_size
        )
    }
  }
  
  # Add sample IDs
  if (add_ids == TRUE) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = sample), max.overlaps = 20, point.padding = 3
    )
  }
  
  # Add ellipses around groups
  if (!is.null(ellipse_by)) {
    p <- p +
      ggforce::geom_mark_ellipse(
        aes(group = .data[[ellipse_by]]),
        expand = unit(2, "mm"),
        color = "grey60", fill = "grey90", linetype = "dashed"
        )
  }
  
  # Plot polishing
  p <- p +
    labs(x = x_lab, y = y_lab, title = title) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, color = "grey20"),
          aspect.ratio = asp)
  
  return(p)
}


# Function to create a PCA biplot
pca_biplot <- function(
    pca_res,                      # PCA results after running `PCA_run()`
    col,                          # Vary point color by this variable from the metadata
    shape,                        # Vary point shape by this variable from the metadata
    x = "PC1",                    # Principal component to plot on the x-axis
    y = "PC2",                    # Principal component to plot on the y-axis
    n_genes = 5,                  # Number of genes to plot loadings for
    add_ids = TRUE,               # Add sample names to points TRUE/FALSE
    title = NULL                  # Add a title as a string; "NULL" is no title
) {
  percent_var <- round(pca_res$variance, 2)
  
  if (add_ids == TRUE) name_idx <- 1:nrow(pca_res$metadata)
  if (add_ids == FALSE) name_idx <- 0
  
  x_nr <- as.integer(sub("PC", "", x))
  y_nr <- as.integer(sub("PC", "", y))
  
  p <- biplot(pca_res,         # `biplot()` function from PCAtools
              x = x,
              y = y,
              selectLab = name_idx,
              showLoadings = TRUE,
              ntopLoadings = n_genes / 2,
              colby = col,
              shape = shape) +
    xlab(paste0(x, " (", percent_var[x_nr], "% of variance)")) +
    ylab(paste0(y, " (", percent_var[y_nr], "% of variance)")) +
    ggtitle(title)
  
  return(p)
}

# Function to create a Venn diagram
pvenn <- function(
    DE_res,               # DE results from extract_DE(): requires 'isDE', 'gene', and 'contrast' columns
    contrasts,            # Select the following contrasts from the 'contrast' column
    cat_names = NULL,     # Give non-default category (contrast) names
    cat_colors = NULL,    # Give non-default (RColorBrewer Set 2) category names
    auto_scale = TRUE     # Scale the circles by category size
    ) {
  
  # Create the list
  venn_list <- map(
    .x = contrasts,
    .f = function(x) DE_res |> filter(contrast == x, isDE == TRUE) |> pull(gene)
  )
  
  # Give the list elements names
  if (is.null(cat_names)) cat_names <- contrasts
  names(venn_list) <- cat_names
  
  # Set the colors
  if (is.null(cat_colors)) {
    cat_colors <- suppressWarnings(
      RColorBrewer::brewer.pal(length(venn_list), name = "Set2")[1:length(venn_list)]
    )
  }
  
  # Make the Venn diagram
  ggvenn::ggvenn(venn_list,
                 fill_color = cat_colors,
                 auto_scale = auto_scale)
}
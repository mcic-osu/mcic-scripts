# Packages
if(!require(janitor)) install.packages("janitor")
if(!require(tidyverse)) install.packages("tidyverse")
library(janitor)
library(tidyverse)

# Run DE analysis
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


# Get DE results
extract_DE <- function(
    comp,                   # Vector of 2 with focal levels of factor 'fac'
    fac,                    # Focal factor (column in metadata with the levels in 'comp')
    dds,                    # DESeq2 object
    sig_only = FALSE,       # Only return significant results
    p_tres = 0.05,          # Adj. p-value threshold for DE significance
    lfc_tres = 0.5,         # Log2-fold change threshold for DE significance
    mean_tres = 0,          # Mean expr. level threshold for DE significance
    count_df = NULL,        # Optional: df with normalized counts
    annot = NULL            # Optional: df with gene annotations
  ) {

  # Get DEseq results
  res <- results(dds,
                 contrast = c(fac, comp),
                 lfcThreshold = lfc_tres,
                 alpha = p_tres,
                 tidy = TRUE) |>
    dplyr::rename(gene = row,
                  lfc = log2FoldChange,
                  mean = baseMean) |>
    mutate(group1 = comp[1],
           group2 = comp[2],
           contrast = paste0(comp, collapse = "_")) |>
    dplyr::select(-lfcSE, -stat) |>
    arrange(padj) |>
    as_tibble()

  # Include mean normalized counts
  if (!is.null(count_df)) {
    fcount_df <- count_df |> dplyr::filter(.data[[fac]] %in% comp)
    
    group_means <- fcount_df |>
      group_by(gene, .data[[fac]]) |>
      summarize(mean = mean(count), .groups = "drop") |>
      pivot_wider(id_cols = gene, values_from = mean, names_from = all_of(fac))
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
  message(comp[1], " vs ", comp[2], " - Nr DEG: ", n_sig,
          " (", n_up, "/", n_down, " up/down in group1)")
  
  return(res)
}

# Run the extract_DE function for all pairwise comparisons of a factor
extract_DE_all <- function(dds, fac, count_df = NULL) {
  lvls <- levels(colData(dds)[[fac]])
  combs <- as.list(as.data.frame(combn(lvls, 2)))
  
  extract_DE2 <- function(lvl1, lvl2) {
    extract_DE(
      comp = c(lvl1, lvl2),
      fac = fac,
      dds = dds,
      count_df = count_df
    )
  } 
  
  map_dfr(.x = combs, .f = ~extract_DE2(.x[1], .x[2]))
}

# Function to provide shrunken LFC estimates
shrink_lfc <- function(
    dds,                # DESeq object
    fac,                # focal factor (column name)
    comp,               # vector of two with focal factor levels (column values)
    lfc_tres = 1,       # If set to 0, original p-values will be used; otherwise, LFC-based s-values
    p_tres = 0.05,      # P-value threshold
    s_tres = 0.005      # S-value threshold
  ) {

  coef <- paste0(fac, "_", paste0(comp, collapse = "_vs_"))
  message("\nCoefficient: ", coef)

  if (!coef %in% resultsNames(dds)) {
    message("This contrast is not in resultsNames, so we need to rerun DESeq()")
    dds[[fac]] <- relevel(dds[[fac]], ref = comp[2])
    design(dds) <- as.formula(paste("~", fac))
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


# Function to normalize counts
norm_counts <- function(
    dds,                     # DESeq object
    transform = "rlog",      # Normalization/transformation method: either 'vst', 'rlog', or 'lib_size'
    annot = NULL,            # Gene ID column should be named 'gene'
    return_matrix = FALSE    # Don't transform to tidy format & don't add metadata
    ) {
  
  # Normalize the counts
  if (transform == "vst") {
    count_mat <- assay(vst(dds, blind = TRUE))
  } else if (transform == "rlog") {
    # Suppress messages to avoid "vst is much faster transformation - message"
    count_mat <- suppressMessages(assay(rlog(dds, blind = TRUE)))
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

  if (!is.null(annot)) {
    count_df <- count_df |>
      dplyr::left_join(annot, by = "gene")
  }
  
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
    interactive = FALSE,
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

  # Interactive text
  if ("gene_name" %in% colnames(deseq_results)) {
      glue_string <- "Gene ID: {gene}
                      Gene name: {gene_name}
                      Description: {gene_description}"
  } else {
      glue_string <- "Gene ID: {gene}"
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
               aes(x = mean,
                   y = lfc,
                   text = glue::glue(glue_string)),
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

  if (interactive) ggplotly(p, tooltip = "text") else return(p)
}


# Volcano plot
pvolc <- function(DE_df,
                  sig_only = TRUE,
                  contrasts = "all",
                  colors = NULL,
                  interactive = FALSE,
                  plot_grid = FALSE,
                  grid_rows = NULL,
                  grid_cols = NULL,
                  facet_scales = "fixed") {

  fcontrasts <- contrasts

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

  # Interactive text
  if ("gene_name" %in% colnames(DE_df)) {
    glue_string <- "Gene ID: {gene}
                    Gene name: {gene_name}
                    Description: {gene_description}"
  } else {
    glue_string <- "Gene ID: {gene}"
  }

  # Labels
  xlab <- expression("Log"[2]*"-fold change")
  ylab <- expression("-Log"[10]*" of adj. P")

  # Make the plot
  p <- ggplot(DE_df) +
    aes(x = lfc,
        y = padj,
        text = glue_string) +
    geom_point(data = dplyr::filter(DE_df, isDE == FALSE),
               fill = "grey80",
               size = 2, shape = 21, color = "grey40", alpha = 0.3) +
    geom_point(data = dplyr::filter(DE_df, isDE == TRUE),
               aes(fill = contrast),
               size = 2, shape = 21, color = "grey20", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = xlab,
         y = ylab) +
    guides(fill = "none") +
    theme(panel.grid.minor = element_blank())

  if (is.null(colors)) {
    p <- p + scale_fill_brewer(palette = "Dark2")
  } else {
    p <- p + scale_fill_manual(values = colors)
  }

  # When no focal contrast is specified, show all with a facet
  if ((length(fcontrasts) > 1 || fcontrasts[1] == "all") & plot_grid == FALSE) {
    p <- p + facet_wrap(vars(contrast),
                        nrow = 1,
                        scales = facet_scales)
  }
  if (plot_grid == TRUE) {
    p <- p + facet_grid(rows = vars(!!sym(grid_rows)),
                        cols = vars(!!sym(grid_cols)),
                        scales = facet_scales)
  }
  
  if (interactive == TRUE) {
    suppressPackageStartupMessages(library(plotly))
    ggplotly(p, tooltip = "text")
  } else {
    return(p)
  }
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
      
      cols <- suppressWarnings(
        RColorBrewer::brewer.pal(n = length(levs), name = palette)[1:length(levs)]
      )
      col_list[[group]] <- cols
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
pheat <- function(genes,
                  count_mat,
                  meta_df,                    # Should have sample IDs as rownames!
                  groups = NULL,
                  show_rownames = TRUE,
                  show_colnames = FALSE,
                  cluster_rows = TRUE,
                  logtrans = FALSE,
                  scale = "none",
                  annotation_colors = NULL,
                  id_labsize = 10,
                  ...) {

  suppressPackageStartupMessages(library(pheatmap))

  # Arrange metadata according to the columns with included factors
  if (!is.null(groups)) {
    meta_df <- meta_df |>
      dplyr::select(all_of(groups)) |>
      arrange(across(all_of(groups)))
  }

  # Select and arrange count matrix
  fcount_mat <- count_mat[match(genes, rownames(count_mat)),
                          match(rownames(meta_df), colnames(count_mat))]
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
                                     width = 20, ellipsis = "")

  # Function to create the plot
  p <- pheatmap(
    fcount_mat,
    annotation_col = meta_df,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
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
pbox <- function(
    gene,
    count_df,
    x_by,
    col_by = NULL,
    annot = NULL,   
    cols = NULL,
    log_scale = FALSE,
    theme_base_size = 13,
    ymin = NA,
    save_plot = FALSE,
    plotdir = "results/figures/geneplots"
    ) {

  fgene <- gene

  if (!is.null(annot)) {
    annot <- annot |> janitor::clean_names()  
    id_col_name <- colnames(annot)[1] 
    
    g_descrip <- annot |>
      dplyr::filter(.data[[id_col_name]] == fgene) |>
      pull(gene_description) |>
      str_trunc(width = 50)

    if ("gene_name" %in% colnames(annot)) {
      g_name <- annot |>
        dplyr::filter(.data[[id_col_name]] == fgene) |>
        pull(gene_name)
      ptitle <- g_name
      psub <- paste0(g_descrip, "\n", fgene)
    } else {
      ptitle <- fgene
      psub <- paste0(g_descrip)
    }
  } else {
    ptitle <- fgene
    psub <- NULL
  }

  # Filter to contain only focal ID
  count_df <- count_df |> dplyr::filter(gene %in% fgene)

  # Make the plot
  p <- ggplot(count_df) +
    aes(y = count, x = .data[[x_by]])
  
  if (!is.null(col_by)) p <- p + aes(color = .data[[col_by]])
  
  # Add boxplot
  p <- p + geom_boxplot(outlier.shape = NA)
  
  # Add points
  if (!is.null(col_by)) {
    # With color-aesthetic, use jitter-doge and smaller points:
    p <- p +
      geom_point(position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
                 size = 1.5)
  } else {
    # Withour color aesthetic
    p <- p +
      geom_point(position = position_jitter(width = 0.1, height = 0),
                 size = 2.5)
  }
  
  # Plot formatting
  p <- p +
    labs(title = ptitle,
         subtitle = psub,
         x = NULL) +
    theme_bw(base_size = theme_base_size) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold.italic"),
          plot.subtitle = element_text(hjust = 0.5, face = "plain"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm"))

  # No legend needed if color-aes is same as x-aes
  if (!is.null(col_by)) if (x_by == col_by) p <- p + guides(color = "none")

  if (log_scale == TRUE) {
    p <- p +
      scale_y_log10(labels = scales::comma) +
      labs(y = "Normalized count (log10-scale)")
  } else {
    p <- p +
      scale_y_continuous(labels = scales::comma,
                         limits = c(ymin, NA),
                         expand = expansion(mult = c(0.03, 0.03))) +
      labs(y = "Normalized count")
  }

  if (is.null(cols)) {
    p <- p + scale_color_brewer(palette = "Set1")
  } else {
    p <- p + scale_color_manual(values = cols)
  }

  if (save_plot == TRUE) {
    dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
    plotfile <- file.path(plotdir, paste0(gene, ".png"))
    ggsave(filename = plotfile, plot = p,
           width = 6, height = 6, dpi = "retina")
  }

  return(p)
}


# Wrapper function to combine 4 single-gene boxplots
p4box <- function(genes, count_df,
                  annot = NULL, cols = NULL,
                  x_by = "tissue", col_by = "tissue", ...) {
  
  library(patchwork)

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
    dds,
    transform = "rlog",
    ntop = 500
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
    dds,                 # DESeq object
    remove_prop = 0.1,   # Remove this proportion of variables (genes) with the lowest variance
    transform = "rlog",  # Type of normalization
    mat_norm = NULL      # Optionally, input a pre-made normalized matrix     
) {
  
  # Install/load the PCAtools package
  if (!require(PCAtools)) {
    BiocManager::install("PCAtools")
    library(PCAtools)
  }
  
  # Normalize the data
  if (!is.null(mat_norm)) {
    if (transform == "vst") mat_norm <- assay(vst(dds, blind = TRUE))
    if (transform == "rlog") mat_norm <- assay(rlog(dds, blind = TRUE))
  }
  
  # Run the PCA
  pca_res <- pca(mat_norm,
                 metadata = colData(dds),
                 removeVar = remove_prop)
  
  return(pca_res)
}

# Function to create a regular PCA plot from PCAtools PCA results
pca_plot <- function(
    pca_res,                      # PCA results after running `PCA_run()`
    x = "PC1",                    # Principal component to plot on the x-axis
    y = "PC2",                    # Principal component to plot on the y-axis
    col = NULL,                   # Vary point color by this variable from the metadata
    fill = NULL,
    shape = NULL,                 # Vary point shape by this variable from the metadata
    pt_size = 5,
    pt_shape = NULL,
    add_ids = FALSE,              # Add sample names to points TRUE/FALSE
    title = NULL,                 # Add a title as a string; "NULL" is no title
    interactive = FALSE,          # Use ggiraph to make the plot interactive
                                  # NOTE: need to call 'girafe(ggobj = p)' on output!
    pc_var = TRUE                 # Indicate % of variation in axis titles
) {
  if (class(pca_res) == "pca") {
    # Extract the data from the PCAtools object
    meta <- as.data.frame(pca_res$metadata)
    meta$sample <- rownames(meta)
    df <- as.data.frame(pca_res$rotated)
    df$sample <- rownames(df)
    df <- dplyr::left_join(df, meta, by = "sample")
    
    percent_var <- round(pca_res$variance, 2)
    
  } else {
    percent_var <- pca_res$percent_var
    
    df <- pca_res$df
    df$sample <- rownames(df)
  }
  
  # Sample names
  if (add_ids == TRUE) names <- pca_res$yvars
  if (add_ids == FALSE) names <- 0
  
  # Axis labels
  if (pc_var == TRUE) {
    x_nr <- as.integer(sub("PC", "", x))
    y_nr <- as.integer(sub("PC", "", y))
    x_lab <- paste0(x, " (", percent_var[x_nr], "% of variance)")
    y_lab <- paste0(y, " (", percent_var[y_nr], "% of variance)")
  } else {
    x_lab <- x
    y_lab <- y
  }
  
  # Create the base plot
  p <- ggplot(data = df) +
    aes(x = .data[[x]], y = .data[[y]])
  
  # Color aesthethic
  if (!is.null(col)) {
    p <- p +
      aes(color = .data[[col]]) +
      scale_color_brewer(palette = "Dark2")
  }
  
  # Fill aesthetic
  if (!is.null(fill)) {
    p <- p +
      aes(fill = .data[[fill]]) +
      scale_fill_brewer(palette = "Dark2")
    if (is.null(shape)) pt_shape <- 21
  }
  
  # Shape aesthetic
  if (!is.null(shape)) {
    n_shapes <- length(unique(df[[shape]]))
    shapes <- c(21:25)[1:n_shapes]
    
    p <- p +
      aes(shape = .data[[shape]]) +
      scale_shape_manual(values = shapes) +
      guides(fill = guide_legend(override.aes = list(shape = 21)))
  }
  
  # Add points
  if (interactive == FALSE) {
    if (is.null(shape)) {
      p <- p + geom_point(size = pt_size, shape = pt_shape)
    } else {
      if (!is.null(shape)) p <- p + geom_point(size = pt_size)
    }
  } else {
    if (!require(ggiraph)) install.packages("ggiraph")
    if (is.null(shape)) {
      p <- p + geom_point_interactive(aes(tooltip = sample),
                                      size = pt_size, shape = pt_shape)
    } else {
      p <- p + geom_point_interactive(aes(tooltip = sample), size = pt_size)
    }
  }
  
  # Plot polishing
  p <- p +
    labs(x = x_lab, y = y_lab, title = title) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5,
                                    color = "grey20"))
  
  # Add sample IDs
  if (add_ids == TRUE) {
    if (!require(ggrepel)) install.packages("ggrepel")
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = sample), max.overlaps = 20, point.padding = 3
        )
  }
  
  return(p)
  #if (interactive == FALSE) return(p) else girafe(ggobj = p)
}


# Function to create a PCA biplot
pca_biplot <- function(
    pca_res,                      # PCA results after running `PCA_run()`
    x = "PC1",                    # Principal component to plot on the x-axis
    y = "PC2",                    # Principal component to plot on the y-axis
    col = "Treatment",            # Vary point color by this variable from the metadata
    shape = "Irrigation",         # Vary point shape by this variable from the metadata
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
    ggtitle(title) +
    theme_bw()
  
  return(p)
}

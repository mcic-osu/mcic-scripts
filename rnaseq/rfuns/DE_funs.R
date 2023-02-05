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
  res <- results(dds, contrast = c(fac, comp), tidy = TRUE) |>
    as_tibble() |>
    rename(gene = row,
           lfc = log2FoldChange,
           mean = baseMean) |>
    mutate(group1 = comp[1],
           group2 = comp[2],
           contrast = paste0(comp, collapse = "_")) |>
    select(-lfcSE, -stat) |>
    arrange(padj)

  # Include mean normalized counts
  if (!is.null(count_df)) {
    fcount_df <- count_df |>
      filter(.data[[fac]] %in% comp) |>
      group_by(gene, tissue) |>
      summarize(mean = mean(count), .groups = "drop") |>
      pivot_wider(id_cols = gene, values_from = mean, names_from = tissue)
    colnames(fcount_df)[2:3] <- c("mean_A", "mean_B")
    fcount_df <- fcount_df |> mutate(mean = (mean_A + mean_B) / 2)
    res <- res |> left_join(fcount_df, by = "gene")

    # Determine whether a gene is DE
    res <- res |>
      mutate(isDE = ifelse(
        padj < p_tres & abs(lfc) > lfc_tres & (mean_A > mean_tres | mean_B > mean_tres),
        TRUE, FALSE)
    )
  }

  # Determine whether a gene is DE
  res <- res |> mutate(
    isDE = ifelse(padj < p_tres & abs(lfc) > lfc_tres & mean > mean_tres, TRUE, FALSE)
    )

  # Report
  nsig <- sum(res$isDE, na.rm = TRUE)
  message(comp[1], " vs ", comp[2], " - Nr DEGs: ", nsig)

  # Only take significant genes
  if (sig_only == TRUE) res <- res |> filter(isDE == TRUE)

  # Add gene annotation
  if (!is.null(annot)) {
    res <- res |>
      left_join(select(annot, gene, gene_name, description), res, by = "gene") |>
      arrange(padj)
  }

  return(res)
}


# Function to provide shrunken LFC estimates
shrink_lfc <- function(
    dds,                # DESeq object
    fac,                # focal factor (column name)
    comp,               # vector of two with focal factor levels (column values)
    lfc_tres = 1,       # If set to 0, original p-values will be used; otherwise, LFC-based s-values
    p_tres = 0.05,       # P-value threshold
    s_tres = 0.005)     # S-value threshold
  {

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
    select(-lfcSE) |>
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
norm_counts <- function(dds,
                        meta_df = NULL,
                        return_matrix = FALSE) {

  # Normalize the counts
  dds <- estimateSizeFactors(dds)
  counts_raw <- dds@assays@data$counts
  count_mat <- sweep(counts_raw, 2, sizeFactors(dds), "/")

  if (return_matrix == TRUE) return(count_mat)

  # Get metadata from dds if not provided
  if (is.null(meta_df)) {
    message("Extracting metadata from dds")
    meta_df <- colData(dds) |>
      as.data.frame() |>
      rownames_to_column("sample")
  }

  # Get a long-format count df with metadata
  count_df <- as.data.frame(count_mat) |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "sample", values_to = "count") |>
    inner_join(meta_df, by = "sample")

  return(count_df)
}



# Function to make an MA plot
pMA <- function(
    deseq_results,               # should be DE results df from extract_DE()
    rm_padj_na = TRUE,
    ptsize_nonsig = 1,
    ptsize_sig = 1,
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
                    Description: {description}"
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
               fill = "grey50",
               color = "grey50",
               size = ptsize_nonsig,
               alpha = alpha_nonsig) +
    geom_point(data = dplyr::filter(deseq_results, isDE == TRUE),
               aes(x = mean,
                   y = lfc,
                   text = glue::glue(glue_string)),
               shape = 21,
               fill = "blue",
               color = "blue",
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
          plot.title = element_text(hjust = 0.5))

  if (!is.null(contrast)) {
    if (length(fcontrast) > 1) p <- p + facet_wrap(vars(contrast))
  }

  if (interactive) ggplotly(p, tooltip = "text") else return(p)
}


# Volcano plot
pvolc <- function(DE_df,
                  sig_only = TRUE,
                  contrasts = NULL,
                  cols = NULL,
                  interactive = FALSE,
                  facet_scales = "fixed") {

  fcontrasts <- contrasts

  # Remove genes with NA as the adj. pvalue
  DE_df <- DE_df |> filter(!is.na(padj))

  # Remove non-significant genes
  if (sig_only == TRUE) DE_df <- DE_df |> filter(isDE == TRUE)

  # Get log10 of padj to remove Infinite values
  DE_df <- DE_df |>
    mutate(padj = -log10(padj)) |>
    filter(!is.infinite(padj))

  # Only take results for the the focal contrast
  if (!is.null(fcontrasts) && fcontrasts[1] != "all") {
    message("Selecting only the focal contrasts...")
    DE_df <- DE_df |> filter(contrast %in% fcontrasts)
  }

  # Interactive text
  if ("gene_name" %in% colnames(DE_df)) {
    glue_string <- "Gene ID: {gene}
                    Gene name: {gene_name}
                    Description: {description}"
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
    geom_point(data = filter(DE_df, isDE == FALSE),
               fill = "grey80",
               size = 3, shape = 21, color = "grey40", alpha = 0.3) +
    geom_point(data = filter(DE_df, isDE == TRUE),
               aes(fill = contrast),
               size = 3, shape = 21, color = "grey20", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = xlab,
         y = ylab) +
    guides(fill = "none") +
    theme(panel.grid.minor = element_blank())

  if (is.null(cols)) {
    p <- p + scale_fill_brewer(palette = "Dark2")
  } else {
    p <- p + scale_fill_manual(values = cols)
  }

  # When no focal contrast is specified, show all with a facet
  if (length(fcontrasts) > 1 || fcontrasts[1] == "all") {
    p <- p + facet_wrap(vars(contrast),
                        nrow = 1,
                        scales = facet_scales)
  }

  if (interactive == TRUE) ggplotly(p, tooltip = "text")

  return(p)
}


# Heatmap plot showing abundances for multiple/many genes
pheat <- function(genes,
                  count_mat,
                  meta_df,
                  groups,
                  show_rownames = TRUE,
                  show_colnames = FALSE,
                  cluster_rows = TRUE,
                  logtrans = TRUE,
                  annotation_colors = NULL,
                  id_labsize = 10,
                  ...) {

  library(pheatmap)

  # Select groups and genes
  fmeta <- meta_df[, groups, drop = FALSE]
  #fmeta <- fmeta |> mutate(across(everything(), as.character))

  ## Arrange metadata according to the columns with included factors
  if (length(groups) == 1) fmeta <- fmeta |> arrange(.data[[groups[1]]])
  if (length(groups) == 2) fmeta <- fmeta |> arrange(.data[[groups[1]]],
                                                      .data[[groups[2]]])
  if (length(groups) == 3) fmeta <- fmeta |> arrange(.data[[groups[1]]],
                                                      .data[[groups[2]]],
                                                      .data[[groups[3]]])

  ## Select and arrange count matrix
  fcount_mat <- count_mat[match(genes, rownames(count_mat)),
                          match(rownames(fmeta), colnames(count_mat))]
  fcount_mat <- as.matrix(fcount_mat)

  ## Log-transform
  if (logtrans == TRUE) {
    fcount_mat <- log10(fcount_mat)
    fcount_mat[fcount_mat == -Inf] <- 0
  }

  ## If few features are included, reduce the cell (row) height
  cellheight <- ifelse(length(genes) > 20, NA, 20)
  id_labsize <- ifelse(length(genes) > 40, 6, id_labsize)

  ## Truncate long taxon names
  row.names(fcount_mat) <- str_trunc(row.names(fcount_mat),
                                     width = 20, ellipsis = "")

  ## Function to create the plot
  pheatmap(fcount_mat, annotation_col = fmeta,
           cluster_rows = cluster_rows, cluster_cols = FALSE,
           show_rownames = show_rownames, show_colnames = show_colnames,
           annotation_colors = annotation_colors,
           cellheight = cellheight,
           fontsize = 9, fontsize_row = id_labsize, cex = 1,
           ...)
}

# Boxplot showing abundances for a single gene
pbox <- function(gene,
                 count_df,
                 x_by,
                 col_by = NULL,
                 annot = NULL,
                 cols = NULL,
                 log_scale = FALSE,
                 theme_base_size = 13,
                 save_plot = FALSE,
                 plotdir = "results/figures/geneplots") {

  fgene <- gene

  if (!is.null(annot)) {
    g_descrip <- annot |>
      filter(gene == fgene) |>
      pull(description) |>
      str_trunc(width = 50)

    if ("gene_name" %in% colnames(annot)) {
      g_name <- annot |> filter(gene == fgene) |> pull(gene_name)
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
  count_df <- count_df |> filter(gene %in% fgene)

  # Make the plot
  p <- ggplot(count_df) +
    aes(y = count, x = .data[[x_by]], color = .data[[col_by]]) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0),
               size = 1.5) +
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

  if (x_by == col_by) p <- p + guides(color = "none")

  if (log_scale == TRUE) {
    p <- p +
      scale_y_log10(labels = scales::comma) +
      labs(y = "Normalized count (log10-scale)")
  } else {
    p <- p +
      scale_y_continuous(labels = scales::comma,
                         limits = c(0, NA),
                         expand = expansion(mult = c(0.01, 0.02))) +
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

  print(p)
  return(p)
}

# Wrapper function to combine 4 single-gene boxplots
p4box <- function(genes, count_df,
                  annot = NULL, cols = NULL,
                  x_by = "tissue", col_by = "tissue", ...) {

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

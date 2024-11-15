# A function to create taxonomic abundance barplots for metabarcoding data
# ...and a few color vectors to use in this function
# Jelmer Poelstra, last updated 2024-11-15

# Color sets
cols1 <- c("#a74bb4", "#62b54f", "#7064d3", "#b5b348", "#dd6fc5",
           "#4db18a", "#ce416e", "#45aecf", "#d55035", "#7784cb",
           "#cc8e3e", "#ac6090", "#647934", "#df837e", "#9c5736")

# This one is from microViz::distinct_palette(pal = "brewerPlus", add = NA)
cols_brewerplus <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                     "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                     "#FFFF99", "#B15928", "#1ff8ff", "#1B9E77", "#D95F02",
                     "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
                     "#666666", "#4b6a53", "#b249d5", "#7edc45", "#5c47b8",
                     "#cfd251", "#ff69b4", "#69c86c", "#cd3e50", "#83d5af",
                     "#da6130", "#5e79b2", "#c29545", "#532a5a", "#5f7b35",
                     "#c497cf", "#773a27", "#7cb9cb", "#594e50", "#d3c4a8",
                     "#c17e7f")
                     
# This one is from microViz::distinct_palette(pal = "kelly", add = NA)
cols_kelly <- c("#f3c300", "#875692", "#f38400", "#a1caf1", "#be0032",
                "#c2b280", "#848482", "#008856", "#e68fac", "#0067a5",
                "#f99379", "#604e97", "#f6a600", "#b3446c", "#dcd300",
                "#882d17", "#8db600", "#654522", "#e25822", "#2b3d26")

# Function to create a barplot showing taxon abundances
pbar <- function(
    ps = NULL,              # Provide a phyloseq object (ps) with RELATIVE abundances as the input
    taxrank = "Phylum",     # Taxonomic rank to summarize abundance by Or 'Family', 'Genus', etc
    x_var = "Sample",       # What to plot along the x-axis
                            #   'Sample' for indiv. samples, or a column name from sample_data(ps)) (quoted string)
    facet_var = NULL,       # Which column in sample_data to facet by (quoted string)
    facet_var2 = NULL,      # Which column in sample_data to also facet by (quoted string)
    xlab = NULL,            # X-axis label
    abund_tres = 0.01,      # Lump taxa with abundances below this threshold into a category 'other (rare)'
                            #   (use 'NA' for no threshold)
    focal_taxa = NULL,      # Instead of filtering taxa by abundance, use the taxa listed in this vector
    na_to_unknown = TRUE,   # Change 'NA' taxa to a category 'unknown' in the graph
    sort_by_abund = TRUE,   # Sort the taxa by abundance in the graph (rather than alphabetically)
    colors = cols_kelly,    # Vector of colors. Presets are 'cols1', 'cols_brewerplus', and 'cols_kelly' (default)
                            # You can also provide your own vector of colors.
    abund_df = NULL         # Alternative to providing a phyloseq object (ps) as input: an abundance df from abund_stats()
    ) {

  # Compute abundance stats if needed
  if (is.null(abund_df)) {
    abund_df <- abund_stats(
      ps = ps,
      taxrank = taxrank,
      abund_tres = abund_tres,
      focal_taxa = focal_taxa,
      na_to_unknown = na_to_unknown,
      sort_by_abund = sort_by_abund,
      groupby = c(x_var, facet_var, facet_var2)
      )
  }
  
  # Set colors
  ntax <- length(unique(na.omit(abund_df[[taxrank]])))
  if (is.null(colors)) {
    # [Package gives too many installation problems - turning this functionality off]
    #colors <- randomcoloR::distinctColorPalette(k = ntax)
    stop("Error: please specify colors")
  } else {
    colors <- colors[1:ntax]
  }
  
  # Set last color ('other' category) to grey:
  if (any(abund_df[[taxrank]] == "other (rare)")) {
    if (any(abund_df[[taxrank]] == "unknown")) {
      colors[length(colors) - 1] <- "grey80"
      colors[length(colors)] <- "grey60"
    } else {
      colors[length(colors)] <- "grey80"
    }
  } else if (any(abund_df[[taxrank]] == "unknown")) {
    colors[length(colors)] <- "grey60"
  }
  
  # Create the plot
  p <- ggplot(abund_df) +
    aes(x = .data[[x_var]],
        y = Abundance,
        fill = .data[[taxrank]]) +
    geom_col(color = "grey20",
             position = position_stack(reverse = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.005)),
                       labels = scales::label_percent()) +
    scale_fill_manual(values = colors,
                      guide = guide_legend(ncol = 1, reverse = TRUE)) +
    labs(x = xlab) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  
  if (!is.null(facet_var)) {
    if (is.null(facet_var2)) {
      p <- p +
        facet_grid(cols = vars(.data[[facet_var]]),
                   scales = "free_x", space = "free")
    } else {
      p <- p +
        facet_grid(cols = vars(.data[[facet_var]]),
                   rows = vars(.data[[facet_var2]]),
                   scales = "free_x", space = "free")
    }
  }
    
  if (x_var == "Sample") p <- p + theme(axis.text.x = element_text(angle = 270))
  
  print(p)
}

# Helper function to compute per-taxon abundances - used for/in the pbar() function
abund_stats <- function(
    ps,
    taxrank,
    groupby = NULL,
    abund_tres = 0.01,        # NA => no abundance threshold
    focal_taxa = NULL,        # Alternative to abundance threshold, use vector of taxa to keep
    na_to_unknown = TRUE,
    sort_by_abund = TRUE
    ) {
  
  # If using a list of focal taxa, don't use an abundance threshold
  if (!is.null(focal_taxa)) {
    abund_tres <- NA
    focal_taxa <- c(focal_taxa, "unknown")
  }
  
  # Agglomerate by a taxrank
  ps <- tax_glom(ps, taxrank = taxrank, NArm = FALSE)
  
  # Turn the phyloseq object into a dataframe
  df <- psmelt(ps)
  meta <- as(sample_data(ps), "data.frame") |> rownames_to_column("Sample")
  
  # Merge different NA taxa
  if (any(is.na(df[[taxrank]]))) {
    NA_df <- df |>
      filter(is.na(.data[[taxrank]])) |>
      group_by(Sample) |>
      summarize(Abundance = sum(Abundance)) |>
      left_join(meta, by = "Sample")
    NA_df[[taxrank]] <- NA
    df <- df |>
      filter(!is.na(.data[[taxrank]])) |>
      bind_rows(NA_df)
  }
  
  # Change NA taxa to "unknown"
  if (na_to_unknown == TRUE) {
    #levels(df[[taxrank]]) <- c(levels(df[[taxrank]]), "unknown")
    df[[taxrank]][is.na(df[[taxrank]])] <- "unknown"
  }
  
  # If not plotting by sample, compute mean by a grouping variable
  if (!is.null(groupby)) {
    df <- df |>
      group_by_at(c(groupby, taxrank)) |> 
      summarize(Abundance = mean(Abundance), .groups = "drop")
    
    if (groupby[1] == "Sample") {
      df <- df |>
        left_join(meta, by = "Sample") |>
        # Remove duplicated columns!
        select(!contains(".y"))
      colnames(df) <- sub("\\.x$", "", colnames(df))
    }
  }
  
  # Change low-abundance or non-focal taxa to "other"
  if (!is.na(abund_tres)) {
    # Get vector of taxa to be lumped
    to_lump <- df |>
      group_by(.data[[taxrank]]) |>
      summarize(mean_abund = mean(Abundance)) |>
      filter(mean_abund <= abund_tres) |>
      pull(.data[[taxrank]])
    message(length(to_lump), " low-abundance taxa will be lumped into a new category 'other'")
    
  } else if (!is.null(focal_taxa)) {
    to_lump <- df |>
      filter(! .data[[taxrank]] %in% focal_taxa) |>
      pull(.data[[taxrank]])
    message(length(to_lump), " non-focal taxa will be lumped into a new category 'other'")
  }
  
  if (!is.na(abund_tres) || !is.null(focal_taxa)) {
    # Create a df for the lumped taxa, with summed abundance
    other_df <- df |>
      filter(.data[[taxrank]] %in% to_lump) |>
      group_by_at(groupby) |> 
      summarize(Abundance = sum(Abundance))
    
    if (groupby[1] == "Sample") {
      other_df <- other_df |>
        left_join(meta, by = "Sample") |>
        select(!contains(".y"))
      colnames(other_df) <- sub("\\.x$", "", colnames(other_df))
    }
    other_df[[taxrank]] <- "other (rare)"
    
    # Filter out original low-abund taxa and add lumped one
    df <- df |>
      filter(! .data[[taxrank]] %in% to_lump) |>
      bind_rows(other_df)
    
    # Sort ASVs by mean overall abundance
    if (sort_by_abund == TRUE) {
      tax_ord <- df |>
        group_by(.data[[taxrank]]) |>
        summarize(abund = mean(Abundance)) |>
        arrange(-abund) |>
        drop_na() |>
        pull(.data[[taxrank]])
      df[[taxrank]] <- factor(df[[taxrank]], levels = tax_ord)
    }
    
    # Make sure the 'other' category is the last factor level
    df[[taxrank]] <- fct_relevel(df[[taxrank]], "other (rare)", after = Inf)
    
  } else if (sort_by_abund == TRUE) {
    # Sort ASVs by mean overall abundance
    tax_ord <- df |>
      group_by(.data[[taxrank]]) |>
      summarize(abund = mean(Abundance)) |>
      arrange(-abund) |>
      drop_na() |>
      pull(.data[[taxrank]])
    df[[taxrank]] <- factor(df[[taxrank]], levels = tax_ord)
  }
  
  # Change NA taxa to "unknown"
  if (na_to_unknown == TRUE) {
    df[[taxrank]] <- fct_relevel(df[[taxrank]], "unknown", after = Inf)
  }
  
  return(df)
}

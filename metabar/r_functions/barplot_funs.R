# Functions to create taxonomic abundance barplots for metabarcoding data

# Nice set of colors for taxon barplots
mycols <- c("#a74bb4", "#62b54f", "#7064d3", "#b5b348", "#dd6fc5",
            "#4db18a", "#ce416e", "#45aecf", "#d55035", "#7784cb",
            "#cc8e3e", "#ac6090", "#647934", "#df837e", "#9c5736")

# Function to create a barplot showing taxon abundances
pbar <- function(
    ps = NULL,              # Either provide a phyloseq object (ps) or an abundance df from abund_stats()
    abund_df = NULL,        # Either provide a phyloseq object (ps) or an abundance df from abund_stats()
    taxrank = "Phylum",     # Or 'Family', 'Genus', etc
    x_var = "Sample",       # What to plot along the x-axis ('Sample' for indiv. samples, or provide a column name from sample_data(ps)) (quoted string)
    facet_var = NULL,       # What to facet by (quoted string)
    facet_var2 = NULL,      # Second variable to facet by (quoted string)
    xlab = NULL,            # X-axis label
    abund_tres = 0.01,      # Lump taxa with abundances below this threshold into a category 'other (rare)' (use 'NA' for no threshold)
    focal_taxa = NULL,      # Instead of filtering taxa by abundance, use the taxa listed in this vector
    na_to_unknown = TRUE,   # Change 'NA' taxa to a category 'unknown' in the graph
    sort_by_abund = TRUE,   # Sort the taxa by abundance in the graph (rather than alphabetically)
    colors = mycols         # Provide a set of colors
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
    colors <- randomcoloR::distinctColorPalette(k = ntax)
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
                      guide = guide_legend(reverse = TRUE)) +
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

# Function to compute per-taxon abundances - used for/in the pbar() function
abund_stats <- function(
    ps,
    taxrank,
    groupby = NULL,
    abund_tres = 0.01,        # NA => no abundance threshold
    focal_taxa = NULL,        # Alternative to abundance threshold, use vector of taxa to keep
    na_to_unknown = TRUE,
    sort_by_abund = TRUE
    ) {
  # ps <- ps_prop_filt; groupby = c("Sample", "INFECTION"); abund_tres = 0.01; focal_taxa = NULL
  # na_to_unknown = TRUE; sort_by_abund = TRUE; taxrank = "Family"
  
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

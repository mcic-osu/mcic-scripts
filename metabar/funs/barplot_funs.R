## Barplot showing taxon abundances
pbar <- function(ps, taxrank,
                 xvar = "Sample", facetvar = NULL,
                 xlab = NULL,
                 abund_tres = 0.01,
                 na_to_unknown = TRUE,
                 sort_by_abund = TRUE,
                 cols = NULL) {
  #taxrank="phylum"; abund_tres=NA; cols=mycols; ps=fCS
  #na_to_unknown = TRUE; sort_by_abund = TRUE;
  #xvar = "travel_time"; facetvar = "Age"
  
  df <- abund_stats(ps = ps, taxrank = taxrank, abund_tres = abund_tres,
                    na_to_unknown = na_to_unknown, sort_by_abund = sort_by_abund,
                    groupby = c(xvar, facetvar))
  
  ## Set colors
  ntax <- length(unique(na.omit(df[[taxrank]])))
  if (is.null(cols)) {
    cols <- palette(rainbow(ntax))
  } else {
    cols <- cols[1:ntax]
  }
  
  ## Set last color ('other' category) to grey:
  if (any(df[[taxrank]] == "other (rare)")) {
    if (any(df[[taxrank]] == "unknown")) {
      cols[length(cols) - 1] <- "grey80"
      cols[length(cols)] <- "grey60"
    } else {
      cols[length(cols)] <- "grey80"
    }
  } else if (any(df[[taxrank]] == "unknown")) {
    cols[length(cols)] <- "grey60"
  }
  
  ## Create the plot
  p <- ggplot(df) +
    aes(x = .data[[xvar]],
        y = Abundance,
        fill = .data[[taxrank]]) +
    geom_col(color = "grey20",
             position = position_stack(reverse = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.005)),
                       labels = scales::label_percent()) +
    scale_fill_manual(values = cols,
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = xlab) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  
  ## If plotting by sample, always facet by a grouping variable
  if (xvar == "Sample" | !is.null(facetvar)) {
    p <- p +
      facet_grid(cols = vars(.data[[facetvar]]),
                 scales = "free_x", space = "free")
  }
  if (xvar == "Sample") p <- p + theme(axis.text.x = element_text(angle = 270))
  
  print(p)
}

## Barplot showing taxon abundances
abund_stats <- function(ps, taxrank,
                        groupby = NULL,
                        abund_tres = 0.01,
                        na_to_unknown = TRUE,
                        sort_by_abund = TRUE) {
  #taxrank="phylum"; abund_tres=NA; cols=mycols; ps=ps_prop
  #na_to_unknown = TRUE; sort_by_abund = TRUE; groupby = "Sample"
  
  ## Agglomerate by a taxrank
  ps <- tax_glom(ps, taxrank = taxrank, NArm = FALSE)
  
  ## Turn the phyloseq object into a dataframe
  df <- psmelt(ps)
  
  ## Merge different NA taxa
  if (any(is.na(df[[taxrank]]))) {
    NAs <- df %>%
      filter(is.na(.data[[taxrank]])) %>%
      group_by(Sample) %>%
      summarize(Abundance = sum(Abundance)) %>%
      left_join(sample_data(ps), by = c("Sample" = "sampleID"))
    NAs[[taxrank]] <- NA
    df <- df %>%
      filter(!is.na(.data[[taxrank]])) %>%
      bind_rows(NAs)
  }
  
  ## If not plotting by sample, compute mean by a grouping variable
  if (!is.null(groupby)) {
    df <- df %>%
      group_by_at(c(groupby, taxrank)) %>% 
      #group_by(.data[[taxrank]], .data[[groupby]]) %>%
      summarize(Abundance = mean(Abundance), .groups = "drop")
    
    if (groupby[1] == "Sample")
      df <- df %>% left_join(sample_data(ps), by = c("Sample" = "sampleID"))
  }
  
  ## Change low-abundance taxa to "other"
  if (!is.na(abund_tres)) {
    
    ## Get vector of taxa to be lumped
    to_lump <- df %>%
      group_by(.data[[taxrank]]) %>%
      summarize(mean_abund = mean(Abundance)) %>%
      filter(mean_abund <= abund_tres) %>%
      pull(.data[[taxrank]])
    message(length(to_lump), " taxa will be lumped into a new category 'other'")
    
    ## Create a df for the lumped taxa, with summed abundance
    other <- df %>%
      filter(.data[[taxrank]] %in% to_lump) %>%
      #group_by(.data[[groupby]]) %>%
      group_by_at(groupby) %>% 
      summarize(Abundance = sum(Abundance))
    
    if (groupby[1] == "Sample") {
      other <- other %>% left_join(sample_data(ps), by = c("Sample" = "sampleID"))
    }
    
    other[[taxrank]] <- "other (rare)"
    
    ## Filter out original low-abund taxa and add lumped one
    df <- df %>%
      filter(! .data[[taxrank]] %in% to_lump) %>%
      bind_rows(other)
    
    ## Sort ASVs by mean overall abundance
    if (sort_by_abund == TRUE) {
      tax_ord <- df %>%
        group_by(.data[[taxrank]]) %>%
        summarize(abund = mean(Abundance)) %>%
        arrange(-abund) %>%
        drop_na() %>%
        pull(.data[[taxrank]])
      df[[taxrank]] <- factor(df[[taxrank]], levels = tax_ord)
    }
    
    ## Make sure the 'other' category is the last factor level
    df[[taxrank]] <- fct_relevel(df[[taxrank]], "other (rare)", after = Inf)
    
  } else if (sort_by_abund == TRUE) {
    ## Sort ASVs by mean overall abundance
    tax_ord <- df %>%
      group_by(.data[[taxrank]]) %>%
      summarize(abund = mean(Abundance)) %>%
      arrange(-abund) %>%
      drop_na() %>%
      pull(.data[[taxrank]])
    df[[taxrank]] <- factor(df[[taxrank]], levels = tax_ord)
  }
  
  ## Change NA taxa to "unknown"
  if (na_to_unknown == TRUE) {
    levels(df[[taxrank]]) <- c(levels(df[[taxrank]]), "unknown")
    df[[taxrank]][is.na(df[[taxrank]])] <- "unknown"
    df[[taxrank]] <- fct_relevel(df[[taxrank]], "unknown", after = Inf)
  }
  
  return(df)
}

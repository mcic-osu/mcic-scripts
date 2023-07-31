# Function to plot a tree with ggtree
plot_tree <- function(tree,
                      alignment = NULL,
                      annot = NULL,
                      label_columns = NULL,
                      color_column = NULL,
                      node_labels = TRUE,
                      label_size = 3,
                      fig_file = "tree.png",
                      print_fig = FALSE,
                      plot_msa = FALSE,
                      msa_offset = "auto") {

  if (!is.null(annot)) {
    p <- ggtree(tree) %<+% annot
    
    if (is.null(label_columns)) {
        p <- p + geom_tiplab()
    } else {
        if (!is.na(label_columns[1])) {
            if(!is.null(color_column)) {
                p <- p + geom_tiplab(aes_string(label = label_columns[1],
                                                color = color_column),
                                     align = TRUE, size = label_size)
                if (color_column == label_columns[1]) p <- p + guides(color = "none")
            } else {
                p <- p + geom_tiplab(aes_string(label = label_columns[1]),
                                     size = label_size)
            }
        }
    
        if (!is.na(label_columns[2]))
            p <- p + geom_tiplab(aes_string(label = label_columns[2]),
                                align = TRUE, linetype = "blank", offset = 3.5,
                                size = label_size)

        if (!is.na(label_columns[3]))
            p <- p + geom_tiplab(aes_string(label = label_columns[3]),
                                align = TRUE, linetype = "blank", offset = 5,
                                size = label_size)

    }
  } else {
    # Tree when there is no 'annot' dataframe
    p <- ggtree(tree) + geom_tiplab()
  }
  
  # Node labels
  if (node_labels == TRUE) {
    p <- p +
      geom_text2(aes(subset = !isTip, label = label),
                 color = "grey50", nudge_y = 0.4, nudge_x = -0.05, size = 3)
  }
  
  # Wide margins etc to show full labels
  p <- p +
    theme(plot.margin = margin(0.2, 8, 0.2, 0.2, "cm")) +
    coord_cartesian(clip = "off")

  if (plot_msa == TRUE) {
    # Determine offset for Multiple Sequence Alignment (MSA) visualization
    if (msa_offset == "textlen") {
        text_len <- max(nchar(annot$tree_label))
        tree_dp <- max(node.depth(tree))
        msa_offset <- (text_len * 0.005) + (tree_dp * 0.02)
        if (msa_offset > 0.5) msa_offset <- 0.5
    }

    # Add Multiple Sequence Alignment (MSA) visualization
    if (msa_offset != "auto") {
        msaplot(p, alignment, offset = msa_offset, width = 2)
    } else {
        msaplot(p, alignment, width = 2)
    }
    
    # Plot dimensions with an MSA
    plot_width <- 14
    plot_height <- 7

  } else {
    # Plot dimensions when there is no MSA
    plot_width <- 8
    plot_height <- 10
  }
  
  # Save the plot to file
  ggsave(fig_file, width = plot_width, height = plot_height)
  if (print_fig == TRUE) print(p)
}

# Function to prep the annotation dataframe associated with the tree
prep_annot_blast <- function(annot_file, tree_tips, show_strain) {
  
  colnames_annot <- c("subject", "query", "p_ident", "align_len", "n_mismatch",
                      "n_gap", "q_start", "q_end", "s_start", "s_end",
                      "evalue", "bitscore",
                      "s_taxon", "s_len", "s_fullname")
  
  annot_raw <- read_tsv(annot_file,
                        col_names = colnames_annot,
                        col_types = cols()) %>%
    select(subject, query, s_taxon, s_fullname)
  
  query_rows <- data.frame(subject = unique(annot_raw$query))
  
  annot <- annot_raw %>%
    distinct(subject, .keep_all = TRUE) %>%
    bind_rows(query_rows) %>%
    filter(subject %in% tree_tips) %>% 
    mutate(is_query = ifelse(subject %in% query_rows$subject, TRUE, FALSE),
           is_genome = ifelse(grepl("complete genome", s_fullname), TRUE, FALSE),
           strain_or_isolate = str_extract(s_fullname, "strain .*|isolate .*"),
           strain_or_isolate = gsub("\\w gene.*|,.*|polyprotei", "", strain_or_isolate),
           tree_label = ifelse(
             subject %in% query_rows$subject,
             subject,
             if (show_strain == TRUE) {
               paste0(s_taxon, " - ", strain_or_isolate, " (", subject,")")
             } else {
               paste0(s_taxon, " (", subject,")")
             }),
           tree_label = sub(" - NA", "", tree_label)) %>%
    select(-query, -s_fullname)
  
  return(annot)
}
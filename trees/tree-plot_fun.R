## Function to plot the tree
plot_tree <- function(tree,
                      alignment = NULL,
                      annot = NULL, label_cols = NA, color_col = NULL,
                      fig_file = "tree.png",
                      plot_msa = FALSE, msa_offset = "auto") {

  if (!is.null(annot)) {
    p <- ggtree(tree) %<+% annot
    
    if (is.na(label_cols[1])) p <- p + geom_tiplab(aes_string(label = annot_col))
    
    if (!is.na(label_cols[1])) {
      if(!is.null(color_col)) {
        p <- p + geom_tiplab(aes_string(label = label_cols[1],
                                        color = color_col),
                             align = TRUE)
        if (color_col == label_cols[1]) p <- p + guides(color = "none")
      } else {
        p <- p + geom_tiplab(aes_string(label = label_cols[1]))
      }
    }
  } else {
    p <- ggtree(tree) +
      geom_tiplab()
  }
  
  if (!is.na(label_cols[2])) {
    p <- p + geom_tiplab(aes_string(label = label_cols[2]),
                         align = TRUE, linetype = "blank", offset = 3.5)
  }
  if (!is.na(label_cols[3])) {
    p <- p + geom_tiplab(aes_string(label = label_cols[3]),
                         align = TRUE, linetype = "blank", offset = 5)
  }

  p <- p +
    theme(plot.margin = margin(0.2, 8, 0.2, 0.2, "cm")) +
    coord_cartesian(clip = "off")

  if (plot_msa == TRUE) {
    ## Determine offset for Multiple Sequence Alignment (MSA) visualization
    if (msa_offset == "textlen") {
        text_len <- max(nchar(annot$tree_label))
        tree_dp <- max(node.depth(tree))
        msa_offset <- (text_len * 0.005) + (tree_dp * 0.02)
        if (msa_offset > 0.5) msa_offset <- 0.5
    }

    ## Add Multiple Sequence Alignment (MSA) visualization
    if (msa_offset != "auto") {
        msaplot(p, alignment, offset = msa_offset, width = 2)
    } else {
        msaplot(p, alignment, width = 2)
    }

    plot_width <- 14
    plot_height <- 7

  } else {
    plot_width <- 8
    plot_height <- 10
  }
  
  ## Save the plot to file
  ggsave(fig_file, width = plot_width, height = plot_height)
}

## Function to prep the annotation dataframe associated with the tree
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
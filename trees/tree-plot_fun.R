## Function to prep the annotation dataframe associated with the tree
prep_annot <- function(annot_file, tree_tips, show_strain) {

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
           s_fullname = sub("Bovine viral diarrhea virus", "BVD", s_fullname),
           s_taxon = sub("Bovine viral diarrhea virus", "BVD", s_taxon),
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

## Function to plot the tree
plot_tree <- function(tree, annot, alignment, fig_file,
                      msa_offset) {

  ## Plot the tree
  p <- ggtree(tree) %<+% annot +
    geom_tiplab(aes(label = tree_label, color = is_query)) +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = c("black", "red"), guide = "none") +
    theme(plot.margin = margin(0.2, 1, 0.2, 0.2, "cm")) +
    guides(size = "none", color = "none")

  ## Add Multiple Sequence Alignment (MSA) visualization
  if (msa_offset == "auto") {
      text_len <- max(nchar(annot$tree_label))
      tree_dp <- max(node.depth(tree))
      msa_offset <- (text_len * 0.005) + (tree_dp * 0.02)
      if (msa_offset > 0.5) msa_offset <- 0.5

      message("## Max nr of characters of tip labels: ", text_len)
      message("## Max node depth of the tree: ", tree_dp)
      message("## MSA_offset: ", msa_offset)
  }

  msaplot(p, alignment, offset = msa_offset, width = 2)

  ggsave(fig_file, width = 14, height = 7)
}

# Function to plot a timetree
timetree <- function(
    tree,                         # Tree or treedata object
    mrsd,                         # Most recent sampling (tip) date
    ci_name = NULL,               # Name of the column/element with the dating CIs to plot
    meta = NULL,                  # Dataframe with metadata
    tiplab_column = NULL,         # Column from 'meta' with tiplabels to use (default: tree's own tiplabels)
    color_column = NULL,          # Column from 'meta' to color tiplabels by 
    tiplab_size = 4,              # Font size of the tiplabels
    xlab_size = 12,               # Font size of the x-axis (time scale) labels
    scale_breaks = NULL           # Breaks along the x-axis (time scale)
    ) {
  
  # Determine breaks (time) along the x-axis
  if (!is.null(scale_breaks)) {
    n_breaks <- NULL
  } else {
    scale_breaks <- waiver()
    n_breaks <- 8
  }
  
  # Determine the total size of the tree
  if (class(tree) == "treedata") {
    tree_size <- sum(tree@phylo$edge.length)
  } else {
    tree_size <- sum(tree$edge.length)
  }
  
  # Base tree
  p <- ggtree(tree, layout = "rectangular", mrsd = mrsd)
  
  # Add annotation dataframe if provided
  if (!is.null(meta)) {
    meta$tiplab <- meta[[tiplab_column]]
    tiplab_column <- "tiplab"
    p <- p %<+% meta
  }
  
  # Add tip labels
  if (!is.null(tiplab_column)) {
    message("# Using custom tip labels...")
    p <- p + suppressWarnings(
      geom_tiplab(aes_string(color = color_column, label = tiplab_column),
                  align = TRUE, linesize = 0,
                  size = tiplab_size, fontface = "bold")
    )
  } else {
    message("# Using default tip labels...")
    p <- p + suppressWarnings(
      geom_tiplab(aes_string(color = color_column),
                  align = TRUE, linesize = 0,
                  size = tiplab_size, fontface = "bold")
    )
  }
  
  # Add dating confidence intervals (CIs)
  if (!is.null(ci_name)) {
    p <- p +
      geom_range(range = ci_name,
                 color = "red", size = 1, alpha = 0.3)
  }
  
  # Formatting - scales etc
  p <- p +
    scale_color_brewer(palette = "Dark2") +
    geom_rootedge(rootedge = tree_size / 50) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(breaks = scale_breaks, n.breaks = n_breaks,
                       labels = scales::comma)
  
  # Formatting - theme
  p <- p +
    theme_tree2() +
    theme(plot.margin = margin(0.2, 3, 0.2, 0.75, "cm"),
          axis.text.x = element_text(size = xlab_size, color = "grey50"),
          axis.line.x.bottom = element_line(color = "grey50"),
          axis.ticks.x.bottom = element_blank(),
          panel.grid.major.x = element_line(linetype = "longdash"),
          legend.position = "top")
  
  suppressWarnings(print(p))
}

# Function to convert a BactDating result object to a treedata object
# for plotting with ggtree.
bactdate2treedata <- function(bactdate_obj) {
  
  tree_list <- as.treedata.resBactDating(bactdate_obj)
  
  tree_data <- methods::new(
    'treedata',
    phylo = tree_list[[1]],
    data = as_tibble(as.data.frame(tree_list[[2]]))
    )
  
  return(tree_data)
}


# Function to plot a tree with ggtree
ptree <- function(
    tree,                         # Tree or treedata object
    meta = NULL,                  # Dataframe with metadata
    tiplab_column = NULL,         # Column from 'meta' with tiplabels to use (default: tree's own tiplabels)
    color_column = NULL,          # Column from 'meta' to color tiplabels by 
    tiplab_size = 4,              # Font size of the tiplabels
    xlab_size = 12               # Font size of the x-axis labels
  ) {
  
  # Determine the total size of the tree
  if (class(tree) == "treedata") {
    tree_size <- sum(tree@phylo$edge.length)
  } else {
    tree_size <- sum(tree$edge.length)
  }
  
  # Base tree
  p <- ggtree(tree, layout = "rectangular")
  
  # Add annotation dataframe if provided
  if (!is.null(meta)) {
    meta$tiplab <- meta[[tiplab_column]]
    tiplab_column <- "tiplab"
    p <- p %<+% meta
  }
  
  # Add tip labels
  if (!is.null(tiplab_column)) {
    message("# Using custom tip labels...")
    p <- p + suppressWarnings(
      geom_tiplab(aes_string(color = color_column, label = tiplab_column),
                  align = TRUE, linesize = 0,
                  size = tiplab_size, fontface = "bold")
    )
  } else {
    message("# Using default tip labels...")
    p <- p + suppressWarnings(
      geom_tiplab(aes_string(color = color_column),
                  align = TRUE, linesize = 0,
                  size = tiplab_size, fontface = "bold")
    )
  }
  
  # Formatting - scales etc
  p <- p +
    scale_color_brewer(palette = "Dark2") +
    geom_rootedge(rootedge = tree_size / 50) +
    coord_cartesian(clip = "off")
  
  # Formatting - theme
  p <- p +
    theme(plot.margin = margin(0.2, 3, 0.2, 0.75, "cm"),
          axis.text.x = element_text(size = xlab_size, color = "grey50"),
          axis.line.x.bottom = element_line(color = "grey50"),
          axis.ticks.x.bottom = element_blank(),
          panel.grid.major.x = element_line(linetype = "longdash"),
          legend.position = "top")
  
  suppressWarnings(print(p))
}


# Function to plot a tree with ggtree, optionally with an MSA
# NOTE 2023-10-25: OUTDATED
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

# Packages
if (!"tidyverse" %in% installed.packages()) install.packages("tidyverse")
if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")
if (!"KEGGREST" %in% installed.packages()) BiocManager::install("KEGGREST")

# KEGG DATABASE FUNCTIONS ------------------------------------------------------
# Get the genes belonging to a certain KEGG pathway
pw2genes <- function(pathway_id) {
  print(pathway_id)
  
  pw <- keggGet(pathway_id)
  if (is.null(pw[[1]]$GENE)) return(NA)
  
  genes_init <- pw[[1]]$GENE[c(TRUE, FALSE)]
  genes <- unlist(lapply(
    X = strsplit(genes_init, split = ";", fixed = TRUE),
    FUN = function(x) x[1]
    ))
  
  return(genes)
}

# Get the KOs belonging to a certain KEGG pathway
pw2kos <- function(pathway_id) {
  cat(pathway_id, "... ")
  
  pw <- keggGet(pathway_id)
  cat(pw[[1]]$NAME[[1]])
  
  if (is.null(pw[[1]]$GENE)) {
    # If no KOs are found
    cat(" - No KOs found!\n")
    ko_df <- tibble(
      pathway = pathway_id,
      pathway_nr = sub("[[:alpha:]]+", "", pathway_id),
      KO = NA
      )
  } else {
    # If KOs are found
    kos <- pw[[1]]$GENE[c(FALSE, TRUE)]
    kos <- unique(sub(".*KO:(K\\d+).*", "\\1", x = kos))
    cat(" - Nr of KOs:", length(kos), "\n")
    
    ko_df <- tibble(
      pathway = pathway_id,
      pathway_nr = sub("[[:alpha:]]+", "", pathway_id),
      KO = kos
      )
  }
  
  return(ko_df)
}

# Get the KO numbers associated with a module from the KEGG database
mod2kos <- function(idx, modules) {
  module <- names(modules)[idx]
  mod_vec <- keggGet(module)[[1]]$ORTHOLOGY
  mod_df <- data.frame(
    KO_id = names(mod_vec),
    KO_descrip = mod_vec,
    module = module, module_descrip = modules[idx],
    row.names = NULL
    )
  return(mod_df)
}

# Function to get a KEGG pathway (ko-ID) associated with a KEGG K-term
ko2pw <- function(ko_id, outdir) {
  cat("# KO ID:", ko_id, " ")
  tryCatch(
    {
      keggres <- keggGet(ko_id)
      
      if (!is.null(keggres[[1]]$PATHWAY)) pw <- names(keggres[[1]]$PATHWAY) else pw <- NA
      if (!is.null(keggres[[1]]$PATHWAY)) pw_descr <- keggres[[1]]$PATHWAY else pw_descr <- NA
      if (!is.null(keggres[[1]]$MODULE)) mod <- names(keggres[[1]]$MODULE) else mod <- NA
      if (!is.null(keggres[[1]]$MODULE)) mod_descr <- keggres[[1]]$MODULE else mod_descr <- NA
      
      pathway_df <- data.frame(
        KO_id = ko_id,
        symbol = ifelse(is.null(keggres[[1]]$SYMBOL), NA, keggres[[1]]$SYMBOL),
        name = ifelse(is.null(keggres[[1]]$NAME), NA, keggres[[1]]$NAME),
        pw_mod = pw,
        pw_mod_descr = pw_descr,
        row.names = NULL
      )
      cat("// n pathways:", nrow(pathway_df))
      
      if (!is.na(mod[1])) {
        module_df <- data.frame(
          KO_id = ko_id,
          symbol = ifelse(is.null(keggres[[1]]$SYMBOL), NA, keggres[[1]]$SYMBOL),
          name = ifelse(is.null(keggres[[1]]$NAME), NA, keggres[[1]]$NAME),
          pw_mod = mod,
          pw_mod_descr = mod_descr,
          row.names = NULL
        )
        kegg_df <- bind_rows(pathway_df, module_df)
        cat(" // n modules:", nrow(module_df), "\n")
      } else {
        kegg_df <- pathway_df
        cat("\n")
      }
      if (nrow(kegg_df) > 0) return(kegg_df) else return(NULL)
    },
    error = function(cond) {
      message("keggGet failure for ", cond)
      return(NULL)
    }
  )
}

# Get NCBI ID for a gene ID
get_NCBI_id <- function(geneID) {
  geneID_NCBI <- entrez_search(db = "gene", term = geneID)$ids
  message(geneID, " - ", geneID_NCBI)
  if (is_empty(geneID_NCBI)) geneID_NCBI <- NA
  return(data.frame(geneID, geneID_NCBI))
}

# Function to get the description (technically: 'Name') of a KEGG pathway,
# given its pathway ID ('ko' or 'map' IDs).
# Needs tryCatch because some IDs fail (example of a failing pathway ID: "ko01130"),
# and the function needs to keep running.
# Use like so: `pw_descrip <- map_dfr(.x = kegg_ids, .f = kegg_descrip)`
pw2descrip <- function(kegg_pathway) {
  message(kegg_pathway)
  tryCatch( {
    description <- keggGet(kegg_pathway)[[1]]$NAME
    return(data.frame(term = kegg_pathway, description))
  }, error = function(cond) {
    message("keggGet failure")
    return(NULL)
  }
  )
}

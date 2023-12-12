# Packages
if (! "tidyverse" %in% installed.packages()) install.packages("tidyverse")
if (! "BiocManager" %in% installed.packages()) install.packages("BiocManager")
if (! "KEGGREST" %in% installed.packages()) BiocManager::install("KEGGREST")

# KEGG DATABASE FUNCTIONS ------------------------------------------------------
# Get the genes belonging to a certain KEGG pathway
p2genes <- function(pathway_id) {
  print(pathway_id)
  
  pw <- keggGet(pathway_id)
  if (is.null(pw[[1]]$GENE)) return(NA)
  pw2 <- pw[[1]]$GENE[c(TRUE, FALSE)]
  pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = TRUE),
                       function(x) x[1]))
  
  return(pw2)
}

# Get the genes belonging to a certain KEGG pathway
pw2kos <- function(pathway_id) {
  cat(pathway_id, "... ")
  
  pw <- keggGet(pathway_id)
  cat(pw[[1]]$NAME[[1]])
  
  if (is.null(pw[[1]]$GENE)) {
    cat(" - No KOs found!\n")
    ko_df <- tibble(pathway = pathway_id,
                    pathway_nr = sub("[[:alpha:]]+", "", pathway_id),
                    KO = NA)
  } else {
    kos <- pw[[1]]$GENE[c(FALSE, TRUE)]
    kos <- unique(sub(".*KO:(K\\d+).*", "\\1", x = kos))
    cat(" - Nr of KOs:", length(kos), "\n")
    
    ko_df <- tibble(pathway = pathway_id,
                    pathway_nr = sub("[[:alpha:]]+", "", pathway_id),
                    KO = kos)
  }
  
  return(ko_df)
}

# Get the KO numbers associated with a module from the KEGG database
mod2kos <- function(idx, modules) {
  module <- names(modules)[idx]
  mod_vec <- keggGet(module)[[1]]$ORTHOLOGY
  mod_df <- data.frame(KO_id = names(mod_vec), KO_descrip = mod_vec,
                       module = module, module_descrip = modules[idx],
                       row.names = NULL)
  return(mod_df)
}

# Function to get a KEGG pathway (ko-ID) associated with a KEGG K-term
# (Example K_term: "K13449")
ko2pw <- function(K_term, outdir) {
  cat("K_term:", K_term, " ")
  
  tryCatch(
    {
      kegg_info <- keggGet(K_term)
      pathway_df <- data.frame(pathway_description = kegg_info[[1]]$PATHWAY) |>
        rownames_to_column("pathway_id") |>
        mutate("K_term" = K_term)
      
      cat(" Nr of pathways:", nrow(pathway_df), "\n")
      
      if(nrow(pathway_df) > 0) {
        pathway_df_file <- file.path(outdir, paste0(K_term, ".txt"))
        write_tsv(pathway_df, pathway_df_file)
        return(pathway_df)
      } else {
        return(NULL)
      }
    },
    error = function(cond) {
      message("keggGet failure")
      message(cond)
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
get_kegg_descrip <- function(kegg_pathway) {
  message(kegg_pathway)
  tryCatch( {
    description <- keggGet(kegg_pathway)[[1]]$NAME
    return(data.frame(kegg_pathway, description))
  }, error = function(cond) {
    message("keggGet failure")
    return(NULL)
  }
  )
}

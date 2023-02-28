# Get mean abundances of KOs that belong to a specified module or pathway
# The `module` argument can also represent a pathway
getmod <- function(module, mat, meta, kegg_map) {
  KO_ids <- filter(kegg_map, pw_mod == module) %>% pull(KO_id)
  KO_ids_present <- KO_ids[KO_ids %in% rownames(mat)]
  
  if (length(KO_ids_present) > 1) {
    # Take the mean count across all KOs in the focal module/pathway
    mod_mean <- round(colSums(mat[KO_ids_present, ]) / length(KO_ids_present))
    mat_foc <- t(as.matrix(mod_mean))
    rownames(mat_foc) <- module
  } else if (length(KO_ids_present) == 1) {
    mat_foc <- t(as.matrix(mat[KO_ids_present, ]))
    rownames(mat_foc) <- module
  } else {
    message("No KOs found for module ", module)
  }
  return(mat_foc)
}

# Function to get the KO numbers associated with a module,
# straight from the KEGG database
get_kos_in_mod <- function(idx, modules) {
  module <- names(modules)[idx]
  mod_vec <- keggGet(module)[[1]]$ORTHOLOGY
  mod_df <- data.frame(KO_id = names(mod_vec), KO_descrip = mod_vec,
                       module = module, module_descrip = modules[idx],
                       row.names = NULL)
  return(mod_df)
}

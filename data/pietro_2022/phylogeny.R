# packages necessary

if (!require(rotl)) {
  install.packages("rotl")
}

if (!require(ape)) {
  install.packages("ape")
}

# phylogeny function ----
phylo.tree <- function(data) {
  
  taxa <- tnrs_match_names(unique(data$species), 
                           context_name = "Animals",
                           do_approximate_matching = T)
  
  # correct errors with Dermestes maculatus and Drosophila melanogaster ----
  if(any(taxa$search_string == "dermestes_maculatus")) {
    taxa <- update(taxa, 
                   taxon_name = "dermestes_maculatus",
                   new_ott_id = inspect(taxa, taxon_name = "dermestes_maculatus")[1, 4])
  }
  
  if(any(taxa$search_string == "drosophila_melanogaster")){
    taxa <- update(taxa, 
                   taxon_name = "drosophila_melanogaster",
                   new_ott_id = inspect(taxa, taxon_name = "drosophila_melanogaster")[1, 4])
  }
  
  taxon_map <- structure(taxa$search_string,
                         names = taxa$unique_name)
  
  tr <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))])
  
  otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
  tr$tip.label <- taxon_map[otl_tips]
  tr$node.label <- NULL
  
  tr <- compute.brlen(tr)
  cor <- vcv(tr, cor = T)
  
  tr_corrected <- tr
  tr_corrected$tip.label <- str_replace(str_to_title(taxon_map[otl_tips]), "_", " ")
  assign("tr", tr_corrected, envir = .GlobalEnv)
  assign("cor", cor, envir = .GlobalEnv)
}

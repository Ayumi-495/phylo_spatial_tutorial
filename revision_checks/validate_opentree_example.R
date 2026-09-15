#!/usr/bin/env Rscript

# Query the current OpenTree service and confirm the tutorial's safe example.

suppressPackageStartupMessages(library(rotl))

all_taxa <- c(
  "Escherichia coli", "Chlamydomonas reinhardtii",
  "Drosophila melanogaster", "Arabidopsis thaliana",
  "Rattus norvegicus", "Mus musculus", "Cavia porcellus",
  "Xenopus laevis", "Saccharomyces cerevisiae", "Danio rerio"
)

full_query <- rotl::tnrs_match_names(
  names = all_taxa, do_approximate_matching = FALSE
)
stopifnot(
  identical(full_query$unique_name[full_query$search_string == "escherichia coli"],
            "Escherichia coli"),
  identical(full_query$ott_id[full_query$search_string == "escherichia coli"],
            474506L)
)
full_tree <- suppressWarnings(rotl::tol_induced_subtree(
  ott_ids = full_query$ott_id, label_format = "name"
))
stopifnot(any(grepl("^mrcaott", full_tree$tip.label)))

example_taxa <- setdiff(all_taxa, "Escherichia coli")
example_query <- rotl::tnrs_match_names(
  names = example_taxa, do_approximate_matching = FALSE
)
example_tree <- suppressWarnings(rotl::tol_induced_subtree(
  ott_ids = example_query$ott_id, label_format = "name"
))
example_labels <- gsub("_", " ", example_tree$tip.label)

stopifnot(
  !any(grepl("^mrcaott", example_tree$tip.label)),
  setequal(example_labels, example_taxa),
  length(example_labels) == length(example_taxa)
)

cat("OPENTREE_EXAMPLE_PASSED\n")

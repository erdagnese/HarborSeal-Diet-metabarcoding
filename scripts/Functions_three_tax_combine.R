library(tidyverse)

#' Read taxonomy file
#' @param file_path Path to taxonomy file
read_taxonomy <- function(file_path) {
  df <- read.table(file_path, 
                   sep = "\t", 
                   header = TRUE, 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)
  return(df)
}

#' Standardize taxonomy columns
#' @param df Data frame with taxonomy data
#' @param source Source of taxonomy ("local" or "universal")
#' Standardize taxonomy columns with flexible naming
#' @param df Data frame with taxonomy data
#' @param prefix Prefix for taxonomy and confidence columns (e.g., "prey", "full", "db1", "db2")
standardize_taxonomy_cols <- function(df, prefix) {
  df <- df %>%
    dplyr::rename(
      featureid = "Feature ID",
      !!paste0(prefix, "_taxonomy") := "Taxon",
      !!paste0(prefix, "_confidence") := "Confidence"
    )
  
  return(df)
}


#' Find the deepest taxonomic match between two taxonomy strings
#' @param tax1 First taxonomy string
#' @param tax2 Second taxonomy string
find_deepest_match <- function(tax1, tax2) {
  # Handle cases where either taxonomy string is NA
  if (is.na(tax1) || is.na(tax2)) {
    return(NA_character_)
  }
  
  # Split the taxonomy strings into parts
  parts1 <- strsplit(tax1, ";")[[1]]
  parts2 <- strsplit(tax2, ";")[[1]]
  
  # Find the minimum length to avoid out-of-bounds errors
  max_len <- min(length(parts1), length(parts2))
  
  # Initialize a vector to store matching parts
  matching_parts <- c()
  
  # Iterate through the ranks and find matches
  for (i in 1:max_len) {
    if (parts1[i] == parts2[i] && !is.na(parts1[i]) && !is.na(parts2[i])) {
      matching_parts <- c(matching_parts, parts1[i])
    } else {
      break
    }
  }
  
  # If there are matching parts, return them (even if only kingdom level)
  if (length(matching_parts) > 0) {
    return(paste(matching_parts, collapse = ";"))
  } else {
    return(NA_character_)
  }
}

#' Count the number of taxonomic ranks in a taxonomy string
#' @param tax_string Taxonomy string separated by semicolons
count_tax_ranks <- function(tax_string) {
  sum(str_count(tax_string, ";")) + 1 - sum(str_count(tax_string, "NA"))
}

#' Determine the rank of a taxonomy string (kingdom, phylum, class, etc.)
#' @param tax_string Taxonomy string separated by semicolons
get_tax_rank <- function(tax_string) {
  # Handle NA values
  if (is.na(tax_string)) {
    return(NA_character_)
  }
  
  # Define the ranks
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # Count the number of ranks in the taxonomy string
  rank_count <- count_tax_ranks(tax_string)
  
  # Return the corresponding rank
  if (rank_count > 0 && rank_count <= length(ranks)) {
    return(ranks[rank_count])
  } else {
    return(NA_character_)
  }
}

#' Combine prey, full_local, and universal taxonomies and calculate LCAs and ranks
combine_taxonomies <- function(prey_df, full_df, universal_df) {
  # Combine all datasets by featureid
  combined_df <- prey_df %>%
    left_join(full_df, by = "featureid") %>%
    left_join(universal_df, by = "featureid")
  
  # Calculate LCAs and ranks
  combined_df <- combined_df %>%
    rowwise() %>%
    mutate(
      # LCA between all three taxonomies
      lca_all = find_deepest_match(
        find_deepest_match(prey_taxonomy, full_taxonomy),
        universal_taxonomy
      ),
      lca_all_rank = get_tax_rank(lca_all),
      
      # LCA between universal and full_local
      lca_uni_full = find_deepest_match(universal_taxonomy, full_taxonomy),
      lca_uni_full_rank = get_tax_rank(lca_uni_full),
      
      # LCA between full_local and prey
      lca_full_prey = find_deepest_match(full_taxonomy, prey_taxonomy),
      lca_full_prey_rank = get_tax_rank(lca_full_prey),
      
      # LCA between uni_prey and prey
      lca_uni_prey = find_deepest_match(universal_taxonomy, prey_taxonomy),
      lca_uni_prey_rank = get_tax_rank(lca_uni_prey),
      
      # ranks for each of the taxonomy assignments
      prey_rank = get_tax_rank(prey_taxonomy),
      full_local_rank = get_tax_rank(full_taxonomy),
      universal_rank = get_tax_rank(universal_taxonomy),
      
      # Highest confidence and its rank
      highest_confidence = case_when(
        prey_confidence >= full_confidence & prey_confidence >= universal_confidence ~ "prey_confidence",
        full_confidence >= prey_confidence & full_confidence >= universal_confidence ~ "full_confidence",
        universal_confidence >= prey_confidence & universal_confidence >= full_confidence ~ "universal_confidence"
      ),
      highest_conf_rank = case_when(
        highest_confidence == "prey_confidence" ~ get_tax_rank(prey_taxonomy),
        highest_confidence == "full_confidence" ~ get_tax_rank(full_taxonomy),
        highest_confidence == "universal_confidence" ~ get_tax_rank(universal_taxonomy)
      ),
      # Calculate the difference between the highest and second highest confidence
      conf_diff_highest_second = case_when(
        highest_confidence == "prey_confidence" ~ prey_confidence - max(full_confidence, universal_confidence),
        highest_confidence == "full_confidence" ~ full_confidence - max(prey_confidence, universal_confidence),
        highest_confidence == "universal_confidence" ~ universal_confidence - max(prey_confidence, full_confidence)
      )
    ) %>%
    ungroup()
  
  return(combined_df)
}


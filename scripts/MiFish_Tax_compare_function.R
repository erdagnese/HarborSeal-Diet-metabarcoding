library(tidyverse)

# Helper function to clean taxonomy strings
clean_taxonomy <- function(tax) {
  if (is.na(tax)) return(NA_character_)
  # Remove trailing NAs and empty levels
  taxa <- unlist(strsplit(tax, ";"))
  taxa <- taxa[taxa != "NA" & taxa != "" & !is.na(taxa)]
  paste(taxa, collapse = ";")
}

# Helper function to check if species level matches (handles genus changes)
check_species_match <- function(tax1, tax2) {
  if (is.na(tax1) || is.na(tax2)) return(FALSE)
  
  # Clean taxonomy strings
  tax1_clean <- clean_taxonomy(tax1)
  tax2_clean <- clean_taxonomy(tax2)
  
  if (is.na(tax1_clean) || is.na(tax2_clean)) return(FALSE)
  
  taxa1 <- unlist(strsplit(tax1_clean, ";"))
  taxa2 <- unlist(strsplit(tax2_clean, ";"))
  
  # Check if both have at least 7 levels (species level)
  if (length(taxa1) >= 7 && length(taxa2) >= 7) {
    # Check if species names match (level 7)
    if (taxa1[7] == taxa2[7]) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}

# Helper function to find LCA (Lowest Common Ancestor) - vectorized
find_lca <- function(tax1, tax2) {
  map2_chr(tax1, tax2, function(t1, t2) {
    if (is.na(t1) || is.na(t2)) return(NA_character_)
    
    # Clean taxonomy strings first
    t1_clean <- clean_taxonomy(t1)
    t2_clean <- clean_taxonomy(t2)
    
    if (is.na(t1_clean) || is.na(t2_clean)) return(NA_character_)
    
    taxa1 <- unlist(strsplit(t1_clean, ";"))
    taxa2 <- unlist(strsplit(t2_clean, ";"))
    
    common_taxa <- c()
    min_length <- min(length(taxa1), length(taxa2))
    
    for (i in 1:min_length) {
      if (taxa1[i] == taxa2[i]) {
        common_taxa <- c(common_taxa, taxa1[i])
      } else {
        break
      }
    }
    
    if (length(common_taxa) == 0) return(NA_character_)
    
    return(paste(common_taxa, collapse = ";"))
  })
}

# Helper function to get taxonomic rank level - vectorized (cleaned)
get_rank_level <- function(taxonomy) {
  map_dbl(taxonomy, function(tax) {
    if (is.na(tax)) return(0)
    tax_clean <- clean_taxonomy(tax)
    if (is.na(tax_clean)) return(0)
    length(unlist(strsplit(tax_clean, ";")))
  })
}

compare_taxonomy_tables <- function(localTAX, globalTAX) {
  # Preprocess tables
  localTAX <- localTAX %>%
    dplyr::rename("local" = "Taxon", "ConfidenceLoc" = "Confidence")
  
  globalTAX <- globalTAX %>%
    dplyr::rename("global" = "Taxon", "ConfidenceGlob" = "Confidence")
  
  # Combine tables and calculate additional fields
  tax_compare <- left_join(localTAX, globalTAX, by = "feature_id") %>%
    mutate(
      # Clean taxonomy strings
      local_clean = map_chr(local, clean_taxonomy),
      global_clean = map_chr(global, clean_taxonomy),
      
      # Check for species-level matches (handles genus changes)
      species_match = map2_lgl(local, global, check_species_match),
      
      # Calculate LCA using cleaned taxonomies
      lca = find_lca(local_clean, global_clean),
      lca_rank = get_rank_level(lca),
      
      # Extract ranks using cleaned taxonomies
      local_rank = get_rank_level(local_clean),
      global_rank = get_rank_level(global_clean),
      
      # Determine highest confidence
      highest_confidence = case_when(
        ConfidenceLoc >= ConfidenceGlob ~ "local",
        TRUE ~ "global"
      ),
      
      # Initialize final taxonomy and rule
      final_taxonomy = NA_character_,
      rule_applied = NA_character_
    )
  
  # Apply decision rules in specified order
  tax_compare <- tax_compare %>%
    rowwise() %>%
    mutate(
      # -------------------------------------------------------------------------
      # Rule 1: Exact match of clean local and global OR species-level match
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy) & (local_clean == global_clean | species_match),
        local_clean,  # Use local when species match (even if genus differs)
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & (local_clean == global_clean | species_match),
        "Rule 01: Exact or species-level match",
        rule_applied
      ),
      
      # -------------------------------------------------------------------------
      # Rule 2: Local has high confidence (>=0.98) and high rank (>=6)
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy) & 
          ConfidenceLoc >= 0.98 & 
          local_rank >= 6,  # Genus level or higher
        local_clean,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & 
          ConfidenceLoc >= 0.98 & 
          local_rank >= 6,
        "Rule 02: High local confidence and high rank",
        rule_applied
      ),
      
      # -------------------------------------------------------------------------
      # Rule 3: Local rank highest AND local confidence highest
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy) & 
          local_rank > global_rank & 
          highest_confidence == "local",
        local_clean,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & 
          local_rank > global_rank & 
          highest_confidence == "local",
        "Rule 03: Local higher rank and confidence",
        rule_applied
      ),
      
      # -------------------------------------------------------------------------
      # Rule 4: Local rank higher, global confidence higher, but LCA rank >= 5
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy) & 
          local_rank > global_rank & 
          highest_confidence == "global" &
          lca_rank >= 5,  # Family level or higher
        local_clean,  # Take LOCAL (not LCA)
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & 
          local_rank > global_rank & 
          highest_confidence == "global" &
          lca_rank >= 5,
        "Rule 04: Local higher rank, global higher confidence, good LCA - take local",
        rule_applied
      ),
      
      # -------------------------------------------------------------------------
      # Rule 5: Global rank higher, LCA rank is 2, local confidence higher, local rank < 4
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy) & 
          global_rank > local_rank & 
          lca_rank == 2 & 
          highest_confidence == "local" &
          local_rank < 4,
        global_clean,  # Take GLOBAL
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & 
          global_rank > local_rank & 
          lca_rank == 2 & 
          highest_confidence == "local" &
          local_rank < 4,
        "Rule 05: Global higher rank, LCA=2, local confidence higher but low rank - take global",
        rule_applied
      ),
      
      # -------------------------------------------------------------------------
      # Rule 6: Both ranks low (<4), use LCA
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy) & 
          local_rank < 4 & 
          global_rank < 4,
        lca,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & 
          local_rank < 4 & 
          global_rank < 4,
        "Rule 06: Both ranks low (<4) - use LCA",
        rule_applied
      ),
      
      # -------------------------------------------------------------------------
      # NEW RULE 6.5: Both confidences low (<0.9) but good LCA (rank >=5)
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy) & 
          ConfidenceLoc < 0.9 & 
          ConfidenceGlob < 0.9 &
          lca_rank >= 5,  # Good LCA at family level or higher
        lca,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & 
          ConfidenceLoc < 0.9 & 
          ConfidenceGlob < 0.9 &
          lca_rank >= 5,
        "Rule 06.5: Both confidences low but good LCA - use LCA",
        rule_applied
      ),
      
      # -------------------------------------------------------------------------
      # Rule 7: Global has higher rank and higher confidence
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy) & 
          global_rank > local_rank & 
          highest_confidence == "global",
        global_clean,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & 
          global_rank > local_rank & 
          highest_confidence == "global",
        "Rule 07: Global higher rank and confidence",
        rule_applied
      ),
      
      # -------------------------------------------------------------------------
      # Rule 8: Fallback - use higher confidence
      # -------------------------------------------------------------------------
      final_taxonomy = if_else(
        is.na(final_taxonomy),
        case_when(
          ConfidenceLoc >= ConfidenceGlob ~ local_clean,
          TRUE ~ global_clean
        ),
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied),
        "Rule 08: Fallback to higher confidence",
        rule_applied
      )
    ) %>%
    ungroup()
  
  return(tax_compare)
}

# Summary function
summarize_rules <- function(results_df) {
  rule_counts <- results_df %>%
    count(rule_applied) %>%
    rename(Rule = rule_applied, Count = n) %>%
    arrange(Rule)
  
  return(rule_counts)
}
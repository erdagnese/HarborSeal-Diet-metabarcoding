#' Find deepest taxonomic match beyond Eukaryota
#' @param tax1 First taxonomy string
#' @param tax2 Second taxonomy string
# Find deepest taxonomic match beyond Eukaryota
find_deepest_match <- function(tax1, tax2) {
  parts1 <- strsplit(tax1, ";")[[1]]
  parts2 <- strsplit(tax2, ";")[[1]]
  
  # Find matching elements up to the shorter length
  max_len <- min(length(parts1), length(parts2))
  matching_parts <- c()
  
  for(i in 1:max_len) {
    if(parts1[i] == parts2[i]) {
      matching_parts <- c(matching_parts, parts1[i])
    } else {
      break
    }
  }
  
  # If we only match at Eukaryota or have no matches, return NA
  if(length(matching_parts) <= 1) {
    return(NA_character_)
  } else {
    return(paste(matching_parts, collapse = ";"))
  }
}

secondary_taxonomic_comparison <- function(local_results, universal_df, config) {
  # Add universal taxonomy and confidence to local results
  comparison_df <- local_results %>%
    left_join(universal_df %>% 
                select(featureid, universal_taxonomy, universal_confidence),
              by = "featureid")
  
  # 1. Process cases with final_local_taxonomy
  results_with_local <- comparison_df %>%
    filter(!is.na(final_local_taxonomy)) %>%
    mutate(
      # First try to find LCA with universal
      lca_with_prey = mapply(find_deepest_match, prey_taxonomy, universal_taxonomy),
      lca_with_full = mapply(find_deepest_match, full_taxonomy, universal_taxonomy),
      
      # Count ranks in local_final_taxonomy to exclude "Eukaryota" only cases
      local_final_ranks = map_dbl(final_local_taxonomy, count_tax_ranks),
      
      final_taxonomy = case_when(
        # 1. Exact match between final_local and universal
        final_local_taxonomy == universal_taxonomy ~ final_local_taxonomy,
        
        # 2. Trust local if it has high confidence AND has more ranks than just "Eukaryota"
        as.numeric(local_final_confidence) >= config$high_confidence_threshold & 
          local_final_ranks > 1 ~ final_local_taxonomy,
        
        # 3. Eukaryota cases
        # If local_final_taxonomy is "Eukaryota" and universal has high confidence, use universal_taxonomy
        final_local_taxonomy == "Eukaryota" & universal_confidence >= config$high_confidence_threshold ~ universal_taxonomy,
        
        # If universal_taxonomy is "Eukaryota" and local has high confidence, use local_final_taxonomy
        universal_taxonomy == "Eukaryota" & as.numeric(local_final_confidence) >= config$high_confidence_threshold ~ final_local_taxonomy,
        
        # If both are "Eukaryota", retain "Eukaryota"
        final_local_taxonomy == "Eukaryota" & universal_taxonomy == "Eukaryota" ~ "Eukaryota",
        
        # 4. Trust universal if it has very high confidence and is not just "Eukaryota"
        universal_confidence >= config$very_high_confidence_threshold & 
          universal_taxonomy != "Eukaryota" & 
          as.numeric(local_final_confidence) < config$high_confidence_threshold ~ universal_taxonomy,
        
        # 5. If local_final_confidence is higher than universal_confidence and above medium_confidence_threshold,
        # AND the taxonomy has more ranks than just "Eukaryota", use local_full_taxonomy
        as.numeric(local_final_confidence) > universal_confidence & 
          as.numeric(local_final_confidence) >= config$medium_confidence_threshold & 
          local_final_ranks > 1 ~ final_local_taxonomy,
        
        # 6. For needs_review_with_universal cases, try LCA
        !is.na(lca_with_prey) & !is.na(lca_with_full) & nchar(lca_with_prey) > nchar(lca_with_full) ~ lca_with_prey,
        !is.na(lca_with_prey) & !is.na(lca_with_full) & nchar(lca_with_prey) <= nchar(lca_with_full) ~ lca_with_full,
        !is.na(lca_with_prey) ~ lca_with_prey,
        !is.na(lca_with_full) ~ lca_with_full,
        
        # 7. Fallback to needs_review
        TRUE ~ NA_character_
      ),
      
      universal_decision = case_when(
        final_local_taxonomy == universal_taxonomy ~ "exact_match_with_universal",
        as.numeric(local_final_confidence) >= config$high_confidence_threshold & 
          local_final_ranks > 1 ~ "trust_local_high",
        final_local_taxonomy == "Eukaryota" & universal_confidence >= config$high_confidence_threshold ~ "eukaryota_universal_high",
        universal_taxonomy == "Eukaryota" & as.numeric(local_final_confidence) >= config$high_confidence_threshold ~ "eukaryota_local_high",
        universal_confidence >= config$very_high_confidence_threshold & 
          universal_taxonomy != "Eukaryota" ~ "trust_universal_very_high",
        as.numeric(local_final_confidence) > universal_confidence & 
          as.numeric(local_final_confidence) >= config$medium_confidence_threshold & 
          local_final_ranks > 1 ~ "trust_local_medium_high",
        !is.na(lca_with_prey) & !is.na(lca_with_full) & nchar(lca_with_prey) > nchar(lca_with_full) ~ "lca_universal_prey",
        !is.na(lca_with_prey) & !is.na(lca_with_full) & nchar(lca_with_prey) <= nchar(lca_with_full) ~ "lca_universal_full_local",
        !is.na(lca_with_prey) ~ "lca_universal_prey",
        !is.na(lca_with_full) ~ "lca_universal_full_local",
        TRUE ~ "needs_review_with_universal"
      ),
      
      universal_deciding_confidence = case_when(
        universal_decision == "exact_match_with_universal" ~ "universal_confidence",
        universal_decision == "trust_local_high" ~ "local_final_confidence",
        universal_decision == "eukaryota_universal_high" ~ "universal_confidence",
        universal_decision == "eukaryota_local_high" ~ "local_final_confidence",
        universal_decision == "trust_universal_very_high" ~ "universal_confidence",
        universal_decision == "trust_local_medium_high" ~ "local_final_confidence",
        universal_decision == "lca_universal_both_local" ~ "consensus",
        universal_decision == "lca_universal_prey" ~ "consensus",
        universal_decision == "lca_universal_full_local" ~ "consensus",
        TRUE ~ NA_character_
      ),
      
      universal_final_confidence = case_when(
        universal_deciding_confidence == "universal_confidence" ~ universal_confidence,
        universal_deciding_confidence == "local_final_confidence" ~ as.numeric(local_final_confidence),
        TRUE ~ NA_real_
      )
    ) %>%
    select(-lca_with_prey, -lca_with_full, -local_final_ranks)  # Remove temporary columns
  
  # 2. Process cases where local comparison gave lca or needs_review
  results_needing_review <- comparison_df %>%
    filter(is.na(final_local_taxonomy)) %>%
    mutate(
      agrees_with_prey = prey_taxonomy == universal_taxonomy,
      agrees_with_full = full_taxonomy == universal_taxonomy,
      final_taxonomy = case_when(
        # If universal agrees with either local, use that
        agrees_with_prey & prey_confidence >= config$high_confidence_threshold ~ prey_taxonomy,
        agrees_with_full & full_confidence >= config$high_confidence_threshold ~ full_taxonomy,
        # If either local confidence is higher than universal, use local_final_taxonomy
        prey_confidence >= universal_confidence & prey_confidence >= config$high_confidence_threshold ~ prey_taxonomy,
        full_confidence >= universal_confidence & full_confidence >= config$high_confidence_threshold ~ full_taxonomy,
        # If universal has very high confidence and locals disagree
        universal_confidence >= config$very_high_confidence_threshold &
          !agrees_with_prey & !agrees_with_full ~ universal_taxonomy,
        TRUE ~ NA_character_
      ),
      universal_decision = case_when(
        agrees_with_prey & prey_confidence >= config$high_confidence_threshold ~ "universal_agrees_with_prey",
        agrees_with_full & full_confidence >= config$high_confidence_threshold ~ "universal_agrees_with_full",
        prey_confidence >= universal_confidence & prey_confidence >= config$high_confidence_threshold ~ "trust_prey_high_confidence",
        full_confidence >= universal_confidence & full_confidence >= config$high_confidence_threshold ~ "trust_full_high_confidence",
        universal_confidence >= config$very_high_confidence_threshold ~ "trust_universal_very_high",
        local_decision == "lca" ~ "maintain_lca",
        TRUE ~ "needs_manual_review"
      ),
      universal_deciding_confidence = case_when(
        universal_decision == "universal_agrees_with_prey" ~ "prey_confidence",
        universal_decision == "universal_agrees_with_full" ~ "full_confidence",
        universal_decision == "trust_prey_high_confidence" ~ "prey_confidence",
        universal_decision == "trust_full_high_confidence" ~ "full_confidence",
        universal_decision == "trust_universal_very_high" ~ "universal_confidence",
        universal_decision == "maintain_lca" ~ "consensus",
        TRUE ~ NA_character_
      ),
      universal_final_confidence = case_when(
        universal_deciding_confidence == "prey_confidence" ~ prey_confidence,
        universal_deciding_confidence == "full_confidence" ~ full_confidence,
        universal_deciding_confidence == "universal_confidence" ~ universal_confidence,
        TRUE ~ NA_real_
      )
    ) %>%
    select(-agrees_with_prey, -agrees_with_full)
  
  # Combine all results and arrange
  final_results <- bind_rows(
    results_with_local,
    results_needing_review
  ) %>%
    arrange(featureid)
  
  return(final_results)
}
  
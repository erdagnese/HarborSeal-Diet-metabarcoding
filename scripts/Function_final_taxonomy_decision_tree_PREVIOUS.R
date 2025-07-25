# Define the rank ordering from highest to lowest
rank_order <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# Function to find the index of the lowest rank
find_lowest_rank_index <- function(ranks, rank_order) {
  # Convert ranks to their positions in the rank_order vector
  rank_positions <- sapply(ranks, function(r) which(rank_order == r))
  
  # Return the index of the lowest rank (highest position)
  which.max(rank_positions)
}

summarize_rules <- function(combined_df) {
  # Count the number of rows determined by each rule
  rule_counts <- combined_df %>%
    count(rule_applied) %>%
    rename(Rule = rule_applied, Count = n)
  
  # Count the number of rows with NA in final_taxonomy
  na_count <- combined_df %>%
    filter(is.na(final_taxonomy)) %>%
    nrow()
  
  # Add the NA count to the summary
  summary_df <- rule_counts %>%
    bind_rows(tibble(Rule = "NA (No rule applied)", Count = na_count))
  
  return(summary_df)
}

decision_tree <- function(combined_df) {
  combined_df <- combined_df %>%
    rowwise() %>%
    mutate(
      # Initialize final_taxonomy and rule_applied as NA
      final_taxonomy = NA_character_,
      rule_applied = NA_character_,
      
      # Rule 1: If all three taxonomies match exactly, use lca_all
      final_taxonomy = if_else(
        is.na(final_taxonomy) & prey_taxonomy == full_taxonomy & prey_taxonomy == universal_taxonomy,
        lca_all,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & prey_taxonomy == full_taxonomy & prey_taxonomy == universal_taxonomy,
        "Rule 1",
        rule_applied
      ),
      
      # Rule 2: If lca_all_rank is genus or family, use the full_local_taxonomy
      final_taxonomy = if_else(
        is.na(final_taxonomy) & lca_all_rank %in% c("genus", "family"),
        full_taxonomy,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & lca_all_rank %in% c("genus", "family"),
        "Rule 2",
        rule_applied
      ),
      
      # Rule 3: if universal and full local match use the lca_uni_full
      final_taxonomy = if_else(
        is.na(final_taxonomy) & full_taxonomy == universal_taxonomy & lca_uni_full_rank != "kingdom",
        lca_uni_full,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & full_taxonomy == universal_taxonomy & lca_uni_full_rank != "kingdom",
        "Rule 3",
        rule_applied
      ),
      
      # Rule 4: If lca_all_rank is order or class
      final_taxonomy = if_else(
        is.na(final_taxonomy) & lca_all_rank %in% c("order", "class"),
        case_when(
          prey_taxonomy == full_taxonomy ~ case_when(
            highest_confidence == "universal_confidence" ~ universal_taxonomy,
            highest_confidence %in% c("full_confidence", "prey_confidence") ~ lca_full_prey
          ),
          lca_full_prey_rank %in% c("species", "genus") ~ case_when(
            highest_confidence == "full_confidence" ~ full_taxonomy,
            highest_confidence == "prey_confidence" ~ lca_full_prey
          )
        ),
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & lca_all_rank %in% c("order", "class"),
        "Rule 4",
        rule_applied
      ),
      
      # Rule 5: lca_all_rank and lca_full_prey_rank is "kingdom", but lca_uni_full isn't
      final_taxonomy = if_else(
        is.na(final_taxonomy) & lca_uni_prey_rank == "kingdom" & lca_full_prey_rank == "kingdom" & lca_uni_full_rank %in% c("species","genus","family","order", "class", "phylum"),
        case_when(
          universal_confidence > 0.99 ~ universal_taxonomy,
          universal_confidence <= 0.99 ~ lca_uni_full
        ),
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & lca_uni_prey_rank == "kingdom" & lca_full_prey_rank == "kingdom" & lca_uni_full_rank %in% c("species","genus","family","order", "class", "phylum"),
        "Rule 5",
        rule_applied
      ),
      
      # Rule 6: if prey_rank and full_rank are "kingdom" and the universal isn't take the universal
      final_taxonomy = if_else(
        is.na(final_taxonomy) & prey_rank == "kingdom" & full_local_rank == "kingdom" & universal_rank != "kingdom",
        universal_taxonomy,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & prey_rank == "kingdom" & full_local_rank == "kingdom" & universal_rank != "kingdom",
        "Rule 6",
        rule_applied
      ),
      
      # Rule 7: if local_full_confidence is the highest confidence and not at kingdom level
      final_taxonomy = if_else(
        is.na(final_taxonomy) & highest_confidence == "full_confidence" & full_local_rank != "kingdom",
        case_when(
          lca_all_rank != "kingdom" & universal_confidence < 0.99 ~ full_taxonomy,
          lca_full_prey_rank == "species" & full_confidence >= 0.99 ~ lca_full_prey,
          universal_rank %in% c("phylum","kingdom") ~ full_taxonomy,
          lca_uni_full_rank %in% c("class","order","family","genus","species") ~ lca_uni_full,
          lca_uni_full_rank %in% c("phylum","kingdom") & universal_rank %in% c("class","order","family","genus","species") & universal_confidence >= 0.99 ~ universal_taxonomy,
          lca_all_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank == "kingdom" & lca_uni_prey_rank == "kingdom" ~ lca_all,
          lca_all_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_uni_prey_rank == "kingdom" & lca_full_prey_rank != "kingdom" ~ lca_full_prey,
          lca_all_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_uni_prey_rank != "kingdom" & lca_full_prey_rank == "kingdom" ~ lca_uni_prey
        ),
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & highest_confidence == "full_confidence" & full_local_rank != "kingdom",
        "Rule 7",
        rule_applied
      ),
      
      # Rule 8: full_confidence is highest but it's only at kingdom level
      final_taxonomy = if_else(
        is.na(final_taxonomy) & highest_confidence == "full_confidence" & full_local_rank == "kingdom",
        case_when(
          universal_rank != "kingdom" ~ universal_taxonomy,
          universal_rank == "kingdom" ~ lca_uni_full
        ),
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & highest_confidence == "full_confidence" & full_local_rank == "kingdom",
        "Rule 8",
        rule_applied
      ),
      
      # Rule 9: prey confidence is highest but only at kingdom level
      final_taxonomy = if_else(
        is.na(final_taxonomy) & highest_confidence == "prey_confidence",
        case_when(
          prey_rank == "kingdom" & lca_uni_full_rank == "kingdom" & universal_confidence > 0.97 ~ universal_taxonomy,
          prey_rank == "kingdom" & lca_uni_full_rank != "kingdom" ~ lca_uni_full,
          prey_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank != "kingdom" ~ lca_full_prey,
          prey_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank == "kingdom" & lca_uni_prey_rank != "kingdom" ~ lca_uni_prey,
          prey_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank == "kingdom" & lca_uni_prey_rank == "kingdom" ~ lca_all,
          prey_rank != "kingdom" & lca_all_rank != "kingdom" & prey_taxonomy == full_taxonomy & full_confidence > universal_confidence ~ lca_full_prey,
          prey_rank != "kingdom" & lca_all_rank != "kingdom" & prey_taxonomy != full_taxonomy & full_confidence > universal_confidence ~ full_taxonomy,
          prey_rank != "kingdom" & lca_all_rank != "kingdom" & prey_taxonomy == full_taxonomy & full_confidence < universal_confidence ~ lca_all
        ),
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & highest_confidence == "prey_confidence" &
          (prey_rank == "kingdom" & lca_uni_full_rank == "kingdom" & universal_confidence > 0.97  |
             prey_rank == "kingdom" & lca_uni_full_rank != "kingdom" |
             prey_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank != "kingdom" |
             prey_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank == "kingdom" & lca_uni_prey_rank != "kingdom" |
             prey_rank == "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank == "kingdom" & lca_uni_prey_rank == "kingdom" |
             prey_rank != "kingdom" & lca_all_rank != "kingdom" & prey_taxonomy == full_taxonomy & full_confidence > universal_confidence |
             prey_rank != "kingdom" & lca_all_rank != "kingdom" & prey_taxonomy != full_taxonomy & full_confidence > universal_confidence |
             prey_rank != "kingdom" & lca_all_rank != "kingdom" & prey_taxonomy == full_taxonomy & full_confidence < universal_confidence),
        "Rule 9",
        rule_applied
      ),
      
      # Rule 10: universal confidence is the highest 
      final_taxonomy = if_else(
        is.na(final_taxonomy) & highest_confidence == "universal_confidence",
        case_when(
          universal_rank == "species" & full_local_rank %in% c("phylum", "kingdom", "class", "order", "family", "genus") ~ universal_taxonomy,
          universal_rank != "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank == "kingdom" ~ universal_taxonomy,
          universal_rank != "kingdom" & lca_uni_full_rank != "kingdom" & lca_full_prey_rank != "kingdom" & lca_uni_prey_rank != "kingdom" ~ lca_uni_prey,
          find_lowest_rank_index(c(universal_rank, full_local_rank, prey_rank), rank_order) == 1 ~ universal_taxonomy
        ),
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & highest_confidence == "universal_confidence" & 
          (universal_rank == "species" & full_local_rank %in% c("phylum", "kingdom", "class", "order", "family", "genus") |
             universal_rank != "kingdom" & lca_uni_full_rank == "kingdom" & lca_full_prey_rank == "kingdom" |
             universal_rank != "kingdom" & lca_uni_full_rank != "kingdom" & lca_full_prey_rank != "kingdom" & lca_uni_prey_rank != "kingdom" |
             find_lowest_rank_index(c(universal_rank, full_local_rank, prey_rank), rank_order) == 1),
        "Rule 10",
        rule_applied
      ),
      
      # Rule 11: universal rank is kingdom
      final_taxonomy = if_else(
        is.na(final_taxonomy) & universal_rank == "kingdom", 
        case_when(
          lca_full_prey_rank != "kingdom" ~ lca_full_prey,
          lca_full_prey_rank == "kingdom" & full_confidence > prey_confidence ~ full_taxonomy,
          lca_full_prey_rank == "kingdom" & prey_confidence > full_confidence ~ lca_all
        ),
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & universal_rank == "kingdom" &
          (lca_full_prey_rank == "kingdom" |  
             lca_full_prey_rank == "kingdom" & full_confidence > prey_confidence |
             lca_full_prey_rank == "kingdom" & prey_confidence > full_confidence),
        "Rule 11",
        rule_applied
      ),
      
      # Rule 12: Compare ranks of lca_uni_full, lca_uni_prey, and lca_all, at this point we don't want to weight the prey at over universal
      final_taxonomy = if_else(
        is.na(final_taxonomy) & 
          lca_all_rank == "kingdom",
        {
          # Create vectors of ranks and corresponding taxonomies
          ranks <- c(lca_uni_full_rank, lca_uni_prey_rank, lca_all_rank)
          taxonomies <- c(lca_uni_full, lca_uni_prey, lca_all)
          
          # Find the index of the lowest rank
          lowest_index <- find_lowest_rank_index(ranks, rank_order)
          
          # Check if there are ties (multiple ranks with the same lowest rank)
          lowest_rank <- ranks[lowest_index]
          tied_indices <- which(ranks == lowest_rank)
          
          if (length(tied_indices) == 1) {
            # No ties: Return the taxonomy corresponding to the lowest rank
            taxonomies[lowest_index]
          } else {
            # Ties: Check if the taxonomies match
            tied_taxonomies <- taxonomies[tied_indices]
            if (length(unique(tied_taxonomies)) == 1) {
              # Taxonomies match: Return the first tied taxonomy
              tied_taxonomies[1]
            } else {
              # Taxonomies don't match: Leave as NA
              NA_character_
            }
          }
        },
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied) & 
          lca_all_rank == "kingdom",
        "Rule 12",
        rule_applied
      ),
      
      # Default: Leave as NA for now
      final_taxonomy = if_else(
        is.na(final_taxonomy),
        NA_character_,
        final_taxonomy
      ),
      rule_applied = if_else(
        is.na(rule_applied),
        "No rule applied",
        rule_applied
      )
    ) %>%
    ungroup()
  
  return(combined_df)
}
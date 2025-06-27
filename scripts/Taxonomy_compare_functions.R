library(tidyverse)

#' Configuration class for different classification methods
#' @param method The classification method ("naive-bayes", "vsearch", or "blast")
#' @param score_type Type of score used ("confidence", "consensus", or "evalue")
#' @param very_high_confidence_threshold Threshold for very high confidence assignments
#' @param high_confidence_threshold Threshold for high confidence assignments
#' @param medium_confidence_threshold Threshold for medium confidence assignments
#' @param low_confidence_threshold Threshold for low confidence assignments
ClassificationConfig <- function(method, 
                                 score_type,
                                 very_high_confidence_threshold,
                                 high_confidence_threshold,
                                 medium_confidence_threshold,
                                 low_confidence_threshold) {
  list(
    method = method,
    score_type = score_type,
    very_high_confidence_threshold = very_high_confidence_threshold,
    high_confidence_threshold = high_confidence_threshold,
    medium_confidence_threshold = medium_confidence_threshold,
    low_confidence_threshold = low_confidence_threshold
  )
}

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

#' Merge two taxonomy datasets
#' @param tax_df1 First taxonomy data frame
#' @param tax_df2 Second taxonomy data frame
#' @param prefix1 Prefix used in first dataset
#' @param prefix2 Prefix used in second dataset
merge_taxonomies <- function(tax_df1, tax_df2, prefix1, prefix2) {
  tax1_col <- paste0(prefix1, "_taxonomy")
  tax2_col <- paste0(prefix2, "_taxonomy")
  conf1_col <- paste0(prefix1, "_confidence")
  conf2_col <- paste0(prefix2, "_confidence")
  
  merged_df <- full_join(
    tax_df1,
    tax_df2,
    by = "featureid"
  ) %>%
    mutate(
      !!tax1_col := if_else(is.na(!!sym(tax1_col)), "Unassigned", !!sym(tax1_col)),
      !!tax2_col := if_else(is.na(!!sym(tax2_col)), "Unassigned", !!sym(tax2_col)),
      !!conf1_col := if_else(is.na(!!sym(conf1_col)), 0, !!sym(conf1_col)),
      !!conf2_col := if_else(is.na(!!sym(conf2_col)), 0, !!sym(conf2_col))
    )
  
  return(merged_df)
}

#' Count number of taxonomic ranks in a taxonomy string
#' @param tax_string Taxonomy string separated by semicolons
count_tax_ranks <- function(tax_string) {
  sum(str_count(tax_string, ";")) + 1 - sum(str_count(tax_string, "NA"))
}

#' Find the lowest common rank between two taxonomies
#' @param tax1 First taxonomy
#' @param tax2 Second taxonomy
find_common_taxonomy <- function(tax1, tax2) {
  ranks1 <- str_split(tax1, ";")[[1]]
  ranks2 <- str_split(tax2, ";")[[1]]
  
  common_ranks <- c()
  for(i in seq_along(ranks1)) {
    if(i <= length(ranks2) && ranks1[i] == ranks2[i] && ranks1[i] != "NA") {
      common_ranks <- c(common_ranks, ranks1[i])
    } else {
      break
    }
  }
  
  if(length(common_ranks) > 0) {
    return(paste(common_ranks, collapse = ";"))
  } else {
    return(NA_character_)
  }
}

#' Compare taxonomic assignments based on confidence levels and decision rules
#' @param merged_df Merged taxonomy data frame
#' @param config Classification configuration
#' @param prefix1 Prefix for first dataset columns
#' @param prefix2 Prefix for second dataset columns
compare_taxonomic_assignments <- function(merged_df, config, prefix1, prefix2) {
  # Construct column names
  tax1_col <- paste0(prefix1, "_taxonomy")
  tax2_col <- paste0(prefix2, "_taxonomy")
  conf1_col <- paste0(prefix1, "_confidence")
  conf2_col <- paste0(prefix2, "_confidence")
  
  result_df <- merged_df %>%
    mutate(
      final_taxonomy = NA_character_,
      decision = NA_character_,
      tax1_ranks = map_dbl(!!sym(tax1_col), count_tax_ranks),
      tax2_ranks = map_dbl(!!sym(tax2_col), count_tax_ranks)
    )
  
  # 1. Find exact matches
  exact_matches <- result_df %>%
    filter(!!sym(tax1_col) == !!sym(tax2_col)) %>%
    mutate(
      final_taxonomy = !!sym(tax1_col),
      decision = "exact_match"
    )
  
  remaining_df <- result_df %>%
    anti_join(exact_matches, by = "featureid")
  
  # 2. Handle Eukaryota cases with differing rank depths
  eukaryota_cases <- remaining_df %>%
    filter(
      (!!sym(tax1_col) == "Eukaryota" & !!sym(conf1_col) >= config$very_high_confidence_threshold) |
        (!!sym(tax2_col) == "Eukaryota" & !!sym(conf2_col) >= config$very_high_confidence_threshold)
    )
  
  eukaryota_processed <- eukaryota_cases %>%
    rowwise() %>%  # This is important to ensure each row is processed independently
    mutate(
      final_taxonomy = case_when(
        # First dataset is Eukaryota, use second if it has more ranks and high confidence
        !!sym(tax1_col) == "Eukaryota" & tax2_ranks > 1 & 
          !!sym(conf2_col) >= config$high_confidence_threshold ~ !!sym(tax2_col),
        # Second dataset (more diverse) is Eukaryota, but only use first if it has more ranks AND high confidence
        !!sym(tax2_col) == "Eukaryota" & tax1_ranks > 1 & 
          !!sym(conf1_col) >= config$high_confidence_threshold ~ !!sym(tax1_col),
        # If Dataset1 is Eukaryota and Dataset2 has medium confidence, take two ranks after Eukaryota
        !!sym(tax1_col) == "Eukaryota" & tax2_ranks > 2 & 
          !!sym(conf2_col) >= config$medium_confidence_threshold ~ {
            tax2_parts <- strsplit(!!sym(tax2_col), ";")[[1]]
            tax2_parts <- tax2_parts[tax2_parts != "NA"]
            if(length(tax2_parts) >= 3) {
              paste(tax2_parts[1:3], collapse = ";")
            } else {
              paste(tax2_parts, collapse = ";")
            }
          },
        # If Dataset2 is Eukaryota and Dataset1 has medium confidence, take one rank after Eukaryota
        !!sym(tax2_col) == "Eukaryota" & tax1_ranks > 1 & 
          !!sym(conf1_col) >= config$medium_confidence_threshold ~ {
            tax1_parts <- strsplit(!!sym(tax1_col), ";")[[1]]
            tax1_parts <- tax1_parts[tax1_parts != "NA"]
            if(length(tax1_parts) >= 2) {
              paste(tax1_parts[1:2], collapse = ";")
            } else {
              paste(tax1_parts, collapse = ";")
            }
          },
        # For any remaining cases where both are Eukaryota or no other conditions met
        TRUE ~ "Eukaryota"
      ),
      decision = case_when(
        !!sym(tax1_col) == "Eukaryota" & tax2_ranks > 1 & 
          !!sym(conf2_col) >= config$high_confidence_threshold ~ "eukaryota_resolution_high_conf",
        !!sym(tax2_col) == "Eukaryota" & tax1_ranks > 1 & 
          !!sym(conf1_col) >= config$high_confidence_threshold ~ "eukaryota_resolution_high_conf",
        !!sym(tax1_col) == "Eukaryota" & tax2_ranks > 2 & 
          !!sym(conf2_col) >= config$medium_confidence_threshold ~ "eukaryota_resolution_medium_conf_two_ranks",
        !!sym(tax2_col) == "Eukaryota" & tax1_ranks > 1 & 
          !!sym(conf1_col) >= config$medium_confidence_threshold ~ "eukaryota_resolution_medium_conf_one_rank",
        TRUE ~ "eukaryota_resolution_default"
      )
    ) %>%
    ungroup()
  
  remaining_df <- remaining_df %>%
    anti_join(eukaryota_cases, by = "featureid")
  
  # 3. Both very high confidence but different
  very_high_both <- remaining_df %>%
    filter(
      !!sym(conf1_col) >= config$very_high_confidence_threshold &
        !!sym(conf2_col) >= config$very_high_confidence_threshold &
        !!sym(tax1_col) != !!sym(tax2_col)
    ) %>%
    mutate(
      decision = "manual_very_high_both"
    )
  
  remaining_df <- remaining_df %>%
    anti_join(very_high_both, by = "featureid")
  
  # 4. High confidence first dataset cases
  high_first_cases <- remaining_df %>%
    filter(
      !!sym(conf1_col) >= config$high_confidence_threshold &
        (!!sym(conf2_col) <= config$medium_confidence_threshold)
    ) %>%
    mutate(
      final_taxonomy = !!sym(tax1_col),
      decision = paste0("high_confidence_", prefix1)
    )
  
  high_both <- remaining_df %>%
    filter(
      (!!sym(conf1_col) >= config$high_confidence_threshold |
         !!sym(conf1_col) >= config$very_high_confidence_threshold) &
        (!!sym(conf2_col) >= config$medium_confidence_threshold) &
        !!sym(tax1_col) != !!sym(tax2_col)
    ) %>%
    mutate(
      decision = "manual_high_both"
    )
  
  remaining_df <- remaining_df %>%
    anti_join(high_first_cases, by = "featureid") %>%
    anti_join(high_both, by = "featureid")
  
  # 5. Low/medium first but high second
  high_second_cases <- remaining_df %>%
    filter(
      !!sym(conf1_col) < config$high_confidence_threshold &
        !!sym(conf2_col) >= config$high_confidence_threshold
    ) %>%
    mutate(
      final_taxonomy = !!sym(tax2_col),
      decision = paste0("high_confidence_", prefix2)
    )
  
  remaining_df <- remaining_df %>%
    anti_join(high_second_cases, by = "featureid")
  
  # 6. Both low/medium - find common taxonomy
  low_medium_both <- remaining_df %>%
    filter(
      !!sym(conf1_col) < config$high_confidence_threshold &
        !!sym(conf2_col) < config$high_confidence_threshold
    ) %>%
    rowwise() %>%
    mutate(
      final_taxonomy = find_common_taxonomy(!!sym(tax1_col), !!sym(tax2_col)),
      decision = if_else(!is.na(final_taxonomy), "lca", "needs_review")
    )
  
  remaining_df <- remaining_df %>%
    anti_join(low_medium_both, by = "featureid") %>%
    mutate(
      decision = "needs_review"
    )
  
  # Combine all results
  final_results <- bind_rows(
    exact_matches,
    eukaryota_processed,
    very_high_both,
    high_first_cases,
    high_both,
    high_second_cases,
    low_medium_both,
    remaining_df
  ) %>%
    select(featureid, 
           !!sym(tax1_col), !!sym(conf1_col),
           !!sym(tax2_col), !!sym(conf2_col),
           final_taxonomy, decision) %>%
    mutate(
      deciding_confidence = case_when(
        decision == "exact_match" ~ paste0(prefix2, "_confidence"),  # returns "full_confidence"
        decision == "eukaryota_resolution_high_conf" & !!sym(tax1_col) == "Eukaryota" ~ paste0(prefix2, "_confidence"),
        decision == "eukaryota_resolution_high_conf" & !!sym(tax2_col) == "Eukaryota" ~ paste0(prefix1, "_confidence"),
        decision == "eukaryota_resolution_medium_conf_two_ranks" ~ paste0(prefix2, "_confidence"),
        decision == "eukaryota_resolution_medium_conf_one_rank" ~ paste0(prefix1, "_confidence"),
        decision == "eukaryota_resolution_default" ~ if_else(!!sym(tax1_col) == "Eukaryota", 
                                                             paste0(prefix1, "_confidence"), 
                                                             paste0(prefix2, "_confidence")),
        decision == paste0("high_confidence_", prefix1) ~ paste0(prefix1, "_confidence"),
        decision == paste0("high_confidence_", prefix2) ~ paste0(prefix2, "_confidence"),
        decision == "manual_very_high_both" ~ if_else(!!sym(conf1_col) > !!sym(conf2_col),
                                                      paste0(prefix1, "_confidence"),
                                                      paste0(prefix2, "_confidence")),
        decision == "manual_high_both" ~ if_else(!!sym(conf1_col) > !!sym(conf2_col),
                                                 paste0(prefix1, "_confidence"),
                                                 paste0(prefix2, "_confidence")),
        decision == "lca" ~ "consensus",
        TRUE ~ NA_character_
      ),
      final_confidence = case_when(
        deciding_confidence == paste0(prefix1, "_confidence") ~ !!sym(conf1_col),
        deciding_confidence == paste0(prefix2, "_confidence") ~ !!sym(conf2_col),
        deciding_confidence == "consensus" ~ NA_real_,
        TRUE ~ NA_real_
      )
    ) %>%
    arrange(featureid)
  
  return(final_results)
}
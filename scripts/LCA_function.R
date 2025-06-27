library(tidyverse)

handle_needs_attention <- function(mismatches) {
  # Helper function to find the lowest common taxonomic rank
  find_lowest_common_rank <- function(local_rank, universal_rank) {
    ranks <- c("k", "p", "c", "o", "f", "g", "s")
    matching_indices <- which(local_rank == universal_rank)
    if (length(matching_indices) == 0) {
      return(NA)
    } else {
      return(ranks[max(matching_indices)])
    }
  }
  
  # Helper function to create final taxonomic columns
  create_final_taxonomic_columns <- function(data) {
    ranks <- c("k", "p", "c", "o", "f", "g", "s")
    final_columns <- setNames(rep(NA, length(ranks)), paste0("final_", ranks))
    hc_index <- match(data$lowest_common_rank, ranks)
    for (i in 1:length(ranks)) {
      ifelse(i <= hc_index, 
             final_columns[paste0("final_", ranks[i])] <- data[[paste0("L_", ranks[i])]],
             final_columns[paste0("final_", ranks[i])] <- NA)
    }
    return(as.list(final_columns))
  }
  
  # Main processing
  mismatch2 <- mismatches %>%
    separate(local, c("L_k","L_p","L_c","L_o","L_f","L_g","L_s"), sep = ";") %>%
    separate(global, c("U_k","U_p","U_c","U_o","U_f","U_g","U_s"), sep = ";")
  
  misLowRank <- mismatch2 %>%
    rowwise() %>%
    mutate(lowest_common_rank = find_lowest_common_rank(c_across(starts_with("L_")), c_across(starts_with("U_")))) %>%
    filter(!is.na(lowest_common_rank)) %>%
    mutate(final_columns = list(create_final_taxonomic_columns(across(c(lowest_common_rank, starts_with("L_")))))) %>%
    unnest_wider(final_columns)
  
  FinalLowCom <- misLowRank %>% 
    unite("LocalTax", L_k:L_s, sep = ";", na.rm = TRUE) %>% 
    unite("UniTax", U_k:U_s, sep = ";", na.rm = TRUE) %>%
    unite("FinalTax", final_k:final_s, sep = ";", na.rm = TRUE) %>%
    select(-lowest_common_rank)
  
  return(FinalLowCom)
}
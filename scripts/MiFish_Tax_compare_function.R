library(tidyverse)

compare_taxonomy_tables <- function(localTAX, globalTAX) {
  # Preprocess local taxonomy table
  localTAX <- localTAX %>%
    dplyr::rename("local" = "Taxon",
                  "ConfidenceLoc" = "Confidence")
  
  # Preprocess global taxonomy table
  globalTAX <- globalTAX %>%
    dplyr::rename("global" = "Taxon",
                  "ConfidencUni" = "Confidence")
  
  # Combine tables
  tax_compare <- left_join(localTAX, globalTAX, by = "Feature ID")
  
  # Find matches and mismatches
  matches <- tax_compare %>%
    filter(local == global) 
  
  mismatches <- anti_join(tax_compare, matches, by = "Feature ID")
  
  # Decision making process
  decision1 <- mismatches %>%
    filter(local %in% c("Eukaryota;Chordata;Actinopteri", "Eukaryota", "Eukaryota;Chordata")) %>%
    mutate(final = global)
  
  mismatches2 <- anti_join(mismatches, decision1, by = "Feature ID")
  
  decision2 <- mismatches2 %>%
    filter(ConfidenceLoc >= 98) %>% 
    mutate(final = local)
  
  mismatches3 <- anti_join(mismatches2, decision2, by = "Feature ID")
  
  decision3 <- bind_rows(decision1, decision2)
  
  # Combine all results
  final_Mifish_tax <- bind_rows(
    matches %>% mutate(final = local),
    decision3
  ) 
  
  # Handle needs attention cases
  needs_attention_resolved <- handle_needs_attention(mismatches3)
  
  # Combine final results
  final_results <- bind_rows(
    final_Mifish_tax,
    needs_attention_resolved %>% mutate(final = FinalTax) %>% select(-FinalTax)
  )
  
  return(final_results)
}
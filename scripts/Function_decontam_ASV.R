library(tidyverse)

posctrl_decontam_data <- function(data, positive_control_samples) {
  # Ensure the data is in a tidy format
  data_long <- data %>%
    tidyr::pivot_longer(cols = -featureid, names_to = "Sample", values_to = "nReads")
  
  # Identify ASVs present in any of the positive control samples
  positive_control_asvs <- data_long %>%
    dplyr::filter(Sample %in% positive_control_samples & nReads > 0) %>%
    dplyr::pull(featureid) %>%
    unique()  # Ensure unique ASVs
  
  # Filter the data to retain only ASVs with reads in any positive control
  filtered_data <- data %>%
    dplyr::filter(featureid %in% positive_control_asvs)
  
  # Calculate metrics for non-positive control samples
  other_samples_summary <- data_long %>%
    dplyr::filter(featureid %in% positive_control_asvs & !Sample %in% positive_control_samples) %>%
    dplyr::group_by(featureid) %>%
    dplyr::summarise(
      Avg_Reads_Samples = mean(nReads, na.rm = TRUE),
      Sum_Reads_Samples = sum(nReads, na.rm = TRUE),
      Max_Reads_Samples = max(nReads, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Subset the data to retain only positive control columns
  positive_control_data <- filtered_data %>%
    dplyr::select(featureid, dplyr::all_of(positive_control_samples))
  
  # Merge the metrics with the positive control data
  result <- positive_control_data %>%
    dplyr::left_join(other_samples_summary, by = "featureid")
  
  return(result)
}

# Example usage:
# Assuming your data is in a data frame called `asv_data` with columns: featureid, Sample1, Sample2, ..., PositiveControl1, PositiveControl2, ...
# asv_data <- data.frame(
#   featureid = c("ASV1", "ASV2", "ASV3"),
#   Sample1 = c(10, 0, 5),
#   Sample2 = c(0, 20, 15),
#   PositiveControl1 = c(5, 0, 10),
#   PositiveControl2 = c(3, 0, 8),
#   ...
# )

# processed_data <- posctrl_decontam_data(asv_data, c("PositiveControl1", "PositiveControl2"))
# print(processed_data)

library(tidyverse)

resolution_18S_tax <- function(tax) {
  tax %>%
    mutate(final_taxonomy = case_when(
      str_starts(final_taxonomy, "Eukaryota;Chordata;Mammalia") ~ "Eukaryota;Chordata;Mammalia",
      str_starts(final_taxonomy, "Eukaryota;Chordata;Actinopteri") ~ "Eukaryota;Chordata;Actinopteri",
      str_starts(final_taxonomy, "Eukaryota;Mollusca;Cephalopoda") ~ "Eukaryota;Mollusca;Cephalopoda",
      str_starts(final_taxonomy, "Eukaryota;Arthropoda;Malacostraca;Decapoda;Crangonidae;Crangon") ~ "Eukaryota;Arthropoda;Malacostraca;Decapoda;Crangonidae;Crangon",
      str_starts(final_taxonomy, "Eukaryota;Arthropoda;Malacostraca;Decapoda;Pandalidae;Pandalus") ~ "Eukaryota;Arthropoda;Malacostraca;Decapoda;Pandalidae;Pandalus",
      str_starts(final_taxonomy, "Eukaryota;Annelida;Polychaeta;Phyllodocida;Nereididae") ~ "Eukaryota;Annelida;Polychaeta;Phyllodocida;Nereididae",
      str_starts(final_taxonomy, "Eukaryota;Arthropoda;Malacostraca;Decapoda;Callianassidae;Neotrypaea") ~ "Eukaryota;Arthropoda;Malacostraca;Decapoda;Callianassidae;Neotrypaea",
      TRUE ~ final_taxonomy  # Keep the original taxonomy if it doesn't match the above conditions
    ))
}


resolution_MiFish_tax <- function(tax) {
  tax %>%
    mutate(final_taxonomy = case_when(
      str_starts(final_taxonomy, "Eukaryota;Chordata;Mammalia;Carnivora;Phocidae") ~ "Eukaryota;Chordata;Mammalia;Carnivora;Phocidae;Phoca",
      str_starts(final_taxonomy, "Eukaryota;Chordata;Mammalia;Carnivora;Otariidae;Zalophus") ~ "Eukaryota;Chordata;Mammalia;Carnivora;Otariidae;Zalophus",
      str_starts(final_taxonomy, "Eukaryota;Chordata;Actinopteri;Perciformes;Sebastidae;Sebastes") ~ "Eukaryota;Chordata;Actinopteri;Perciformes;Sebastidae;Sebastes",
      str_starts(final_taxonomy, "EEukaryota;Chordata;Hyperoartia;Petromyzontiformes;Petromyzontidae;Lampetra") ~ "Eukaryota;Chordata;Hyperoartia;Petromyzontiformes;Petromyzontidae;Lampetra",
      TRUE ~ final_taxonomy  # Keep the original taxonomy if it doesn't match the above conditions
    ))
}
# Example usage:
# Assuming your taxonomy data is in a data frame called `tax` with columns: featureid and final_taxonomy
# tax <- data.frame(
#   featureid = c("ASV1", "ASV2", "ASV3", "ASV4"),
#   final_taxonomy = c(
#     "Eukaryota;Chordata;Mammalia;Primates;Hominidae;Homo",
#     "Eukaryota;Chordata;Actinopteri;Clupeiformes;Clupeidae;Clupea",
#     "Eukaryota;Arthropoda;Insecta;Diptera;Drosophilidae;Drosophila",
#     "Eukaryota;Chordata;Mammalia;Carnivora;Felidae;Felis"
#   )
# )

# modified_tax <- resolution_18S_tax(tax)
# print(modified_tax)
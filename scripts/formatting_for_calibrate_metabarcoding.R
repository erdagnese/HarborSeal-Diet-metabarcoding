### WEH Aug 29, 2025 - Erin D
## This sweet of functions help format the data for model input

################DEFINE SUB-FUNCTIONS
library(tidyverse)
library(stringr)

# Function to extract final taxon
extract_final_taxon <- function(full_taxonomy) {
  # Split by semicolon and remove empty elements
  taxa_levels <- str_split(full_taxonomy, ";")[[1]] %>% 
    str_trim() %>% 
    .[. != ""]
  
  # Work backwards through the taxonomy levels
  for (i in length(taxa_levels):1) {
    current_taxon <- taxa_levels[i]
    # Skip if the taxon is NA or empty
    if (!is.na(current_taxon) && current_taxon != "NA" && current_taxon != "") {
      return(current_taxon)
    }
  }
  
  # If all levels are NA/empty, return NA
  return(NA_character_)
}
# Function to parse sample names based on region
library(tidyverse)
library(stringr)

# Enhanced parser that handles complex station names
parse_sample_name <- function(samplename, region) {
  parts <- str_split(samplename, "_")[[1]]
  
  if (region == "ControlFeeding") {
    # Format: PvCF_D1_01_TR1 (exactly 4 parts)
    creek <- parts[1]
    station <- parts[2]
    biol <- as.integer(parts[3])
    tech <- as.integer(str_remove(parts[4], "TR"))
  } else if (region == "TissueMix") {
    creek <- "TM"
    biol <- as.integer(parts[2])
    
    # Check if last part is a tech rep (digits or R+digits)
    last_part <- tail(parts, 1)
    if (str_detect(last_part, "^(R?\\d+)$")) {
      tech <- as.integer(str_remove(last_part, "R"))
      # Station is everything between biol and tech
      station <- paste(parts[3:(length(parts)-1)], collapse = "_")
    } else {
      tech <- 1
      # Station is everything after biol
      station <- paste(parts[3:length(parts)], collapse = "_")
    }
  } else {
    return(list(creek = NA, station = NA, biol = NA, tech = NA))
  }
  
  return(list(creek = creek, station = station, biol = biol, tech = tech))
}


# Final transformation function
modified_sample_data <- function(asv_data) {
  asv_data %>%
    filter(region %in% c("ControlFeeding", "TissueMix")) %>%
    rowwise() %>%
    mutate(
      parsed = list(parse_sample_name(samplename, region)),
      creek = parsed$creek,
      station = parsed$station,
      biol = parsed$biol,
      tech = parsed$tech,
      time = 1,
      species = final_taxonomy,
      Nreads = count
    ) %>%
    ungroup() %>%
    # Keep samplename in the output for verification
    select(samplename, time, creek, station, biol, tech, species, Nreads) %>%
    filter(!is.na(creek))  # Remove any unparseable samples
}


# Run the transformation
#sample_data <- modified_sample_data(samples18S_asv_data)

library(tidyverse)

analyze_mock_species_in_samples <- function(sample_df, combined_mock_data) {
  # Get unique species from combined mock data
  mock_species <- unique(combined_mock_data$species)
  
  # First, get the creek information for each sample (assuming it's consistent per samplename)
  sample_creek_info <- sample_df %>%
    distinct(samplename, creek)
  
  # Sum reads by sample and species for mock species only
  mock_reads_by_sample <- sample_df %>%
    filter(species %in% mock_species) %>%
    group_by(samplename, species) %>%
    summarise(species_reads = sum(Nreads), .groups = 'drop')
  
  # Pivot wider to get one column per mock species
  mock_species_wide <- mock_reads_by_sample %>%
    pivot_wider(
      names_from = species, 
      values_from = species_reads, 
      values_fill = 0
    )
  
  # Calculate total metrics for each sample
  sample_totals <- sample_df %>%
    group_by(samplename) %>%
    summarise(
      total_reads = sum(Nreads),
      total_mock_sp_reads = sum(Nreads[species %in% mock_species]),
      prop_mock_reads = total_mock_sp_reads / total_reads,
      .groups = 'drop'
    )
  
  # Combine the wide format mock species with the total metrics and creek info
  sample_metrics <- sample_totals %>%
    left_join(mock_species_wide, by = "samplename") %>%
    left_join(sample_creek_info, by = "samplename") %>%
    # Replace NA with 0 for any mock species columns that might be missing
    mutate(across(any_of(mock_species), ~ replace_na(., 0)))
  
  # Calculate n_mock_species based on columns with reads > 0
  sample_metrics <- sample_metrics %>%
    rowwise() %>%
    mutate(
      n_mock_species = sum(c_across(any_of(mock_species)) > 0)
    ) %>%
    ungroup() %>%
    # Reorder columns to put creek in a useful position
    select(samplename, creek, total_reads, n_mock_species, total_mock_sp_reads, prop_mock_reads, everything()) %>%
    arrange(desc(prop_mock_reads))
  
  return(sample_metrics)
}

# You can also create a summary to help with filtering decisions
summary_stats <- function(metrics_df) {
  metrics_df %>%
    summarise(
      n_samples = n(),
      avg_prop_mock = mean(prop_mock_reads),
      median_prop_mock = median(prop_mock_reads),
      min_prop_mock = min(prop_mock_reads),
      max_prop_mock = max(prop_mock_reads),
      avg_mock_species = mean(n_mock_species),
      median_mock_species = median(n_mock_species)
    )
}

filter_samples_by_mock_content <- function(metrics_df, creek_filter = NULL,
                                           min_prop_mock = 0.5, min_mock_species = 2,
                                           min_total_reads = 1000) {
  
  filtered_metrics <- metrics_df %>%
    {if(!is.null(creek_filter)) filter(., creek == creek_filter) else .} %>%
    filter(prop_mock_reads >= min_prop_mock,        # Corrected column name
           n_mock_species >= min_mock_species,      # Corrected column name
           total_reads >= min_total_reads)          # Correct column name
  
  return(filtered_metrics)
}

# Your main function stays the same
filter_sample_data_by_metrics <- function(sample_df, metrics_df, creek_filter = NULL,
                                          min_prop_mock = 0.5, min_mock_species = 2,
                                          min_total_reads = 1000) {
  
  # Apply the same filtering criteria to get samples to keep
  filtered_metrics <- filter_samples_by_mock_content(
    metrics_df, 
    creek_filter = creek_filter,
    min_prop_mock = min_prop_mock,
    min_mock_species = min_mock_species,
    min_total_reads = min_total_reads
  )
  
  # Get the samplenames to keep
  samples_to_keep <- filtered_metrics$samplename
  
  # Filter the sample data
  filtered_data <- sample_df %>%
    filter(samplename %in% samples_to_keep)
  
  return(filtered_data)
}


# Usage example:
#sample_metrics <- analyze_mock_species_in_samples(sample_data, combined_mock2)
#wild_sample_metrics <- analyze_mock_species_in_samples(wild_sample_data, combined_mock2)
library(ggplot2)
library(dplyr)
library(tidyr)

plot_taxonomy_proportions <- function(data, 
                                      taxonomy_col = "ID_18S",
                                      count_col = "count",
                                      group_cols = c("community", "tech"),
                                      top_n_taxa = 12,
                                      plot_title = "Taxonomic Proportions by Group") {
  
  # Calculate proportions
  prop_data <- data %>%
    group_by(across(all_of(c(group_cols, taxonomy_col)))) %>%
    summarise(total_count = sum(.data[[count_col]]), .groups = "drop_last") %>%
    mutate(proportion = total_count / sum(total_count)) %>%
    ungroup()
  
  # Get top taxa for each group
  top_taxa <- prop_data %>%
    group_by(.data[[taxonomy_col]]) %>%
    summarise(avg_prop = mean(proportion)) %>%
    arrange(desc(avg_prop)) %>%
    slice_head(n = top_n_taxa) %>%
    pull(.data[[taxonomy_col]])
  
  # Filter and order data
  plot_data <- prop_data %>%
    mutate(!!sym(taxonomy_col) := ifelse(.data[[taxonomy_col]] %in% top_taxa, 
                                         .data[[taxonomy_col]], "Other"),
           !!sym(taxonomy_col) := factor(.data[[taxonomy_col]], 
                                         levels = c(top_taxa, "Other"))) %>%
    group_by(across(all_of(c(group_cols, taxonomy_col)))) %>%
    summarise(proportion = sum(proportion), .groups = "drop")
  
  # Create combined grouping variable for plotting
  plot_data <- plot_data %>%
    unite("group", all_of(group_cols), sep = " | ")
  
  # Generate plot
  p <- ggplot(plot_data, aes(x = group, y = proportion, fill = .data[[taxonomy_col]])) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = plot_title,
         x = paste(group_cols, collapse = " | "),
         y = "Proportion of Counts",
         fill = taxonomy_col) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(size = 7)) +
    scale_y_continuous(labels = scales::percent)
  
  if (length(top_taxa) > 5) {
    p <- p + guides(fill = guide_legend(ncol = 2))
  }
  
  return(p)
}

plot_proportion_comparison <- function(data, 
                                       taxonomy_col = "species",
                                       count_col = "nReads",
                                       existing_prop_col = "propReads",
                                       group_cols = c("Community", "tech"),
                                       top_n_taxa = 12,
                                       plot_title = "Proportion Comparison") {
  
  # Calculate proportions from counts
  calc_data <- data %>%
    group_by(across(all_of(c(group_cols, taxonomy_col)))) %>%
    summarise(total_count = sum(.data[[count_col]]), .groups = "drop_last") %>%
    mutate(calculated_prop = total_count / sum(total_count)) %>%
    ungroup()
  
  # Get existing proportions if column exists
  if (existing_prop_col %in% names(data)) {
    existing_data <- data %>%
      group_by(across(all_of(c(group_cols, taxonomy_col)))) %>%
      summarise(existing_prop = sum(.data[[existing_prop_col]]), .groups = "drop")
    
    # Combine both datasets
    plot_data <- full_join(calc_data, existing_data, 
                           by = c(group_cols, taxonomy_col))
  } else {
    message("No existing proportion column found - plotting only calculated proportions")
    plot_data <- calc_data %>%
      mutate(existing_prop = NA_real_)
  }
  
  # Get top taxa based on calculated proportions
  top_taxa <- calc_data %>%
    group_by(.data[[taxonomy_col]]) %>%
    summarise(avg_prop = mean(calculated_prop)) %>%
    arrange(desc(avg_prop)) %>%
    slice_head(n = top_n_taxa) %>%
    pull(.data[[taxonomy_col]])
  
  # Prepare data for plotting
  plot_data <- plot_data %>%
    mutate(!!sym(taxonomy_col) := ifelse(.data[[taxonomy_col]] %in% top_taxa, 
                                         .data[[taxonomy_col]], "Other"),
           !!sym(taxonomy_col) := factor(.data[[taxonomy_col]], 
                                         levels = c(top_taxa, "Other"))) %>%
    group_by(across(all_of(c(group_cols, taxonomy_col)))) %>%
    summarise(calculated_prop = sum(calculated_prop, na.rm = TRUE),
              existing_prop = sum(existing_prop, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_longer(cols = c(calculated_prop, existing_prop),
                 names_to = "proportion_type",
                 values_to = "proportion") %>%
    unite("group", all_of(group_cols), sep = " | ") %>%
    mutate(proportion_type = factor(proportion_type,
                                    levels = c("calculated_prop", "existing_prop"),
                                    labels = c("From Counts", "Existing Proportions")))
  
  # Generate plot
  p <- ggplot(plot_data, aes(x = group, y = proportion, fill = .data[[taxonomy_col]])) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ proportion_type, ncol = 2) +
    labs(title = plot_title,
         x = paste(group_cols, collapse = " | "),
         y = "Proportion",
         fill = taxonomy_col) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(size = 7),
          strip.text = element_text(face = "bold")) +
    scale_y_continuous(labels = scales::percent)
  
  if (length(top_taxa) > 5) {
    p <- p + guides(fill = guide_legend(ncol = 2))
  }
  
  return(p)
}


# Modified formatting function with external station mapping
format_metabarcoding_data <- function(input_metabarcoding_RDS, input_mock_comm_RDS, 
                                      station_mapping, reference_species = "Alosa sapidissima") {
  require(tidyverse)
  
  Observation <- input_metabarcoding_RDS
  Mock <- input_mock_comm_RDS
  
  # Process Mock data FIRST - join station_idx before using it
  Mock <- Mock %>% 
    group_by(Community) %>% 
    mutate(station_idx = cur_group_id()) %>% 
    ungroup() %>%
    rename(Nreads = nReads) %>%
    mutate(
      time = 1,
      creek = 1,
      biol = 1,
      station = station_idx,
      stationname = Community
    ) %>% 
    replace(is.na(.), 0) %>% 
    filter(start_prop > 0 & Nreads > 0 & !is.na(Nreads)) %>% 
    unite(c(time, creek, station, biol, tech), col = "Sample", sep = ".", remove = FALSE)
  
  # Keep species present in both datasets
  keepSpecies <- intersect(Mock$species, Observation$species)
  Observation <- Observation %>% filter(species %in% keepSpecies)
  Mock <- Mock %>% filter(species %in% keepSpecies)
  
  # Create species index with reference species as the last index
  all_species <- unique(c(Mock$species, Observation$species))
  
  # Remove reference species from the main list
  main_species <- setdiff(all_species, reference_species)
  
  # Create ordered species list: main species first, reference species last
  ordered_species <- c(main_species, reference_species)
  
  # Create species index mapping
  sp_list <- data.frame(
    species = ordered_species,
    species_idx = 1:length(ordered_species),
    stringsAsFactors = FALSE
  )
  
  # Process Mock data with the new species index
  Mock <- Mock %>% 
    left_join(sp_list, by = "species") %>% 
    mutate(
      speciesname = species,
      species = species_idx
    ) %>% 
    group_by(station, tech, time, creek, biol) %>% 
    mutate(b_proportion = start_prop/sum(start_prop)) %>% 
    ungroup()
  
  # Process Observation data with the provided station mapping
  Observation <- Observation %>% 
    left_join(sp_list, by = "species") %>% 
    mutate(
      speciesname = species,
      species = species_idx
    ) %>% 
    left_join(station_mapping, by = "stationname") %>% 
    group_by(station, time, creek, biol) %>% 
    mutate(
      tech = match(tech, unique(tech)),
      tot_reads = sum(Nreads)
    ) %>% 
    filter(tot_reads > 0) %>% 
    ungroup() %>% 
    dplyr::select(-tot_reads)
  
  # Return final object
  list(
    Observation = Observation,
    Mock = Mock,
    N_pcr_mock = unique(input_mock_comm_RDS$N_pcr_mock),
    NSpecies = nrow(sp_list),
    station_list = station_mapping,
    sp_list = sp_list
  )
}

# example
# library(here)
# qmdata <- format_metabarcoding_data(input_metabarcoding_RDS,input_mock_comm_RDS, station_mapping = station_mapping, reference_species = "Actinopteri")

# Metabarcoding calibration model function
# WEH - ERD 2025

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

library(rstan)
library(MCMCpack) #for rdirichelet function
library(gridExtra)
library(unikn)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


N_pcr_cycles <- 39
MOCK <- qmdata$Mock  

##### FIXED VERSION WITH spread() ####
alrTransform <- function(MOCK){
  require(tidyverse)
  require(compositions)
  
  p_mock <- MOCK %>% 
    dplyr::select(species, station, tech, b_proportion) %>% 
    rename(tech_rep = tech) %>% 
    spread(key = species, value = b_proportion, fill = 1e-9) %>%  # Use spread() instead
    ungroup() 
  
  colnames(p_mock)[3:(length(unique(MOCK$species))+2)] <- paste0("alr_", 1:length(unique(MOCK$species)))
  
  p_mock <- alr(p_mock[,3:ncol(p_mock)]) %>% as.matrix() %>% as.data.frame()
  p_mock[,length(unique(MOCK$species))] <- 0  #add reference zero expressly
  names(p_mock)[length(unique(MOCK$species))] <- paste0("alr_", length(unique(MOCK$species)))
  
  p_mock <-  cbind(MOCK %>% dplyr::select(tech, station) %>% distinct(),
                   p_mock) %>% ungroup()
  names(p_mock)[1] <- "tech_rep"
  
  return(p_mock)
}
  
  
makeDesign <- function(obs, #obs is a named list with elements Observation, Mock, N_pcr_mock, sp_list
                       N_pcr_cycles){ #N_pcr_cycles is the number of PCR cycles in your experimental/enviro samples; currently a single value, could be made into a vector if this number varies
  #library(tidyverse)
  library(MCMCpack)
  library(compositions)
  library(rstan)
  library(dplyr)
  
  
  mock <- obs$Mock
  observed <- obs$Observation
  
  p_mock_all <- alrTransform(mock)
  
  mock_3 <- mock %>% 
    dplyr::select(species, station, tech, Nreads) %>% 
    ungroup() %>% 
    mutate(species = paste0("sp_", species)) %>% 
    pivot_wider(names_from = species, values_from = Nreads, 
                values_fill = 0)  # FIXED: named list
  
  N_pcr_mock <- rep(obs$N_pcr_mock, nrow(p_mock_all)) #assumes all have the same Npcr
  
  
  p_samp_all <- observed %>% 
    ungroup() %>% 
    unite(time, creek, station, biol, col = "station") %>% 
    dplyr::select(station, tech, species, Nreads) %>% 
    rename(tech_rep = tech) %>% 
    mutate(species = paste0("sp_", species)) %>% 
    arrange(species) %>% 
    group_by(station, tech_rep,species) %>% 
    dplyr::summarise(Nreads= sum(Nreads)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = species, values_from = Nreads, 
                values_fill = 0)
  N_pcr_samp <- rep(N_pcr_cycles, nrow(p_samp_all))
  
  ########################################################################
  #### Create data frames that can be read into Stan model
  ########################################################################
  
  NOM <- as.name(colnames(p_mock_all)[1])
  formula_a <- eval(NOM) ~ N_pcr_mock -1
  model_frame <- model.frame(formula_a, p_mock_all)
  model_vector_a_mock <- model.matrix(formula_a, model_frame) %>% as.numeric()
  N_pcr_mock_small <- cbind(N_pcr_mock, p_mock_all) %>%  filter(tech_rep == 1) %>% pull(N_pcr_mock)
  formula_b <- eval(NOM) ~ N_pcr_mock_small -1
  model_frame <- model.frame(formula_b, p_mock_all%>% filter(tech_rep==1))
  model_vector_a_mock_small <- model.matrix(formula_b, model_frame) %>% as.numeric()
  
  N_obs_mock       <- nrow(p_mock_all)
  
  # unknown communities second
  # species compositions (betas)
  
  NOM <- as.name(colnames(p_samp_all)[1])    
  
  p_samp_all$station <- as.factor(p_samp_all$station)
  N_station = length(unique(p_samp_all$station))
  p_samp_all$tech_rep <- as.factor(p_samp_all$tech_rep)
  if(N_station == 1){
    formula_b <- eval(NOM) ~ 1  
  } else {
    formula_b <- eval(NOM) ~ station
  }
  
  model_frame <- model.frame(formula_b, p_samp_all)
  model_matrix_b_samp <- model.matrix(formula_b, model_frame)
  
  # choose a single representative for each station to make predictions to
  model_frame <- model.frame(formula_b, p_samp_all %>% filter(tech_rep==1))
  model_matrix_b_samp_small <- model.matrix(formula_b, model_frame)
  
  # efficiencies (alpha)
  formula_a <- eval(NOM) ~ N_pcr_samp -1
  model_frame <- model.frame(formula_a, p_samp_all)
  model_vector_a_samp <- model.matrix(formula_a, model_frame) %>% as.numeric()
  N_pcr_samp_small <- cbind(N_pcr_samp, p_samp_all) %>% filter(tech_rep == 1) %>% pull(N_pcr_samp)
  formula_b <- eval(NOM) ~ N_pcr_samp_small -1
  
  model_frame <- model.frame(formula_b, p_samp_all %>% filter(tech_rep==1))
  model_vector_a_samp_small <- model.matrix(formula_b, model_frame) %>% as.numeric()
  
  #counters 
  N_obs_samp_small <- nrow(model_matrix_b_samp_small)
  N_obs_samp <- nrow(p_samp_all)
  N_b_samp_col <- ncol(model_matrix_b_samp)
  
  
  #### Make Stan objects
  
  stan_data <- list(
    N_species = ncol(p_samp_all)-2,   # Number of species in data
    N_obs_samp = nrow(p_samp_all), # Number of observed community samples and tech replicates ; this will be Ncreek * Nt * Nbiol * Ntech * 2 [for upstream/downstream observations]
    N_obs_mock = nrow(p_mock_all), # Number of observed mock samples, including tech replicates
    N_obs_samp_small = nrow(p_samp_all %>% filter(tech_rep == 1)), # Number of unique observed community samples ; this will be Ncreek * Nt * Nbiol * 2 [for upstream/downstream observations]
    
    # Observed data of community matrices
    sample_data = p_samp_all %>% dplyr::select(contains("sp")),
    sample_vector = unique(p_samp_all$station),
    mock_data   = mock_3 %>% dplyr::select(contains("sp")),
    sp_list = obs$sp_list,
    
    # True proportions for mock community
    #mock_true_prop = p_mock_all %>% dplyr::select(contains("sp")),
    alr_mock_true_prop = p_mock_all %>% dplyr::select(contains("alr")),
    
    # vectors of PCR numbers
    N_pcr_samp = N_pcr_samp,
    N_pcr_mock = N_pcr_mock,
    
    # Design matrices: field samples
    N_b_samp_col = N_b_samp_col,
    model_matrix_b_samp = model_matrix_b_samp,
    model_matrix_b_samp_small = as.array(model_matrix_b_samp_small),
    model_vector_a_samp = model_vector_a_samp,
    model_vector_a_samp_small = as.array(model_vector_a_samp_small),
    
    # Design matrices: mock community samples
    model_vector_a_mock = as.array(model_vector_a_mock),
    
    # Priors
    alpha_prior = c(0,0.5),  # normal prior
    beta_prior = c(0,5),    # normal prior
    tau_prior = c(1,2)   # gamma prior
  )
  
  return(stan_data)
  
}

  #example
  #stan_metabarcoding_data <- makeDesign(metabarcoding_data, N_pcr_cycles = 43)    
  
  
###########################################
###########################################
QM_likelihood <- function(stanmodelname, stan_metabarcoding_data){
    M <- stan_model(stanmodelname)
    
    stanOpt <- optimizing(M, data=stan_metabarcoding_data, iter=30000,draws=0,
                          verbose=T,  
                          tol_param=1e-40,
                          algorithm="LBFGS",
                          hessian = TRUE)
    
    MLest <- stanOpt$par[grep("int_samp_small", names(stanOpt$par))] %>%
      matrix(ncol = stan_metabarcoding_data$N_species) %>% 
      as.data.frame()
    names(MLest) <- stan_metabarcoding_data$sp_list$species
    rownames(MLest) <- stan_metabarcoding_data$sample_vector
    ML_a <- stanOpt$par[grep("alpha\\[", names(stanOpt$par))]  
    ML_a <- data.frame("alpha_est" = ML_a, 
                       "species" = stan_metabarcoding_data$sp_list$species)
    
    
    return(list(
      ML_modelfit = stanOpt,
      ML_estimates = MLest,
      ML_alpha_est = ML_a
    ))
    
  }
  
  #example; note this can fail stochastically; run it several times if need be
  # ML_out <- QM_likelihood(here("quant_metabar_rosetta_noSampleEta.stan"), stan_metabarcoding_data)
###########################################
###########################################
  
###########################################
###########################################
  
  
  
  
  
  
  QM_bayes <- function(stanmodelname, stan_metabarcoding_data, NCHAINS = 3, WARMUP = 500, ITER = 1500){
    require(tidyverse)
    require(rstan)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    
    
    stan_pars <- c( 
      "alpha",
      "beta",
      "eta_mock",
      "tau",
      "mu_samp",
      "mu_mock",
      "int_samp_small"
    )
    
    stanMod = stan(file = stanmodelname ,data = stan_metabarcoding_data,
                   verbose = FALSE, chains = NCHAINS, thin = 1,
                   warmup = WARMUP, iter = ITER,
                   control = list(adapt_init_buffer = 175,
                                  max_treedepth=12,
                                  stepsize=0.01,
                                  adapt_delta=0.7,
                                  metric="diag_e"),
                   pars = stan_pars,
                   refresh = 10,
                   boost_lib = NULL 
    )  
    
    
    mean_est <- summary(stanMod, par = "int_samp_small")$summary[,1] %>%
      matrix(ncol = stan_metabarcoding_data$N_species, byrow = TRUE) %>% 
      as.data.frame()
    names(mean_est) <- stan_metabarcoding_data$sp_list$species
    rownames(mean_est) <- stan_metabarcoding_data$sample_vector
    ci25_est <- summary(stanMod, par = "int_samp_small")$summary[,5] %>%
      matrix(ncol = stan_metabarcoding_data$N_species, byrow = TRUE) %>% 
      as.data.frame()
    names(ci25_est) <- stan_metabarcoding_data$sp_list$species
    rownames(ci25_est) <- stan_metabarcoding_data$sample_vector
    
    ci75_est <- summary(stanMod, par = "int_samp_small")$summary[,7] %>%
      matrix(ncol = stan_metabarcoding_data$N_species, byrow = TRUE) %>% 
      as.data.frame()
    names(ci75_est) <- stan_metabarcoding_data$sp_list$species
    rownames(ci75_est) <- stan_metabarcoding_data$sample_vector
    
    mean_a_est <- summary(stanMod, par = "alpha")$summary[,1]
    mean_a_est <- data.frame("alpha_est" = mean_a_est, 
                             "species" = stan_metabarcoding_data$sp_list$species)
    
    #unique(meta.samples$sample)
    
    return(list(
      Bayes_modelfit = stanMod,
      Bayes_estimates = mean_est,
      Bayes_25ci = ci25_est,
      Bayes_75ci = ci75_est,
      Bayes_alpha_est = mean_a_est
    ))
    
  }
  
##example
  #QM_bayes_out <- QM_bayes(here("quant_metabar_rosetta_noSampleEta.stan"), stan_metabarcoding_data)
  #QM_bayes_out$Bayes_estimates
  ###########################################
  ###########################################
  
  
  
  ########################      
  
  #all-in-one function that calls those above:
  #inputs == input_csv, stanmodelname
  #output == fitted stan model object
  
  run_QM_model <- function(input_metabarcoding_RDS, 
                                      input_mock_comm_RDS, 
                                      N_pcr_cycles = 39,
                                      station_mapping,
                                      reference_species = "Actinopteri",
                                      stanmodel, 
                                      method = "ML"  # "ML" or "Bayes"
                                      ){
    require(dplyr)
    require(here)
    
    
    metabarcoding_data <- format_metabarcoding_data(input_metabarcoding_RDS = input_metabarcoding_RDS,
                                                    input_mock_comm_RDS = input_mock_comm_RDS, 
                                                    station_mapping = station_mapping, 
                                                    reference_species = reference_species)
    
    stan_metabarcoding_data <- makeDesign(metabarcoding_data, N_pcr_cycles = N_pcr_cycles) 
    
    if (method == "ML"){
      ML_out <- QM_likelihood(stanmodel, stan_metabarcoding_data)
    } 
    
    if (method == "Bayes"){
      QM_bayes_out <- QM_bayes(stanmodel, stan_metabarcoding_data)
    } 
    
    if (!method %in% c("ML", "Bayes")) {
      print("Pick a valid method: `ML' or `Bayes'")
    }
    
    
    if (method == "ML"){
      ML_out <- QM_likelihood(stanmodel, stan_metabarcoding_data)
    } 
    
    if (method == "Bayes"){
      QM_bayes_out <- QM_bayes(stanmodel, stan_metabarcoding_data)
    } 
    
    if (!method %in% c("ML", "Bayes")) {
      print("Pick a valid method: `ML' or `Bayes'")
    }
    
    
    if (method == "ML"){
      fit_out <- list(metabarcoding_data = metabarcoding_data,
                      ML_out = ML_out)
    } 
    
    if (method == "Bayes"){
      fit_out <- list(metabarcoding_data = metabarcoding_data,
                      QM_bayes_out = QM_bayes_out)
    }
    
    return(fit_out)
  }
  
  ## ML example
  # library(here)
  # ML_out <- run_QM_model(input_metabarcoding_RDS,
  #                        input_mock_comm_RDS,
  #                        station_mapping = station_mapping,
  #                        reference_species = "Actinopteri",
  #                        stanmodel = here("scripts","models","quant_metabar_rosetta_noSampleEta.stan"),
  #                        method = "ML"
  #                )

      
  ##bayesian example
  # library(here)
  # Bayes_out <- run_QM_model(input_metabarcoding_RDS,
 #                        input_mock_comm_RDS,
 #                        station_mapping = station_mapping,
 #                        reference_species = "Actinopteri",
 #                        stanmodel = here("scripts","models","quant_metabar_rosetta_noSampleEta.stan"),
  #                       method = "Bayes"
  # )
  #
  # 
  # 
  # 
 
 library(dplyr)
 library(tidyr)
 library(stringr)
 library(purrr)
 
 summarize_bayes_output <- function(bayes_out, stan_data, station_mapping, sample_mapping) {
   # Function to summarize Bayesian output with comprehensive mapping
   
   # Check if we have the actual stanfit object or just summaries
   if (!"Bayes_estimates" %in% names(bayes_out)) {
     stop("bayes_out does not contain Bayes_estimates. Please check the structure.")
   }
   
   # 1. Extract the mean estimates from Bayes_estimates
   mean_est <- bayes_out$Bayes_estimates
   
   # Convert to tibble for easier manipulation
   mean_est_tibble <- mean_est %>%
     as.data.frame() %>%
     rownames_to_column("model_sample_name") %>%
     pivot_longer(-model_sample_name, names_to = "species", values_to = "mean_proportion")
   
   # 2. Remove duplicates from sample_mapping to avoid many-to-many relationships
   unique_sample_mapping <- sample_mapping %>%
     distinct(model_sample_name, .keep_all = TRUE)
   
   # 3. Join with sample mapping to get all metadata
   full_data <- mean_est_tibble %>%
     left_join(unique_sample_mapping, by = "model_sample_name")
   
   # 4. Create summary statistics for each sample-species combination
   summary_stats <- full_data %>%
     group_by(original_sample_name, model_sample_name, time, creek, 
              station_num, station_original, station_type, biol, species) %>%
     summarise(
       mean_prop = mean(mean_proportion),
       median_prop = median(mean_proportion),
       .groups = 'drop'
     )
   
   # 5. Create species-level summaries
   species_summary <- full_data %>%
     group_by(species) %>%
     summarise(
       max_proportion = max(mean_proportion),
       mean_proportion = mean(mean_proportion),
       .groups = 'drop'
     )
   
   # 6. Prepare data for plotting - convert time to character first for separation
   plotting_data <- summary_stats %>%
     mutate(time_char = as.character(time)) %>%
     separate(time_char, into = c("month", "year"), sep = 2, remove = FALSE) %>%
     unite(newtime, c(year, month), sep = "-", remove = FALSE)
   
   # Return all the useful data frames in a list
   return(list(
     mean_estimates = mean_est_tibble,
     full_posterior = full_data,
     summary_stats = summary_stats,
     species_summary = species_summary,
     plotting_data = plotting_data,
     sample_mapping = unique_sample_mapping  # Return for verification
   ))
 }
 
library(ggridges)
library(lubridate)
 # 1. Fixed ridgeline plot function (handles numeric time values)
plot_bayes_ridges <- function(posterior_data, 
                               target_sampletype = NULL, 
                               target_station = NULL,
                               min_proportion = 0.001,
                               time_as_factor = TRUE) {
   # Prepare data - your tibble already has the separated components
   plot_data <- posterior_data %>%
     filter(proportion > min_proportion)
   
   # Convert time to appropriate format (numeric or factor)
   if (time_as_factor) {
     plot_data <- plot_data %>% mutate(time = factor(time))
   } else {
     plot_data <- plot_data %>% mutate(time = as.numeric(time))
   }
   
   # Apply filters if specified
   if (!is.null(target_sampletype)) {
     plot_data <- plot_data %>% filter(sampletype %in% target_sampletype)
   }
   if (!is.null(target_station)) {
     plot_data <- plot_data %>% filter(stationname %in% target_station)
   }
   
   # Check if we have data after filtering
   if (nrow(plot_data) == 0) {
     warning("No data remaining after filtering. Check your target_sampletype and target_station parameters.")
     return(ggplot() + geom_blank() + labs(title = "No data available after filtering"))
   }
   
   # Create the plot
   p <- plot_data %>%
     ggplot(aes(x = proportion, y = time, fill = time)) + 
     geom_density_ridges(rel_min_height = 0.01) +
     labs(x = "Proportion of DNA after QM correction", 
          y = "Time Point", 
          fill = "Time Point",
          title = paste("Posterior Distributions:", 
                        ifelse(!is.null(target_sampletype), paste("Sample type:", target_sampletype), ""),
                        ifelse(!is.null(target_station), paste("Station:", target_station), ""))) +
     theme_bw() +
     facet_wrap(~species_name, scales = "free_x")
   
   return(p)
 }
 
 # 2. Flexible bar plot function (for summary data)
 plot_bayes_bars <- function(summary_data,
                             target_stations = NULL, 
                             target_sampletypes = NULL,
                             min_proportion = 0.001) {
   
   # Prepare data - assuming summary_data has mean_prop, species_name, station_name, etc.
   plot_data <- summary_data %>%
     filter(mean_prop > min_proportion)
   
   # Apply filters if specified
   if (!is.null(target_stations)) {
     plot_data <- plot_data %>% filter(stationname %in% target_stations)
   }
   if (!is.null(target_sampletypes)) {
     plot_data <- plot_data %>% filter(sampletype %in% target_sampletypes)
   }
   
   # Create the plot
   p <- plot_data %>%
     ggplot(aes(x = stationname, fill = species_name, y = mean_prop)) +
     geom_col(position = "stack") +
     labs(x = "Station", 
          y = "Mean Proportion of DNA after QM correction", 
          fill = "Species") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   # Add faceting if you want to separate by sampletype or time
   if ("sampletype" %in% names(plot_data) && length(unique(plot_data$sampletype)) > 1) {
     p <- p + facet_wrap(~sampletype, scales = "free_x")
   }
   
   return(p)
 }
 
 # Also need to fix the timeseries function since it has the same issue:
 plot_bayes_timeseries <- function(posterior_data,
                                   target_station = NULL,
                                   target_species = NULL,
                                   min_proportion = 0.001,
                                   highlight_period = NULL,
                                   time_as_numeric = TRUE) {
   
   # Prepare data
   plot_data <- posterior_data %>%
     filter(proportion > min_proportion)
   
   # Convert time to appropriate format
   if (time_as_numeric) {
     plot_data <- plot_data %>% mutate(time_num = as.numeric(time))
   } else {
     plot_data <- plot_data %>% mutate(time_num = factor(time))
   }
   
   # Apply filters if specified
   if (!is.null(target_station)) {
     plot_data <- plot_data %>% filter(station_name %in% target_station)
   }
   if (!is.null(target_species)) {
     plot_data <- plot_data %>% filter(species_name %in% target_species)
   }
   
   # Calculate summary statistics for plotting
   summary_data <- plot_data %>%
     group_by(time_num, station_name, species_name, sampletype) %>%
     summarise(
       mean_prop = mean(proportion),
       median_prop = median(proportion),
       lower_95 = quantile(proportion, 0.025),
       upper_95 = quantile(proportion, 0.975),
       .groups = 'drop'
     )
   
   # Create the plot
   p <- summary_data %>%
     ggplot(aes(x = time_num, y = mean_prop, color = species_name)) +
     geom_point() +
     geom_errorbar(aes(ymin = lower_95, ymax = upper_95), alpha = 0.5) +
     labs(x = "Time Point", 
          y = "Proportion of DNA after QM correction", 
          color = "Species") +
     theme_bw()
   
   # Facet by species if multiple species are shown
   if (length(unique(summary_data$species_name)) > 1) {
     p <- p + facet_wrap(~species_name, scales = "free_y")
   }
   
   # Add highlight period if specified (only works if time is numeric)
   if (!is.null(highlight_period) && time_as_numeric) {
     p <- p + 
       annotate("rect", 
                xmin = highlight_period$start, 
                xmax = highlight_period$end, 
                ymin = 0, ymax = Inf, 
                alpha = 0.2, fill = "darkred")
   }
   
   return(p)
 }
 
 # 4. Sample composition plot by sampletype and station
 plot_sample_composition <- function(posterior_data, 
                                     target_sampletypes = NULL,
                                     target_stations = NULL,
                                     min_proportion = 0.01,
                                     facet_by_station = TRUE,
                                     show_biol_reps = TRUE) {
   
   # Calculate mean proportions for each sample-species combination
   summary_data <- posterior_data %>%
     group_by(sample_name, time, sampletype, stationname, biol, species_name) %>%
     summarise(mean_prop = mean(proportion), .groups = 'drop') %>%
     filter(mean_prop > min_proportion)
   
   # Apply filters if specified
   if (!is.null(target_sampletypes)) {
     summary_data <- summary_data %>% filter(sampletype %in% target_sampletypes)
   }
   if (!is.null(target_stations)) {
     summary_data <- summary_data %>% filter(stationname %in% target_stations)
   }
   
   # Create a better x-axis label that shows biological replicates clearly
   if (show_biol_reps) {
     summary_data <- summary_data %>%
       mutate(sample_label = paste("BioRep", biol))
   } else {
     summary_data <- summary_data %>%
       mutate(sample_label = sample_name)
   }
   
   # Create the plot
   p <- summary_data %>%
     ggplot(aes(x = sample_label, fill = species_name, y = mean_prop)) +
     geom_col(position = "stack") +
     labs(x = "Biological Replicate", 
          y = "Mean Proportion", 
          fill = "Species",
          title = paste("Sample Composition:",
                        ifelse(!is.null(target_sampletypes), 
                               paste("Type:", paste(target_sampletypes, collapse = ", ")), ""),
                        ifelse(!is.null(target_stations),
                               paste("Stations:", paste(target_stations, collapse = ", ")), ""))) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   # Add faceting by station if requested
   if (facet_by_station && "station_name" %in% names(summary_data)) {
     p <- p + facet_wrap(~station_name, scales = "free_x")
   }
   
   # Alternatively, facet by sampletype if not faceting by station
   else if (!facet_by_station && "sampletype" %in% names(summary_data)) {
     p <- p + facet_wrap(~sampletype, scales = "free_x")
   }
   
   return(p)
 }
 
 # 5. Enhanced version with both sampletype and station faceting
 plot_sample_composition_grid <- function(posterior_data, 
                                          target_sampletypes = NULL,
                                          target_stations = NULL,
                                          min_proportion = 0.01) {
   
   # Calculate mean proportions
   summary_data <- posterior_data %>%
     group_by(sample_name, time, sampletype, stationname, biol, species_name) %>%
     summarise(mean_prop = mean(proportion), .groups = 'drop') %>%
     filter(mean_prop > min_proportion)
   
   # Apply filters
   if (!is.null(target_sampletypes)) {
     summary_data <- summary_data %>% filter(sampletype %in% target_sampletypes)
   }
   if (!is.null(target_stations)) {
     summary_data <- summary_data %>% filter(stationname %in% target_stations)
   }
   
   # Create the plot with grid faceting
   p <- summary_data %>%
     mutate(biol_label = paste("BioRep", biol)) %>%
     ggplot(aes(x = biol_label, fill = species_name, y = mean_prop)) +
     geom_col(position = "stack") +
     labs(x = "Biological Replicate", 
          y = "Mean Proportion", 
          fill = "Species",
          title = "Sample Composition by Station and Sample Type") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
           strip.text = element_text(size = 9)) +
     facet_grid(sampletype ~ stationname, scales = "free_x", space = "free_x")
   
   return(p)
 }
 
 
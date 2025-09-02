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


## NGN Run QM Model
# Author: Eily Allan - but really all from Ryan Kelly and Ole Shelton 
# Person running: Eily
# Last modified: 10/17/22 by Eily
# Date of run: 10/17/22 by Eily 

# Overview 
# 

# Inputs: 
# 1)  

# Outputs: 
# 1) 

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(rstan)
library(MCMCpack) #for rdirichelet function
library(here)
library(gridExtra)
library(unikn)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here("Scripts","calibrate_metabarcoding.R"))

# read in taxa table 
taxa_table_2  <- readRDS( here("analysis","combined_data_long_reads_c.RDS"))
colnames(taxa_table_2)

mock_2 <- readRDS("/Users/zgold/Documents/GitHub/Metabarcodings_Signal_from_Noise/signal2noise/data/mifish_mock_community_data.RDS")

enviro_2 <- taxa_table_2 %>% 
  filter(., River %in% c('Twentymile','Kenai')) %>% 
  unite(col="station",c(Site,Depth), sep=".") %>% 
  dplyr::select(time =Sample_date,creek=River, station =station,biol=Bio_rep,tech=Tech_rep, species=Species, Nreads= nReads) %>% 
  dplyr::mutate(tech = replace_na(tech, 1))

# only focus on overlapping species

mock_2 <- mock_2 %>% 
  filter(., Cycles ==39) %>% 
  filter(., str_detect(community,"North")) %>% 
  mutate(., CommType = case_when(str_detect(community,"Even")~"Even",
                                 str_detect(community,"Skew_1")~"Skew_1",
                                 str_detect(community,"Skew_2")~"Skew_2")) %>%
  filter(., CommType !="Skew_1") %>% 
 filter(., ID_mifish %in% c("Oncorhynchus gorbuscha" ,  "Oncorhynchus keta"   ,     "Oncorhynchus kisutch"   ,  "Oncorhynchus mykiss"   ,   "Oncorhynchus nerka" ,     
                                "Oncorhynchus tshawytscha")) %>% 
  group_by(community, tech_rep) %>% 
  mutate(., totReads=sum(nReads),
         propReads = nReads/totReads,
         tot_start = sum(start_conc_ng),
         start_prop=start_conc_ng/tot_start) %>% 
  dplyr::select(Community=community, CommType,Tech_Rep=tech_rep,species=ID_mifish, nReads,totReads,propReads,start_prop, N_pcr_mock=Cycles) %>% 
  mutate(Tech_Rep = as.numeric(Tech_Rep),
         N_pcr_mock = as.numeric(N_pcr_mock)) 
  
# Prepare for stan model 
input_metabarcoding_RDS <- enviro_2
input_mock_comm_RDS <- mock_2
format_metabarcoding_data <- function(input_metabarcoding_RDS, input_mock_comm_RDS){
  require(tidyverse)
  
  Observation <- input_metabarcoding_RDS 
  Mock <- input_mock_comm_RDS %>% 
    mutate(time = 1,    #for formatting compatibility
           creek = 1,
           biol = 1) %>% 
    rename(Nreads = nReads) %>% 
    #unite(c(Genus, Species), col = "species", sep = " ") %>% 
    unite(c(Community, CommType), col = "Community", sep = "_") %>% 
    rename(station = Community) %>% 
    replace(is.na(.), 0) %>% 
    filter(start_prop > 0) %>% 
    filter(Nreads > 0) %>% #omit things that are absent from the mocks
    filter(! is.na(Nreads)) %>% 
    rename(tech = Tech_Rep,
           b_proportion = start_prop) %>% 
    #mutate(stationname = station) %>% 
    unite(c(time, creek, station,biol, tech), col = "Sample", sep = ".", remove = F) 
  #unite(c(time, creek, stationname, biol), col = "Sample", sep = ".", remove = F) %>%
  #rename(station = stationname) 
  
  # only keep species present in both mocks and observations
  keepSpecies <- intersect(Mock$species, Observation$species)
  Observation <- Observation %>% 
    filter(species %in% keepSpecies) 
  Mock <- Mock %>% filter(species %in% keepSpecies)

  # index species to a common standard 
  sp_list <- data.frame(
    species = c(Mock$species, Observation$species) %>% unique(),
    species_idx = NA)
  sp_list$species_idx <- match(sp_list$species, unique(sp_list$species)) 
  
  #reindex and renormalize to deal with omitted species
  Mock <- Mock %>% 
    left_join(sp_list) %>% 
    mutate(station = match(station, unique(station)),
           speciesname = species,
           species = species_idx) %>% 
    group_by(station, tech, time, creek, biol) %>% 
    mutate(b_proportion = b_proportion/sum(b_proportion)) %>% 
    ungroup()
  
  Observation <- Observation %>% 
    left_join(sp_list) %>% 
    mutate(stationname = station,
           speciesname = species,
           species = species_idx) %>% 
    mutate(station = case_when(stationname == "Knik Arm.Surface" ~ 1,
                               stationname == "Little Susitna.Surface" ~ 2,
                               stationname == "Susitna.Surface" ~ 3,
                               stationname == "Twentymile.Surface" ~ 4,
                               stationname == "Bridge.Surface" ~ 5,
                               stationname == "Dock.Surface" ~ 6,
                               stationname == "Bridge.5m" ~ 7,
                               stationname == "Dock.5m" ~ 8)) %>% 
    group_by(station, time, creek,biol) %>% 
    mutate(tech = match(tech, unique(tech)),
           tot_reads=sum(Nreads)) %>% 
    filter(., tot_reads >0) %>% 
    ungroup() %>% 
    dplyr::select(-tot_reads)
  
  station_list <- data.frame(
    station = Observation$stationname %>% unique(),
    station_idx = Observation$station %>% unique())
  
  #list object containing named elements Observation, Mock, and Npcr
  
  return(
    metabarcoding_data <- list(
      Observation = Observation,
      Mock = Mock,
      N_pcr_mock = 39, 
      NSpecies = nrow(sp_list),
      station_list = station_list,
      sp_list = sp_list
    ))
}



qmdata <- format_metabarcoding_data(enviro_2, mock_2)
obs <- qmdata

N_pcr_cycles <- 39
MOCK <- qmdata$Mock   

alrTransform <- function(MOCK){
  require(tidyverse)
  require(compositions)
  
  p_mock <- MOCK %>% 
    dplyr::select(species, station, tech, b_proportion) %>% 
    rename(tech_rep = tech) %>% 
    pivot_wider(names_from = species, values_from = b_proportion, values_fill = 1e-9) %>% 
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
    pivot_wider(names_from = species, values_from = Nreads, values_fill = 0)
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
    pivot_wider(names_from = species, values_from = Nreads, values_fill = 0)
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


stan_metabarcoding_data <- makeDesign(qmdata, N_pcr_cycles = 39)
stan_metabarcoding_data$alr_mock_true_prop
stan_metabarcoding_data$mock_data
stan_metabarcoding_data$sample_vector

####################################################################
# Run ML QM model
####################################################################
stanmodelname <- "/Users/zgold/Documents/GitHub/quantitative_salmon_culverts/Scripts/functions/quant_metabar_rosetta_noSampleEta.stan"
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

ML_out <- QM_likelihood("/Users/zgold/Documents/GitHub/quantitative_salmon_culverts/Scripts/functions/quant_metabar_rosetta_noSampleEta.stan", stan_metabarcoding_data)

#write_rds(ML_out, here("Output","metabarcoding","ML_out.RDS"))

ML_out$ML_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE, sep="_") %>% 
  mutate(station = case_when(station ==1 ~"Knik Arm.Surface",
                             station ==2 ~"Little Susitna.Surface",
                             station ==3 ~"Susitna.Surface",
                             station ==4 ~"Twentymile.Surface",
                             station ==5 ~"Bridge.Surface",
                             station ==6 ~"Dock.Surface",
                             station ==7 ~"Bridge.5m",
                             station ==8 ~"Dock.5m")) %>% 
  filter(value > 0.001) %>%
  filter(., creek %in% c('Twentymile','Kenai')) %>% 
  ggplot(aes(x = biol, fill = species, y = value)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~creek ~station ~time) 

ggsave(here("Figures","20230225_proportions_after_qm_ML.png"))

####################################################################
# Run Bayesian QM model
####################################################################

bayes_out <- QM_bayes("/Users/zgold/Documents/GitHub/quantitative_salmon_culverts/Scripts/functions/quant_metabar_rosetta_noSampleEta.stan", stan_metabarcoding_data)

saveRDS(bayes_out, here("analysis","20231210_bayesout_salmon_only.RDS"))

summaryout <- summary(bayes_out$Bayes_modelfit)$summary
write.csv(summaryout, here("analysis","20231211_bayesout_salmon_eul_only_summary.csv"))

bayes_out <- readRDS( here("analysis","20231210_bayesout_salmon_only.RDS"))

bayes_out$Bayes_alpha_est
bayes_out$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") -> tibble_bayes

bayes_out$Bayes_25ci %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") -> tibble_bayes_25

bayes_out$Bayes_75ci %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") -> tibble_bayes_75

tibble_bayes$Bayes_25ci <- tibble_bayes_25$value
tibble_bayes$Bayes_75ci <- tibble_bayes_75$value

bayes_out$Bayes_modelfit
pars2 <- rstan::extract(bayes_out$Bayes_modelfit, par = "int_samp_small")

species_mapping <- stan_metabarcoding_data$sp_list
Stations <-  stan_metabarcoding_data$sample_vector %>% as.tibble()
Stations %>% 
  mutate(., station_idx =row_number())-> Stations

pars2$int_samp_small %>% 
    as_tibble() %>% 
    pivot_longer(., cols = `1.1`:`197.6`, names_to="station.species", values_to = "est") %>%
    separate(station.species, into = c("station_idx","species_idx"), remove = F, sep="\\.") %>%
    mutate(., station_idx = as.numeric(station_idx),
           species_idx = as.numeric(species_idx)) %>%
     left_join(species_mapping) %>%
       left_join(Stations) -> int_samp_small_tibble

saveRDS(int_samp_small_tibble, here("analysis","int_samp_small_tibble.RDS"))

summary(bayes_out$Bayes_modelfit, par = "int_samp_small")$summary[,1] %>%
  matrix(ncol = stan_metabarcoding_data$N_species, byrow = TRUE) %>% 
  as.data.frame()
names(mean_est) <- stan_metabarcoding_data$sp_list$species
rownames(mean_est) <- stan_metabarcoding_data$sample_vector

tibble_bayes %>% 
  group_by(species) %>% 
  dplyr::summarise(max(value))

int_samp_small_tibble %>% 
  group_by(species) %>% 
  dplyr::summarise(mean(est))

tibble_bayes %>% 
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE,sep = "_") %>% 
  mutate(station = case_when(station == "Knik Arm" ~ 1,
                             station == "Little Susitna" ~ 2,
                             station == "Susitna" ~ 3,
                             station == "Twentymile" ~ 4,
                             station == "Bridge" ~ 5,
                             station == "Dock" ~ 6)) %>% 
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE) -> qm_data_plotting

saveRDS(qm_data_plotting, here("analysis","qm_data_plotting"))

library(ggridges)

int_samp_small_tibble %>% 
  mutate(., sample=value) %>% 
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE,sep = "_") %>% 
  mutate(station = case_when(station ==  1~ "Knik Arm",
                             station == 2 ~ "Little Susitna",
                             station == 3~"Susitna",
                             station == 4~"Twentymile",
                             station == 5~"Bridge",
                             station == 6~"Dock")) %>% 
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  filter(., creek %in% c('Twentymile')) %>% 
  filter(., station=="Twentymile") %>% 
  ungroup() %>% 
  mutate(time=as.Date(time)) %>% 
  ggplot(., aes(x = est, y = time,fill = time)) + 
  geom_density_ridges( rel_min_height = 0.01) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  labs(x="Proportion of DNA after QM correction", fill="Species", color="") +
  theme_bw() +
  facet_wrap(~species)


qm_data_plotting %>% 
  filter(value > 0.001) %>% 
  filter(., creek %in% c('Twentymile','Kenai')) %>% 
  mutate(facetorder = factor(creek, levels=c('Twentymile','Kenai'))) %>% 
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
  group_by(creekstntime) %>% 
  mutate(sumprop = sum(value)) %>% 
  mutate(avgprop = value/sumprop) %>%
  ggplot(aes(x = newtime, fill = species, y = avgprop)) +
  geom_col() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_grid(~facetorder) +
  ylab("Proportion of DNA after QM correction") +
  labs(x="Date (YY-MM)", y="Proportion of DNA after QM correction", fill="Species", color="") %>% 
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  theme_bw()

ggsave(here("Figures","20230225_proportions_after_qm.png"), width = 12, height = 8)

qm_data_plotting$newtime %>% class()
qm_data_plotting %>% 
  filter(value > 0.001) %>% 
  filter(., creek %in% c('Twentymile')) %>% 
  mutate(facetorder = factor(creek, levels=c('Twentymile','Kenai'))) %>% 
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
  group_by(creekstntime) %>% 
  mutate(sumprop = sum(value)) %>% 
  mutate(avgprop = value/sumprop) %>%
  mutate(., Sample_date=as.Date(time)) %>% 
  ggplot(aes(x = Sample_date, fill = species, y = avgprop)) +
  geom_point() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_grid(~species) +
  ylab("Proportion of DNA after QM correction") +
  labs(x="Date (YY-MM)", y="Proportion of DNA after QM correction", fill="Species", color="") +
  theme_bw() + annotate("rect", xmin = as.Date("2021-08-24"), xmax = as.Date("2021-10-21"), ymin = 0, ymax = 0.4, alpha = .2,fill ="darkred") 

ggsave(here("Figures","20230225_proportions_after_qm.png"), width = 12, height = 8)


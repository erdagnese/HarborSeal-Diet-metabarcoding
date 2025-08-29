### script to make the CB wild_sample metadata

library(here)

wild_md <- read.table(here("data", "Wild_CB_Samples", "tourmaline-18S_CB2","00-data","metadata.tsv"), sep = "\t", header = T)
scatlog <- read.csv(here("input-files","ScatLogMetadata.csv"), header=T,check.names = F)
scatrunlist <- read.csv(here("input-files","run_list_metadata.csv"), header=T,check.names = F)

library(tidyverse)
library(stringr)

scatlog_md <- left_join(scatlog,scatrunlist, by = "Sample_ID")
scatlog_md %>% 
  mutate(sampleid = ifelse(is.na(sampleid),Sample_ID, sampleid)) %>% 
  select(!year) %>% 
  separate(Collection_Date, into = c("day","month","year"),
           sep = "-", remove = TRUE) %>% 
  mutate(year = ifelse(str_length(year) == 2, paste0("20", year), year)) -> scatlog_md

# Create a clean version of sampleid in wild_md that matches scatlog_md format
wild_md2 <- wild_md %>%
  mutate(
    # Standardize case
    samplename = ifelse(!str_starts(sampleid, "Kang"), 
                        str_replace(sampleid, "^PV", "Pv"), 
                        sampleid),
    
    # Extract tech number
    tech = case_when(
      str_starts(samplename, "Kang") ~ "1",
      str_detect(samplename, "_\\d+$") ~ str_extract(samplename, "\\d+$"),
      str_detect(samplename, "-\\d+$") ~ str_extract(samplename, "\\d+$"),
      TRUE ~ "1"
    ),
    
    # Remove tech numbers and _A suffixes
    clean_sampleid = ifelse(str_starts(samplename, "Kang"),
                            samplename,
                            str_replace(samplename, "[_\\-]\\d+$", "")),
    
    clean_sampleid = ifelse(str_detect(clean_sampleid, "^Pv18.*_A"),
                            str_replace(clean_sampleid, "_A.*", ""),
                            clean_sampleid)
  ) %>%
  # Handle Pv18 number padding separately
  mutate(
    pv18_num = ifelse(str_detect(clean_sampleid, "^Pv18[_\\-]"),
                      str_extract(clean_sampleid, "(?<=[_\\-])\\d+"),
                      NA_character_),
    
    clean_sampleid = ifelse(!is.na(pv18_num),
                            paste0("Pv18-", str_pad(pv18_num, width = 4, pad = "0", side = "left")),
                            clean_sampleid)
  ) %>%
  select(-pv18_num)

#let's save this because we have done this so many times it's frankly dumb not to have it for use
write.csv(scatlog_md, here("input-files","SPS_scatlog_runlist_metadata.csv"), row.names = F)

# Now you can join using the clean_sampleid
wild_md3 <- wild_md2 %>%
  left_join(scatlog_md, by = c("clean_sampleid" = "sampleid")) %>% 
  mutate(
    month = case_when(
      month == "Jan" ~ "01",
      month == "Feb" ~ "02", 
      month == "Mar" ~ "03",
      month == "Apr" ~ "04",
      month == "May" ~ "05",
      month == "Jun" ~ "06",
      month == "Jul" ~ "07",
      month == "Aug" ~ "08",
      month == "Sep" ~ "09",
      month == "Oct" ~ "10",
      month == "Nov" ~ "11",
      month == "Dec" ~ "12",
      TRUE ~ month  # fallback
    )) %>% 
  mutate(
    Sample_ID = ifelse(is.na(Sample_ID), samplename, Sample_ID)) %>% 
  mutate(
    Collection_Location = case_when(
      str_detect(sampleid, "Kang") ~ "ControlSample",
      TRUE ~ Collection_Location  # keep original value if not Kang
    )
  ) %>% 
  mutate(
    # Ensure day, month_number, and year are all character for pasting
    across(c(day, month, year), as.character),
    
    # Create date column (YYYY-MM-DD format)
    date = ymd(paste(year, month, day, sep = "-")),
    
    # Alternative if you want numeric date (days since epoch)
    date_numeric = as.numeric(date)
  ) %>% 
  mutate(
    Collection_Location = str_remove_all(Collection_Location, "\\s+")
  ) %>% 
  group_by(Collection_Location, date) %>%
  mutate(
    biol = row_number()  # running number within each location-date group
  ) %>%
  ungroup()

write.csv(wild_md3, here("input-files","FULL_CBrun2_metadata.csv"),row.names = F)
#this will map back to the sample names and metadata


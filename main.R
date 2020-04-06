# Extract DHS treatment parameters and ACT/non-ACT drug coverage, then prepare malaria model treatment coverage inputs
# Author: AB Hogan
# Date: 6 April 2020

# METHODOLOGY

# For each country, extract the proportion of fevers that sought medical treatment, and the proportion that received anti-malarial that is ACT. Multiply by 80% (estimated scaling factor for receiving appropriate treatment). DHS weightings are used to aggregate individual-level data at the country level.

# Packages
library(rdhs)
library(tidyverse)
library(readr)

# Selects desired surveys and years
survs <- dhs_surveys(countryIds = c("AF", "AO", "BD", "BJ", "BF", "BO", "BR", "BU", "CF", "CM", "CI", "CG", "CD", "CO", "DR", "EC", "EK", "ER", "ET", "GA", "GH", "GM", "GN", "GU", "GY", "HT", "HN", "IA", "ID", "KE", "KH", "KM", "LA", "LB", "MD", "MW", "ML", "MM", "MN", "MR", "MZ", "NG", "NM", "NC", "NI", "PE", "PH", "PK", "RW", "SD", "SN", "SL", "TD", "TZ", "TG", "TH", "TL", "VN", "YE", "UG", "ZM", "ZW"), surveyYear = c(2010:2019))

# Selects the desired datasets
datasets <- dhs_datasets(
  surveyIds = survs$SurveyId,
  fileFormat = "FL",
  fileType = "KR")

# Downloads the chosen datasets and selects variables
downloads <- get_datasets(datasets$FileName)

question_labels <- search_variable_labels(datasets$FileName, search_terms = c("artemisinin"))
vars <- c("v005", "b5", "b8", "h22", "h32z", "ml13e")
questions <- search_variables(datasets$FileName, variables = vars)

# Print the list of questions (note: assumes questions/codes are the same across all surveys, but this is not always the case)
print(questions[1:20, 1:3])

# Extracts the data
extract <- extract_dhs(questions, add_geo = TRUE)
extract_bound <-
  rbind_labelled(
    extract$AOKR71FL,
    extract$BJKR71FL,
    extract$BFKR7AFL,
    extract$BUKR70FL,
    extract$CMKR61FL,
    extract$CIKR62FL,
    extract$CDKR61FL,
    extract$GNKR62FL,
    extract$GHKR7BFL,
    extract$KEKR72FL,
    extract$LBKR71FL,
    extract$MWKR7IFL,
    extract$MLKR7HFL,
    extract$MZKR7AFL,
    extract$NGKR7AFL,
    extract$SLKR72FL,
    extract$TZKR7BFL,
    extract$TGKR71FL,
    extract$UGKR7BFL,
    extract$ZMKR71FL,
    extract$CMKR61FL,
    extract$CGKR61FL,
    extract$ETKR71FL,
    extract$TDKR71FL,
    extract$GAKR61FL,
    extract$GMKR61FL,
    extract$MDKR71FL,
    extract$SNKR7ZFL,
    extract$ZWKR72FL
  )

# save raw extracted data
saveRDS(extract_bound,  file = "DHS_treatment.RDS")

# read in
extract_bound <- readRDS("DHS_treatment.RDS")
countrylist <- read_csv("countrylist.csv")

# aggregate at country level using DHS weightings
dat_had_fever <- extract_bound %>%
  group_by(SurveyId) %>%
  filter(h22 == 1, b5 == 1, b8 < 5, b8 >=1) %>%
  summarise(num_children_fever = sum(v005/1e6))

dat_received_treatment <- extract_bound %>%
  group_by(SurveyId) %>%
  filter(h22==1, b5 ==1, h32z==1, b8 < 5, b8 >=1) %>%
  summarise(received_treatment = sum(v005/1e6))

dat_received_act <- extract_bound %>%
  group_by(SurveyId) %>%
  filter(h22 == 1, b5 == 1, b8 < 5, b8 >=1, ml13e == 1, h32z == 1) %>%
  summarise(received_act = sum(v005/1e6))

# calculate coverages and rearrange so in appropriate format
dat_all <- left_join(dat_had_fever, dat_received_treatment, by = "SurveyId") %>%
  left_join(dat_received_act, by = "SurveyId") %>%
  mutate(drug_cov = received_treatment/num_children_fever,
         drug_cov_1_0 = received_act/num_children_fever,
         drug_cov_0_0 = drug_cov - drug_cov_1_0) %>%
  dplyr::select(-c(num_children_fever, received_treatment, received_act, drug_cov)) %>%
  mutate(COUNTRYCODE = substr(SurveyId, 1, 2))

dat_out <- countrylist %>%
  left_join(dat_all, by = "COUNTRYCODE") %>%
  dplyr::select(ISO, NAME_0, CONTINENT, drug_cov_0_0, drug_cov_1_0) %>%
  mutate(drug_0_efficacy = 0.75,
         drug_1_efficacy = 0.95)

# read in hardcoded values for countries where no data available
#...

# save final data for model inputs
write_csv(dat_out, "drug_coverage_2020.csv")


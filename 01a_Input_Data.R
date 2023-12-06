##########################################################
# Name of file: 01a_Input_Data.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@phs.scot
# Original date: 09 November 2021
# Latest update author (if not using version control) - Chris Robertson chrisobertson@nhs.net
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats

# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the cohort and merges in the Q Covid risk groups
#                         reads in vaccination and testing data
#                         selects only those records belonging to children and young people under 18
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(tidyverse)
library(lubridate)
library(survival)
library(dplyr)
library(labelled)
library(readr)
library(survminer)
#Load data

min_age = 0 # Minimum age of the cohort we are interested in
max_age = 100 # Maximum age of the cohort we are interested in

urti_codes = c("J00", "J02", "J03", "J04", "J05", "J06")
flu_codes = c("J09", "J10", "J11", "J12",
               "J13", "J14", "J15", "J16", "J17", "J18")
bronch_codes = c("J20", "J21", "J40")
asthma_codes = c("J45", "J46")
covid_three_char_codes =  c("U08", "U09","U10")
unspecified_codes = c("J22")
three_char_admittance_codes = c(urti_codes, flu_codes, bronch_codes, unspecified_codes)
                      # 3 character ICD-10 codes

covid_four_char_codes = c("U071", "U072")
rsv_codes = c("J121", "J125", "J210") #, "B794") 
strep_codes = c( "B950", "J020", "J029", "J030", "J039", "A38X")
wheeze_codes = c("R062")

four_char_admittance_codes = c(covid_four_char_codes, rsv_codes)

# These are 4 characters so need a different check

Location <- "/conf/"  # Server
setwd("/conf/EAVE/GPanalysis/progs/TM/winter_pressure_hospitalisation_description/")
#Location <- "//isdsf00d03/"  # Desktop
a_begin = as.Date("2022-09-01")
a_begin_minus_6_months = as.Date("2022-03-01")
a_end = as.Date("2023-01-31")
a_analysis_date <- a_end #change this to try and recreate historical analyses
Vaccinations = readRDS(paste0(Location, "EAVE/GPanalysis/data/temp/vaccine_cleaned.rds")) %>% select(-health_board)

EAVE_endpoints <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/severe_endpoints2022-06-23.rds")) #n=770,429rows

#smr01 <- readRDS(paste0(Location,"EAVE/GPanalysis/data/SMR01_allstays.rds")) #SMR hospital data
#smr01 = readRDS("/conf/EAVE/GPanalysis/data/smr01_2023_03_30.rds")
#smr01 <- smr01 %>% mutate(ADMISSION_DATE = as_date(ADMISSION_DATE), DISCHARGE_DATE = as_date(DISCHARGE_DATE))

# We need to handle CIS stuff
# We do this by taking the admission and discharge dates from the 1st and last entries
# and then the condition codes from the first entry
#smr01 = smr01 %>% arrange(EAVE_LINKNO, CIS_MARKER) %>% group_by(EAVE_LINKNO, CIS_MARKER) %>%
#  summarise(ADMISSION_DATE = min(ADMISSION_DATE), DISCHARGE_DATE = max(DISCHARGE_DATE), 
#  MAIN_CONDITION = first(MAIN_CONDITION), OTHER_CONDITION_1 = first(OTHER_CONDITION_1),
#  OTHER_CONDITION_2 = first(OTHER_CONDITION_2), OTHER_CONDITION_3 = first(OTHER_CONDITION_3),
#  OTHER_CONDITION_4 = first(OTHER_CONDITION_4), OTHER_CONDITION_5 = first(OTHER_CONDITION_5),
#  ADMISSION_TYPE = first(ADMISSION_TYPE)) %>% ungroup()

#write_rds(smr01, "smr01_cis")
smr01 = readRDS("smr01_cis")

print("Number of CIS rows in SMR01")
print(nrow(smr01))

# How many emergency hospitalisations during this time?
z_num_emergency_hosps = smr01 %>% filter(ADMISSION_DATE >= a_begin & ADMISSION_DATE <= a_end) %>%
  filter(ADMISSION_TYPE >= 30 & ADMISSION_TYPE <= 40)
print("Number of emergency hospitalisations during time period")
print(nrow(z_num_emergency_hosps))

EAVE_cohort_refresh = readRDS(paste0(Location, "EAVE/GPanalysis/data/EAVE_LINKNO_refresh.rds")) %>%
  arrange(EAVE_LINKNO, desc(DATE_TRANSFER_OUT)) %>%
  filter(!duplicated(EAVE_LINKNO))

# Select people currently living in Scotland
EAVE_cohort_refresh = EAVE_cohort_refresh %>%
  filter(unvalidatedCHI_flag != 1) %>%
  filter(is.na(NRS.Date.Death) | NRS.Date.Death > a_begin) %>%
  filter(is.na(DATE_TRANSFER_IN) | DATE_TRANSFER_IN < a_begin) %>%
  filter(is.na(DATE_TRANSFER_OUT) | DATE_TRANSFER_OUT > a_end) %>%
  select(-unvalidatedCHI_flag, -chili_checked)

#EAVE_demographics <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_demographics_SK.rds")
EAVE_demographics = readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Dates2021-07-28.rds")) %>%
  mutate(ur6_2016 = as.numeric(substring(ur6_2016_name, 1, 2))) %>%
  select(-NRS.Date.Death)
#  select(EAVE_LINKNO, Sex, ageYear, simd2020_sc_quintile, )

z_healthboard_lookup <- read_csv("/conf/EAVE/GPanalysis/data/lookups/Datazone2011lookup.csv") %>%
  select(DZ2011_Code, HB_Name) %>%
  rename(DataZone = DZ2011_Code, health_board = HB_Name)# %>%
 # filter(!duplicated(datazone2011))

EAVE_cohort = EAVE_cohort_refresh %>% left_join(EAVE_demographics, by=join_by(EAVE_LINKNO_old == EAVE_LINKNO)) %>%
  left_join(z_healthboard_lookup, by="DataZone") %>%
  filter(!is.na(ur6_2016))

EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))
EAVE_cohort  <- EAVE_cohort %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
EAVE_cohort$eave_weight[is.na(EAVE_cohort$eave_weight)] <- mean(EAVE_cohort$eave_weight, na.rm=T)

# Get the people we're interested in
EAVE_cohort = EAVE_cohort %>%
  mutate(ageYear = ageYear + 3) %>% # This is the age at march 2020
  filter(!is.na(ageYear) & ageYear >= min_age & ageYear <= max_age)

# Sort out the EAVE weights - we look in our datasets to see if we have a record of a person
# and if so we give them a weight of 1
bnf <- readRDS(paste0(Location,"EAVE/GPanalysis/data/BNF_paragraphs.rds"))
#those with a covid test
cdw_full  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/CDW_full.rds"))
cdw_full <- cdw_full %>% mutate(date_ecoss_specimen = as_date(date_ecoss_specimen)) %>% 
  filter(date_ecoss_specimen <= a_analysis_date)
#all deaths
all_hospitalisations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/automated_any_hospitalisation_post_01022020.rds"))
pis_asthma <- readRDS(paste0(Location,"EAVE/GPanalysis/data/PIS_ASTHMA_2021-09-03.rds"))
#there's an error with the mutate function below so ignore lines 193-194 and just directly remove IDs below
pis_asthma <- pis_asthma %>% mutate(dispensed_full_date = as_date(dispensed_full_date, format="%Y%m%d")) %>% 
  filter(dispensed_full_date >= as_date("2019-03-01"))
all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds")) %>% 
  filter(!duplicated(EAVE_LINKNO))
#all known to exist - give a weight of 1 and downweight the rest
z_ids <- c(Vaccinations$EAVE_LINKNO, all_deaths$EAVE_LINKNO, 
           cdw_full$EAVE_LINKNO, all_hospitalisations$EAVE_LINKNO, bnf$EAVE_LINKNO, pis_asthma$EAVE_LINKNO) %>% unique()
#summary(filter(EAVE_cohort, !(EAVE_LINKNO %in% z_ids))$eave_weight)
z_N <- round(sum(EAVE_cohort$eave_weight) )
z_k <- sum(EAVE_cohort$EAVE_LINKNO %in% z_ids)
z_m <- round(sum(filter(EAVE_cohort, (EAVE_LINKNO %in% z_ids))$eave_weight))
z <- EAVE_cohort %>% mutate(ew = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight*(z_N - z_k)/(z_N - z_m)) )
EAVE_cohort <- z %>% dplyr::select(-eave_weight) %>% dplyr::rename(eave_weight=ew)

# See who was admitted in the past few months
z_hospitalised = smr01 %>% filter(ADMISSION_DATE >= a_begin & ADMISSION_DATE <= a_end) %>%
  filter(as.numeric(difftime(DISCHARGE_DATE, ADMISSION_DATE, units="days")) > 0) %>%
  filter(ADMISSION_TYPE >= 30 & ADMISSION_TYPE <= 40)

print("Number of rows within the time period")
print(nrow(z_hospitalised))

# Tag people who were admitted for certain things
z_hospitalised = z_hospitalised %>% mutate(acute_resp_admission = if_else(
  substr(MAIN_CONDITION, 0, 3) %in% three_char_admittance_codes |
  substr(MAIN_CONDITION, 0, 4) %in% four_char_admittance_codes, 1, 0))

z_hospitalised = z_hospitalised %>% mutate(
  acute_resp_admission_1 = if_else(
  substr(OTHER_CONDITION_1, 0, 3) %in% three_char_admittance_codes |
  substr(OTHER_CONDITION_1, 0, 4) %in% four_char_admittance_codes, 1, 0),
  acute_resp_admission_2 = if_else(
  substr(OTHER_CONDITION_2, 0, 3) %in% three_char_admittance_codes |
  substr(OTHER_CONDITION_2, 0, 4) %in% four_char_admittance_codes, 1, 0),
  acute_resp_admission_3 = if_else(
  substr(OTHER_CONDITION_3, 0, 3) %in% three_char_admittance_codes |
  substr(OTHER_CONDITION_3, 0, 4) %in% four_char_admittance_codes, 1, 0),
  acute_resp_admission_4 = if_else(
  substr(OTHER_CONDITION_4, 0, 3) %in% three_char_admittance_codes |
  substr(OTHER_CONDITION_4, 0, 4) %in% four_char_admittance_codes, 1, 0),
  acute_resp_admission_5 = if_else(
  substr(OTHER_CONDITION_5, 0, 3) %in% three_char_admittance_codes |
  substr(OTHER_CONDITION_5, 0, 4) %in% four_char_admittance_codes, 1, 0)) %>%
  mutate(acute_resp_secondary = if_else(
    acute_resp_admission_1 == 1 | acute_resp_admission_2 == 1 |
    acute_resp_admission_3 == 1 | acute_resp_admission_4 == 1 | acute_resp_admission_5 == 1, 1, 0)
  ) %>%
  select(-acute_resp_admission_1, -acute_resp_admission_2, -acute_resp_admission_3, -acute_resp_admission_4,
         -acute_resp_admission_5)

z_flu_tests <- readRDS("/conf/EAVE/GPanalysis/data/Ecoss_Flu_Apr2023.RDS") 

z_flu_tests = z_flu_tests %>% mutate(date_ecoss_specimen = as.Date(datespec),
                                     pos = if_else(denom == 1 & neg == 0, 1, 0)) %>%
  filter(pos == 1) %>%
  select(EAVE_LINKNO, date_ecoss_specimen)

z_covid_positive_tests = cdw_full %>% select(EAVE_LINKNO, date_ecoss_specimen, test_result) %>%
  filter(test_result == "POSITIVE") %>%
  filter(date_ecoss_specimen > a_begin & date_ecoss_specimen < a_end)

# Now tag anyone who was admitted with a positive covid or flu test
z_hospitalised_covid = z_hospitalised %>% 
  mutate(id = row_number()) %>%
  left_join(z_covid_positive_tests, by="EAVE_LINKNO") %>%
  mutate(test_diff = ADMISSION_DATE - date_ecoss_specimen) %>%
  filter(!is.na(test_diff) & test_diff > -2 & test_diff < 14) %>%
  arrange(id, test_diff) %>%
  filter(!duplicated(id)) %>% # remove any duplicate rows as we only really need one
  mutate(covid_admit_with_test = 1) %>%
  select(id, covid_admit_with_test)

z_hospitalised = z_hospitalised %>%
  mutate(id = row_number()) %>%
  left_join(z_hospitalised_covid, by="id") %>%
  select(-id) %>%
  mutate(covid_admit_with_test = if_else(is.na(covid_admit_with_test), 0, 1))

z_hospitalised_flu = z_hospitalised %>% 
  mutate(id = row_number()) %>%
  left_join(z_flu_tests, by="EAVE_LINKNO") %>%
  mutate(test_diff = ADMISSION_DATE - date_ecoss_specimen) %>%
  filter(!is.na(test_diff) & test_diff > -2 & test_diff < 14) %>%
  arrange(id, test_diff) %>%
  filter(!duplicated(id)) %>% # remove any duplicate rows as we only really need one
  mutate(flu_admit_with_test = 1) %>%
  select(id, flu_admit_with_test)

z_hospitalised = z_hospitalised %>%
  mutate(id = row_number()) %>%
  left_join(z_hospitalised_flu, by="id") %>%
  select(-id) %>%
  mutate(flu_admit_with_test = if_else(is.na(flu_admit_with_test), 0, 1))

# Remove anyone who didn't have a respiratory condition as main or secondary cause of admission
z_hospitalised = z_hospitalised %>% filter(acute_resp_admission == 1 | 
                                           acute_resp_secondary == 1 | 
                                           covid_admit_with_test == 1 |
                                           flu_admit_with_test == 1)

print("Number of rows containing relevant conditions")
print(nrow(z_hospitalised))

## Now set up the death stuff
all_deaths = readRDS(paste0(Location, "EAVE/GPanalysis/data/all_deaths.rds"))

# Look at primary cause of deaths
z_deaths = all_deaths %>% filter(substr(UNDERLYING_CAUSE_OF_DEATH, 0, 3) %in% three_char_admittance_codes |
                                 substr(UNDERLYING_CAUSE_OF_DEATH, 0, 4) %in% four_char_admittance_codes) %>%
  select(EAVE_LINKNO, NRS.Date.Death, UNDERLYING_CAUSE_OF_DEATH) %>%
  rename(dod_cause = NRS.Date.Death, cause_of_death = UNDERLYING_CAUSE_OF_DEATH) %>% mutate(primary_resp_cod = 1)

# Look at secondary cause of deaths - we look in the 10 secondary death cause fields
z_deaths_secondary = all_deaths %>% filter(
                                   substr(CAUSE_OF_DEATH_CODE_0, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_0, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_1, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_1, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_2, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_2, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_3, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_3, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_4, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_4, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_5, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_5, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_6, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_6, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_7, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_7, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_8, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_8, 0, 4) %in% four_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_9, 0, 3) %in% three_char_admittance_codes |
                                   substr(CAUSE_OF_DEATH_CODE_9, 0, 4) %in% four_char_admittance_codes ) %>%
  select(EAVE_LINKNO, NRS.Date.Death, UNDERLYING_CAUSE_OF_DEATH,
         CAUSE_OF_DEATH_CODE_0, CAUSE_OF_DEATH_CODE_1, CAUSE_OF_DEATH_CODE_2,
         CAUSE_OF_DEATH_CODE_3, CAUSE_OF_DEATH_CODE_4, CAUSE_OF_DEATH_CODE_5,
         CAUSE_OF_DEATH_CODE_6, CAUSE_OF_DEATH_CODE_7, CAUSE_OF_DEATH_CODE_8,
         CAUSE_OF_DEATH_CODE_9) %>%
  rename(dod_secondary_cause = NRS.Date.Death) %>% mutate(secondary_resp_cod = 1)

# Death cohort becomes our primary dataframe for our interest in whod died in the community
death_cohort = full_join(z_deaths, 
                         z_deaths_secondary %>% select(EAVE_LINKNO, UNDERLYING_CAUSE_OF_DEATH, secondary_resp_cod), 
                         by="EAVE_LINKNO") %>%
  mutate(primary_resp_cod = if_else(is.na(primary_resp_cod), 0, 1),
         secondary_resp_cod = if_else(is.na(secondary_resp_cod), 0, 1),
         cause_of_death = if_else(is.na(cause_of_death), UNDERLYING_CAUSE_OF_DEATH, cause_of_death)) %>%
  filter(dod_cause > a_begin & dod_cause < a_end)

# Tag the specific pathogen
z_hospitalised = z_hospitalised %>% mutate(flu_admit = if_else(
  substr(MAIN_CONDITION, 0, 3) %in% flu_codes, 1, 0))

z_hospitalised = z_hospitalised %>% mutate(rsv_admit = if_else(
  substr(MAIN_CONDITION, 0, 4) %in% rsv_codes, 1, 0))

z_hospitalised = z_hospitalised %>% mutate(covid_admit = if_else(
  substr(MAIN_CONDITION, 0, 3) %in% covid_three_char_codes |
  substr(MAIN_CONDITION, 0, 4) %in% covid_four_char_codes, 1, 0))

z_hospitalised = z_hospitalised %>% mutate(urti_admit = if_else(
  substr(MAIN_CONDITION, 0, 3) %in% urti_codes, 1, 0))

z_hospitalised = z_hospitalised %>% mutate(bronch_admit = if_else(
  substr(MAIN_CONDITION, 0, 3) %in% bronch_codes, 1, 0))

z_hospitalised = z_hospitalised %>% mutate(lrti_admit = if_else(
  substr(MAIN_CONDITION, 0, 3) %in% unspecified_codes, 1, 0))

#z_hospitalised = z_hospitalised %>% mutate(strep_admit = if_else(
#    substr(MAIN_CONDITION, 0, 4) %in% strep_codes, 1, 0))

hosp_cohort = z_hospitalised %>% left_join(EAVE_cohort, by="EAVE_LINKNO") #%>%
  #filter(!is.na(ageYear) & !is.na(ur6_2016)) 

print("Number of rows contained within the EAVE cohort")
print(nrow(hosp_cohort))

hosp_cohort = hosp_cohort %>% # Remove anyone who didn't match
  filter(ADMISSION_TYPE >= 30 & ADMISSION_TYPE < 40) # keep only emergency admissions

print("Number of non emergency admissions")
print(nrow(hosp_cohort))

# Now handle the secondary cause of admissions
z_urti_secondary = check_for_condition_in_secondary_codes(hosp_cohort, urti_codes, c())
z_urti_secondary = z_urti_secondary %>% rename(urti_admit_secondary = secondary)
z_bronch_secondary = check_for_condition_in_secondary_codes(hosp_cohort, bronch_codes, c())
z_bronch_secondary = z_bronch_secondary %>% rename(bronch_admit_secondary = secondary)
z_lrti_unspecified_secondary = check_for_condition_in_secondary_codes(hosp_cohort, unspecified_codes, c())
z_lrti_unspecified_secondary = z_lrti_unspecified_secondary %>% rename(lrti_admit_secondary = secondary)
z_flu_secondary = check_for_condition_in_secondary_codes(hosp_cohort, flu_codes, c())
z_flu_secondary = z_flu_secondary %>% rename(flu_admit_secondary = secondary)
z_rsv_secondary = check_for_condition_in_secondary_codes(hosp_cohort, rsv_codes, c())
z_rsv_secondary = z_rsv_secondary %>% rename(rsv_admit_secondary = secondary)

z_covid_secondary = check_for_condition_in_secondary_codes(hosp_cohort, covid_three_char_codes,
                                                           covid_four_char_codes)
z_covid_secondary = z_covid_secondary %>% rename(covid_admit_secondary = secondary)


#z_strep_secondary = check_for_condition_in_secondary_codes(z_hospitalised, c(),
#                                                           strep_codes)
#z_strep_secondary = z_strep_secondary %>% rename(strep_admit_secondary = secondary)

hosp_cohort = hosp_cohort %>% left_join(z_flu_secondary, by=c("EAVE_LINKNO", "CIS_MARKER"))
hosp_cohort = hosp_cohort %>% left_join(z_rsv_secondary, by=c("EAVE_LINKNO", "CIS_MARKER"))
hosp_cohort = hosp_cohort %>% left_join(z_covid_secondary, by=c("EAVE_LINKNO", "CIS_MARKER"))
hosp_cohort = hosp_cohort %>% left_join(z_urti_secondary, by=c("EAVE_LINKNO", "CIS_MARKER"))
hosp_cohort = hosp_cohort %>% left_join(z_bronch_secondary, by=c("EAVE_LINKNO", "CIS_MARKER"))
hosp_cohort = hosp_cohort %>% left_join(z_lrti_unspecified_secondary, by=c("EAVE_LINKNO", "CIS_MARKER"))

#hosp_cohort = hosp_cohort %>% left_join(z_strep_secondary, by=c("EAVE_LINKNO", "CIS_MARKER"))

hosp_cohort = hosp_cohort %>% left_join(Vaccinations, by="EAVE_LINKNO")

#risk groups
rg <- readRDS( "/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid_all.rds")
rg <- filter(rg,!duplicated(EAVE_LINKNO))
rg = rg %>% mutate(rg_bmi = bmi_impute) # Store the original bmi values so we know how many are missing
# These QCovid risk groups are not relevant to CYP
z_vars_to_remove = c("Q_BMI", "Q_DIAG_COPD", "Q_DIAG_CHD", "Q_DIAG_DEMENTIA", "Q_DIAG_PARKINSONS")

#individuals with no values in rg have no risk conditions
z <- hosp_cohort %>% 
  left_join(dplyr::select(rg,-(Sex:simd2020_sc_quintile), -DataZone, -ur6_2016_name) , by=c("EAVE_LINKNO_old" = "EAVE_LINKNO"))
z <- z %>% mutate_at(vars(Q_DIAG_AF:Q_DIAG_CKD_LEVEL), ~replace(., is.na(.), 0))
z <- z %>% mutate_at(vars(Q_DIAG_AF:Q_DIAG_CKD_LEVEL), ~as.numeric(.))
z <- z %>% mutate(n_risk_gps_old = fct_explicit_na(n_risk_gps, na_level="0"))
hosp_cohort <- z

hosp_cohort <- hosp_cohort %>% dplyr::select(-bmi_impute)

z_cyp_cohort = hosp_cohort %>% filter(ageYear < 18) %>% 
  mutate_at(z_vars_to_remove, function(x, na.rm = FALSE) ( 0 ) )

z_adult_cohort = hosp_cohort %>% filter(ageYear > 17)

hosp_cohort = rbind(z_adult_cohort, z_cyp_cohort)

qcovid_col_idx = grepl("Q_", colnames(hosp_cohort)) & !grepl("Q_BMI", colnames(hosp_cohort))

# Recalculate the number of risk groups
hosp_cohort$n_risk_gps = rowSums(hosp_cohort[,qcovid_col_idx] != 0)
hosp_cohort$n_risk_gps = cut(hosp_cohort$n_risk_gps, breaks=c(-1, 0, 1, 2, 3, 4, max(hosp_cohort$n_risk_gps)))

hosp_cohort$n_risk_gps = factor(hosp_cohort$n_risk_gps, levels(hosp_cohort$n_risk_gps), 
                                labels=c("0", "1", "2", "3", "4", "5+"))
# Produce a summary of the QCovid conditions

print("Entire Cohort")
print(colSums(hosp_cohort[,qcovid_col_idx]))

qcovid_col_idx = grep("Q_", colnames(z_cyp_cohort))
print("CYP")
print(colSums(z_cyp_cohort[,qcovid_col_idx] != 0))

hosp_cohort = hosp_cohort %>% mutate(age_gp = cut(ageYear, breaks=c(-1, 0, 1, 2, 5, 17, seq(24,80, by=5),max(ageYear))))
hosp_cohort$age_gp = factor(hosp_cohort$age_gp, levels(hosp_cohort$age_gp), labels = c("0", "1", "2", "3-5", "6-17", "18-24", "25-29",
                                                           "30-34","35-39","40-44","45-49","50-54","55-59",
                                                           "60-64","65-69","70-74","75-79","80+"))

# Sort out the covid vaccine status
hosp_cohort = hosp_cohort %>%
  mutate(vacc_1_diff = ADMISSION_DATE - date_vacc_1,
         vacc_2_diff = ADMISSION_DATE - date_vacc_2,
         vacc_3_diff = ADMISSION_DATE - date_vacc_3,
         vacc_4_diff = ADMISSION_DATE - date_vacc_4,
         vacc_5_diff = ADMISSION_DATE - date_vacc_5) %>%
  mutate(covid_vs = case_when(
                        !is.na(vacc_5_diff) & vacc_5_diff > 14 ~ "v5_2+",
                        #!is.na(vacc_5_diff) & vacc_5_diff > 0 ~ "v5_0:2",
                        !is.na(vacc_4_diff) & vacc_4_diff > 14 ~ "v4_2+",
                        #!is.na(vacc_4_diff) & vacc_4_diff > 0 ~ "v4_0:2",
                        !is.na(vacc_3_diff) & vacc_3_diff > 14 ~ "v3_2+",
                        #!is.na(vacc_3_diff) & vacc_3_diff > 0 ~ "v3_0:2",
                        !is.na(vacc_2_diff) & vacc_2_diff > 14 ~ "v2_2+",
                        #!is.na(vacc_2_diff) & vacc_2_diff > 0 ~ "v2_0:2",
                        !is.na(vacc_1_diff) & vacc_1_diff > 14 ~ "v1_2+",
                        #!is.na(vacc_1_diff) & vacc_1_diff > 0 ~ "v1_0:2",
                        is.na(vacc_1_diff) | vacc_1_diff < 15 ~ "uv"))

hosp_cohort$covid_vs = factor(hosp_cohort$covid_vs,  labels = c("Unvaccinated", 
                                                                #"1st Dose 0 - 14 days",
                                                                "1st Dose 14+ days",
                                                                "2nd Dose 14+ days",
                                                                #"3rd Dose 0 - 14 days",
                                                                "3rd Dose 14+ days",
                                                                #"4th Dose 0 - 14 days",
                                                                "4th Dose 14+ days",
                                                                "5th Dose 14+ days"))

#all_hospitalisations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/automated_any_hospitalisation_post_01022020.rds"))

z_flu_vacc = readRDS(paste0(Location, "EAVE/GPanalysis/data/cleaned_data/FLUvaccine_dvprod.rds")) %>%
  mutate(occurrence_time = as.Date(occurrence_time))

# Only take the 2022 - 2023 flu season
z_flu_vacc = z_flu_vacc %>% filter(occurrence_time > as.Date("2022-09-01")) %>%
  select(EAVE_LINKNO, occurrence_time, dose_number) %>%
  filter(dose_number == 1) %>%
  arrange(EAVE_LINKNO, occurrence_time) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  rename(date_flu_vacc_1 = occurrence_time) %>%
  select(-dose_number)


hosp_cohort = hosp_cohort %>% left_join(z_flu_vacc, by="EAVE_LINKNO")

hosp_cohort = hosp_cohort %>%
  mutate(vacc_1_diff = ADMISSION_DATE - date_flu_vacc_1) %>%
  mutate(flu_vs = case_when(!is.na(vacc_1_diff) & vacc_1_diff > 14 ~ "v1_2+",
                            !is.na(vacc_1_diff) & vacc_1_diff > 0 ~ "v1_0:2",
                            is.na(vacc_1_diff) | vacc_1_diff < 1 ~ "uv"))

hosp_cohort$flu_vs = factor(hosp_cohort$flu_vs, labels = c("Unvaccinated", "0 - 14 days", "14+ days"))

# Handle ICU data
icu_raw = readRDS("/conf/EAVE/GPanalysis/data/SICSAG_episode_level_.rds") %>%
  mutate(AdmitUnit = as.Date(AdmitUnit))

z_icu_hosp = hosp_cohort %>% left_join(icu_raw, by="EAVE_LINKNO") %>%
  mutate(row_num = row_number()) %>%
  filter(!is.na(id)) %>% # Remove anyone who doesn't have an ICU admission
  arrange(EAVE_LINKNO, ADMISSION_DATE) %>%
  mutate(admit_date_diff = abs(ADMISSION_DATE - AdmitHosp)) %>%
  filter(admit_date_diff <= 0) %>% # Check whether the ICU admission is close enough to the hospital admission
  filter(DiscDate <= DISCHARGE_DATE) %>% # Check this ICU admission corresponds to that hospitalisation
  mutate(icu_admit = 1) %>%
  select(row_num, icu_admit)

print(nrow(hosp_cohort))

# ICU admissions will correspond to multiple hospitalisations, so we make sure we don't increase
# the number of rows in the cohort by taking the closest ICU admittance
hosp_cohort = hosp_cohort %>% mutate(row_num = row_number()) %>%
  left_join(z_icu_hosp, by="row_num") %>%
  mutate(icu_admit = if_else(is.na(icu_admit), 0, icu_admit)) %>%
  select(-row_num)

print(nrow(hosp_cohort))

# Now set the difference between adults and CYP for ICU
hosp_cohort = hosp_cohort %>% mutate(icu_admit_age = case_when(
  ageYear > 17 & icu_admit == 1 ~ "Adult ICU Admission",
  ageYear < 18 & icu_admit == 1 ~ "CYP ICU Admission",
  TRUE ~ "No ICU Admission"
))

hosp_cohort = hosp_cohort %>% mutate(death_28 = (!is.na(NRS.Date.Death) & (as.numeric(NRS.Date.Death - ADMISSION_DATE) < 28)))

# Now set the difference between adults and CYP for death
hosp_cohort = hosp_cohort %>% mutate(death_28_age = case_when(
  ageYear > 17 & death_28 == 1 ~ "Adult Death",
  ageYear < 18 & death_28 == 1 ~ "CYP Death",
  TRUE ~ "Lived"
))


# Sort out of the deaths according to cause of death

# This sets up two binary variables - death_cause and death_hosp_cause
# death_cause indicates whether a person died of a respiratory disease 
# (that we are interested in)
# And death_hosp_cause indicates whether a person died of a respiratory disease
# that is the same as their main cause of admission
hosp_cohort = hosp_cohort %>% left_join(z_deaths, by="EAVE_LINKNO") %>%
  mutate(cause_of_death_substr = case_when(
    substr(cause_of_death, 0, 3) %in% three_char_admittance_codes ~ substr(cause_of_death, 0, 3),
    substr(cause_of_death, 0, 4) %in% four_char_admittance_codes ~ substr(cause_of_death, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(Main_diag_admit_substr = case_when(
    substr(MAIN_CONDITION, 0, 3) %in% three_char_admittance_codes ~ substr(MAIN_CONDITION, 0, 3),
    substr(MAIN_CONDITION, 0, 4) %in% four_char_admittance_codes ~ substr(MAIN_CONDITION, 0, 4),
    TRUE ~ ""
    
  )) %>%
  mutate(death_hosp_cause = if_else(!is.na(dod_cause) & Main_diag_admit_substr == cause_of_death_substr, 1, 0)) %>%
  mutate(death_cause = if_else(!is.na(dod_cause) & cause_of_death_substr != "", 1, 0)) %>%
  select(-cause_of_death_substr)

# Now we handle situations where the respiratory disease was mentioned as a secondary cause of death
hosp_cohort = hosp_cohort %>% left_join(z_deaths_secondary, by="EAVE_LINKNO") %>%
  mutate(cause_of_death_0_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_0, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_0, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_0, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_0, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(cause_of_death_1_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_1, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_1, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_1, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_1, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(cause_of_death_2_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_2, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_2, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_2, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_2, 0, 4),
    TRUE ~ ""
  )) %>%  mutate(cause_of_death_3_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_3, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_3, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_3, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_3, 0, 4),
    TRUE ~ ""
  )) %>%  mutate(cause_of_death_4_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_4, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_4, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_4, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_4, 0, 4),
    TRUE ~ ""
  )) %>%  mutate(cause_of_death_5_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_5, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_5, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_5, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_5, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(cause_of_death_6_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_6, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_6, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_6, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_6, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(cause_of_death_7_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_7, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_7, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_7, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_7, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(cause_of_death_8_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_8, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_8, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_8, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_8, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(cause_of_death_9_substr = case_when(
    substr(CAUSE_OF_DEATH_CODE_9, 0, 3) %in% three_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_9, 0, 3),
    substr(CAUSE_OF_DEATH_CODE_9, 0, 4) %in% four_char_admittance_codes ~ substr(CAUSE_OF_DEATH_CODE_9, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(Main_diag_admit_substr = case_when(
    substr(MAIN_CONDITION, 0, 3) %in% three_char_admittance_codes ~ substr(MAIN_CONDITION, 0, 3),
    substr(MAIN_CONDITION, 0, 4) %in% four_char_admittance_codes ~ substr(MAIN_CONDITION, 0, 4),
    TRUE ~ ""
  )) %>%
  mutate(death_hosp_secondary_cause = case_when
         (!is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_0_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_1_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_2_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_3_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_4_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_5_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_6_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_7_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_8_substr ~ 1,
          !is.na(dod_secondary_cause) & Main_diag_admit_substr == cause_of_death_9_substr ~ 1,
          TRUE ~ 0)) %>%
  mutate(death_secondary_cause = case_when
          (!is.na(dod_secondary_cause) & cause_of_death_0_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_1_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_2_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_3_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_4_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_5_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_6_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_7_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_8_substr != "" ~ 1,
          !is.na(dod_secondary_cause) & cause_of_death_9_substr != "" ~ 1,
           TRUE ~ 0)) %>%
  select(-cause_of_death_0_substr, -cause_of_death_1_substr, -cause_of_death_2_substr, 
         -cause_of_death_3_substr, -cause_of_death_4_substr, -cause_of_death_5_substr,
         -cause_of_death_6_substr, -cause_of_death_7_substr, -cause_of_death_8_substr,
         -cause_of_death_9_substr, -Main_diag_admit_substr)
  


hosp_cohort  = hosp_cohort %>%
  mutate(hosp_los = as.numeric(difftime(DISCHARGE_DATE, ADMISSION_DATE, units="days"))) %>%
  mutate(hosp_los_gp = cut(hosp_los, breaks=c(-1, 0, 1, 2, 5, 10, 20, max(hosp_los))))

hosp_cohort$hosp_los_gp = factor(hosp_cohort$hosp_los_gp, levels(hosp_cohort$hosp_los_gp), labels=c("0","1", "2", "3-5", "6-9", "10-19", "20+"))

# Add in ethnicity
ethnicity = readRDS("/conf/EAVE/GPanalysis/data/lookups/EAVE_Ethnicity_2022.rds") %>% 
  filter(!duplicated(EAVE_LINKNO))

ethnicity = ethnicity %>% mutate(ethnic_gp = case_when(
  substr(ethnic_code, 0, 1) == "1" ~ "White",
  substr(ethnic_code, 0, 1) == "2" ~ "Mixed",
  substr(ethnic_code, 0, 1) == "3" ~ "Asian",
  substr(ethnic_code, 0, 1) == "4" ~ "Black",
  substr(ethnic_code, 0, 1) == "5" ~ "Black",
  substr(ethnic_code, 0, 1) == "6" ~ "Other",
  TRUE ~ "Unknown"
  ))

hosp_cohort = hosp_cohort %>% left_join(ethnicity, by="EAVE_LINKNO") %>%
  mutate(ethnic_gp = if_else(is.na(ethnic_gp), "Unknown", ethnic_gp))

# Add in the number of emeregency admissions in the 6 months
# leading up to the admission date
z_emergency = smr01 %>% select(EAVE_LINKNO, ADMISSION_TYPE, ADMISSION_DATE) %>%
  filter(ADMISSION_TYPE >= 30 & ADMISSION_TYPE < 40) %>% # keep only emergency admissions
  filter(ADMISSION_DATE >= a_begin_minus_6_months & ADMISSION_DATE < a_begin) %>%
  mutate(prev_admit_date = ADMISSION_DATE) %>%
  select(-ADMISSION_TYPE) %>%
  group_by(EAVE_LINKNO) %>%
  summarise(n = n()) %>%
  rename(num_prev_admission = n)

hosp_cohort = hosp_cohort %>% left_join(z_emergency, by="EAVE_LINKNO") %>%
  mutate(num_prev_admission = ifelse(is.na(num_prev_admission), 0, num_prev_admission)) %>%
  mutate(num_prev_admission_gp = cut(num_prev_admission, breaks=c(-1, 0, 1, 2, 3, 4, 5, max(num_prev_admission))))

hosp_cohort$num_prev_admission_gp = factor(hosp_cohort$num_prev_admission_gp, levels(hosp_cohort$num_prev_admission_gp),
                                           labels=c("0", "1", "2", "3", "4", "5", "6+"))

# Sort out the urban/rural classification
hosp_cohort = hosp_cohort %>% mutate(urban_rural_classification = case_when(
  ur6_2016 == 1 | ur6_2016 == 2 | ur6_2016 == 3  ~ "Urban",
  ur6_2016 == 4 | ur6_2016 == 5 | ur6_2016 == 6 ~ "Rural",
  TRUE ~ "Unknown"
))

# Sort out the extended hospital stay variable
hosp_cohort = hosp_cohort %>% mutate(extended_los = if_else(hosp_los > 5, 1, 0))

# We know that everyone in this cohort exists so we can set their weight to 1
hosp_cohort = hosp_cohort %>% mutate(eave_weight = 1)

# Tidy up healthboard data
hosp_cohort = hosp_cohort %>% mutate(health_board = if_else(is.na(health_board), "Unknown", health_board))

# Add in smoking
z_smoking <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds") %>% 
  select(EAVE_LINKNO, EAVE_Smoking_Status_Best, EAVE_Smoking_Status_Worst) %>%
  filter(!duplicated(EAVE_LINKNO))

hosp_cohort = hosp_cohort %>% left_join(z_smoking, by=join_by(EAVE_LINKNO_old == EAVE_LINKNO))

#remove data sets not needed
rm(bnf, pis_asthma, cdw_full, all_deaths)
#rm(Vaccinations, smr01, rg, icu_raw, EAVE_demographics, all_hospitalisations,
#   EAVE_cohort_refresh, EAVE_endpoints)
remove(list=ls(pa="^z"))

##########################################
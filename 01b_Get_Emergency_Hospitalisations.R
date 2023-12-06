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

cardiac_codes = c("I")
cancer_codes = c("C")
trauma_codes = c("S", "T")

one_char_admittance_codes = c("C", "I", "S", "T")

print("Number of CIS rows in SMR01")
print(nrow(smr01))

# How many emergency hospitalisations during this time?
z_num_emergency_hosps = smr01 %>% filter(ADMISSION_DATE >= a_begin & ADMISSION_DATE <= a_end) %>%
  filter(ADMISSION_TYPE >= 30 & ADMISSION_TYPE <= 40)
print("Number of emergency hospitalisations during time period")
print(nrow(z_num_emergency_hosps))

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
  filter(as.numeric(difftime(DISCHARGE_DATE, ADMISSION_DATE, units="days")) > 0)

print("Number of rows within the time period")
print(nrow(z_hospitalised))

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

print("Number of rows containing relevant conditions")
print(nrow(z_hospitalised))

emergency_hosp_cohort = z_hospitalised %>% left_join(EAVE_cohort, by="EAVE_LINKNO") %>%
  filter(!is.na(ageYear)) 

print("Number of rows contained within the EAVE cohort")
print(nrow(emergency_hosp_cohort))

emergency_hosp_cohort = emergency_hosp_cohort %>% # Remove anyone who didn't match
  filter(ADMISSION_TYPE >= 30 & ADMISSION_TYPE < 40) %>% # keep only emergency admissions
  arrange(EAVE_LINKNO, ADMISSION_DATE) %>% 
  filter(!duplicated(EAVE_LINKNO)) # Keep only the first admission

print("Number of non emergency admissions")
print(nrow(emergency_hosp_cohort))

emergency_hosp_cohort = emergency_hosp_cohort %>% left_join(Vaccinations, by="EAVE_LINKNO")

# These QCovid risk groups are not relevant to CYP
z_vars_to_remove = c("Q_BMI", "Q_DIAG_COPD", "Q_DIAG_CHD", "Q_DIAG_DEMENTIA", "Q_DIAG_PARKINSONS")

#individuals with no values in rg have no risk conditions
z <- emergency_hosp_cohort %>% 
  left_join(dplyr::select(rg,-(Sex:simd2020_sc_quintile), -DataZone, -ur6_2016_name) , by=c("EAVE_LINKNO_old" = "EAVE_LINKNO"))
z <- z %>% mutate_at(vars(Q_DIAG_AF:Q_DIAG_CKD_LEVEL), ~replace(., is.na(.), 0))
z <- z %>% mutate_at(vars(Q_DIAG_AF:Q_DIAG_CKD_LEVEL), ~as.numeric(.))
z <- z %>% mutate(n_risk_gps_old = fct_explicit_na(n_risk_gps, na_level="0"))
emergency_hosp_cohort <- z

emergency_hosp_cohort <- emergency_hosp_cohort %>% dplyr::select(-bmi_impute)

z_cyp_cohort = emergency_hosp_cohort %>% filter(ageYear < 18) %>% 
  mutate_at(z_vars_to_remove, function(x, na.rm = FALSE) ( 0 ) )

z_adult_cohort = emergency_hosp_cohort %>% filter(ageYear > 17)

emergency_hosp_cohort = rbind(z_adult_cohort, z_cyp_cohort)

qcovid_col_idx = grepl("Q_", colnames(emergency_hosp_cohort)) & !grepl("Q_BMI", colnames(emergency_hosp_cohort))

# Recalculate the number of risk groups
emergency_hosp_cohort$n_risk_gps = rowSums(emergency_hosp_cohort[,qcovid_col_idx] != 0)
emergency_hosp_cohort$n_risk_gps = cut(emergency_hosp_cohort$n_risk_gps, breaks=c(-1, 0, 1, 2, 3, 4, max(emergency_hosp_cohort$n_risk_gps)))

emergency_hosp_cohort$n_risk_gps = factor(emergency_hosp_cohort$n_risk_gps, levels(emergency_hosp_cohort$n_risk_gps), 
                                labels=c("0", "1", "2", "3", "4", "5+"))
# Produce a summary of the QCovid conditions

print("Entire Cohort")
print(colSums(emergency_hosp_cohort[,qcovid_col_idx]))

qcovid_col_idx = grep("Q_", colnames(z_cyp_cohort))
print("CYP")
print(colSums(z_cyp_cohort[,qcovid_col_idx] != 0))

emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(age_gp = cut(ageYear, breaks=c(-1, 0, 1, 2, 5, 17, seq(24,80, by=5),max(ageYear))))
emergency_hosp_cohort$age_gp = factor(emergency_hosp_cohort$age_gp, levels(emergency_hosp_cohort$age_gp), labels = c("0", "1", "2", "3-5", "6-17", "18-24", "25-29",
                                                           "30-34","35-39","40-44","45-49","50-54","55-59",
                                                           "60-64","65-69","70-74","75-79","80+"))

# Sort out the covid vaccine status
emergency_hosp_cohort = emergency_hosp_cohort %>%
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
                        #!is.na(vacc_2_diff) & vacc_2_diff > 0 ~ "v3_0:2",
                        !is.na(vacc_1_diff) & vacc_1_diff > 14 ~ "v1_2+",
                        #!is.na(vacc_1_diff) & vacc_1_diff > 0 ~ "v1_0:2",
                        is.na(vacc_1_diff) | vacc_1_diff < 15 ~ "uv"))

emergency_hosp_cohort$covid_vs = factor(emergency_hosp_cohort$covid_vs,  labels = c("Unvaccinated", 
                                                                #"1st Dose 0 - 14 days",
                                                                "1st Dose 14+ days",
                                                                "2nd Dose 14+ days",
                                                                #"3rd Dose 0 - 14 days",
                                                                "3rd Dose 14+ days",
                                                                #"4th Dose 0 - 14 days",
                                                                "4th Dose 14+ days",
                                                                #"5th Dose 0 - 14 days",
                                                                "5th Dose 14+ days"))


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


emergency_hosp_cohort = emergency_hosp_cohort %>% left_join(z_flu_vacc, by="EAVE_LINKNO")

emergency_hosp_cohort = emergency_hosp_cohort %>%
  mutate(vacc_1_diff = ADMISSION_DATE - date_flu_vacc_1) %>%
  mutate(flu_vs = case_when(!is.na(vacc_1_diff) & vacc_1_diff > 14 ~ "v1_2+",
                            !is.na(vacc_1_diff) & vacc_1_diff > 0 ~ "v1_0:2",
                            is.na(vacc_1_diff) | vacc_1_diff < 1 ~ "uv"))

emergency_hosp_cohort$flu_vs = factor(emergency_hosp_cohort$flu_vs, labels = c("Unvaccinated", "0 - 14 days", "14+ days"))

# Handle ICU data
z_icu_hosp = emergency_hosp_cohort %>% left_join(icu_raw, by="EAVE_LINKNO") %>%
  mutate(row_num = row_number()) %>%
  filter(!is.na(id)) %>% # Remove anyone who doesn't have an ICU admission
  arrange(EAVE_LINKNO, ADMISSION_DATE) %>%
  mutate(admit_date_diff = abs(ADMISSION_DATE - AdmitHosp)) %>%
  filter(admit_date_diff <= 0) %>% # Check whether the ICU admission is close enough to the hospital admission
  filter(DiscDate <= DISCHARGE_DATE) %>% # Check this ICU admission corresponds to that hospitalisation
  mutate(icu_admit = 1) %>%
  select(row_num, icu_admit)

print(nrow(emergency_hosp_cohort))

# ICU admissions will correspond to multiple hospitalisations, so we make sure we don't increase
# the number of rows in the cohort by taking the closest ICU admittance
emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(row_num = row_number()) %>%
  left_join(z_icu_hosp, by="row_num") %>%
  mutate(icu_admit = if_else(is.na(icu_admit), 0, icu_admit)) %>%
  select(-row_num)

print(nrow(emergency_hosp_cohort))

# Now set the difference between adults and CYP for ICU
emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(icu_admit_age = case_when(
  ageYear > 17 & icu_admit == 1 ~ "Adult ICU Admission",
  ageYear < 18 & icu_admit == 1 ~ "CYP ICU Admission",
  TRUE ~ "No ICU Admission"
))

emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(death_28 = (!is.na(NRS.Date.Death) & (as.numeric(NRS.Date.Death - ADMISSION_DATE) < 28)))

# Now set the difference between adults and CYP for death
emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(death_28_age = case_when(
  ageYear > 17 & death_28 == 1 ~ "Adult Death",
  ageYear < 18 & death_28 == 1 ~ "CYP Death",
  TRUE ~ "Lived"
))


# Sort out of the deaths according to cause of death
emergency_hosp_cohort  = emergency_hosp_cohort %>%
  mutate(hosp_los = as.numeric(difftime(DISCHARGE_DATE, ADMISSION_DATE, units="days"))) %>%
  mutate(hosp_los_gp = cut(hosp_los, breaks=c(-1, 0, 1, 2, 5, 10, 20, max(hosp_los))))

emergency_hosp_cohort$hosp_los_gp = factor(emergency_hosp_cohort$hosp_los_gp, levels(emergency_hosp_cohort$hosp_los_gp), labels=c("0","1", "2", "3-5", "6-9", "10-19", "20+"))

# Add in ethnicity
ethnicity = readRDS("/conf/EAVE/GPanalysis/data/lookups/EAVE_Ethnicity_2022.rds")

ethnicity = ethnicity %>% mutate(ethnic_gp = case_when(
  substr(ethnic_code, 0, 1) == "1" ~ "White",
  substr(ethnic_code, 0, 1) == "2" ~ "Mixed",
  substr(ethnic_code, 0, 1) == "3" ~ "Asian",
  substr(ethnic_code, 0, 1) == "4" ~ "Black",
  substr(ethnic_code, 0, 1) == "5" ~ "Black",
  substr(ethnic_code, 0, 1) == "6" ~ "Other",
  TRUE ~ "Unknown"
  ))

emergency_hosp_cohort = emergency_hosp_cohort %>% left_join(ethnicity, by="EAVE_LINKNO") %>%
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

emergency_hosp_cohort = emergency_hosp_cohort %>% left_join(z_emergency, by="EAVE_LINKNO") %>%
  mutate(num_prev_admission = ifelse(is.na(num_prev_admission), 0, num_prev_admission)) %>%
  mutate(num_prev_admission_gp = cut(num_prev_admission, breaks=c(-1, 0, 1, 2, 3, 4, 5, max(num_prev_admission))))

emergency_hosp_cohort$num_prev_admission_gp = factor(emergency_hosp_cohort$num_prev_admission_gp, levels(emergency_hosp_cohort$num_prev_admission_gp),
                                           labels=c("0", "1", "2", "3", "4", "5", "6+"))

# Sort out the urban/rural classification
emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(urban_rural_classification = case_when(
  ur6_2016 == 1 | ur6_2016 == 2 | ur6_2016 == 3  ~ "Urban",
  ur6_2016 == 4 | ur6_2016 == 5 | ur6_2016 == 6 ~ "Rural",
  TRUE ~ "Unknown"
))

# Sort out the extended hospital stay variable
emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(extended_los = if_else(hosp_los > 5, 1, 0))

# We know that everyone in this cohort exists so we can set their weight to 1
emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(eave_weight = 1)

# Tidy up healthboard data
emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(health_board = if_else(is.na(health_board), "Unknown", health_board))

# Handle causes of admission
# Tag people who were admitted for certain things
emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(acute_resp_admission = if_else(
  substr(MAIN_CONDITION, 0, 3) %in% three_char_admittance_codes |
    substr(MAIN_CONDITION, 0, 4) %in% four_char_admittance_codes, 1, 0)) %>%
  mutate(acute_trauma_admission = if_else(
    substr(MAIN_CONDITION, 0, 1) %in% trauma_codes, 1, 0
  ), 
  acute_cardiac_admission = if_else(
    substr(MAIN_CONDITION, 0, 1) %in% cardiac_codes, 1, 0
  ),
  acute_cancer_admission = if_else(
    substr(MAIN_CONDITION, 0, 1) %in% cancer_codes, 1, 0
  ))

emergency_hosp_cohort = emergency_hosp_cohort %>% mutate(
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
         -acute_resp_admission_5) %>% 
  mutate(
           acute_cancer_admission_1 = if_else(
             substr(OTHER_CONDITION_1, 0, 1) %in% cancer_codes, 1, 0),
           acute_cancer_admission_2 = if_else(
             substr(OTHER_CONDITION_2, 0, 1) %in% cancer_codes, 1, 0),
           acute_cancer_admission_3 = if_else(
             substr(OTHER_CONDITION_3, 0, 1) %in% cancer_codes, 1, 0),
           acute_cancer_admission_4 = if_else(
             substr(OTHER_CONDITION_4, 0, 1) %in% cancer_codes, 1, 0),
           acute_cancer_admission_5 = if_else(
             substr(OTHER_CONDITION_5, 0, 1) %in% cancer_codes, 1, 0)) %>%
  mutate(acute_cancer_secondary = if_else(
    acute_cancer_admission_1 == 1 | acute_cancer_admission_2 == 1 |
      acute_cancer_admission_3 == 1 | acute_cancer_admission_4 == 1 | acute_cancer_admission_5 == 1, 1, 0)
  ) %>%
  select(-acute_cancer_admission_1, -acute_cancer_admission_2, -acute_cancer_admission_3, -acute_cancer_admission_4,
         -acute_cancer_admission_5) %>% 
  mutate(
    acute_trauma_admission_1 = if_else(
      substr(OTHER_CONDITION_1, 0, 1) %in% trauma_codes, 1, 0),
    acute_trauma_admission_2 = if_else(
      substr(OTHER_CONDITION_2, 0, 1) %in% trauma_codes, 1, 0),
    acute_trauma_admission_3 = if_else(
      substr(OTHER_CONDITION_3, 0, 1) %in% trauma_codes, 1, 0),
    acute_trauma_admission_4 = if_else(
      substr(OTHER_CONDITION_4, 0, 1) %in% trauma_codes, 1, 0),
    acute_trauma_admission_5 = if_else(
      substr(OTHER_CONDITION_5, 0, 1) %in% trauma_codes, 1, 0)) %>%
  mutate(acute_trauma_secondary = if_else(
    acute_trauma_admission_1 == 1 | acute_trauma_admission_2 == 1 |
      acute_trauma_admission_3 == 1 | acute_trauma_admission_4 == 1 | acute_trauma_admission_5 == 1, 1, 0)
  ) %>%
  select(-acute_trauma_admission_1, -acute_trauma_admission_2, -acute_trauma_admission_3, -acute_trauma_admission_4,
         -acute_trauma_admission_5) %>% 
  mutate(
    acute_cardiac_admission_1 = if_else(
      substr(OTHER_CONDITION_1, 0, 1) %in% cardiac_codes, 1, 0),
    acute_cardiac_admission_2 = if_else(
      substr(OTHER_CONDITION_2, 0, 1) %in% cardiac_codes, 1, 0),
    acute_cardiac_admission_3 = if_else(
      substr(OTHER_CONDITION_3, 0, 1) %in% cardiac_codes, 1, 0),
    acute_cardiac_admission_4 = if_else(
      substr(OTHER_CONDITION_4, 0, 1) %in% cardiac_codes, 1, 0),
    acute_cardiac_admission_5 = if_else(
      substr(OTHER_CONDITION_5, 0, 1) %in% cardiac_codes, 1, 0)) %>%
  mutate(acute_cardiac_secondary = if_else(
    acute_cardiac_admission_1 == 1 | acute_cardiac_admission_2 == 1 |
      acute_cardiac_admission_3 == 1 | acute_cardiac_admission_4 == 1 | acute_cardiac_admission_5 == 1, 1, 0)
  ) %>%
  select(-acute_cardiac_admission_1, -acute_cardiac_admission_2, -acute_cardiac_admission_3, -acute_cardiac_admission_4,
         -acute_cardiac_admission_5) 



#remove data sets not needed
rm(bnf, pis_asthma, cdw_full, all_deaths)
#rm(Vaccinations, smr01, rg, icu_raw, EAVE_demographics, all_hospitalisations,
#   EAVE_cohort_refresh, EAVE_endpoints)
remove(list=ls(pa="^z"))

##########################################
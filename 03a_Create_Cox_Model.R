library(survival)
library(ggsurvfit)
library(tidycmprsk)
library(forcats)
library(lubridate)

set.seed(54)

# Which definition to use - set this to TRUE to include secondary conditions
broad_defintion = 1

# Set this to 1 to use a covid test as a definition of a covid hospitalisation
# This will override the other definiton above for covid specificially!
use_tests = 0

# Whether to use vaccination as an adjustment
use_vaccination = 0 

# Use smoking as an adjustment
use_smoking = 0

if (broad_defintion == 1) {
  df_events = hosp_cohort %>% filter(acute_resp_admission == 1 | acute_resp_secondary == 1)
} else if (use_tests == 1) {
  df_events = hosp_cohort %>% filter(acute_resp_admission == 1 | flu_admit_with_test == 1 | covid_admit_with_test == 1)
} else {
  df_events = hosp_cohort %>% filter(acute_resp_admission == 1)
}

df_events = df_events %>%
  select(EAVE_LINKNO, EAVE_LINKNO_old, ageYear, age_gp, Sex, ur6_2016, simd2020_sc_quintile, ethnic_gp, ADMISSION_DATE,
         flu_admit, covid_admit, rsv_admit, flu_admit_secondary, covid_admit_secondary, rsv_admit_secondary, 
         covid_admit_with_test, flu_admit_with_test, num_prev_admission_gp, eave_weight, health_board,
         date_vacc_1, date_vacc_2, date_vacc_3, date_vacc_4, date_vacc_5, date_flu_vacc_1,
         icu_admit_age, hosp_los, death_28, EAVE_Smoking_Status_Worst, EAVE_Smoking_Status_Best, NRS.Date.Death) %>%
  rename(event_date = ADMISSION_DATE) %>%
  mutate(event = 1)

print("Number of hospitalisations at start of cox model")
print(nrow(df_events))
df_events = df_events %>%
  arrange(EAVE_LINKNO, event_date) %>% 
  filter(!duplicated(EAVE_LINKNO)) # We only want the first admission to hospital

print("Number of first hospitalisations")
print(nrow(df_events))

n_events = nrow(df_events)
n_samples_per_event = 10
n_controls = n_events * n_samples_per_event

# Find out who has an event so we don't select them as a control
# THIS MIGHT LEAD TO BIAS - CHECK WITH CHRIS
event_LINKNOs = df_events %>% select(EAVE_LINKNO) %>% unique() %>% pull()

z_eligible = EAVE_cohort %>% filter(!(EAVE_LINKNO %in% event_LINKNOs))
print("Number of eligible controls")
print(sum(z_eligible$eave_weight))

# Select our controls
df_controls = EAVE_cohort %>%
  select(EAVE_LINKNO, EAVE_LINKNO_old, ageYear, Sex, ur6_2016, simd2020_sc_quintile, NRS.Date.Death, eave_weight, health_board) %>%
  filter(!(EAVE_LINKNO %in% event_LINKNOs)) %>%
  mutate(event = 0, event_date = a_end, flu_admit = 0, covid_admit = 0, rsv_admit = 0,
         flu_admit_secondary = 0, covid_admit_secondary = 0, rsv_admit_secondary = 0,
         covid_admit_with_test = 0, flu_admit_with_test = 0)  %>%
  mutate(event_date = if_else(is.na(NRS.Date.Death), a_end, NRS.Date.Death)) %>% # censor anyone who died
  mutate(age_gp = cut(ageYear, breaks=c(-1, 0, 1, 2, 5, 17, seq(24,80, by=5),max(ageYear))))

z_ethnicity = ethnicity %>% select(EAVE_LINKNO, ethnic_gp)

# Add in ethnicity
df_controls = df_controls %>% left_join(z_ethnicity, by="EAVE_LINKNO") %>%
  mutate(ethnic_gp = if_else(is.na(ethnic_gp), "Unknown", ethnic_gp))

# Previous emergency hospital admissions
z_emergency = smr01 %>% select(EAVE_LINKNO, ADMISSION_TYPE, ADMISSION_DATE) %>%
  filter(ADMISSION_TYPE >= 30 & ADMISSION_TYPE < 40) %>% # keep only emergency admissions
  filter(ADMISSION_DATE >= a_begin_minus_6_months & ADMISSION_DATE < a_begin) %>%
  mutate(prev_admit_date = ADMISSION_DATE) %>%
  select(-ADMISSION_TYPE) %>%
  group_by(EAVE_LINKNO) %>%
  summarise(n = n()) %>%
  rename(num_prev_admission = n)

df_controls = df_controls %>% left_join(z_emergency, by="EAVE_LINKNO") %>%
  mutate(num_prev_admission = ifelse(is.na(num_prev_admission), 0, num_prev_admission)) %>%
  mutate(num_prev_admission_gp = cut(num_prev_admission, breaks=c(-1, 0, 1, 2, 3, 4, 5, 21))) %>%
  select(-num_prev_admission)

z_vacc = Vaccinations %>% select(EAVE_LINKNO, date_vacc_1, date_vacc_2, date_vacc_3, date_vacc_4, date_vacc_5)

# Add in Vaccination statuses
df_controls = df_controls %>% left_join(z_vacc, by="EAVE_LINKNO") 

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



df_controls = df_controls %>% left_join(z_flu_vacc, by="EAVE_LINKNO")

# Set up the factor levels to match those in the hosp cohort
df_controls$num_prev_admission_gp = 
  factor(df_controls$num_prev_admission_gp, levels(df_controls$num_prev_admission_gp),
                                           labels=c("0", "1", "2", "3", "4", "5", "6+"))

df_controls$age_gp = 
  factor(df_controls$age_gp, 
         levels(df_controls$age_gp), 
         labels = c("0", "1", "2", "3-5", "6-17", "18-24", "25-29",
                    "30-34","35-39","40-44","45-49","50-54","55-59",
                    "60-64","65-69","70-74","75-79","80+"))

# Add in smoking
z_smoking <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds") %>% 
  select(EAVE_LINKNO, EAVE_Smoking_Status_Best, EAVE_Smoking_Status_Worst) %>%
  filter(!duplicated(EAVE_LINKNO))

df_controls = df_controls %>% left_join(z_smoking, by=join_by(EAVE_LINKNO_old == EAVE_LINKNO))

# Controls can have an ICU admission as some of them might have been admitted to hospital with
# something else
z_icu_hosp = df_controls %>% left_join(icu_raw, by="EAVE_LINKNO") %>%
  filter(!is.na(id)) %>% # Remove anyone who doesn't have an ICU admission
  filter(!duplicated(EAVE_LINKNO)) %>% # Only keep one row per person
  filter(AdmitUnit > a_begin & AdmitUnit < a_end) %>% # Ensure the ICU admission is within our period of interest
  mutate(icu_admit = 1)  %>% 
  select(EAVE_LINKNO, icu_admit) 

df_controls = df_controls %>%
  left_join(z_icu_hosp, by="EAVE_LINKNO") %>%
  mutate(icu_admit_age = case_when(
    ageYear > 17 & icu_admit == 1 ~ "Adult ICU Admission",
    ageYear < 18 & icu_admit == 1 ~ "CYP ICU Admission",
    TRUE ~ "No ICU Admission"
  )) %>%
  select(-icu_admit)

# This is a hack to display the number of deaths in the control cohort
# Our definition of death for the cases is different (28 days after a hospitalisation)
# But this doesn't make sense in the control context
df_controls = df_controls %>%
  mutate(death_28 = !is.na(NRS.Date.Death))


## Add in hospital length of stay for the first hospitalisation for the controls
## Some of them will obviously not have any stays


z_smr01 = smr01 %>% 
  filter(ADMISSION_DATE >= a_begin & ADMISSION_DATE <= a_end) %>%
  filter(ADMISSION_TYPE >= 30 & ADMISSION_TYPE <= 40) %>%
  mutate(hosp_los = as.numeric(difftime(DISCHARGE_DATE, ADMISSION_DATE, units = "days"))) %>% # We're only interested in stays of over a day
  filter(hosp_los > 0) %>%
  arrange(EAVE_LINKNO, ADMISSION_DATE) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  select(EAVE_LINKNO, hosp_los)
df_controls = df_controls %>% left_join(z_smr01, by="EAVE_LINKNO") %>%
  mutate(hosp_los = if_else(is.na(hosp_los), 0, hosp_los))

z_controls = df_controls %>%
  slice_sample(n = n_controls)

print("Number of selected controls")
print(sum(z_controls$eave_weight))

df_all = rbind(df_events, z_controls)
print("Number of combined controls and events")
print(sum(df_all$eave_weight))

# Set up the risk groups
rg <- readRDS( "/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid_all.rds")
rg <- filter(rg,!duplicated(EAVE_LINKNO))

#individuals with no values in rg have no risk conditions
z <- df_all %>% 
  left_join(dplyr::select(rg,-(Sex:simd2020_sc_quintile), -DataZone, -ur6_2016_name) , by=c("EAVE_LINKNO_old" = "EAVE_LINKNO"))
z <- z %>% mutate_at(vars(Q_DIAG_AF:Q_DIAG_CKD_LEVEL), ~replace(., is.na(.), 0))
z <- z %>% mutate_at(vars(Q_DIAG_AF:Q_DIAG_CKD_LEVEL), ~as.numeric(.))
z <- z %>% mutate(n_risk_gps = fct_explicit_na(n_risk_gps, na_level="0"))

# Remove the QCovid risk groups that don't make sense for cyp
z_vars_to_remove = c("Q_BMI", "Q_DIAG_COPD", "Q_DIAG_CHD", "Q_DIAG_DEMENTIA", "Q_DIAG_PARKINSONS")

z_cyp_cohort = z %>% filter(ageYear < 18) %>% 
  mutate_at(z_vars_to_remove, function(x, na.rm = FALSE) ( 0 ) )

z_adult_cohort = z %>% filter(ageYear > 17)

z = rbind(z_adult_cohort, z_cyp_cohort)

df_all <- z

z <- df_controls %>% 
  left_join(dplyr::select(rg,-(Sex:simd2020_sc_quintile), -DataZone, -ur6_2016_name) , by=c("EAVE_LINKNO_old" = "EAVE_LINKNO"))
z <- z %>% mutate_at(vars(Q_DIAG_AF:Q_DIAG_CKD_LEVEL), ~replace(., is.na(.), 0))
z <- z %>% mutate_at(vars(Q_DIAG_AF:Q_DIAG_CKD_LEVEL), ~as.numeric(.))
z <- z %>% mutate(n_risk_gps = fct_explicit_na(n_risk_gps, na_level="0"))

# Remove the QCovid risk groups that don't make sense for cyp
z_vars_to_remove = c("Q_BMI", "Q_DIAG_COPD", "Q_DIAG_CHD", "Q_DIAG_DEMENTIA", "Q_DIAG_PARKINSONS")

z_cyp_cohort = z %>% filter(ageYear < 18) %>% 
  mutate_at(z_vars_to_remove, function(x, na.rm = FALSE) ( 0 ) )

z_adult_cohort = z %>% filter(ageYear > 17)

df_controls = rbind(z_adult_cohort, z_cyp_cohort)

df_all = df_all %>% rename(bmi = bmi_impute, simd = simd2020_sc_quintile)
df_controls = df_controls %>% rename(bmi = bmi_impute, simd = simd2020_sc_quintile)

# Use Steven's BMI imputation
z_bmi_imp = readRDS("/conf/EAVE/GPanalysis/analyses/imputation/data/df_imp.rds") %>%
  select(EAVE_LINKNO, Q_BMI) %>%
  rename(bmi_imp = Q_BMI)

df_all = df_all %>% left_join(z_bmi_imp, by="EAVE_LINKNO") 

df_all = df_all %>% mutate(bmi_gp = cut(bmi_imp, breaks=c(0, 18.5, 25, 30, 35, 40, max(bmi_imp, na.rm=TRUE))))


df_controls = df_controls %>% left_join(z_bmi_imp, by="EAVE_LINKNO") 

df_controls = df_controls %>% mutate(bmi_gp = cut(bmi_imp, breaks=c(0, 18.5, 25, 30, 35, 40, max(bmi_imp, na.rm=TRUE))))

# Sort out the urban/rural classification
df_all = df_all %>% mutate(urban_rural_classification = case_when(
  ur6_2016 == 1 | ur6_2016 == 2 | ur6_2016 == 3 ~ "Urban",
  ur6_2016 == 4 | ur6_2016 == 5 | ur6_2016 == 6 ~ "Rural",
  TRUE ~ "Unknown"
))

df_controls = df_controls %>% mutate(urban_rural_classification = case_when(
  ur6_2016 == 1 | ur6_2016 == 2 | ur6_2016 == 3 ~ "Urban",
  ur6_2016 == 4 | ur6_2016 == 5 | ur6_2016 == 6 ~ "Rural",
  TRUE ~ "Unknown"
))

# Set up the reference levels
df_all$age_gp = relevel(df_all$age_gp, ref = "25-29")
df_all$simd = factor(df_all$simd)
df_all$simd = relevel(df_all$simd, ref=5)
df_all$urban_rural_classification = factor(df_all$urban_rural_classification)
df_all$urban_rural_classification = relevel(df_all$urban_rural_classification, ref="Rural")

df_all$bmi_gp = relevel(df_all$bmi_gp, ref = "(18.5,25]")
df_all$ethnic_gp = factor(df_all$ethnic_gp)
df_all$ethnic_gp = relevel(df_all$ethnic_gp, ref = "White")

# Remove anyone who has missing SIMD or Urban rural classification data
df_all = df_all %>% filter(!is.na(simd) & !is.na(urban_rural_classification))
print("Number of rows after removing missing information")
print(sum(df_all$eave_weight))

df_all = df_all %>% mutate(time = event_date - a_begin)
df_all = df_all %>% mutate(health_board = if_else(is.na(health_board), "Unknown", health_board))

df_all = df_all %>%
  mutate(vacc_1_diff = event_date - date_vacc_1,
         vacc_2_diff = event_date - date_vacc_2,
         vacc_3_diff = event_date - date_vacc_3,
         vacc_4_diff = event_date - date_vacc_4,
         vacc_5_diff = event_date - date_vacc_5) %>%
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

df_all$covid_vs = factor(df_all$covid_vs,  labels = c("Unvaccinated", 
                                                                #"1st Dose 0 - 14 days",
                                                                "1st Dose 14+ days",
                                                                "2nd Dose 14+ days",
                                                                #"3rd Dose 0 - 14 days",
                                                                "3rd Dose 14+ days",
                                                                #"4th Dose 0 - 14 days",
                                                                "4th Dose 14+ days",
                                                                "5th Dose 14+ days"))

df_all = df_all %>%
  mutate(vacc_1_diff = event_date - date_flu_vacc_1) %>%
  mutate(flu_vs = case_when(!is.na(vacc_1_diff) & vacc_1_diff > 14 ~ "v1_2+",
                            !is.na(vacc_1_diff) & vacc_1_diff > 0 ~ "v1_0:2",
                            is.na(vacc_1_diff) | vacc_1_diff < 1 ~ "uv"))

df_all$flu_vs = factor(df_all$flu_vs, labels = c("Unvaccinated", "0 - 14 days", "14+ days"))
df_all = df_all %>% mutate(hosp_los_gp = cut(hosp_los, breaks=c(-1, 0, 1, 2, 5, 10, 20, max(hosp_los))))

# Need to set this for displaying the cox model input
df_controls = df_controls %>% mutate(hosp_los_gp = cut(hosp_los, breaks=c(-1, 0, 1, 2, 5, 10, 20, max(hosp_los))))

df_all = df_all %>% mutate(EAVE_Smoking_Status_Best = as.factor(EAVE_Smoking_Status_Best),
                           EAVE_Smoking_Status_Worst = as.factor(EAVE_Smoking_Status_Worst))

# Do the same for the controls
df_controls = df_controls %>% mutate(health_board = if_else(is.na(health_board), "Unknown", health_board))

df_controls = df_controls %>%
  mutate(vacc_1_diff = event_date - date_vacc_1,
         vacc_2_diff = event_date - date_vacc_2,
         vacc_3_diff = event_date - date_vacc_3,
         vacc_4_diff = event_date - date_vacc_4,
         vacc_5_diff = event_date - date_vacc_5) %>%
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

df_controls$covid_vs = factor(df_controls$covid_vs,  labels = c("Unvaccinated", 
                                                                #"1st Dose 0 - 14 days",
                                                                "1st Dose 14+ days",
                                                                "2nd Dose 14+ days",
                                                                #"3rd Dose 0 - 14 days",
                                                                "3rd Dose 14+ days",
                                                                #"4th Dose 0 - 14 days",
                                                                "4th Dose 14+ days",
                                                                "5th Dose 14+ days"))

df_controls = df_controls %>%
  mutate(vacc_1_diff = event_date - date_flu_vacc_1) %>%
  mutate(flu_vs = case_when(!is.na(vacc_1_diff) & vacc_1_diff > 14 ~ "v1_2+",
                            !is.na(vacc_1_diff) & vacc_1_diff > 0 ~ "v1_0:2",
                            is.na(vacc_1_diff) | vacc_1_diff < 1 ~ "uv"))

df_controls$flu_vs = factor(df_controls$flu_vs, labels = c("Unvaccinated", "0 - 14 days", "14+ days"))

df_controls = df_controls %>% mutate(EAVE_Smoking_Status_Best = as.factor(EAVE_Smoking_Status_Best),
                           EAVE_Smoking_Status_Worst = as.factor(EAVE_Smoking_Status_Worst))

#df_controls = df_controls %>% filter(!is.na(simd) & !is.na(urban_rural_classification))

remove(list=ls(pa="^z"))

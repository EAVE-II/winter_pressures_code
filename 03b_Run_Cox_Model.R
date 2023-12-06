# 
adult_vars = c("Sex", "age_gp", 
  "simd","n_risk_gps", "bmi_gp", 
  "num_prev_admission_gp", "urban_rural_classification")

cyp_vars = c("Sex", "age_gp", 
               "simd","n_risk_gps",  
               "num_prev_admission_gp", "urban_rural_classification")

adult_vars_flu = c("Sex", "age_gp", 
               "simd","n_risk_gps", "bmi_gp", 
               "num_prev_admission_gp", "urban_rural_classification")

adult_vars_covid = c("Sex", "age_gp", 
                   "simd","n_risk_gps", "bmi_gp", 
                   "num_prev_admission_gp", "urban_rural_classification")

if (use_vaccination == 1) {
  adult_vars_flu = append(adult_vars_flu, "flu_vs") 
  adult_vars_covid = append(adult_vars_covid, "covid_vs")
}

if (use_smoking == 1) {
  #adult_vars = append(adult_vars_covid, "EAVE_Smoking_Status_Best")
  #adult_vars_flu = append(adult_vars_flu, "EAVE_Smoking_Status_Best") 
  #adult_vars_covid = append(adult_vars_covid, "EAVE_Smoking_Status_Best")
  adult_vars = append(adult_vars_covid, "EAVE_Smoking_Status_Worst")
  adult_vars_flu = append(adult_vars_flu, "EAVE_Smoking_Status_Worst") 
  adult_vars_covid = append(adult_vars_covid, "EAVE_Smoking_Status_Worst")
}

adult_expression = as.formula(paste0("Surv(time, event) ~ ", paste(adult_vars, collapse=" + ")))
cyp_expression = as.formula(paste0("Surv(time, event) ~ ", paste(cyp_vars, collapse=" + ")))
adult_expression_flu = as.formula(paste0("Surv(tstart, tstop, flu_event) ~ ", paste(adult_vars_flu, collapse=" + ")))
adult_expression_covid = as.formula(paste0("Surv(tstart, tstop, covid_event) ~ ", paste(adult_vars_covid, collapse=" + ")))


# Set up the overall factors
df_all$Sex = factor(df_all$Sex)

# Model for entire population
model_all = coxph(adult_expression,data=df_all, weights=eave_weight)
surv_all = survfit(adult_expression,data=df_all, weights=eave_weight)

# Age specific models
df_adults = df_all %>% filter(ageYear > 17)
df_cyp = df_all %>% filter(ageYear < 18)

print("Number of adult rows")
print(sum(df_adults$eave_weight))
print("Adult control/event split")
print(df_adults %>% group_by(event) %>% summarise(n = sum(eave_weight)))

print("Number of CYP rows")
print(sum(df_cyp$eave_weight))
print("CYP control/event split")
print(df_cyp %>% group_by(event) %>% summarise(n = sum(eave_weight)))

# Reset the age levels
df_adults$age_gp = droplevels(df_adults$age_gp)
df_cyp$age_gp = droplevels(df_cyp$age_gp)
df_cyp$ageYear = factor(df_cyp$ageYear)
df_cyp$ageYear = relevel(df_cyp$ageYear, ref = "17")
df_cyp$age_gp = relevel(df_cyp$age_gp, ref = "6-17")

z_df_adults = df_adults %>% mutate(event = as.factor(event))
z_df_cyp = df_cyp %>% mutate(event = as.factor(event))

cuminc_adults = cuminc(Surv(time, event) ~ 1, data = z_df_adults)
cuminc_cyp = cuminc(Surv(time, event) ~ 1, data = z_df_cyp)

cuminc_plot_adults = cuminc_adults %>% ggcuminc() + add_confidence_interval() + ylim(0, 0.1)

cuminc_plot_cyp = cuminc_cyp %>% ggcuminc() + add_confidence_interval() + ylim(0, 0.1) 

ggsave("cuminc_plot_adults.png", cuminc_plot_adults)
ggsave("cuminc_plot_cyp.png", cuminc_plot_cyp)

model_adults = coxph(adult_expression,data=df_adults,weights=eave_weight)
model_fit_adults = cox.zph(model_adults)
coxzph_plot_adults = ggcoxzph(model_fit_adults)

for (i in 1:length(coxzph_plot_adults)) {
  ggsave(paste0(paste0("coxzph_plot_adults_", i), ".png"), coxzph_plot_adults[[i]])
}
model_cyp = coxph(cyp_expression,data=df_cyp,weights=eave_weight)
model_fit_cyp = cox.zph(model_cyp)

# Pathogen specific models
# This is a bit more complex because we need to pick out the controls too

if (use_tests == 1) {
  df_adults = df_adults %>% mutate(covid_event = if_else(covid_admit_with_test == 1, 1, 0))
  df_adults = df_adults %>% mutate(flu_event = if_else(flu_admit_with_test == 1, 1, 0))
} else if (broad_defintion == 1) {
  df_adults = df_adults %>% mutate(flu_event = if_else(flu_admit_secondary == 1 | flu_admit == 1, 1, 0))
  df_cyp = df_cyp %>% mutate(flu_event = if_else(flu_admit_secondary == 1 | flu_admit == 1, 1, 0))
  df_adults = df_adults %>% mutate(covid_event = if_else(covid_admit_secondary == 1 | covid_admit == 1, 1, 0))
  df_cyp = df_cyp %>% mutate(covid_event = if_else(covid_admit_secondary == 1 | covid_admit == 1, 1, 0))
} else {
  df_adults = df_adults %>% mutate(flu_event = if_else(flu_admit == 1, 1, 0))
  df_cyp = df_cyp %>% mutate(flu_event = if_else(flu_admit == 1, 1, 0))
  df_adults = df_adults %>% mutate(covid_event = if_else(covid_admit == 1, 1, 0))
  df_cyp = df_cyp %>% mutate(covid_event = if_else(covid_admit == 1, 1, 0))
}

# Handle the vaccine status 
# Firstly handle this for the flu vaccine as a person only receives one of them
# Set this to the day before our actual start to ensure time is positive
a_begin_cox = as.Date("2022-08-31")
z_flu = df_adults %>% select(-event)
z_flu <- tmerge(z_flu, z_flu, id=EAVE_LINKNO, flu_event = event(event_date, flu_event), tstart=a_begin_cox, tstop=event_date)
z_flu <- tmerge(z_flu,df_adults, id=EAVE_LINKNO, flu_vs=tdc(date_flu_vacc_1))
z_flu$flu_vs = factor(z_flu$flu_vs, labels=c("Unvaccinated", "Vaccinated"))
df_adults_flu = z_flu

df_adults_flu = df_adults_flu %>% mutate(tstart = difftime(tstart, a_begin_cox, unit="days")) %>%
  mutate(tstop = difftime(tstop, a_begin_cox, unit="days"))

# Now we do this for COVID - there is a lot more to do due to the 5 doses!
z_covid = df_adults %>% select(-event)
z_covid <- tmerge(z_covid, z_covid, id=EAVE_LINKNO, covid_event = event(event_date, covid_event), tstart=a_begin_cox, tstop=event_date)
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_1))
names(z_covid)[names(z_covid)=="per1"] <- "pv_uv"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_1 + 14))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v1_0:2"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_2))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v1_2+"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_2 + 14))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v2_0:2"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_3))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v2_2+"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_3 + 14))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v3_0:2"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_4))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v3_2+"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_4 + 14))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v4_0:2"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_5))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v4_2+"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(date_vacc_5 + 14))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v5_0:14"
z_covid <- tmerge(z_covid,df_adults, id=EAVE_LINKNO, per1=tdc(event_date))
names(z_covid)[names(z_covid)=="per1"] <- "pv_v5_2+"

z_names <- names(z_covid)[grepl("pv_", names(z_covid))]
z_covid <- z_covid %>% mutate(pv_period = apply(z_covid[,z_names], 1, sum)) 
z_covid <- z_covid %>% dplyr::select(-all_of(z_names)) 

z_covid <- z_covid %>%
  mutate(pv_period_f = case_when(pv_period==0 ~ "uv",
                                 pv_period==1 ~ "v1_0:2",
                                 pv_period==2 ~ "v1_2+",
                                 pv_period==3 ~ "v2_0:1",
                                 pv_period==4 ~ "v2_2+",
                                 pv_period==5 ~ "v3_0:1",
                                 pv_period==6 ~ "v3_2+",
                                 pv_period==7 ~ "v4_0:1+",
                                 pv_period==8 ~ "v4_2+",
                                 pv_period==9 ~ "v5_0:2",
                                 TRUE ~ "v5_2+")) %>% 
  mutate(pv_period_f = factor(pv_period_f, levels=c("uv","v1_0:2","v1_2+", "v2_0:1","v2_2+", "v3_0:1", "v3_2+",
                                                    "v4_0:2", "v4_2+", "v5_0:2", "v5_2+"))) %>%
  mutate(covid_vs = pv_period_f)

# This is a bit ugly, but basically we combine vaccination groups to increase the 
# number of events. I haven't done this above as we might want to revert it in the future
z_covid = z_covid %>% mutate(covid_vs = case_when(covid_vs == "uv" ~ 0,
                                                           covid_vs == "v1_0:2" ~ 0,
                                                           covid_vs == "v1_2+" ~ 0,
                                                           covid_vs == "v2_0:1" ~ 0,
                                                           covid_vs == "v2_2+" ~ 0,
                                                           covid_vs == "v3_0:1" ~ 1,
                                                           covid_vs == "v3_2+" ~ 1,
                                                           covid_vs == "v4_0:1+" ~ 2,
                                                           covid_vs == "v4_2+" ~ 2,
                                                           covid_vs == "v5_0:2" ~ 2,
                                                           covid_vs == "v5_2+" ~ 2))

z_covid$covid_vs = factor(z_covid$covid_vs)

df_adults_covid = z_covid

df_adults_covid = df_adults_covid %>% mutate(tstart = difftime(tstart, a_begin_cox, unit="days")) %>%
  mutate(tstop = difftime(tstop, a_begin_cox, unit="days"))

df_all = df_all %>% mutate(time = difftime(event_date, a_begin, unit="days"))
print("Adult flu control/event split")
z_df_adults_flu = df_adults_flu %>% group_by(EAVE_LINKNO) %>% summarize(case = sum(flu_event), eave_weight = first(eave_weight))
print(z_df_adults_flu %>% group_by(case) %>% summarise(n = sum(eave_weight)))
print("Adult COVID-19 control/event split")
z_df_adults_covid= df_adults %>% group_by(EAVE_LINKNO) %>% summarize(case = sum(covid_event), eave_weight = first(eave_weight))
print(z_df_adults_covid %>% group_by(case) %>% summarise(n = sum(eave_weight)))

model_flu_adults = coxph(adult_expression_flu,data=df_adults_flu,weights=eave_weight)
model_covid_adults = coxph(adult_expression_covid ,data=df_adults_covid,weights=eave_weight)

# Turn the model coefficients into tibbles for display
model_all_coefs = cox_model_to_tibble(model_all, df_all, "event", adult_vars)
model_adults_coefs = cox_model_to_tibble(model_adults, df_adults, "event", adult_vars)
model_cyp_coefs =  cox_model_to_tibble(model_cyp, df_cyp, "event", cyp_vars)

model_flu_adults_coefs = cox_model_to_tibble(model_flu_adults, df_adults_flu, "flu_event", adult_vars_flu)
model_covid_adults_coefs = cox_model_to_tibble(model_covid_adults, df_adults_covid, "covid_event", adult_vars_covid)

remove(list=ls(pa="^z"))


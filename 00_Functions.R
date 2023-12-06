
#### Summary table using weights ####
# Creates a table of cohort summaries using weights
# For categorical variables, the sum of the weights and the % is calculated
# For numerical variables, the weighed mean and weighted sd are calculated, as well as
# the weighted median and weighted IQR

# Input:
# - data = the dataset (must have weights named eave_weight)
# - dependent = a character of the dependent variables name
# - explanatory = a string of characters of the explanatory variables

# Output:
# A table with each explanatory variable as a row (multiple rows for each category if categorical)
# with two columns of the weighted summaries for the levels in the dependent variable

summary_factorlist_wt <- function(data, dependent, explanatory){
  # Create list to put in summaries into each element
  summary_tbl_list <- list()
  
  for(i in 1:length(explanatory)){
    # Extract variable
    n <- data %>%
      pull(!!sym(explanatory[i]))
    
    # If the variable has a nice name, then we use that
    variable_label = data %>% select(!!sym(explanatory[i]))
    variable_label = variable_label %>% var_label(unlist = TRUE)
    variable_label = as.vector(unlist(variable_label))
    variable_label = ifelse(variable_label == "", explanatory[i], variable_label)

    # If numeric then make weighted mean
    if(is.numeric(n)) {
      z_mean <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(mean = round(weighted.mean(!!sym(explanatory[i]), w = eave_weight, na.rm = TRUE),1),
                  sd = round(sqrt(spatstat.geom::weighted.var(!!sym(explanatory[i]), w = eave_weight)),1)) %>%
        mutate(mean.sd = paste0(mean, " (",sd,")")) %>%
        select(-mean, -sd) %>%
        mutate("characteristic" = explanatory[i]) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = mean.sd) %>%
        relocate(characteristic) %>%
        mutate(levels = "mean.sd")
      
      
      z_median <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(median = spatstat.geom::weighted.median(!!sym(explanatory[i]), w = eave_weight),
                  q1 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = eave_weight, probs = 0.25),
                  q3 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = eave_weight, probs = 0.75)) %>%
        mutate("characteristic" = explanatory[i]) %>%
        mutate(iqr = paste0(paste0(round(q1,1), ","), round(q3, 1))) %>%
        mutate(median.iqr = paste0(median, " (",iqr,")")) %>%
        select(-q1, -q3, -median, -iqr) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = median.iqr) %>%
        relocate(characteristic) %>%
        mutate(levels = "median.iqr")
      
      # Combine!!
      summary_tbl_list[[i]] <- full_join(z_mean, z_median)
      
      
      # Else get sum of weights of each level
    } else if (length(unique(data %>% pull(!!sym(dependent) ) ) ) ==1) {
      # This is for when there is only one level in the dependent variable
      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i])) %>%
        dplyr::summarise(n = sum(eave_weight)) %>%
        ungroup() %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate(n = ifelse(n < 5, "< 5", paste0(round(n, 0)))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        dplyr::mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        dplyr::rename("levels"= explanatory[i], !!dependent := n_perc) %>%
        mutate(characteristic = variable_label) %>%
        relocate(characteristic)
      
    } else {
      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i]), !!sym(dependent)) %>%
        dplyr::summarise(n = sum(eave_weight)) %>%
        ungroup() %>%
        group_by(!!sym(dependent)) %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate(n = ifelse(n < 5, "< 5", paste0(round(n,0)))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = n_perc) %>%
        dplyr::rename("levels"=explanatory[i]) %>%
        mutate(characteristic = variable_label) %>%
        relocate(characteristic)
      
    }
  }

  # Combine list together to make dataset
  summary_tbl_wt <- summary_tbl_list %>%
    reduce(full_join)
  
  summary_tbl_wt
}


summary_factorlist_tm <- function(data, dependent, explanatory){
  # A modification of Rachel's code which deals only with 1 dependent level
  # Create list to put in summaries into each element
  summary_tbl_list <- list()
  
  for(i in 1:length(explanatory)){
    # Extract variable
    n <- data %>%
      pull(!!sym(explanatory[i]))
    
    # If the variable has a nice name, then we use that
    variable_label = data %>% select(!!sym(explanatory[i]))
    variable_label = variable_label %>% var_label(unlist = TRUE)
    variable_label = as.vector(unlist(variable_label))
    variable_label = ifelse(variable_label == "", explanatory[i], variable_label)
    
      # This is for when there is only one level in the dependent variable
      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i])) %>%
        dplyr::summarise(n = sum(eave_weight)) %>%
        ungroup() %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        dplyr::mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        dplyr::rename("levels"= explanatory[i], !!dependent := n_perc) %>%
        mutate(characteristic = variable_label) %>%
        relocate(characteristic)
      
  }
  
  #summary_tbl_wt = summary_tbl_list[[1]]
  
  #for (i in 2:length(summary_tbl_list)) {
  #  summary_tbl_wt = summary_tbl_wt %>% add_row(summary_tbl_list[[i]])
  #}
  # Combine list together to make dataset
  #summary_tbl_wt <- summary_tbl_list %>%
  #  reduce(full_join)
  
  summary_tbl_list
}

redact_low_counts = function(df_count) {
  # Redacts any values in the dataframe that have a count of less than 5
  z_df = df_count %>% mutate(n = if_else(n < 5, "< 5", paste0(n)))
  return(z_df)
}

get_reference_levels = function(z_df, vars) {
  ret = list()
  
  for (i in 1:length(vars)) {
    v = vars[i]
    ref = levels(z_df %>% select(!!sym(v)) %>% pull)[1]
    
    tmp = tibble(coef_name = paste0(vars[i], ref), coef=0, HR = 1, LCL = 1, UCL = 1,
                 p_value = 0, `se(coef)` = NA, z = NA)
    ret[[i]] = tmp
  }
  
  return(ret %>% reduce(full_join))
}

get_reference_levels_lr = function(z_df, vars) {
  ret = list()
  
  for (i in 1:length(vars)) {
    v = vars[i]
    ref = levels(z_df %>% select(!!sym(v)) %>% pull)[1]
    
    tmp = tibble(coef_name = paste0(vars[i], ref), OR = 1, LCL = 1, UCL = 1,
                 p_value = 0, `Std. Error` = NA, `z value` = NA)
    ret[[i]] = tmp
  }
  
  return(ret %>% reduce(full_join))
}

cox_model_to_tibble = function(model, z_df, event_name, vars) {
  vals = as_tibble(summary(model)$coefficients, rownames="coef_name") %>%
    rename(HR = `exp(coef)`, p_value = `Pr(>|z|)`)
  con = as_tibble(confint(model))
  colnames(con) = c("LCL", "UCL")
  con$LCL = exp(con$LCL)
  con$UCL = exp(con$UCL)
  
  # Remove the robust se if it exists
  vals = vals %>% select(-one_of("robust se"))
  
  ret = vals %>% cbind(con)
  
  ref_levels = get_reference_levels(z_df, vars)
  print(ref_levels)
  print(ret)
  ret = ret %>% rbind(ref_levels)
  
  event_count = create_event_count(z_df, vars, event_name)
  
  ret = ret %>% left_join(event_count, by="coef_name")
  
  return(ret %>% arrange(coef_name))
}

lr_model_to_tibble = function(model, z_df, event_name, vars) {
  vals = as_tibble(summary(model)$coefficients, rownames="coef_name") %>%
    rename(OR = Estimate, p_value = `Pr(>|z|)`)
  
  vals$OR = exp(vals$OR)
  con = as_tibble(confint(model))
  colnames(con) = c("LCL", "UCL")
  con$LCL = exp(con$LCL)
  con$UCL = exp(con$UCL)
  
  ref_levels = get_reference_levels_lr(z_df, vars)
  print(ref_levels)
  
  ret = vals %>% cbind(con)
  print(ret)
  ret = ret %>% rbind(ref_levels)
  
  event_count = create_event_count(z_df, vars, event_name)
  
  ret = ret %>% left_join(event_count, by="coef_name")
  
  return(ret %>% arrange(coef_name))
}

check_for_condition_in_secondary_codes = function(hosp_df, three_char_codes, four_char_codes = c()) {
  hosp_df = hosp_df %>% mutate(
    admission_1 = if_else(
      substr(OTHER_CONDITION_1, 0, 3) %in% three_char_codes |
        substr(OTHER_CONDITION_1, 0, 4) %in% four_char_codes, 1, 0),
    admission_2 = if_else(
      substr(OTHER_CONDITION_2, 0, 3) %in% three_char_codes |
        substr(OTHER_CONDITION_2, 0, 4) %in% four_char_codes, 1, 0),
    admission_3 = if_else(
      substr(OTHER_CONDITION_3, 0, 3) %in% three_char_codes |
        substr(OTHER_CONDITION_3, 0, 4) %in% four_char_codes, 1, 0),
    admission_4 = if_else(
      substr(OTHER_CONDITION_4, 0, 3) %in% three_char_codes |
        substr(OTHER_CONDITION_4, 0, 4) %in% four_char_codes, 1, 0),
    admission_5 = if_else(
      substr(OTHER_CONDITION_5, 0, 3) %in% three_char_codes |
        substr(OTHER_CONDITION_5, 0, 4) %in% four_char_codes, 1, 0)) %>%
    mutate(secondary = if_else(
      admission_1 == 1 | admission_2 == 1 |
        admission_3 == 1 | admission_4 == 1 | admission_5 == 1, 1, 0)
    ) %>%
    select(-admission_1, -admission_2, -admission_3, -admission_4,
           -admission_5) %>%
    select(EAVE_LINKNO, CIS_MARKER, secondary)

  return(hosp_df)
}

create_event_count = function(dataset, vars, event_name) {
  summary_lst = list()
  
  for (i in 1:length(vars)) {
    n <- dataset %>% filter(!!sym(event_name) == 1) %>% count(!!sym(vars[i])) %>%
      rename(levels = !!sym(vars[i]))
    
    var_sum = n %>% mutate(coef_name = paste0(vars[i], levels)) %>% select(-levels)
    summary_lst[[i]] = var_sum
  }
  
  event_count = summary_lst %>% reduce(full_join)
  return(event_count)
}
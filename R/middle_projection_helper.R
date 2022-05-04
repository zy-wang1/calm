
#' @export
ipw_middle <- function(task, lik, ipw_args, fold_number){

  cf_likelihood_control = ipw_args$cf_likelihood_control
  cf_likelihood_treatment = ipw_args$cf_likelihood_treatment
  intervention_list_treatment <- ipw_args$intervention_list_treatment
  intervention_list_control <- ipw_args$intervention_list_control
  # cf_task_treatment <- ipw_args$cf_task_treatment
  # cf_task_control <- ipw_args$cf_task_control
  static_likelihood <- ipw_args$static_likelihood
  # # todo: extend for stochastic
  # cf_task_treatment <- cf_likelihood_treatment$enumerate_cf_tasks(task)[[1]]
  # cf_task_control <- cf_likelihood_control$enumerate_cf_tasks(task)[[1]]

  cf_task_treatment <- cf_likelihood_treatment$enumerate_cf_tasks(task)[[1]]
  cf_task_control <- cf_likelihood_control$enumerate_cf_tasks(task)[[1]]


  intervention_nodes <-
    # intersect(
    union(names(intervention_list_treatment), names(intervention_list_control))
  #   ,
  #   names(task$npsem)
  # )

  temp_node_names <- names(task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

  Y <- task$get_tmle_node(last(temp_node_names), format = T)[[1]]

  # get list of all possible predicted lkds
  obs_data <- task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
  obs_variable_names <- colnames(obs_data)
  # ZW todo: to handle long format and wide format
  # ZW todo: see if observed_likelihood needs to change to targeted likelihood

  intervention_variables <- map_chr(task$npsem[intervention_nodes], ~.x$variables)
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  intervention_levels_treat <- map_dbl(intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
  intervention_levels_control <- map_dbl(intervention_list_control, ~.x$value %>% as.character %>% as.numeric)

  list_H <- get_obs_H_full(task, obs_data,
                           current_likelihood = static_likelihood,
                           # current_likelihood = lik,
                           cf_task_treatment, cf_task_control,
                           intervention_variables, intervention_levels_treat, intervention_levels_control,
                           fold_number = fold_number
                           # , bound = 0.05
                           )


  list_newH <- list()
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      if (ind_var %in% loc_Z) list_newH[[ind_var]] <- (list_H[[ind_var]] * Y) %>% as.matrix
      if (ind_var %in% loc_RLY) list_newH[[ind_var]] <- (list_H[[ind_var]] * Y) %>% as.matrix
    }
  }
  names(list_newH) <- temp_node_names

  return(list_newH)

  # ZW todo: for categorical variables
}


#' @export
gradient_generator_middle <- function(tmle_task, lik,  node, include_outcome = T, ipw_args = NULL, fold_number){

  task <- tmle_task$get_regression_task(node)
  if (include_outcome) {
    IC <- ipw_middle(tmle_task, lik,  ipw_args, fold_number)[[node]] %>% as.vector
    cols <- task$add_columns(data.table(IC = IC))
  } else {
    cols <- task$add_columns(NULL)
  }
  task <- task$clone()
  nodes <- task$nodes
  if (include_outcome) {
    nodes$outcome <- "IC"
  }
  nodes$covariates <- c(nodes$covariates, tmle_task$npsem[[node]]$variables)

  task$initialize(
    task$internal_data,
    nodes = nodes,
    folds = task$folds,
    column_names = cols,
    row_index = task$row_index,
    outcome_type = "continuous"
  )
  return(task)
}












#' @export
ipw_middle_survival <- function(task, lik, ipw_args, fold_number){

  cf_likelihood_control = ipw_args$cf_likelihood_control
  cf_likelihood_treatment = ipw_args$cf_likelihood_treatment
  intervention_list_treatment <- ipw_args$intervention_list_treatment
  intervention_list_control <- ipw_args$intervention_list_control
  static_likelihood <- ipw_args$static_likelihood

  # # todo: extend for stochastic
  cf_task_treatment <- cf_likelihood_treatment$enumerate_cf_tasks(task)[[1]]
  cf_task_control <- cf_likelihood_control$enumerate_cf_tasks(task)[[1]]

  intervention_nodes <- union(names(intervention_list_treatment), names(intervention_list_control))

  temp_node_names <- names(task$npsem)
  loc_delta_nodes <- grep("delta_", temp_node_names)
  if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
  loc_A_E <- grep("^A_E", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  # loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
  if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

  # use correct outcome and correct full covariates
  outcome_node <- ipw_args$outcome_node
  tau <- ipw_args$tau
  if (is.null(outcome_node)) outcome_node <- last(temp_node_names)
  if (is.null(tau)) tau <- strsplit(self$outcome_node, "_")[[1]] %>% last %>% as.numeric  # tau is the last time point involved in the outcome
  loc_last_node <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]] %>% last %>% as.numeric) <= tau) %>% last
  last_node <- temp_node_names[loc_last_node]

  Y <- task$get_tmle_node(outcome_node)
  Y[is.na(Y)] <- 0

  # get list of all possible predicted lkds
  obs_data <- task$data %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))
  obs_variable_names <- colnames(obs_data)
  # ZW todo: to handle long format and wide format
  # ZW todo: see if observed_likelihood needs to change to targeted likelihood

  intervention_variables <- map_chr(task$npsem[intervention_nodes], ~.x$variables)
  intervention_variables_loc <- map_dbl(paste0(intervention_variables, "$"), ~grep(.x, obs_variable_names))
  intervention_levels_treat <- map_dbl(intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
  intervention_levels_control <- map_dbl(intervention_list_control, ~.x$value %>% as.character %>% as.numeric)

  list_H <- get_obs_H_full_survival(task, obs_data,
                           current_likelihood = static_likelihood,
                           # current_likelihood = lik,
                           cf_task_treatment, cf_task_control,
                           intervention_variables, intervention_levels_treat, intervention_levels_control,
                           fold_number = fold_number
                           # , bound = 0.05
                           , loc_last_node = loc_last_node
  )


  list_newH <- list()
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      # if there is missingness, match the get_regression_task structure
      if (!is.null(task$npsem[[ind_var]]$censoring_node$variables)) {
        if_observed <- task$get_tmle_node(task$npsem[[ind_var]]$censoring_node$variables)  # force Y 0 where censored due to bivariate trt nodes
        if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
        if_observed[is.na(if_observed)] <- F
      }
      if (ind_var %in% loc_Z) list_newH[[ind_var]] <- (list_H[[ind_var]] * Y[if_observed]) %>% as.matrix
      if (ind_var %in% loc_RLY) list_newH[[ind_var]] <- (list_H[[ind_var]] * Y[if_observed]) %>% as.matrix
    }
  }
  names(list_newH) <- temp_node_names

  return(list_newH)

  # ZW todo: for categorical variables
}


#' @export
gradient_generator_middle_survival <- function(tmle_task, lik,  node, include_outcome = T, ipw_args = NULL, fold_number){

  task <- tmle_task$get_regression_task(node)
  if (include_outcome) {
    IC <- ipw_middle_survival(tmle_task, lik,  ipw_args, fold_number)[[node]] %>% as.vector
    cols <- task$add_columns(data.table(IC = IC))
  } else {
    cols <- task$add_columns(data.table(NA))
  }
  task <- task$clone()
  nodes <- task$nodes
  if (include_outcome) {
    nodes$outcome <- "IC"
  }
  nodes$covariates <- c(nodes$covariates, tmle_task$npsem[[node]]$variables)

  task$initialize(
    task$internal_data,
    nodes = nodes,
    folds = task$folds,
    column_names = cols,
    row_index = task$row_index,
    outcome_type = "continuous"
  )
  return(task)
}


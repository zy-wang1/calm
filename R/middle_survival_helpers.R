#' @export
get_obs_H_full_survival <- function(tmle_task, obs_data, current_likelihood,
                           cf_task_treatment, cf_task_control,
                           intervention_variables, intervention_levels_treat, intervention_levels_control,
                           fold_number = "full",
                           bound = NULL,
                           loc_last_node = NULL
) {
  # when called in ipw helper, cf_tasks are created based on the task (can be observed task or integral expanded task)
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_delta_nodes <- grep("delta_", temp_node_names)
  if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
  loc_A_E <- grep("A_E", temp_node_names)
  loc_A_C <- grep("A_C", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  # loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
  intervention_variables_loc <- map_dbl(paste0(intervention_variables, "$"), ~grep(.x, obs_variable_names))

  if (!is.null(loc_last_node)) {
      loc_RLY <- loc_RLY[loc_RLY <= loc_last_node]
      loc_A_E <- loc_A_E[loc_A_E <= loc_last_node]
      loc_A_C <- loc_A_C[loc_A_C <= loc_last_node]
      loc_Z <- loc_Z[loc_Z <= loc_last_node]
      intervention_variables_loc <- intervention_variables_loc[intervention_variables_loc <= loc_last_node]


  }

  # get a list of corresponding H covariates; ordered by nodes, not variables
  list_H <- list()
  # calculate RLY nodes
  for (temp_ind in loc_RLY) {
    loc_A_needed <- loc_A_E
    loc_Z_needed <- loc_Z
    # this is the At indicators for H_RLY; now
    A_ind <- apply(sapply(loc_A_needed, function(k) {
      obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_treat[which(loc_A_needed == k)]
    }), 1, prod) == 1

    # these A probs will be taken as product
    part_A <- lapply(loc_A_needed, function(k) {
      temp_p_A_E <- current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)  # A_E | A_C=1
      if (!is.null(tmle_task$npsem[[k]]$censoring_node$variables)) {  # if there is missingness in A_E
        temp_full <- if_A_E_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k]]$censoring_node$name)
        if (!is.logical(if_A_E_observed)) if_A_E_observed <- if_A_E_observed == 1  # in case it is a binary node
        if_A_E_observed[is.na(if_A_E_observed)] <- F
        temp_full[if_A_E_observed] <- temp_p_A_E
        temp_full[!if_A_E_observed] <- NA
        temp_p_A_E <- temp_full
      }
      k_A_C <- loc_A_C[loc_A_C < k]
      k_A_C <- k_A_C[which.min(abs(k_A_C - k))]  # always let censoring node to lead each time point
      temp_p_A_C <- current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k_A_C], fold_number)  # A_C=1; to be aligned
      if (!is.null(tmle_task$npsem[[k_A_C]]$censoring_node$variables)) {  # if there is missingness in A_E
        temp_full <- if_A_C_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$censoring_node$name)
        if (!is.logical(if_A_C_observed)) if_A_C_observed <- if_A_C_observed == 1  # in case it is a binary node
        if_A_C_observed[is.na(if_A_C_observed)] <- F
        observed_A_C <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$name)
        if_A_C_observed[is.na(observed_A_C)] <- F
        temp_full[if_A_C_observed] <- temp_p_A_C
        temp_full[!if_A_C_observed] <- NA
        temp_p_A_C <- temp_full
      }
      return(temp_p_A_C * temp_p_A_E)
    }) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      temp_p <- current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)
      if (!is.null(tmle_task$npsem[[k]]$censoring_node$variables)) {  # if there is missingness
        temp_full <- if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k]]$censoring_node$variables)
        if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
        if_observed[is.na(if_observed)] <- F
        temp_full[if_observed] <- temp_p
        temp_full[!if_observed] <- NA
        temp_p <- temp_full
      }
      return(temp_p)
    }) %>% pmap_dbl(prod)
    if(length(part_Z) == 0) part_Z <- 1

    if (!is.null(bound)) part_A[part_A < bound] <- bound
    temp_vec <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
    temp_vec[is.na(temp_vec)] <- 0  # due to bivariate trt nodes or g-comp
    if(!is.null(tmle_task$npsem[[temp_ind]]$censoring_node$variables)) {
      if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[temp_ind]]$censoring_node$variables)
      if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
      if_observed[is.na(if_observed)] <- F
      temp_vec <- temp_vec[if_observed]
    }
    list_H[[temp_ind]] <- temp_vec
  }
  # calculate Z nodes
  for (temp_ind in loc_Z) {
    loc_A_needed <- loc_A_E
    loc_RLY_needed <- loc_RLY
    A_ind <-
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_control[which(loc_A_needed == k)]
      }), 1, prod) == 1
    part_A <- lapply(loc_A_needed, function(k) {
      temp_p_A_E <- current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)  # A_E | A_C=1
      if (!is.null(tmle_task$npsem[[k]]$censoring_node$variables)) {  # if there is missingness in A_E
        temp_full <- if_A_E_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k]]$censoring_node$variables)
        if (!is.logical(if_A_E_observed)) if_A_E_observed <- if_A_E_observed == 1  # in case it is a binary node
        if_A_E_observed[is.na(if_A_E_observed)] <- F
        temp_full[if_A_E_observed] <- temp_p_A_E
        temp_full[!if_A_E_observed] <- NA
        temp_p_A_E <- temp_full
      }
      k_A_C <- loc_A_C[loc_A_C < k]
      k_A_C <- k_A_C[which.min(abs(k_A_C - k))]  # always let censoring node to lead each time point
      temp_p_A_C <- current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k_A_C], fold_number)  # A_C=1; to be aligned
      if (!is.null(tmle_task$npsem[[k_A_C]]$censoring_node$variables)) {  # if there is missingness in A_E
        temp_full <- if_A_C_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$censoring_node$variables)
        if (!is.logical(if_A_C_observed)) if_A_C_observed <- if_A_C_observed == 1  # in case it is a binary node
        if_A_C_observed[is.na(if_A_C_observed)] <- F
        observed_A_C <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$name)
        if_A_C_observed[is.na(observed_A_C)] <- F
        temp_full[if_A_C_observed] <- temp_p_A_C
        temp_full[!if_A_C_observed] <- NA
        temp_p_A_C <- temp_full
      }
      return(temp_p_A_C * temp_p_A_E)
    }) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_RLY <- lapply(loc_RLY_needed, function(k) {
      temp_p <- current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)
      if (!is.null(tmle_task$npsem[[k]]$censoring_node$variables)) {  # if there is missingness
        temp_full <- if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k]]$censoring_node$variables)
        if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
        if_observed[is.na(if_observed)] <- F
        temp_full[if_observed] <- temp_p
        temp_full[!if_observed] <- NA
        temp_p <- temp_full
      }
      return(temp_p)
    }) %>% pmap_dbl(prod)
    if (!is.null(bound)) part_A[part_A < bound] <- bound
    temp_vec <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
    temp_vec[is.na(temp_vec)] <- 0  # due to bivariate trt nodes or g-comp
    if(!is.null(tmle_task$npsem[[temp_ind]]$censoring_node$variables)) {
      if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[temp_ind]]$censoring_node$variables)
      if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
      if_observed[is.na(if_observed)] <- F
      temp_vec <- temp_vec[if_observed]
    }
    list_H[[temp_ind]] <- temp_vec
  }
  return(list_H)
}

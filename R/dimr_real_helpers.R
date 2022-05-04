# get_obs_H_full_survival

# # not sure if absolutely needed
# temp_node_names <- names(dimr_initial_likelihood$training_task$npsem)
# temp_node_names <- temp_node_names[-grep("^pred_|^initial_gradient_|^id$|^t$|^delta_", temp_node_names)]

# # if (is.null(outcome)) private$.outcome_node <- last(temp_node_names) else private$.outcome_node <- outcome
# outcome_node <- last(temp_node_names)
# # if (is.null(tau)) tau <- strsplit(self$outcome_node, "_")[[1]] %>% last %>% as.numeric  # tau is the last time point involved in the outcome
# if (is.null(tau)) tau <- strsplit(outcome_node, "_")[[1]] %>% last %>% as.numeric  # tau is the last time point involved in the outcome
# loc_last_node <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]] %>% last %>% as.numeric) <= tau) %>% last  # among the node names; use loc_last_node carefully
# last_node <- temp_node_names[loc_last_node]


# intervention_list_treatment <- treatment
# intervention_list_control <- control
# which_intervention_keep <- sapply(intervention_list_treatment, function(u) which(u$name == temp_node_names)) <= loc_last_node


# tmle_task <- dimr_task

# intervention_nodes <- union(names(intervention_list_treatment), names(intervention_list_control))
# intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
# intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
# intervention_levels_treat <- map_dbl(intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
# intervention_levels_control <- map_dbl(intervention_list_control, ~.x$value %>% as.character %>% as.numeric)



# obs_variable_names <- colnames(obs_data)

# current_likelihood <- dimr_initial_likelihood

# obs_data <- dimr_task$data %>% copy
# cols_to_delete <- names(obs_data)[grep("^pred_|^initial_gradient_|^id$|^t$|^delta_", names(obs_data))]
# if (length(cols_to_delete) > 0) obs_data[, (cols_to_delete) := NULL]

# static_likelihood <- dimr_initial_likelihood

# ipw_args <- list(
#     # cf_likelihood_treatment = self$cf_likelihood_treatment,
#     #              cf_likelihood_control = self$cf_likelihood_control,
#     which_intervention_keep = which_intervention_keep,
#     intervention_list_treatment = intervention_list_treatment[which_intervention_keep],
#     intervention_list_control = intervention_list_control[which_intervention_keep],
#     static_likelihood = static_likelihood,
#     outcome_node = outcome_node,
#     tau = tau,
#     original_likelihood = original_likelihood,
#     initial_gradient_likelihood = initial_gradient_likelihood
# )



get_obs_H_full_survival_dimr <- function(tmle_task,   # can be the dimr task in training, can be partial task for prediction
                                         obs_data,   # might have additional dimr variables
                                         current_likelihood,   # static/delayed dimr task
         # cf_task_treatment, cf_task_control,   # need to be replaced; using intervention
         cf_likelihood_treatment_ori, cf_likelihood_control_ori, cf_likelihood_treatment_ig, cf_likelihood_control_ig,   # use cf ori/ig lkds
         which_intervention_keep, intervention_list_treatment, intervention_list_control,  # use these with ori and ig tasks instead
         intervention_variables, intervention_levels_treat, intervention_levels_control,  # not sure if absolutely needed
         original_likelihood, initial_gradient_likelihood,   # ori and ig lkds are needed now
         fold_number = "validation",
         bound = NULL
) {
    temp_ori_task <- tmle3_Task_drop_censored$new(data = obs_data, npsem = original_likelihood$training_task$npsem)
    temp_ig_data <- obs_data %>% copy
    ig_cols <- names(initial_gradient_likelihood$training_task$data)
    ig_cols <- ig_cols[grep("^G_", ig_cols)]
    temp_ig_data[, (ig_cols) := 0]
    temp_ig_task <- tmle3_Task_drop_censored$new(data = temp_ig_data, npsem = initial_gradient_likelihood$training_task$npsem, id = "id", time = "t")

    # use intervened ori and ig lkds for dimr covariates
    cf_task_treatment_ori <- cf_likelihood_treatment_ori$enumerate_cf_tasks(temp_ori_task)[[1]]
    cf_task_control_ori <- cf_likelihood_control_ori$enumerate_cf_tasks(temp_ori_task)[[1]]
    cf_task_treatment_ig <- cf_likelihood_treatment_ig$enumerate_cf_tasks(temp_ig_task)[[1]]
    cf_task_control_ig <- cf_likelihood_control_ig$enumerate_cf_tasks(temp_ig_task)[[1]]

    # make intervened tasks for dimr models
    cf_dimr_data_treatment <- cf_task_treatment_ori$data %>% copy
    cf_dimr_data_control <- cf_task_control_ori$data %>% copy
    cf_dimr_data_treatment[, (c("id", "t")) := NULL]
    cf_dimr_data_control[, (c("id", "t")) := NULL]

    cols_all <- names(tmle_task$data)
    cols_pred <- cols_all[grep("^pred_", cols_all)]
    cols_ig <- cols_all[grep("^initial_gradient_", cols_all)]

    for (each_col in cols_pred) {
        actual_col <- str_remove(each_col, "^pred_[0-9]_")
        temp_pred <- original_likelihood$get_likelihood(cf_task_treatment_ori, node = actual_col, fold_number = "validation")
        if (is(original_likelihood$training_task$npsem[[actual_col]]$censoring_node, "tmle3_Node")) {
            if_obs <- cf_task_treatment_ori$get_tmle_node(
                original_likelihood$training_task$npsem[[actual_col]]$censoring_node$name
            )
            if_obs[is.na(if_obs)] <- 0
            cur_obs <- cf_task_treatment_ori$get_tmle_node(actual_col)
            if_obs <- if_obs == 1 & !is.na(cur_obs)
            temp_val <- rep(NA, length(if_obs))
            temp_val[if_obs] <- temp_pred
        } else temp_val <- temp_pred
        cf_dimr_data_treatment[, (each_col) := temp_val]
        temp_pred <- original_likelihood$get_likelihood(cf_task_control_ori, node = actual_col, fold_number = "validation")
        if (is(original_likelihood$training_task$npsem[[actual_col]]$censoring_node, "tmle3_Node")) {
            if_obs <- cf_task_control_ori$get_tmle_node(
                original_likelihood$training_task$npsem[[actual_col]]$censoring_node$name
            )
            if_obs[is.na(if_obs)] <- 0
            cur_obs <- cf_task_treatment_ori$get_tmle_node(actual_col)
            if_obs <- if_obs == 1 & !is.na(cur_obs)
            temp_val <- rep(NA, length(if_obs))
            temp_val[if_obs] <- temp_pred
        } else temp_val <- temp_pred
        cf_dimr_data_control[, (each_col) := temp_val]
    }
    for (each_col in cols_ig) {
        temp_pred <- initial_gradient_likelihood$get_likelihood(cf_task_treatment_ig, node = each_col, fold_number = "validation")
        cf_likelihood_treatment_ig$get_likelihood(cf_task_treatment_ig, node = each_col)

        initial_gradient_likelihood$get_likelihood(tmle_task = initial_gradient_likelihood$training_task, each_col, "validation")
        cf_likelihood_treatment_ig$get_likelihood(cf_task_treatment_ig, each_col, fold_number = "validation")
        cf_likelihood_treatment_ig$get_likelihood(cf_task_treatment_ig, each_col)
        cf_likelihood_treatment_ig$factor_list$initial_gradient_tc_R_1$get_mean(cf_task_treatment_ig, each_col, fold_number = "validation")

        if (is(initial_gradient_likelihood$training_task$npsem[[each_col]]$censoring_node, "tmle3_Node")) {
            if_obs <- cf_task_treatment_ig$get_tmle_node(
                initial_gradient_likelihood$training_task$npsem[[each_col]]$censoring_node$name
            )
            if_obs[is.na(if_obs)] <- 0
            temp_val <- rep(0, length(if_obs))  # any censoring just sets gradient as 0
            temp_val[if_obs == 1] <- temp_pred
        } else temp_val <- temp_pred
        cf_dimr_data_treatment[, (each_col) := temp_val]

        temp_pred <- initial_gradient_likelihood$get_likelihood(cf_task_control_ig, node = each_col, fold_number = "validation")
        if (is(initial_gradient_likelihood$training_task$npsem[[each_col]]$censoring_node, "tmle3_Node")) {
            if_obs <- cf_task_control_ig$get_tmle_node(
                initial_gradient_likelihood$training_task$npsem[[each_col]]$censoring_node$name
            )
            if_obs[is.na(if_obs)] <- 0
            temp_val <- rep(0, length(if_obs))  # any censoring just sets gradient as 0
            temp_val[if_obs == 1] <- temp_pred
        } else temp_val <- temp_pred
        cf_dimr_data_control[, (each_col) := temp_val]
    }
    setDT(cf_dimr_data_treatment)
    setDT(cf_dimr_data_control)
    cf_task_dimr_control <- tmle3_Task_drop_censored$new(data = cf_dimr_data_control, npsem = tmle_task$npsem)
    cf_task_dimr_treatment <- tmle3_Task_drop_censored$new(data = cf_dimr_data_treatment, npsem = tmle_task$npsem)


    # when called in ipw helper, cf_tasks are created based on the task (can be observed task or integral expanded task)
    obs_variable_names <- colnames(obs_data)
    temp_node_names <- names(tmle_task$npsem)
    to_remove <- grep("^pred_|^initial_gradient_|^id$|^t$|^delta_", temp_node_names)
    if (length(to_remove) > 0) temp_node_names <- temp_node_names[-to_remove]
    loc_A_E <- grep("^A_E", temp_node_names)
    loc_A_C <- grep("^A_C", temp_node_names)
    loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
    # loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
    loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
    intervention_variables_loc <- map_dbl(paste0(intervention_variables, "$"), ~grep(.x, obs_variable_names))

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
            temp_p_A_E <- current_likelihood$get_likelihood(cf_task_dimr_treatment, temp_node_names[k], fold_number)  # A_E | A_C=1
            if (!is.null(tmle_task$npsem[[k]]$censoring_node$variables)) {  # if there is missingness in A_E
                temp_full <- if_A_E_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k]]$censoring_node$name)
                if (!is.logical(if_A_E_observed)) if_A_E_observed <- if_A_E_observed == 1  # in case it is a binary node
                if_A_E_observed[is.na(if_A_E_observed)] <- F
                if_A_E_observed <- if_A_E_observed & !is.na(cf_task_dimr_treatment$get_tmle_node(temp_node_names[k]))  # only AC nodes need this
                temp_full[if_A_E_observed] <- temp_p_A_E
                temp_full[!if_A_E_observed] <- NA
                temp_p_A_E <- temp_full
            }
            k_A_C <- loc_A_C[loc_A_C < k]
            k_A_C <- k_A_C[which.min(abs(k_A_C - k))]  # always let censoring node to lead each time point
            temp_p_A_C <- current_likelihood$get_likelihoods(cf_task_dimr_treatment, temp_node_names[k_A_C], fold_number)  # A_C=1; to be aligned
            if (!is.null(tmle_task$npsem[[k_A_C]]$censoring_node$name)) {  # if there is missingness in A_E
                temp_full <- if_A_C_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$censoring_node$name)
                if (!is.logical(if_A_C_observed)) if_A_C_observed <- if_A_C_observed == 1  # in case it is a binary node
                if_A_C_observed[is.na(if_A_C_observed)] <- F
                if_A_C_observed <- if_A_C_observed & !is.na(cf_task_dimr_treatment$get_tmle_node(k_A_C))  # only AC nodes need this
                temp_full[if_A_C_observed] <- temp_p_A_C
                temp_full[!if_A_C_observed] <- NA
                temp_p_A_C <- temp_full
            }
            return(temp_p_A_C * temp_p_A_E)
        }) %>% pmap_dbl(prod)  # this is the likelihood of being 1
        part_Z <- lapply(loc_Z_needed, function(k) {
            temp_p <- current_likelihood$get_likelihoods(cf_task_dimr_control, temp_node_names[k], fold_number) /
                current_likelihood$get_likelihoods(cf_task_dimr_control, temp_node_names[k], fold_number)
            if (!is.null(tmle_task$npsem[[k]]$censoring_node$name)) {  # if there is missingness
                temp_full <- if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k]]$censoring_node$name)
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
        if(!is.null(tmle_task$npsem[[temp_ind]]$censoring_node$name)) {
            if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[temp_ind]]$censoring_node$name)
            if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
            if_observed[is.na(if_observed)] <- F
            if_observed <- if_observed & !is.na(tmle_task$get_tmle_node(temp_ind))
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
            temp_p_A_E <- current_likelihood$get_likelihoods(cf_task_dimr_control, temp_node_names[k], fold_number)  # A_E | A_C=1
            if (!is.null(tmle_task$npsem[[k]]$censoring_node$name)) {  # if there is missingness in A_E
                temp_full <- if_A_E_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k]]$censoring_node$name)
                if (!is.logical(if_A_E_observed)) if_A_E_observed <- if_A_E_observed == 1  # in case it is a binary node
                if_A_E_observed[is.na(if_A_E_observed)] <- F
                if_A_E_observed <- if_A_E_observed & !is.na(cf_task_dimr_control$get_tmle_node(temp_node_names[k]))  # only AC nodes need this
                temp_full[if_A_E_observed] <- temp_p_A_E
                temp_full[!if_A_E_observed] <- NA
                temp_p_A_E <- temp_full
            }
            k_A_C <- loc_A_C[loc_A_C < k]
            k_A_C <- k_A_C[which.min(abs(k_A_C - k))]  # always let censoring node to lead each time point
            temp_p_A_C <- current_likelihood$get_likelihoods(cf_task_dimr_treatment, temp_node_names[k_A_C], fold_number)  # A_C=1; to be aligned
            if (!is.null(tmle_task$npsem[[k_A_C]]$censoring_node$name)) {  # if there is missingness in A_E
                temp_full <- if_A_C_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$censoring_node$name)
                if (!is.logical(if_A_C_observed)) if_A_C_observed <- if_A_C_observed == 1  # in case it is a binary node
                if_A_C_observed[is.na(if_A_C_observed)] <- F
                if_A_C_observed <- if_A_C_observed & !is.na(cf_task_dimr_control$get_tmle_node(k_A_C))  # only AC nodes need this
                temp_full[if_A_C_observed] <- temp_p_A_C
                temp_full[!if_A_C_observed] <- NA
                temp_p_A_C <- temp_full
            }
            return(temp_p_A_C * temp_p_A_E)
        }) %>% pmap_dbl(prod)  # this is the likelihood of being 1
        part_RLY <- lapply(loc_RLY_needed, function(k) {
            temp_p <- current_likelihood$get_likelihoods(cf_task_dimr_treatment, temp_node_names[k], fold_number) /
                current_likelihood$get_likelihoods(cf_task_dimr_control, temp_node_names[k], fold_number)
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
            if_observed <- if_observed & !is.na(tmle_task$get_tmle_node(temp_ind))
            temp_vec <- temp_vec[if_observed]
        }
        list_H[[temp_ind]] <- temp_vec
    }
    return(list_H)
}





ipw_middle_survival_dimr <- function(task, lik, ipw_args, fold_number){

    # cf_likelihood_control = ipw_args$cf_likelihood_control
    # cf_likelihood_treatment = ipw_args$cf_likelihood_treatment

    intervention_list_treatment <- ipw_args$intervention_list_treatment
    intervention_list_control <- ipw_args$intervention_list_control

    # # # todo: extend for stochastic
    # cf_task_treatment <- cf_likelihood_treatment$enumerate_cf_tasks(task)[[1]]
    # cf_task_control <- cf_likelihood_control$enumerate_cf_tasks(task)[[1]]

    intervention_nodes <- union(names(intervention_list_treatment), names(intervention_list_control))

    temp_node_names <- names(task$npsem)
    to_remove <- grep("^pred_|^initial_gradient_|^id$|^t$|^delta_", temp_node_names)
    if (length(to_remove) > 0) temp_node_names <- temp_node_names[-to_remove]

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

    # # get list of all possible predicted lkds
    # obs_data <- task$data %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))

    # dimr version: can be either original dimr task or any partial task
    obs_data <- task$data %>% copy
    cols_to_delete <- names(obs_data)[grep("^pred_|^initial_gradient_|^id$|^t$|^delta_", names(obs_data))]
    if (length(cols_to_delete) > 0) obs_data[, (cols_to_delete) := NULL]

    obs_variable_names <- colnames(obs_data)
    # ZW todo: to handle long format and wide format
    # ZW todo: see if observed_likelihood needs to change to targeted likelihood

    intervention_variables <- map_chr(task$npsem[intervention_nodes], ~.x$variables)
    intervention_variables_loc <- map_dbl(paste0(intervention_variables, "$"), ~grep(.x, obs_variable_names))
    intervention_levels_treat <- map_dbl(intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
    intervention_levels_control <- map_dbl(intervention_list_control, ~.x$value %>% as.character %>% as.numeric)

    list_H <- get_obs_H_full_survival_dimr(tmle_task = task, obs_data = obs_data, current_likelihood = ipw_args$static_likelihood,
                                      # current_likelihood = lik,
                                      # cf_task_treatment, cf_task_control,
                                      ipw_args$cf_likelihood_treatment_ori, ipw_args$cf_likelihood_control_ori, ipw_args$cf_likelihood_treatment_ig, ipw_args$cf_likelihood_control_ig,
                                      which_intervention_keep = ipw_args$which_intervention_keep,
                                      intervention_list_treatment, intervention_list_control,
                                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                                      original_likelihood = ipw_args$original_likelihood, initial_gradient_likelihood = ipw_args$initial_gradient_likelihood,   # ori and ig lkds are needed now
                                      fold_number = fold_number
                                      # , bound = 0.05
    )



    list_newH <- list()
    for (ind_var in 1:length(list_H)) {
        if(!is.null(list_H[[ind_var]])) {
            # if there is missingness, match the get_regression_task structure
            if (!is.null(task$npsem[[ind_var]]$censoring_node$variables)) {
                if_observed <- task$get_tmle_node(task$npsem[[ind_var]]$censoring_node$name)  # force Y 0 where censored due to bivariate trt nodes
                if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
                if_observed[is.na(if_observed)] <- F
                # temp_obs <- task$get_tmle_node(ind_var)
                # if_observed <- if_observed & !is.na(temp_obs)
            }
            if (ind_var %in% loc_Z) list_newH[[ind_var]] <- (list_H[[ind_var]] * Y[if_observed]) %>% as.matrix
            if (ind_var %in% loc_RLY) list_newH[[ind_var]] <- (list_H[[ind_var]] * Y[if_observed]) %>% as.matrix
        }
    }
    names(list_newH) <- temp_node_names

    return(list_newH)

    # ZW todo: for categorical variables
}









gradient_generator_middle_survival_dimr <- function(tmle_task, lik,  node, include_outcome = T, ipw_args = NULL, fold_number){

    task <- tmle_task$get_regression_task(node)
    if (include_outcome) {
        IC <- ipw_middle_survival_dimr(tmle_task, lik,  ipw_args, fold_number)[[node]] %>% as.vector
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






# data_from_task <- all_possible_RZLY_1
transform_original_data_to_augmented_data <- function(data_from_task, original_likelihood, initial_gradient_likelihood) {
    # get original node names
    temp_node_names <- names(original_likelihood$training_task$npsem)
    to_remove <- grep("^pred_|^initial_gradient_|^id$|^t$|^delta_", temp_node_names)
    if (length(to_remove) > 0) temp_node_names <- temp_node_names[-to_remove]

    loc_A_C <- grep("^A_C", temp_node_names)
    loc_A_E <- grep("^A_E", temp_node_names)
    loc_Z <- which(sapply(temp_node_names, function(s) paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") == "Z"))
    loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
    if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
    loc_LR <- seq_along(temp_node_names)[-c(loc_A_E, loc_A_C, loc_Z)] %>% tail(-1)  # except for baseline

    A_E_nodes <- grep("^A_E_", temp_node_names, value = T)
    Z_nodes <- grep("^Z_", temp_node_names, value = T)
    RLY_nodes <- temp_node_names[loc_RLY]


    node_name_list <- original_likelihood$training_task$npsem %>% names
    original_npsem <- original_likelihood$training_task$npsem
    # make conditional prob predictions with dummy task
    dummy_original_task <- tmle3_Task_drop_censored$new(data_from_task, original_npsem,
                                                        id = "id", t = "t",
                                                        long_format = original_likelihood$training_task$long_format)
    # make gradient prediction with dummy task
    dummy_gradient_data <- data_from_task %>% copy
    dummy_variable_names <- names(initial_gradient_likelihood$training_task$data)
    dummy_variable_names <- dummy_variable_names[!(dummy_variable_names %in% names(data_from_task))]
    for (variable in dummy_variable_names) {
        dummy_gradient_data[, (variable) := 0]
        dummy_gradient_data[, (variable) := 0]
    }
    dummy_gradient_task <- tmle3_Task_drop_censored$new(data = dummy_gradient_data, npsem = initial_gradient_likelihood$training_task$npsem, time = "t", id = "id")

    augmented_data <- data_from_task %>% copy

    # for each node, add predicted value at 1
    for (node in node_name_list[-1]) {
        temp_pred <- original_likelihood$get_likelihood(dummy_original_task, node = node, fold_number = "validation")
        if (is.null(dummy_original_task$npsem[[node]]$censoring_node$name)) {
            if_observed <- rep(T, dummy_original_task$nrow)
        } else {
            if_observed <- dummy_original_task$get_tmle_node(dummy_original_task$npsem[[node]]$censoring_node$name) == 1
        }
        if_observed[is.na(if_observed)] <- F
        observed <- dummy_original_task$get_tmle_node(node)
        if_observed <- if_observed & !is.na(observed)
        temp_pred[observed[if_observed] == 0] <- 1 - temp_pred[observed[if_observed] == 0]  # make pred for 1
        temp_vec <- rep(NA, dummy_original_task$nrow)
        temp_vec[if_observed] <- temp_pred
        augmented_data[, (paste0("pred_", 1, "_", node)) := temp_vec]
    }

    # for each node, add needed gradient prediction
    option <- strsplit(dummy_variable_names[1], "_")[[1]] %>% last

    # for (option in mediation_options) {
        for (s in loc_Z) {
            temp_pred <- initial_gradient_likelihood$get_likelihood(tmle_task = dummy_gradient_task, node = paste0("initial_gradient_", option, "_", node_name_list[s]), fold_number = "validation")
            node <- node_name_list[s]
            if (is.null(dummy_original_task$npsem[[node]]$censoring_node$name)) {
                if_observed <- rep(T, dummy_original_task$nrow)
            } else {
                if_observed <- dummy_original_task$get_tmle_node(dummy_original_task$npsem[[node]]$censoring_node$name) == 1
            }
            if_observed[is.na(if_observed)] <- F
            observed <- dummy_original_task$get_tmle_node(node)
            if_observed <- if_observed & !is.na(observed)
            # temp_pred[observed[if_observed] == 0] <- 1 - temp_pred[observed[if_observed] == 0]  # gradient pred is already for 1
            temp_vec <- rep(NA, dummy_original_task$nrow)
            temp_vec[if_observed] <- temp_pred
            augmented_data[, (paste0("initial_gradient_", option, "_", node)) := temp_vec]
        }
        for (s in loc_LR) {
            temp_pred <- initial_gradient_likelihood$get_likelihood(tmle_task = dummy_gradient_task, node = paste0("initial_gradient_", option, "_", node_name_list[s]), fold_number = "validation")
            node <- node_name_list[s]
            if (is.null(dummy_original_task$npsem[[node]]$censoring_node$name)) {
                if_observed <- rep(T, dummy_original_task$nrow)
            } else {
                if_observed <- dummy_original_task$get_tmle_node(dummy_original_task$npsem[[node]]$censoring_node$name) == 1
            }
            if_observed[is.na(if_observed)] <- F
            observed <- dummy_original_task$get_tmle_node(node)
            if_observed <- if_observed & !is.na(observed)
            # temp_pred[observed[if_observed] == 0] <- 1 - temp_pred[observed[if_observed] == 0]  # gradient pred is already for 1
            temp_vec <- rep(NA, dummy_original_task$nrow)
            temp_vec[if_observed] <- temp_pred
            augmented_data[, (paste0("initial_gradient_", option, "_", node)) := temp_vec]
        }
    # }

    return(augmented_data)
}



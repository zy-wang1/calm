#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_mediation_exact_survival <- R6Class(
    classname = "Param_mediation",
    portable = TRUE,
    class = TRUE,
    inherit = Param_base,
    public = list(
        initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome = NULL) {
            super$initialize(observed_likelihood, list())
            
            # outcome is used to decide gcomp and clever covariate; default is the last node
            temp_node_names <- names(observed_likelihood$training_task$npsem)
            loc_delta_nodes <- grep("^delta_", temp_node_names)
            if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
            if (is.null(outcome)) private$.outcome_node <- last(temp_node_names) else private$.outcome_node <- outcome
            
            private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
            private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
            # observed_likelihood$get_likelihoods(observed_likelihood$training_task)
        },
        clever_covariates = function(tmle_task = NULL, fold_number = "full", update = T, node = NULL, submodel_type = "EIC") {
            if (is.null(tmle_task)) {  # calculate for obs data task if not specified
                tmle_task <- self$observed_likelihood$training_task
            }
            intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
            
            if (fold_number == "full") {  # tmle
                list_EIC <- private$.list_EIC
            } else if (fold_number == "validation") {  # cvtmle
                list_EIC <- private$.list_EIC_val
            }  # load cached obs task clever covariates in case its for convergence check
            
            if (!is.null(list_EIC) & update == F & identical(tmle_task, self$observed_likelihood$training_task)) {  # for faster convergence check
                if (!is.null(node)) {  # return partial list of covariates if requested
                    return(list_EIC[node])
                } else {
                    return(list_EIC)
                }
            } else {  # note submodel_type; only calculate when i) no cached EIC, ii) forced to update after tlik is updated; or iii) not obs task, such as cf tasks
                rm(list_EIC)
                
                # load full_p list first
                full_task <- self$observed_likelihood$training_task
                full_node_names <- names(full_task$npsem)
                loc_delta_nodes <- grep("delta_", full_node_names)
                if (length(loc_delta_nodes) != 0) full_node_names <- full_node_names[-grep("delta_", full_node_names)]  # remove delta nodes for wide format fitting
                full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))  # exactly the obs data
                full_variable_names <- colnames(full_data)
                list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
                    if (loc_node > 1) {
                        current_variable <- full_task$npsem[[loc_node]]$variables
                        loc_impute <- grep("Y_|A_C_", full_node_names)  # remain alive and uncensored before current variable
                        loc_impute <- loc_impute[loc_impute < loc_node]
                        if (length(loc_impute) == 0) {  # no subject can drop out/die yet
                            temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
                        } else {
                            temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)],
                                                        rule_variables = sapply(loc_impute, function(s) full_task$npsem[[s]]$variables),
                                                        rule_values = rep(1, length(loc_impute))
                            )  # all possible inputs
                        }
                        delta_vars <- names(sapply(paste0("delta_", full_node_names[1:loc_node]), function(x) grep(x, names(full_task$npsem))) %>% compact %>% unlist)
                        if (length(delta_vars) > 0) {
                            temp_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
                            colnames(temp_input)[(ncol(temp_input) - length(delta_vars) + 1):ncol(temp_input)] <- delta_vars
                        }
                        
                        temp_task <- tmle3_Task$new(temp_input, full_task$npsem[c(1:loc_node,
                                                                                  sapply(paste0("delta_", full_node_names[1:loc_node]), function(x) grep(x, names(full_task$npsem))) %>% compact %>% unlist
                        )])
                        temp_target_node <- intersect(self$update_nodes, full_node_names[loc_node])
                        if (length(temp_target_node) == 1) {
                            # for each short task, only the last node (if it is an update_node) needs to be updated
                            setattr(temp_task, "target_nodes", temp_target_node)
                            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
                        } else {
                            # A nodes won't get updated
                            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
                        }
                        data.frame(temp_input[1:which(full_variable_names == current_variable)], output = temp_output) %>% return
                    }
                })
                names(list_all_predicted_lkd) <- full_node_names
                
                if (all(tmle_task$nrow == self$observed_likelihood$training_task$nrow,
                        identical(tmle_task$data[[1]], self$observed_likelihood$training_task$data[[1]])
                )) {  # for cf or obs tasks
                    # ZW todo: extend for dynamic treatments
                    cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
                    cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
                    
                    temp_node_names <- names(tmle_task$npsem)
                    loc_delta_nodes <- grep("delta_", temp_node_names)
                    if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
                    
                    loc_A_E <- grep("A_E", temp_node_names)
                    loc_A_C <- grep("A_C", temp_node_names)
                    loc_Z <- which(sapply(temp_node_names, function(s) paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") == "Z"))
                    loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "A", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
                    
                    obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))  # note this is compatible if tmle_task is a cf task
                    obs_variable_names <- colnames(obs_data)
                    # ZW todo: to handle long format and wide format
                    
                    intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
                    intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
                    intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
                    intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
                    names(intervention_levels_treat) <- names(self$intervention_list_treatment)
                    names(intervention_levels_control) <- names(self$intervention_list_control)
                    
                    list_H <- get_obs_H_list(tmle_task, obs_data, current_likelihood = self$observed_likelihood,
                                             cf_task_treatment, cf_task_control,
                                             intervention_variables, intervention_levels_treat, intervention_levels_control,
                                             fold_number)
                    list_Q <- get_obs_Q_list(tmle_task, obs_data,
                                             intervention_variables, intervention_levels_treat, intervention_levels_control,
                                             list_all_predicted_lkd  # val version decided above for fold_number == "validation"
                    )
                    temp_vec <- tmle_task$get_tmle_node(length(list_Q))
                    temp_vec <- temp_vec[!is.na(temp_vec)]
                    list_Q[[length(list_Q)+1]] <- temp_vec
                    list_delta_Q <- lapply(1:length(list_H), function(i) {
                        if (is.null(list_Q[[i]]))
                            return(NULL)
                        else {
                            temp_i_plus <- first(which(!sapply(list_Q[(i+1):length(list_Q)], is.null)))  # search for the first non-null loc after i
                            if (length(list_Q[[i+temp_i_plus]]) != length(list_Q[[i]])) {  # length not equal means there is a new censoring node
                                loc_A_C_used <- loc_A_C[loc_A_C<i]
                                temp_ind <- tmle_task$get_tmle_node(last(loc_A_C_used))[
                                    tmle_task$get_tmle_node(loc_A_C_used[length(loc_A_C_used) - 1]) == 1
                                ] == 1
                                temp_ind[is.na(temp_ind)] <- F
                                temp_delta_Q <- list_Q[[i+temp_i_plus]] - list_Q[[i]][temp_ind]
                            } else {
                                temp_delta_Q <- list_Q[[i+temp_i_plus]] - list_Q[[i]]
                            }
                            return(temp_delta_Q)
                        }
                    })
                    list_EIC <- lapply(1:length(list_H), function(i) {
                        if (is.null(list_H[[i]])) return(NULL) else
                            return(list_H[[i]]*list_delta_Q[[i]])
                    })
                    names(list_EIC) <- temp_node_names
                    
                    list_EIC_inserted <- lapply(1:length(list_EIC), function(i) {
                        if (!is.null(list_EIC[[i]])) {
                            if (length(list_EIC[[i]] != nrow(obs_data))) {  # fill back 0 indicators in EIC
                                temp_vec <- vec_if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[i]]$censoring_node$name)
                                if (!is.logical(vec_if_observed)) vec_if_observed <- vec_if_observed == 1  # in case it is a binary node
                                vec_if_observed[is.na(vec_if_observed)] <- F
                                temp_vec[vec_if_observed] <- list_EIC[[i]]
                                temp_vec[!vec_if_observed] <- 0
                                return(temp_vec)
                            } else {
                                return(list_EIC[[i]])
                            }
                        }
                    })
                    
                    # last column might be needed for some tmle update functions
                    list_EIC[[length(list_EIC) + 1]] <- do.call(cbind, list_EIC_inserted)
                    names(list_EIC)[length(list_EIC)] <- "IC"  # to use in by dimension convergence
                    
                    if (identical(tmle_task, self$observed_likelihood$training_task)) {  # cache for obs task
                        if (fold_number == "full") {
                            private$.list_EIC <- list_EIC
                        } else if (fold_number == "validation") {
                            private$.list_EIC_val <- list_EIC
                        }
                    }
                    
                    if (!is.null(node)) {  # return partial list of covariates if requested
                        return(list_EIC[node])
                    } else
                        return(list_EIC)
                } else {  # for library tasks; it's only needed in tlik updates, with single node
                    if (is.null(node)) stop("Please specify single update node for library tasks")
                    
                    tmle_task_backup <- tmle_task
                    tmle_task <- self$observed_likelihood$training_task  # let tmle_task be obs task when calculating for library tasks
                    loc_node <- which(names(tmle_task$npsem) == node)
                    obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))  # note this is compatible if tmle_task is a cf task
                    obs_variable_names <- names(obs_data)
                    intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
                    intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
                    intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
                    intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
                    names(intervention_levels_treat) <- names(self$intervention_list_treatment)
                    names(intervention_levels_control) <- names(self$intervention_list_control)
                    
                    current_H <- get_current_H(loc_node,
                                               tmle_task, obs_variable_names,
                                               intervention_variables, intervention_levels_treat, intervention_levels_control,
                                               list_all_predicted_lkd  # this is decided above by fold_number
                    )  # this is what we need for logistic submodel
                    current_Q_next <- get_current_Q(loc_node, which_Q = 1,
                                                    tmle_task, obs_variable_names,
                                                    intervention_variables, intervention_levels_treat, intervention_levels_control,
                                                    list_all_predicted_lkd,  # this is decided above by fold_number
                                                    if_survival = T
                    )
                    current_Q <- get_current_Q(loc_node, which_Q = 0,
                                               tmle_task, obs_variable_names,
                                               intervention_variables, intervention_levels_treat, intervention_levels_control,
                                               list_all_predicted_lkd,  # this is decided above by fold_number
                                               if_survival = T
                    )
                    current_delta_Q <- current_Q_next - current_Q
                    current_EIC <- current_H*current_delta_Q
                    
                    current_EIC <- list(current_EIC)
                    names(current_EIC) <- node
                    
                    return(current_EIC)
                }
                
            }
        },
        estimates = function(tmle_task = NULL, fold_number = "full", update = T) {
            if (is.null(tmle_task)) {
                tmle_task <- self$observed_likelihood$training_task
            }
            intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
            
            temp_node_names <- names(tmle_task$npsem)
            loc_delta_nodes <- grep("delta_", temp_node_names)
            if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
            loc_A <- grep("^A_", temp_node_names)  # not used here; it can include both A_E and A_C
            loc_A_C <- grep("^A_C_", temp_node_names)
            loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
            loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))
            loc_Y <- grep("^Y_", temp_node_names)
            if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
            
            loc_outcome <- which(temp_node_names == self$outcome_node)
            var_outcome <- tmle_task$npsem[[loc_outcome]]$variables
            
            obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))  # note this is compatible if tmle_task is a cf task
            obs_variable_names <- colnames(obs_data)
            # ZW todo: to handle long format and wide format
            
            if (fold_number == "full") {
                list_EIC <- private$.list_EIC
                result <- private$.result
            } else if (fold_number == "validation") {
                list_EIC <- private$.list_EIC_val
                result <- private$.result_val
            }
            
            # load full_p list for survival
            full_task <- self$observed_likelihood$training_task
            full_node_names <- names(full_task$npsem)
            loc_delta_nodes <- grep("delta_", full_node_names)
            if (length(loc_delta_nodes) != 0) full_node_names <- full_node_names[-grep("delta_", full_node_names)]  # remove delta nodes for wide format fitting
            full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))  # note this is compatible if tmle_task is a cf task
            full_variable_names <- colnames(full_data)
            list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
                if (loc_node > 1) {
                    current_variable <- tmle_task$npsem[[loc_node]]$variables
                    loc_impute <- grep("^Y_|^A_C_", full_node_names)  # remain alive and uncensored before current variable
                    if (length(loc_A_C) > 0) loc_impute <- loc_impute[loc_impute <= last(loc_A_C)]  # the event after the last censoring node can be alive/dead; do not impute
                    loc_impute <- loc_impute[loc_impute < loc_node]
                    if (length(loc_impute) == 0) {  # no subject can drop out/die yet
                        temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
                    } else {
                        temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)],
                                                    rule_variables = sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables),
                                                    rule_values = rep(1, length(loc_impute))
                        )  # all possible inputs
                    }
                    delta_vars <- names(sapply(paste0("delta_", full_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
                    if (length(delta_vars) > 0) {
                        temp_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
                        colnames(temp_input)[(ncol(temp_input) - length(delta_vars) + 1):ncol(temp_input)] <- delta_vars
                    }
                    
                    temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[c(1:loc_node,
                                                                              sapply(paste0("delta_", full_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
                    )])
                    temp_target_node <- intersect(self$update_nodes, full_node_names[loc_node])
                    if (length(temp_target_node) == 1) {
                        # for each short task, only the last node (if it is an update_node) needs to be updated
                        setattr(temp_task, "target_nodes", temp_target_node)
                        temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
                    } else {
                        # A nodes won't get updated
                        temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
                    }
                    data.frame(temp_input[1:which(full_variable_names == current_variable)], output = temp_output) %>% return
                }
            })
            names(list_all_predicted_lkd) <- full_node_names
            
            # make sure we force update this after each updating step; this helps speed up updater$check_convergence
            if (!is.null(list_EIC) & !is.null(result) & update == F) {
                return(result)
            } else {
                intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
                intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
                intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
                intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
                
                loc_impute <- grep("^Y_|^A_C_", temp_node_names)  # in integral, Y and A_C always 1
                if (length(loc_A_C) > 0) loc_impute <- loc_impute[loc_impute <= last(loc_A_C)]  # the event after the last censoring node can be alive/dead; do not impute
                
                # nodes to integrate out in the target identification
                # only support univaraite node for now; assume treatment level is one
                all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                                     rule_variables = c(intervention_variables,
                                                                        sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                                     rule_values = c(intervention_levels_treat,
                                                                     rep(1, length(loc_impute))))
                all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                                     rule_variables = c(intervention_variables,
                                                                        sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                                     rule_values = c(intervention_levels_control, 1,
                                                                     rep(1, length(loc_impute))))
                
                # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
                unique_L0 <- obs_data[, tmle_task$npsem[[1]]$variables] %>% unique
                library_L0 <- data.frame(unique_L0, output =
                                             map_dbl(1:nrow(unique_L0), function(which_row) {
                                                 temp_all_comb_0 <- cbind(unique_L0[which_row, ], all_possible_RZLY_0, row.names = NULL) %>% suppressWarnings
                                                 temp_all_comb_1 <- cbind(unique_L0[which_row, ], all_possible_RZLY_1, row.names = NULL) %>% suppressWarnings
                                                 # for all non-A, non-0 variables, calculate the variable by rule
                                                 # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                                 # note that list_all_predicted_lkd is ordered by node
                                                 temp_list_0 <- lapply(loc_Z,
                                                                       function(each_t) {
                                                                           left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output %>% suppressMessages
                                                                       })
                                                 temp_list_1 <- lapply(loc_RLY,
                                                                       function(each_t) {
                                                                           left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output %>% suppressMessages
                                                                       })
                                                 temp_list <- c(temp_list_0, temp_list_1)
                                                 pmap_dbl(temp_list, prod)[temp_all_comb_1[[var_outcome]] == 1] %>% sum %>% return  # only sum up survival term
                                             })
                )
                # substitution estimator
                vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output %>% suppressMessages
                psi <- mean(vec_est)
                
                list_EIC <- self$clever_covariates(tmle_task, fold_number, submodel_type = "EIC")$"IC"
                EIC <- cbind(list_EIC, vec_est - psi)
                IC <- rowSums(EIC)
                
                result <- list(psi = psi, IC = IC
                               # , full_IC = list_EIC
                )
                
                # these are cached; unless likelihood is updated, or we force it to update, they shouldn't be changed
                if (fold_number == "full") {
                    private$.result <- result
                } else if (fold_number == "validation") {
                    private$.result_val <- result
                }
                
                return(result)
            }
        }
    ),
    active = list(
        name = function() {
            param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
            return(param_form)
        },
        cf_likelihood_treatment = function() {
            return(private$.cf_likelihood_treatment)
        },
        cf_likelihood_control = function() {
            return(private$.cf_likelihood_control)
        },
        intervention_list_treatment = function() {
            return(self$cf_likelihood_treatment$intervention_list)
        },
        intervention_list_control = function() {
            return(self$cf_likelihood_control$intervention_list)
        },
        update_nodes = function() {
            # if (is.null(tmle_task)) {
            tmle_task <- self$observed_likelihood$training_task
            # }
            temp_node_names <- names(tmle_task$npsem)
            loc_delta_nodes <- grep("delta_", temp_node_names)
            if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
            loc_A <- grep("A", sapply(strsplit(temp_node_names, "_"), function(x) x[1]))  # A_E or A_C
            if_not_0 <- sapply(strsplit(temp_node_names, "_"), function(x) last(x) != 0)
            nodes_to_update <- temp_node_names[if_not_0 & !((1:length(temp_node_names)) %in% loc_A)]
            # nodes_to_update <- rev(nodes_to_update)
            # nodes_to_update <- nodes_to_update[-length(nodes_to_update)]
            return(nodes_to_update)
        },
        list_EIC = function() {
            return(private$.list_EIC)
        },
        list_EIC_val = function() {
            return(private$.list_EIC_val)
        }, 
        outcome_node = function() {
            return(private$.outcome_node)
        }
    ),
    private = list(
        .type = "mediation_survival",
        .cf_likelihood_treatment = NULL,
        .cf_likelihood_control = NULL,
        .list_EIC = NULL,  # the clever covariates as the EIC
        .list_EIC_val = NULL,
        .result = NULL,
        .result_val = NULL,
        .submodel_type_supported = c("EIC"), 
        .outcome_node
    )
)
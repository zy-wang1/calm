#' Longitudinal Mediation via HAL-EIC TMLE

#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'







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
Param_mediation_hal_survival <- R6Class(
    classname = "Param_mediation_projection",
    portable = TRUE,
    class = TRUE,
    inherit = Param_base,
    public = list(
        initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, static_likelihood = NULL, n_resampling = NULL,
                              outcome = NULL, tau = NULL) {
            if(inherits(observed_likelihood, "Targeted_Likelihood")){
                fold_number <- observed_likelihood$updater$update_fold
            } else {
                fold_number <- "full"
            }
            
            super$initialize(observed_likelihood, list())
            
            temp_node_names <- names(observed_likelihood$training_task$npsem)
            loc_delta_nodes <- grep("^delta_", temp_node_names)
            if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
            loc_A_E <- grep("^A_E", temp_node_names)
            loc_Z <- which(sapply(temp_node_names, function(s) paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") == "Z"))
            loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
            if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
            # tau <- as.numeric(last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2])))
            
            
            A_E_nodes <- grep("^A_E_", temp_node_names, value = T)
            Z_nodes <- grep("^Z_", temp_node_names, value = T)
            RLY_nodes <- temp_node_names[loc_RLY]
            # grep("(R|L|Y).[1-9]$", temp_node_names, value = T)
            
            if (is.null(outcome)) private$.outcome_node <- last(temp_node_names) else private$.outcome_node <- outcome
            if (is.null(tau)) tau <- strsplit(self$outcome_node, "_")[[1]] %>% last %>% as.numeric  # tau is the last time point involved in the outcome
            loc_last_node <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]] %>% last %>% as.numeric) <= tau) %>% last
            last_node <- temp_node_names[loc_last_node]
            
            private$.static_likelihood <- static_likelihood
            update_nodes <- c(Z_nodes, RLY_nodes)
            update_nodes <- update_nodes[sapply(strsplit(update_nodes, "_"), function(u) as.numeric(u[[2]]) <= tau)]  # nodes after tau won't be used
            private$.update_nodes <- update_nodes
            
            private$.supports_outcome_censoring <- T
            # super$initialize(observed_likelihood, list())
            
            
            which_intervention_keep <- sapply(intervention_list_treatment, function(u) which(u$name == temp_node_names)) <= loc_last_node
            
            private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment[which_intervention_keep])
            private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control[which_intervention_keep])
            
            # private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
            # private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
            
            
            tmle_task <- observed_likelihood$training_task
            obs_data <- tmle_task$data %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))
            obs_variable_names <- colnames(obs_data)
            obs_variable_values <- lapply(obs_data, function(eachCol) sort(unique(eachCol[!is.na(eachCol)])))
            loc_A_C_or_Y <- grep("^A_C_|^Y_", temp_node_names)
            loc_Y <- grep("^Y_", temp_node_names)
            
            if (!is.null(n_resampling)) {  # use expanded Monte Carlo samples to train the HAL projection
                
                # for (u in 1:2) {
                #   observed_task <- tmle_task
                #   tlik <- static_likelihood
                #   num_samples = n_resampling * 2^(u-1)
                #   baseline_node = temp_node_names[1]
                #   time_ordering = temp_node_names
                #
                #   mc_data <- as.data.table(tlik$factor_list[[baseline_node]]$sample(, num_samples))
                #   baseline_var <- observed_task$npsem[[baseline_node]]$variables
                #   set(mc_data, , setdiff(names(observed_task$data), baseline_var), (0))
                #   mc_data$t <- 0
                #   mc_data$id <- 1:nrow(mc_data)
                #   mc_task <- tmle3_Task$new(mc_data, observed_task$npsem, id = "id", t = "t", long_format  = observed_task$long_format)
                #   time_ordering <- setdiff(time_ordering, baseline_node)
                #   for(node in time_ordering) {
                #     if (length(grep("^A_C_", node)) == 0) {
                #       if (is(mc_task$npsem[[node]]$censoring_node, "tmle3_Node")) {
                #         sampled <- sample_from_node(mc_task, tlik, node, observed_task, if_only_return_sampled = T)
                #         if_censored <- as.numeric(mc_task$get_tmle_node(mc_task$npsem[[node]]$censoring_node$name)) == 0
                #         mc_task <- mc_task$generate_counterfactual_task(UUIDgenerate(), sampled, force_at_risk = F)
                #         # new_data <- mc_task$data
                #         # set(new_data, which(if_censored), (node), rep(NA, sum(if_censored)))
                #         # mc_task <- tmle3_Task$new(new_data, mc_task$npsem, time = "t", id = "id")
                #       }
                #     } else {
                #       mc_task <- sample_from_node(mc_task, tlik, node, observed_task)
                #     }
                #   }
                #   new_data <- mc_task$data
                #   for (node in time_ordering) {
                #     if (length(grep("^A_C_", node)) == 0) {
                #       if (is(mc_task$npsem[[node]]$censoring_node, "tmle3_Node")) {
                #         if_censored <- as.numeric(mc_task$get_tmle_node(mc_task$npsem[[node]]$censoring_node$name)) == 0
                #         set(new_data, which(if_censored), (node), rep(NA, sum(if_censored)))
                #       }
                #     }
                #   }
                #
                #   if (length(table(new_data[[last(temp_node_names)]])) > 1) break
                # }
                #
                #
                # mc_task <- tmle3_Task_drop_censored$new(new_data, mc_task$npsem, time = "t", id = "id")
                # temp_task <- mc_task
                
                # temp_input <- tmle_task$get_tmle_node(temp_node_names[1])
                # temp_input <- temp_input[sample(nrow(temp_input), abs(round(n_resampling)), replace = T), ]
                # temp_input <- temp_input[, c("t", "id"):=list(0, 1:nrow(temp_input))] %>% setDT
                # for (i in 2:loc_last_node) {
                #   # all_levels <- tmle_task$get_tmle_node(temp_node_names[i]) %>% unique()
                #   # all_levels <- all_levels[!is.na(all_levels)] %>% sort
                #   all_levels <- c(1)
                #   prob_list <- sapply(all_levels, function(l) {
                #     temp_input <- temp_input[, temp_node_names[i]:=l]
                #     temp_task <- tmle3_Task_drop_censored$new(temp_input, tmle_task$npsem[1:i], time = "t", id = "id")
                #     # temp_input <- temp_input[, temp_node_names[i]:=NULL]
                #     return(static_likelihood$get_likelihood(temp_task, node = temp_node_names[i], fold_number = fold_number))
                #   })
                #   if (length(all_levels) == 1) {
                #     new_sample <- rbinom(length(prob_list), 1, prob_list)
                #   } else {
                #     new_sample <- prob_list %>% apply(1, function(u) {
                #       sample(all_levels, 1, replace = F, prob = u)
                #     })
                #   }
                #   if (is(tmle_task$npsem[[temp_node_names[i]]]$censoring_node, "tmle3_Node")) {
                #     temp_sample <- rep(NA, nrow(temp_input))
                #     temp_sample[temp_input[[temp_node_names[which(temp_node_names == tmle_task$npsem[[temp_node_names[i]]]$censoring_node$name)]]] == 1] <- new_sample
                #   } else {
                #     temp_sample <- new_sample
                #   }
                #   if (length(grep("^A_C_", tmle_task$npsem[[i]]$name)) == 1) temp_sample[is.na(temp_sample)] <- 0  # use 0 not NA for censored in A_C_ nodes
                #   temp_input <- temp_input[, temp_node_names[i]:=temp_sample]
                # }
                # # align event process; non increasing; only for survival targets
                # for (k in 1:length(loc_Y)) {
                #   if (k > 1) {
                #     temp_input[temp_input[[temp_node_names[loc_Y[k-1]]]] < temp_input[[temp_node_names[loc_Y[k]]]], temp_node_names[loc_Y[k]] := 0]
                #   }
                # }
                # temp_task <- tmle3_Task_drop_censored$new(temp_input, tmle_task$npsem[1:loc_last_node], id = "id", time = "t")
                
                
                temp_input <- tmle_task$get_tmle_node(temp_node_names[1])
                temp_input <- temp_input[sample(nrow(temp_input), abs(round(n_resampling)), replace = T), ]
                temp_input <- temp_input[, c("t", "id"):=list(0, 1:nrow(temp_input))] %>% setDT
                for (i in 2:loc_last_node) {
                    # all_levels <- tmle_task$get_tmle_node(temp_node_names[i]) %>% unique()
                    # all_levels <- all_levels[!is.na(all_levels)] %>% sort
                    all_levels <- c(1)
                    prob_list <- sapply(all_levels, function(l) {
                        temp_input <- temp_input[, temp_node_names[i]:=l]
                        temp_task <- tmle3_Task_drop_censored$new(temp_input, tmle_task$npsem[1:i], time = "t", id = "id")
                        # temp_input <- temp_input[, temp_node_names[i]:=NULL]
                        return(static_likelihood$get_likelihood(temp_task, node = temp_node_names[i], fold_number = fold_number))
                    })
                    if (length(all_levels) == 1) {
                        new_sample <- rbinom(length(prob_list), 1, prob_list)
                    } else {
                        new_sample <- prob_list %>% apply(1, function(u) {
                            sample(all_levels, 1, replace = F, prob = u)
                        })
                    }
                    if (is(tmle_task$npsem[[temp_node_names[i]]]$censoring_node, "tmle3_Node")) {
                        temp_sample <- rep(NA, nrow(temp_input))
                        loc_node_csr <- which(temp_node_names == tmle_task$npsem[[temp_node_names[i]]]$censoring_node$name)
                        # not currently or previously censored
                        temp_sample[!(temp_input[[temp_node_names[loc_node_csr]]] == 0 | is.na(temp_input[[temp_node_names[loc_node_csr]]]))] <- new_sample
                    } else {
                        temp_sample <- new_sample
                    }
                    # if (length(grep("^A_C_", tmle_task$npsem[[i]]$name)) == 1) temp_sample[is.na(temp_sample)] <- 0  # use 0 not NA for censored in A_C_ nodes
                    temp_input <- temp_input[, temp_node_names[i]:=temp_sample]
                }
                # align event process; non increasing; only for survival targets
                for (k in 1:length(loc_Y)) {
                    if (k > 1) {
                        # 0 is event; if previously dead, then make it dead now
                        # temp_input[temp_input[[temp_node_names[loc_Y[k-1]]]] < temp_input[[temp_node_names[loc_Y[k]]]], temp_node_names[loc_Y[k]] := 0]
                        temp_rule <- temp_input[[temp_node_names[loc_Y[k-1]]]] == 0
                        temp_rule[is.na(temp_rule)] <- F
                        temp_input[temp_rule, (temp_node_names[loc_Y[k]]) := 0]
                        for (ss in temp_node_names[(loc_Y[k-1]+1) : (loc_Y[k]-1)]) temp_input[temp_rule, (ss) := NA]
                    }
                }
                temp_task <- tmle3_Task_drop_censored$new(temp_input, tmle_task$npsem[1:loc_last_node], id = "id", time = "t")
                
                
                
                # Train the gradient
                private$.gradient <- Gradient_wide$new(observed_likelihood,
                                                       ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                                       cf_likelihood_control = self$cf_likelihood_control,
                                                                       intervention_list_treatment = self$intervention_list_treatment[which_intervention_keep],
                                                                       intervention_list_control = self$intervention_list_control[which_intervention_keep],
                                                                       static_likelihood = self$static_likelihood,
                                                                       outcome_node = self$outcome_node,
                                                                       tau = tau
                                                       ),
                                                       projection_task_generator = gradient_generator_middle_survival,
                                                       target_nodes =  self$update_nodes)
                
                private$.gradient$train_projections(temp_task, fold_number = fold_number)
            } else {
                # todo: extend for stochastic
                # private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
                # private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
                # Train the gradient
                private$.gradient <- Gradient_wide$new(observed_likelihood,
                                                       ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                                       cf_likelihood_control = self$cf_likelihood_control,
                                                                       intervention_list_treatment = self$intervention_list_treatment[which_intervention_keep],
                                                                       intervention_list_control = self$intervention_list_control[which_intervention_keep],
                                                                       # cf_task_treatment = self$cf_task_treatment,
                                                                       # cf_task_control = self$cf_task_control,
                                                                       static_likelihood = self$static_likelihood
                                                       ),
                                                       projection_task_generator = gradient_generator_middle_survival,
                                                       target_nodes =  self$update_nodes)
                
                private$.gradient$train_projections(self$observed_likelihood$training_task, fold_number = fold_number)
            }
            
            setattr(self$observed_likelihood, "target_nodes", self$update_nodes)
            for (node in temp_node_names) {
                self$observed_likelihood$get_likelihood(self$observed_likelihood$training_task, node = node, fold_number = fold_number)
            }
            # for (node in self$update_nodes) {
            #   temp_long_task <- private$.gradient$expand_task(observed_likelihood$training_task, node)
            #   self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
            #   self$observed_likelihood$get_likelihood(observed_likelihood$training_task, node, fold_number)
            #   # private$.gradient$expand_task(private$.cf_task_treatment, node)
            #   # private$.gradient$expand_task(private$.cf_task_control, node)
            # }
            
            # list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
            #   if (loc_node > 1) {
            #     # currently only support univariate node for t>0
            #     current_variable <- tmle_task$npsem[[loc_node]]$variables
            #     loc_impute <- grep("Y_|A_C_", temp_node_names)  # remain alive and uncensored before current variable
            #     loc_impute <- loc_impute[loc_impute < loc_node]
            #     if (length(loc_impute) == 0) {
            #       temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)],
            #                                   values = obs_variable_values
            #                                   )  # all possible inputs
            #     } else {
            #       temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)],
            #                                   values = obs_variable_values,
            #                                   rule_variables = sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables),
            #                                   rule_values = rep(1, length(loc_impute))
            #       )  # all possible inputs
            #     }
            #     delta_vars <- names(sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
            #     if (length(delta_vars) > 0) {
            #       temp_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
            #       colnames(temp_input)[(ncol(temp_input) - length(delta_vars) + 1):ncol(temp_input)] <- delta_vars
            #     }
            #
            #     temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[c(1:loc_node,
            #                                                               sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
            #     )])
            #     temp_target_node <- intersect(self$update_nodes, temp_node_names[loc_node])
            #     if (length(temp_target_node) == 1) {
            #       # for each short task, only the last node (if it is an update_node) needs to be updated
            #       setattr(temp_task, "target_nodes", temp_target_node)
            #       for (node in attr(temp_task, "target_nodes")) {
            #         temp_long_task <- private$.gradient$expand_task(temp_task, node)
            #         self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
            #       }
            #       temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
            #     } else {
            #       # A nodes won't get updated
            #       temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
            #     }
            #     data.frame(temp_input, output = temp_output) %>% return
            #   }
            # })
            
        },
        clever_covariates = function(tmle_task = NULL, fold_number = "full", node = NULL, for_fitting = NULL, submodel_type = "EIC") {
            if (is.null(tmle_task)) {
                tmle_task <- self$observed_likelihood$training_task
            }
            update_nodes <- intersect(self$update_nodes, attr(tmle_task, "target_nodes"))
            if(!is.null(node)){
                update_nodes <- c(node)
            }
            islong = F
            if(is.null(update_nodes)){
                update_nodes <- self$update_nodes
            } else {
                islong= T
            }
            EICs <- lapply(update_nodes, function(node){
                temp <- self$gradient$compute_component(tmle_task, node, fold_number = fold_number)$EIC
                if (length(temp) > 0) return(temp) else {
                    if (is(tmle_task$npsem[[node]]$censoring_node, "tmle3_Node")) {
                        temp <- sum(tmle_task$get_tmle_node(tmle_task$npsem[[node]]$censoring_node$name) == 1)
                        return(rep(0, temp) %>% as.matrix(ncol = 1))
                    } else {
                        return(rep(0, tmle_task$nrow) %>% as.matrix(ncol = 1))
                    }
                }
            })
            names(EICs) <- update_nodes
            
            EICs_impute <- lapply(update_nodes, function(node){
                temp_vec <- self$gradient$compute_component(tmle_task, node, fold_number = fold_number)$EIC %>% as.vector
                if (length(temp_vec) == 0) {
                    temp_vec <- rep(0, tmle_task$nrow)
                } else {
                    if (!is.null(tmle_task$npsem[[node]]$censoring_node)) {
                        full_vec <- if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[node]]$censoring_node$name) %>% as.vector
                        if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
                        if_observed <- if_observed & !is.na(if_observed)  # there might be NA in censoring node
                        if_observed <- if_observed & !is.na(tmle_task$get_tmle_node(node) %>% as.vector)  # added rule for NA censoring node value due to death
                        full_vec[if_observed] <- temp_vec
                        full_vec[!if_observed] <- 0
                        temp_vec <- full_vec
                    }
                }
                return(temp_vec)
            })
            EICs[[length(EICs) + 1]] <- do.call(cbind, EICs_impute)
            EICs[[length(EICs)]]
            colnames(EICs[[length(EICs)]]) <- update_nodes
            EICs[[length(EICs)]] <- as.data.table(EICs[[length(EICs)]])
            
            names(EICs)[length(EICs)] <- "IC"
            
            return(EICs)
        },
        estimates = function(tmle_task = NULL, fold_number = "full", est_n_samples = 5000, final = F, seed = NULL) {
            if (is.null(tmle_task)) {
                tmle_task <- self$observed_likelihood$training_task
            }
            
            intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
            
            # clever_covariates happen here (for this param) only, but this is repeated computation
            EIC <- (do.call(cbind, self$clever_covariates(tmle_task, fold_number)$IC))
            
            #TODO need to montecarlo simulate from likleihood to eval parameter.
            
            temp_node_names <- names(self$observed_likelihood$training_task$npsem)
            loc_delta_nodes <- grep("^delta_", temp_node_names)
            if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("^delta_", temp_node_names)]  # remove delta nodes for wide format fitting
            loc_A <- grep("^A_", temp_node_names)
            loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
            loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))
            if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
            tau <- last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2]))
            nodes_Z <- temp_node_names[loc_Z]
            nodes_A <- temp_node_names[loc_A]
            
            obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))
            obs_variable_names <- colnames(obs_data)
            obs_variable_values <- lapply(obs_data, function(eachCol) sort(unique(eachCol[!is.na(eachCol)])))
            # list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
            #   if (loc_node > 1) {
            #     # currently only support univariate node for t>0
            #     current_variable <- tmle_task$npsem[[loc_node]]$variables
            #     loc_impute <- grep("^Y_|^A_C_", temp_node_names)  # remain alive and uncensored before current variable
            #     loc_impute <- loc_impute[loc_impute < loc_node]
            #     if (length(loc_impute) == 0) {
            #       temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)],
            #                                   values = obs_variable_values)  # all possible inputs
            #     } else {
            #       temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)],
            #                                   values = obs_variable_values,
            #                                   rule_variables = sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables),
            #                                   rule_values = rep(1, length(loc_impute))
            #       )  # all possible inputs
            #     }
            #     delta_vars <- names(sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
            #     if (length(delta_vars) > 0) {
            #       temp_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
            #       colnames(temp_input)[(ncol(temp_input) - length(delta_vars) + 1):ncol(temp_input)] <- delta_vars
            #       temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[c(1:loc_node,
            #                                                                 sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
            #       )])
            #     } else {
            #       temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[c(1:loc_node
            #       )])
            #     }
            #
            #
            #
            #     temp_target_node <- intersect(self$update_nodes, temp_node_names[loc_node])
            #     if (length(temp_target_node) == 1) {
            #       temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
            #     } else {
            #       # A nodes won't get updated
            #       temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
            #     }
            #     data.frame(temp_input[1:which(obs_variable_names == current_variable)], output = temp_output) %>% return
            #   }
            # })
            
            intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
            intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
            intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
            intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
            # nodes to integrate out in the target identification
            # only support univaraite node for now; assume treatment level is one
            loc_impute <- grep("Y_|A_C_", temp_node_names)  # in integral, Y and A_C always 1
            # all_possible_RZLY_1 <- expand_values(variables = obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
            #                                      values = obs_variable_values,
            #                                      rule_variables = c(intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
            #                                      rule_values = c(intervention_levels_treat, rep(1, length(loc_impute)))
            # )
            # all_possible_RZLY_0 <- expand_values(variables = obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
            #                                      values = obs_variable_values,
            #                                      rule_variables = c(intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
            #                                      rule_values = c(intervention_levels_control, rep(1, length(loc_impute)))
            # )
            # # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
            # unique_L0 <- obs_data[, tmle_task$npsem[[1]]$variables] %>% unique
            # library_L0 <- data.frame(unique_L0, output =
            #                            map_dbl(1:nrow(unique_L0), function(which_row) {
            #                              temp_all_comb_0 <- cbind(unique_L0[which_row, ], all_possible_RZLY_0, row.names = NULL)
            #                              temp_all_comb_1 <- cbind(unique_L0[which_row, ], all_possible_RZLY_1, row.names = NULL)
            #                              # for all non-A, non-0 variables, calculate the variable by rule
            #                              # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
            #                              # note that list_all_predicted_lkd is ordered by node
            #                              temp_list_0 <- lapply(loc_Z,
            #                                                    function(each_t) {
            #                                                      left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output %>% suppressMessages
            #                                                    })
            #                              temp_list_1 <- lapply(loc_RLY,
            #                                                    function(each_t) {
            #                                                      left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output %>% suppressMessages
            #                                                    })
            #                              temp_list <- c(temp_list_0, temp_list_1)
            #                              pmap_dbl(temp_list, prod) %>% sum %>% return
            #                            })
            # )
            # vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output %>% suppressMessages
            # psi <- mean(vec_est)
            
            psi <- 0
            
            if(final) {
                loc_Y <- grep("^Y_", temp_node_names)  # in integral, Y always 1 (alive)
                loc_A_C <- grep("^A_C_", temp_node_names)  # in integral, A_C always 1
                
                if (!is.null(seed)) set.seed(seed)
                temp <- expand_values_size(variables = obs_variable_names,
                                           # to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                           values = obs_variable_values,
                                           rule_variables = c(intervention_variables, sapply(loc_Y, function(s) tmle_task$npsem[[s]]$variables), sapply(loc_A_C, function(s) tmle_task$npsem[[s]]$variables)),
                                           rule_values = c(intervention_levels_treat, rep(1, length(loc_Y)), rep(1, length(loc_A_C))),
                                           size = est_n_samples
                )
                all_possible_RZLY_1 <- temp[[1]] %>% as.data.table
                prep_list <- temp[[2]]
                rm(temp)
                all_possible_RZLY_1 <- rbindlist(list(all_possible_RZLY_1,
                                                      # rep(0, ncol(all_possible_RZLY_1)) %>% as.list
                                                      all_possible_RZLY_1[nrow(all_possible_RZLY_1), ]
                ))
                all_possible_RZLY_1[nrow(all_possible_RZLY_1), grep("^Y_", colnames(all_possible_RZLY_1)) := 0]
                
                # all_possible_RZLY_0 <- expand_values_size(variables = obs_variable_names,
                #                                           # to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                #                                      values = obs_variable_values,
                #                                      rule_variables = c(intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                #                                      rule_values = c(intervention_levels_control, rep(1, length(loc_impute))),
                #                                      size = est_n_samples
                # )
                all_possible_RZLY_1$t <- 0
                all_possible_RZLY_1$id <- 1:nrow(all_possible_RZLY_1)
                all_possible_RZLY_0 <- all_possible_RZLY_1 %>% copy
                for (i in 1:length(intervention_variables_loc))  all_possible_RZLY_0[, colnames(all_possible_RZLY_0)[intervention_variables_loc[i]] := intervention_levels_control[i]]
                
                mc_task_trt <- tmle3_Task_drop_censored$new(all_possible_RZLY_1, tmle_task$npsem, id = "id", t = "t", long_format  = tmle_task$long_format)
                mc_task_ctrl <- tmle3_Task_drop_censored$new(all_possible_RZLY_0, tmle_task$npsem, id = "id", t = "t", long_format  = tmle_task$long_format)
                
                
                # all_update_nodes <- self$observed_likelihood$updater$tmle_params %>% sapply(function(param) param$update_nodes) %>% unlist %>% unique
                all_update_nodes <- self$update_nodes
                # self$observed_likelihood$updater$tmle_params %>% lapply(function(param) for (node in all_update_nodes) {
                self$observed_likelihood$updater$tmle_params %>% lapply(function(param) for (node in param$update_nodes) {
                    # param$gradient$expand_task(mc_task_trt, node)
                    if (is.null(self$observed_likelihood$cache$get_update_step(self$observed_likelihood$factor_list[[node]], mc_task_trt, fold_number, node))) {
                        param$clever_covariates(mc_task_trt, fold_number, node)    
                    }
                    if (is.null(self$observed_likelihood$cache$get_update_step(self$observed_likelihood$factor_list[[node]], mc_task_ctrl, fold_number, node))) {
                        param$clever_covariates(mc_task_ctrl, fold_number, node)
                    }
                })
                # for (param in self$observed_likelihood$updater$tmle_params) for (node in all_update_nodes) {
                #     # param$gradient$expand_task(mc_task_trt, node)
                #     param$clever_covariates(mc_task_trt, fold_number, node)
                #     param$clever_covariates(mc_task_ctrl, fold_number, node)
                # }
                
                # self$observed_likelihood$updater$tmle_params %>% lapply(function(param) for (node in all_update_nodes) param$gradient$expand_task(mc_task_ctrl, node))
                # for (update_node in self$update_nodes) {
                #     self$clever_covariates(tmle_task = mc_task_trt, fold_number = fold_number, node = update_node, for_fitting = F)
                #     self$clever_covariates(tmle_task = mc_task_ctrl, fold_number = fold_number, node = update_node, for_fitting = F)
                # }
                
                
                self$observed_likelihood$sync_task(mc_task_trt,
                                                   fold_number = fold_number,
                                                   check = T)
                
                self$observed_likelihood$sync_task(mc_task_ctrl,
                                                   fold_number = fold_number,
                                                   check = T)
                
                loc_outcome <- which(temp_node_names == self$outcome_node)
                # for all non-A, non-0 variables, calculate the variable by rule
                # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                # note that list_all_predicted_lkd is ordered by node
                temp_list_0 <- lapply(loc_Z[loc_Z <= loc_outcome],
                                      function(each_t) {
                                          self$observed_likelihood$get_likelihood(mc_task_ctrl, node = temp_node_names[each_t], fold_number = fold_number, check_sync = F)
                                      })
                temp_list_1 <- lapply(loc_RLY[loc_RLY <= loc_outcome],
                                      function(each_t) {
                                          self$observed_likelihood$get_likelihood(mc_task_trt, node = temp_node_names[each_t], fold_number = fold_number, check_sync = F)
                                      })
                temp_list_b <- list(
                    self$observed_likelihood$get_likelihood(mc_task_trt, node = temp_node_names[1], fold_number = fold_number, check_sync = F)
                )
                temp_list <- c(temp_list_0, temp_list_1, temp_list_b)
                psi <- pmap_dbl(temp_list, prod) %>% head(-1) %>% mean
                
                actual_V <- prep_list[1:grep(paste0("^", tmle_task$npsem[[self$outcome_node]]$variables, "$"), names(prep_list))] %>% sapply(length) %>% prod
                psi <- actual_V * psi
            }
            
            
            
            
            
            if (final) {
                # # ZW: allow other outcome node later
                # last_node <- self$outcome_node
                # loc_last_node <- which(temp_node_names == last_node)
                # which_intervention_keep <- sapply(self$intervention_list_treatment, function(u) which(u$name == temp_node_names)) <= loc_last_node
                #
                #
                #
                # # get vec_est and psi by sequentially resampling
                # observed_task <- tmle_task
                # cf_likelihood_treatment <- CF_Likelihood$new(self$observed_likelihood, self$intervention_list_treatment[which_intervention_keep])
                # cf_likelihood_control <- CF_Likelihood$new(self$observed_likelihood, self$intervention_list_control[which_intervention_keep])
                #
                # num_samples = est_n_samples
                # baseline_node = names(observed_task$npsem)[1]
                # tlik <- self$observed_likelihood
                # time_ordering <- names(tmle_task$npsem)
                # time_ordering <- time_ordering[1:loc_last_node]
                #
                # mc_data <- as.data.table(self$observed_likelihood$factor_list[[baseline_node]]$sample(, num_samples))
                # baseline_var <- observed_task$npsem[[baseline_node]]$variables
                # set(mc_data, , setdiff(names(observed_task$data), baseline_var), (0))
                # mc_data$t <- 0
                # mc_data$id <- 1:nrow(mc_data)
                # mc_task_trt <- tmle3_Task$new(mc_data, observed_task$npsem, id = "id", t = "t", long_format  = observed_task$long_format)
                # mc_task_ctrl <- tmle3_Task$new(mc_data, observed_task$npsem, id = "id", t = "t", long_format  = observed_task$long_format)
                # time_ordering <- setdiff(time_ordering, baseline_node)
                # for(node in time_ordering) {
                #   if (node %in% nodes_A) {
                #     if (length(grep("^A_C_", node))>0) {
                #       sampled <- mc_data[, .(id, t)]
                #       sampled[, (node) := 1]
                #       mc_task_trt <- mc_task_trt$generate_counterfactual_task(UUIDgenerate(), sampled)
                #       mc_task_ctrl <- mc_task_ctrl$generate_counterfactual_task(UUIDgenerate(), sampled)
                #     } else {
                #       mc_task_trt <- sample_from_node(mc_task_trt, cf_likelihood_treatment, node, observed_task)
                #       mc_task_ctrl <- sample_from_node(mc_task_ctrl, cf_likelihood_control, node, observed_task)
                #     }
                #   } else if (node %in% nodes_Z) {
                #     sampled <- sample_from_node(mc_task_ctrl, tlik, node, observed_task, if_only_return_sampled = T)
                #     mc_task_trt <- mc_task_trt$generate_counterfactual_task(UUIDgenerate(), sampled)
                #     mc_task_ctrl <- mc_task_ctrl$generate_counterfactual_task(UUIDgenerate(), sampled)
                #   } else {
                #     sampled <- sample_from_node(mc_task_trt, tlik, node, observed_task, if_only_return_sampled = T)
                #     mc_task_trt <- mc_task_trt$generate_counterfactual_task(UUIDgenerate(), sampled)
                #     mc_task_ctrl <- mc_task_ctrl$generate_counterfactual_task(UUIDgenerate(), sampled)
                #   }
                # }
                #
                # psi <- mean(mc_task_trt$data[[self$outcome_node]])
                # psi <- mean(mc_task_trt$data[[last(temp_node_names)]])
            }
            
            
            EIC <- cbind(EIC, rep(0, nrow(EIC))
                         # vec_est - psi
            )
            
            IC <- rowSums(EIC)
            result <- list(psi =
                               psi
                           # list_all_predicted_lkd
                           ,
                           IC = IC, EIC = colMeans(EIC)
                           # , full_EIC = EIC
            )
            return(result)
        }
    ),
    active = list(
        name = function() {
            param_form <- sprintf("E[%s_{%s; %s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$cf_likelihood_control$name)
            return(param_form)
        },
        cf_likelihood_treatment = function() {
            return(private$.cf_likelihood_treatment)
        },
        cf_likelihood_control = function() {
            return(private$.cf_likelihood_control)
        },
        cf_task_treatment = function() {
            return(private$.cf_task_treatment)
        },
        cf_task_control = function() {
            return(private$.cf_task_control)
        },
        intervention_list_treatment = function() {
            return(self$cf_likelihood_treatment$intervention_list)
        },
        intervention_list_control = function() {
            return(self$cf_likelihood_control$intervention_list)
        },
        update_nodes = function() {
            return(c(private$.update_nodes))
        },
        outcome_node = function() {
            return(private$.outcome_node)
        },
        gradient = function(){
            private$.gradient
        },
        static_likelihood = function(){
            private$.static_likelihood
        }
    ),
    private = list(
        .type = "mediation_survival",
        .cf_likelihood_treatment = NULL,
        .cf_likelihood_control = NULL,
        .cf_task_treatment = NULL,
        .cf_task_control = NULL,
        .supports_outcome_censoring = FALSE,
        .gradient = NULL,
        .submodel_type_supported = c("EIC"),
        .update_nodes = NULL,
        .outcome_node = NULL,
        .static_likelihood = NULL
    )
)

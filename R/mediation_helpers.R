#' Helper Functions for Longitudinal Mediation
#'
#' Handles the ARZL(Y) structure
#'
#' @param data a \code{data.frame}, or \code{data.table} containing data for use in estimation
#' @param node_list a list of character vectors, listing the variables that comprise each node
#' @param variable_types a list of variable types, one for each node. If missing, variable types will be guessed
#' @param tmle_task a \code{\link{tmle3_Task}} as constructed via \code{point_tx_task}
#' @param learner_list a list of sl3 learners, one for A and one for Y to be used for likelihood estimation
#' @param ... extra arguments.
#' @export
#' @rdname point_tx
mediation_npsem <- function(node_list, variable_types = NULL) {
  # make tmle_task; NULL variable_types will be learned by sl3 function
  npsem <- c(define_node("L_0", node_list$L_0, variable_type = variable_types$L_0),
             lapply(2:length(node_list), function(k) {
                 define_node(names(node_list)[k],
                             node_list[[k]],
                             names(node_list)[1:(k-1)],
                             variable_type = variable_types[[ names(node_list)[k] ]])
             })
  )

  return(npsem)
}



#' @export
#' @rdname point_tx
mediation_task <- function(data, node_list, variable_types = NULL, ...) {
  setDT(data)

  npsem <- mediation_npsem(node_list, variable_types)

  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
  }

  return(tmle_task)
}

#' @export
#' @rdname point_tx
mediation_task_drop_censored <- function(data, node_list, variable_types = NULL, ...) {
  setDT(data)

  npsem <- mediation_npsem(node_list, variable_types)

  # add censoring nodes
  names(npsem) <- sapply(npsem, `[[`, "name")
  for (k in grep("^A_C_", names(npsem))) {
    current_A_C <- names(npsem)[k]
    current_time <- strsplit(current_A_C, "_") %>% first %>% last %>% as.numeric
    for (s in grep(paste0("_", current_time, "$"), names(npsem))) {
      temp_node <- names(npsem)[s]
      if (temp_node != current_A_C) npsem[[temp_node]][["censoring_node"]] <- npsem[[current_A_C]] else {
        list_loc_last_AC <- paste0("^A_C_", 1:current_time, "$") %>% head(-1) %>% lapply(function(x) grep(x, names(npsem)))
        which_last_AC <- which(sapply(list_loc_last_AC, length) != 0) %>% last
        if (length(which_last_AC) != 0) if (!is.na(which_last_AC)) npsem[[temp_node]][["censoring_node"]] <- npsem[[list_loc_last_AC[[which_last_AC]]]]
      }
    }
  }

  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task_drop_censored$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task_drop_censored$new(data, npsem = npsem, ...)
  }

  return(tmle_task)
}

#' @export
#' @rdname point_tx
mediation_likelihood <- function(tmle_task, learner_list) {
  factor_list <- list()

  # covariates
  L_0_factor <- define_lf(LF_emp, "L_0")
  factor_list[[1]] <- L_0_factor

  # treatment (bound likelihood away from 0 (and 1 if binary))
  # A_type <- tmle_task$npsem[["A"]]$variable_type
  A_type <- tmle_task$npsem[[ grep("^A", names(tmle_task$npsem))[1] ]]$variable_type  # evaluate the first A node for now
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  temp_names <- names(tmle_task$npsem)
  loc_A <- grep("^A", temp_names)
  factor_list[loc_A] <- lapply(loc_A, function(k) {
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], bound = A_bound, type = "density")
  })

  # others
  loc_others <- (1:length(temp_names))[-c(grep("^A", temp_names), 1)]
  factor_list[loc_others] <- lapply(loc_others, function(k) {
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], type = "density")
  })

  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}


#' @export
get_current_Q <- function(loc_node, which_Q,
                          tmle_task, obs_variable_names,
                          intervention_variables, intervention_levels_treat, intervention_levels_control,
                          list_all_predicted_lkd, if_survival = F
) {
  temp_node_names <- names(tmle_task$npsem)  # this has to be the obs task, but only the node names are needed
  loc_delta <- grep("delta_", temp_node_names)
  if (length(loc_delta) > 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
  if (length(unlist(lapply(temp_node_names, function(x) tmle_task$npsem[[x]]$censoring_node$name))) == 0) {  # no missingness node
    loc_A <- grep("A", sapply(temp_node_names, function(x) substr(x, 1, 1)))
  } else {
    loc_A <- grep("A_E", temp_node_names)
    loc_A_C <- grep("A_C", temp_node_names)
  }
  loc_Z <- which(sapply(temp_node_names, function(s) paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "A", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # ind_var is the order among variables
  ind_var <- loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
  # decide which Q to get; this variable will be integrated out; note which_Q=1 just makes sure input include current variable
  ind_var <- ind_var + which_Q

  intervention_variables_loc_needed <- intervention_variables_loc[intervention_variables_loc < ind_var]

  # only update t!=0 and non-A nodes
  if (loc_node %in% c(loc_Z, loc_RLY)) {
    df_to_update <- temp_current <- list_all_predicted_lkd[[loc_node]] # decide input matrix; the last column might not be fully used

    # get Q current
    if (ind_var <= length(obs_variable_names)) {  # Q_X, except for Q_R_tao+1
      data_temp <- df_to_update[1:(ind_var-1)]

      if (!if_survival) {  # regular outcome; just need non-zero last outcome node
        # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
        all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var-1) ),
                                            A = 1,
                                            rule_variables = c(last(obs_variable_names)), rule_values = c(1))
        all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var-1) ),
                                            A = 0,
                                            rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      } else {  # for decreasing survival events, all Y and A_C need to be 1
        loc_impute <- grep("^Y_|^A_C_", temp_node_names)  # remain alive and uncensored before current variable
        # loc_impute <- loc_impute[loc_impute < loc_node]
        all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var-1) ),
                                             rule_variables = c(last(obs_variable_names),
                                                                intervention_variables,
                                                                sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                             rule_values = c(1, intervention_levels_treat, rep(1, length(loc_impute))))
        all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var-1) ),
                                             rule_variables = c(last(obs_variable_names), intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                             rule_values = c(1, intervention_levels_control, rep(1, length(loc_impute))))
      }

      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_temp[1:(ind_var-1)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- temp_all_comb_1 <- unique_input[which_row, ]
                                       if (length(all_possible_rlz_0) > 0) {
                                         temp_all_comb_0 <- suppressWarnings(cbind(temp_all_comb_0, all_possible_rlz_0))
                                         temp_all_comb_1 <- suppressWarnings(cbind(temp_all_comb_1, all_possible_rlz_1))
                                       }
                                       if (length(intervention_variables_loc_needed) > 0) {
                                         for (i in 1:length(intervention_variables_loc_needed)) {
                                           temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                           temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                         }
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_need <- loc_Z[loc_Z >= loc_node + which_Q]  # only integrate out future variables
                                       temp_list_0 <- lapply(loc_Z_need,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output %>% suppressMessages
                                                             })
                                       loc_other_needed <- loc_RLY[loc_RLY >= loc_node + which_Q]
                                       temp_list_1 <- lapply(loc_other_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output %>% suppressMessages
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      Q_current <- left_join(data_temp, library_output)$output %>% suppressMessages
    } else {
      Q_current <- df_to_update[, ind_var-1]  # Q_R_tao+1 is just the output
    }

    # in case survival param create NA in Q integrals
    Q_current[is.na(Q_current)] <- 0
    return(Q_current)
  } else {
    return(NULL)
  }
}


#' @export
get_current_H <- function(loc_node,
                          tmle_task, obs_variable_names,
                          intervention_variables, intervention_levels_treat, intervention_levels_control,
                          list_all_predicted_lkd
) {
  temp_node_names <- names(tmle_task$npsem)  # this has to be the obs task, but only the node names are needed
  loc_delta <- grep("delta_", temp_node_names)
  if (length(loc_delta) > 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting

  if (length(unlist(lapply(temp_node_names, function(x) tmle_task$npsem[[x]]$censoring_node$name))) == 0) {  # no missingness node
    loc_A <- grep("A", sapply(temp_node_names, function(x) substr(x, 1, 1)))
  } else {
    loc_A <- grep("A_E", temp_node_names)
    loc_A_C <- grep("A_C", temp_node_names)
  }
  loc_Z <- which(sapply(temp_node_names, function(s) paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "A", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # ind_var is the order among variables
  ind_var <- loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
  intervention_variables_loc_needed <- intervention_variables_loc[intervention_variables_loc < ind_var]

  temp_current <- list_all_predicted_lkd[[loc_node]] # decide input matrix; the last column might not be fully used

  # get H current
  {
    data_temp <- temp_current[1:ind_var]  # take =1 probs
    temp_ind <- ind_var  # the variable loc

    all_observed_1 <- all_observed_0 <- data_temp
    if (length(intervention_variables_loc_needed) > 0) {
      for (i in 1:length(intervention_variables_loc_needed)) {
        all_observed_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
        all_observed_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
      }
    }

    # R and L (including Y before t = tau) (or not Z not A not t=0, not Y_tau) can be calculated in the same way
    if (loc_node %in% loc_RLY) {
      loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
      loc_Z_needed <- loc_Z[loc_Z < loc_node]  # all needed Z nodes

      # all predicted probs at these locs are needed
      # needed inputs are in either all_observed_1 or all_observed_0

      # these A probs will be taken as product
      part_A <- lapply(loc_A_needed, function(k) {
        if (length(unlist(lapply(temp_node_names, function(x) tmle_task$npsem[[x]]$censoring_node$name))) == 0) {  # no missingness, no censoring probs
          left_join(all_observed_1, list_all_predicted_lkd[[k]])$output %>% suppressMessages
        } else {
          temp_p_A_E <- left_join(all_observed_1, list_all_predicted_lkd[[k]])$output %>% suppressMessages
          k_A_C <- loc_A_C[loc_A_C < k]
          k_A_C <- k_A_C[which.min(abs(k_A_C - k))]  # choose the closest censoring node
          temp_p_A_C <- left_join(all_observed_1, list_all_predicted_lkd[[k_A_C]])$output %>% suppressMessages
          return(temp_p_A_C * temp_p_A_E)
        }
      }) %>% pmap_dbl(prod)

      # these ratios Z probs will be taken as product
      part_Z <- ifelse_vec(length(loc_Z_needed) == 0, rep(1, length(part_A)),
                           lapply(loc_Z_needed, function(k) {
                             left_join(all_observed_0, list_all_predicted_lkd[[k]])$output /
                               left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
                           }) %>% pmap_dbl(prod)) %>% suppressMessages

      H_current <- ifelse(
        apply(data_temp[intervention_variables_loc_needed] == intervention_levels_treat, 1, all)  # this is the A_ind decided by the library tasks
        , 1/part_A*part_Z, 0) %>% as.vector
    }
    # Z nodes
    if (loc_node %in% loc_Z) {
      loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
      loc_RLY_needed <- loc_RLY[loc_RLY < loc_node]

      # these A probs will be taken as product
      part_A <- lapply(loc_A_needed, function(k) {
        if (length(unlist(lapply(temp_node_names, function(x) tmle_task$npsem[[x]]$censoring_node$name))) == 0) {  # no missingness, no censoring probs
          left_join(all_observed_0, list_all_predicted_lkd[[k]])$output %>% suppressMessages
        } else {
          temp_p_A_E <- left_join(all_observed_0, list_all_predicted_lkd[[k]])$output %>% suppressMessages
          k_A_C <- loc_A_C[loc_A_C < k]
          k_A_C <- k_A_C[which.min(abs(k_A_C - k))]  # choose the closest censoring node
          temp_p_A_C <- left_join(all_observed_0, list_all_predicted_lkd[[k_A_C]])$output %>% suppressMessages
          return(temp_p_A_C * temp_p_A_E)
        }
      }) %>% pmap_dbl(prod)

      # these ratios Z probs will be taken as product
      part_LR <- lapply(loc_RLY_needed, function(k) {
        left_join(all_observed_1, list_all_predicted_lkd[[k]])$output /
          left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
      }) %>% pmap_dbl(prod)

      H_current <- ifelse(
        # ZW todo: stochastic intervention
        apply(data_temp[intervention_variables_loc_needed] == intervention_levels_control, 1, all)
        , 1/part_A*part_LR, 0)
    }
  }
  # H_current[H_current > 50] <- 50
  return(H_current)
}


















get_obs_Q_list <- function(tmle_task, obs_data,
                           intervention_variables, intervention_levels_treat, intervention_levels_control,
                           list_all_predicted_lkd, loc_outcome = NULL) {
    # decide validation or full by list_all_predicted_lkd
    obs_variable_names <- colnames(obs_data)
    temp_node_names <- names(tmle_task$npsem)
    loc_delta <- grep("delta_", temp_node_names)
    loc_delta_nodes <- grep("delta_", temp_node_names)
    if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
    
    if (length(unlist(lapply(temp_node_names, function(x) tmle_task$npsem[[x]]$censoring_node$name))) == 0) {  # no missingness node for all nodes
        loc_A <- grep("^A", sapply(temp_node_names, function(x) substr(x, 1, 1)))
        loc_A_C <- loc_var_A_C <- NULL
    } else {
        loc_A <- grep("^A_E", temp_node_names)
        loc_A_C <- grep("^A_C", temp_node_names)
        loc_var_A_C <- sapply(loc_A_C, function(x) which(obs_variable_names == tmle_task$npsem[[x]]$variables))
    }
    loc_Z <- which(sapply(temp_node_names, function(s) paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") == "Z"))
    loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "A", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
    intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
    
    list_Q <- list()
    for (loc_node in 1:length(temp_node_names)) {
        if (loc_node %in% c(loc_Z, loc_RLY)) {
            # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
            loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
            
            if (length(unlist(lapply(temp_node_names, function(x) tmle_task$npsem[[x]]$censoring_node$name))) == 0) {  # no missingness/censoring
                loc_var_last_not_integrated <- loc_current_var-1
                all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:(loc_var_last_not_integrated) ),
                                                     rule_variables = c(last(obs_variable_names), intervention_variables),
                                                     rule_values = c(1, intervention_levels_treat))
                all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:(loc_var_last_not_integrated) ),
                                                     rule_variables = c(last(obs_variable_names), intervention_variables),
                                                     rule_values = c(1, intervention_levels_control))
            } else {  # missingness/right-censoring
                loc_impute <- grep("Y_|A_C_", temp_node_names)  # in integral, Y and A_C always 1
                if (length(loc_A_C) > 0) loc_impute <- loc_impute[loc_impute <= last(loc_A_C)]  # the event after the last censoring node can be alive/dead; do not impute
                loc_previous <- c(loc_Z, loc_RLY)[c(loc_Z, loc_RLY) < loc_node]
                if (length(loc_previous) == 0) {
                    loc_var_last_not_integrated <- loc_current_var - 1
                } else {
                    loc_var_last_not_integrated <- which(obs_variable_names == tmle_task$npsem[[max(loc_previous)]]$variables)
                }
                all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:(loc_var_last_not_integrated) ),
                                                     rule_variables = c(last(obs_variable_names), intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                                     rule_values = c(1, intervention_levels_treat, rep(1, length(loc_impute))))
                all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:(loc_var_last_not_integrated) ),
                                                     rule_variables = c(last(obs_variable_names), intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                                     rule_values = c(1, intervention_levels_control, rep(1, length(loc_impute))))
            }
            
            # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
            if (!is.null(loc_var_A_C)) {
                loc_var_A_C_used <- loc_var_A_C[loc_var_A_C < loc_var_last_not_integrated]
                if (length(loc_var_A_C_used) > 0) {
                    ind_full_observed <- apply(matrix(obs_data[, loc_var_A_C_used] == 1, ncol = length(loc_var_A_C_used)), 1, all)
                    ind_full_observed[is.na(ind_full_observed)] <- F
                    # T, NA usually means a death event happens; but we only care uncensored subjects in D cal with right-censoring
                } else {
                    ind_full_observed <- rep(T, nrow(obs_data))
                }
            } else {
                ind_full_observed <- rep(T, nrow(obs_data))
            }
            unique_input <- obs_data[ind_full_observed, 1:(loc_var_last_not_integrated)] %>% unique
            library_output <- data.frame(unique_input, output =
                                             map_dbl(1:nrow(unique_input), function(which_row) {
                                                 # probs in the integrals, A=1 or A=0 is inserted
                                                 temp_all_comb_0 <- data.frame(unique_input[which_row, ], all_possible_RZLY_0) %>% suppressWarnings
                                                 temp_all_comb_1 <- data.frame(unique_input[which_row, ], all_possible_RZLY_1) %>% suppressWarnings
                                                 for (i in 1:length(intervention_variables_loc)) {
                                                     temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                                     temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                                 }
                                                 # for all non-A, non-0 variables, calculate the variable by rule
                                                 # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                                 loc_Z_needed <- loc_Z[loc_Z >= loc_node]  # integrate out current variable
                                                 temp_list_0 <- lapply(loc_Z_needed,
                                                                       function(each_t) {
                                                                           left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output %>% suppressMessages
                                                                       })
                                                 loc_RLY_needed <- loc_RLY[loc_RLY >= loc_node]
                                                 temp_list_1 <- lapply(loc_RLY_needed,
                                                                       function(each_t) {
                                                                           left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output %>% suppressMessages
                                                                       })
                                                 temp_list <- c(temp_list_0, temp_list_1)
                                                 if (!is.null(loc_outcome)) where_sum <- temp_all_comb_1[[tmle_task$npsem[[loc_outcome]]$variables]] == 1 else where_sum <- rep(T, unique(sapply(temp_list_1, length)))
                                                 pmap_dbl(temp_list, prod)[where_sum] %>% sum %>% return
                                             })
            )
            if (any(is.na(library_output[[ncol(library_output)]]))) {  # only reason to have NA here is Y=0
                if (!all(library_output[is.na(library_output[[ncol(library_output)]]), ncol(library_output)-1] == 0)) stop(
                    "list_all_predicted_lkd does not contain observed inputs for reasons other than Yt=0"
                ) else {
                    library_output[is.na(library_output[[ncol(library_output)]]), ncol(library_output)] <- 0  # zero prob for Yt=0
                }
            }
            temp_vec <- left_join(obs_data[, 1:(loc_var_last_not_integrated)], library_output)$output %>% suppressMessages
            temp_vec <- temp_vec[!is.na(temp_vec)]
            list_Q[[loc_node]] <- temp_vec
        }
    }
    return(list_Q)
}






get_obs_H_list <- function(tmle_task, obs_data, current_likelihood,
                           cf_task_treatment, cf_task_control,
                           intervention_variables, intervention_levels_treat, intervention_levels_control,
                           fold_number = "full"
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_delta <- grep("delta_", temp_node_names)
  loc_delta_nodes <- grep("delta_", temp_node_names)
  if (length(loc_delta_nodes) != 0) temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting

  if (length(unlist(lapply(temp_node_names, function(x) tmle_task$npsem[[x]]$censoring_node$name))) == 0) {  # no missingness node
    loc_A <- grep("A", sapply(temp_node_names, function(x) substr(x, 1, 1)))
  } else {
    loc_A <- grep("A_E", temp_node_names)
    loc_A_C <- grep("A_C", temp_node_names)
  }
  loc_Z <- which(sapply(temp_node_names, function(s) paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) !(paste0(head(strsplit(s, "_")[[1]], -1), collapse = "_") %in% c("A_C", "A_E", "A", "Z")) & tail(strsplit(s, "_")[[1]], 1) != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  list_H <- list()  # get a list of corresponding H covariates; ordered by nodes, not variables
  for (temp_ind in loc_RLY) {  # calculate RLY nodes
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_Z_needed <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
    A_ind <- apply(sapply(loc_A_needed, function(k) {
      obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_treat[tmle_task$npsem[[k]]$name]
    }), 1, all)
    part_A <- lapply(loc_A_needed, function(k) {
      if (length(unlist(lapply(temp_node_names, function(x) tmle_task$npsem[[x]]$censoring_node$name))) == 0) {  # no missingness
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)
      } else {
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
        k_A_C <- k_A_C[which.min(abs(k_A_C - k))]  # choose the closest censoring node
        temp_p_A_C <- current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k_A_C], fold_number)  # A_C=1; to be aligned
        if (!is.null(tmle_task$npsem[[k_A_C]]$censoring_node$variables)) {  # if there is missingness in A_E
          temp_full <- if_A_C_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$censoring_node$name)
          if (!is.logical(if_A_C_observed)) if_A_C_observed <- if_A_C_observed == 1  # in case it is a binary node
          if_A_C_observed[is.na(if_A_C_observed)] <- F
          observed_A_C <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$name)  # maybe due to Y death
          if_A_C_observed[is.na(observed_A_C)] <- F
          temp_full[if_A_C_observed] <- temp_p_A_C
          temp_full[!if_A_C_observed] <- NA
          temp_p_A_C <- temp_full
        }
        return(temp_p_A_C * temp_p_A_E)
      }
    }) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      temp_p <- current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)
      if (!is.null(tmle_task$npsem[[k]]$censoring_node$variables)) {  # if there is missingness
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
    temp_vec <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
    temp_vec[is.na(temp_vec)] <- 0  # due to bivariate trt nodes or g-comp
    if(!is.null(tmle_task$npsem[[temp_ind]]$censoring_node$variables)) {
      if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[temp_ind]]$censoring_node$name)
      if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
      if_observed[is.na(if_observed)] <- F
      temp_vec <- temp_vec[if_observed]
    }
    list_H[[temp_ind]] <- temp_vec
  }
  for (temp_ind in loc_Z) {    # calculate Z nodes
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_RLY_needed <- loc_RLY[loc_RLY < temp_ind]
    A_ind <- apply(sapply(loc_A_needed, function(k) {
      obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_control[tmle_task$npsem[[k]]$name]
    }), 1, all)
    part_A <- lapply(loc_A_needed, function(k) {  # ctrl propensity under A_C=1
      temp_p_A_E <- current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)  # A_E = 0 | A_C=1
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
      temp_p_A_C <- current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k_A_C], fold_number)  # A_C=1; to be aligned
      if (!is.null(tmle_task$npsem[[k_A_C]]$censoring_node$variables)) {  # if there is missingness in A_E
        temp_full <- if_A_C_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$censoring_node$variables)
        if (!is.logical(if_A_C_observed)) if_A_C_observed <- if_A_C_observed == 1  # in case it is a binary node
        if_A_C_observed[is.na(if_A_C_observed)] <- F
        observed_A_C <- tmle_task$get_tmle_node(tmle_task$npsem[[k_A_C]]$name)  # maybe due to Y death
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
        if_observed[is.na(if_observed)] <- F
        temp_full[if_observed] <- temp_p
        temp_full[!if_observed] <- NA
        temp_p <- temp_full
      }
      return(temp_p)
    }) %>% pmap_dbl(prod)
    if(length(part_RLY) == 0) part_RLY <- 1

    temp_vec <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
    temp_vec[is.na(temp_vec)] <- 0  # due to bivariate trt nodes or g-comp
    if(!is.null(tmle_task$npsem[[temp_ind]]$censoring_node$variables)) {
      if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[temp_ind]]$censoring_node$name)
      if (!is.logical(if_observed)) if_observed <- if_observed == 1  # in case it is a binary node
      if_observed[is.na(if_observed)] <- F
      temp_vec <- temp_vec[if_observed]
    }
    list_H[[temp_ind]] <- temp_vec
  }
  # return(list_H)

  list_H_bounded <- list_H
  # list_H_bounded <- lapply(list_H, function(x) {
  #   if (!is.null(x)) {
  #     y <- x
  #     y[y>50] <- 50
  #     return(y)
  #   }
  # })
  return(list_H_bounded)
}

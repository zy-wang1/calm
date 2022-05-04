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
middle_npsem <- function(node_list, variable_types = NULL) {
  # make tmle_task
  npsem <- c(define_node("L_0", node_list$L_0, variable_type = variable_types$L_0),
             lapply(2:length(node_list), function(k) {
               if (k < length(node_list)) {
                 define_node(names(node_list)[k],
                             node_list[[k]],
                             names(node_list)[1:(k-1)],
                             variable_type = variable_types[[ names(node_list)[k] ]])
               } else {
                 define_node(names(node_list)[k],
                             node_list[[k]],
                             names(node_list)[1:(k-1)],
                             variable_type = variable_types[[ names(node_list)[k] ]],
                             scale = TRUE)
               }
             })
  )

  # npsem <- list(
  #   define_node("W", node_list$W, variable_type = variable_types$W),
  #   define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
  #   define_node("Y", node_list$Y, c("A", "W"), variable_type = variable_types$Y, scale = TRUE)
  # )

  return(npsem)
}


#' @export
#' @rdname point_tx
gamma_npsem <- function(node_list, variable_types = NULL, type = 1, j = 1) {
  temp_node_names <- names(node_list)
  which_A_nodes <- grep("A", temp_node_names)
  A_nodes <- temp_node_names[which_A_nodes]
  last_RZLY_node <- if (type == 1) paste0("Z_", j) else temp_node_names[which(temp_node_names == paste0("Z_", j)) - 1]
  which_last_RZLY <- which(temp_node_names == last_RZLY_node)
  # intersect(1:which_last_RZLY, )
  npsem <- lapply(which_A_nodes, function(each_loc_A) {
    which_A_nodes_not_needed <- which_A_nodes[which_A_nodes >= each_loc_A]
    which_predictor_nodes <- setdiff(1:which_last_RZLY, which_A_nodes_not_needed)
    define_node(temp_node_names[each_loc_A],
                node_list[[each_loc_A]],
                temp_node_names[which_predictor_nodes],
                variable_type = variable_types[[each_loc_A]])
  })

  # npsem <- list(
  #   define_node("W", node_list$W, variable_type = variable_types$W),
  #   define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
  #   define_node("Y", node_list$Y, c("A", "W"), variable_type = variable_types$Y, scale = TRUE)
  # )

  return(npsem)
}



#' @export
#' @rdname point_tx
middle_task <- function(data, node_list, variable_types = NULL, ...) {
  setDT(data)

  npsem <- middle_npsem(node_list, variable_types)

  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
  }

  return(tmle_task)
}

#' @export
#' @rdname point_tx
middle_task_drop_censored <- function(data, node_list, variable_types = NULL, ...) {
  setDT(data)

  npsem <- middle_npsem(node_list, variable_types)

  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task_drop_censored$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task_drop_censored$new(data, npsem = npsem, ...)
  }

  return(tmle_task)
}

#' @export
#' @rdname point_tx
middle_likelihood <- function(tmle_task, learner_list) {
  factor_list <- list()

  # covariates
  L_0_factor <- define_lf(LF_emp, "L_0")
  factor_list[[1]] <- L_0_factor

  # treatment (bound likelihood away from 0 (and 1 if binary))
  # A_type <- tmle_task$npsem[["A"]]$variable_type
  A_type <- tmle_task$npsem[[ grep("A", names(tmle_task$npsem))[1] ]]$variable_type  # evaluate the first A node
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  temp_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_names)
  factor_list[loc_A] <- lapply(loc_A, function(k) {
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], bound = A_bound, type = "density")
  })
  # A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = A_bound)

  # others
  loc_others <- (1:length(temp_names))[-c(grep("A", temp_names), 1)]
  factor_list[loc_others] <- lapply(loc_others, function(k) {
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], type = "density")
  })

  # # outcome
  # Y_factor <- define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")
  #
  # # construct and train likelihood
  # factor_list <- list(W_factor, A_factor, Y_factor)

  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  # temp <- likelihood$list_all_predicted_lkd
  return(likelihood)
}








#' @export
get_current_newH <- function(loc_node,
                             tmle_task, obs_data,
                             intervention_variables, intervention_levels_treat, intervention_levels_control,
                             list_all_predicted_lkd
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)  # this has to be the obs task, but only the node names are needed
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # ind_var is the order among variables
  ind_var <- loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
  intervention_variables_loc_needed <- intervention_variables_loc[intervention_variables_loc < ind_var]

  # only update t!=0 and non-A nodes
  if (loc_node %in% c(loc_Z, loc_RLY)) {
    temp_current <- list_all_predicted_lkd[[loc_node]]
    Xt_input <- temp_current[ind_var]
    # loc_to_update <- Xt_input == 1
    # df_to_update <- temp_current[loc_to_update, ]  # the possible input where L_t = 1
    df_to_update <- df_to_update_0 <- temp_current
    df_to_update[[ind_var]] <- 1
    df_to_update_0[[ind_var]] <- 0

    # df_to_update_0 <- df_to_update
    # df_to_update_0[ind_var] <- 0  # replacing the above input with Lt = 0


    # get Q1 current
    {
      data_temp <- df_to_update[1:ind_var]
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 1,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 0,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_temp[1:(ind_var)] %>% unique
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
                                       loc_Z_need <- loc_Z[loc_Z > loc_node]  # only integrate out future variables
                                       temp_list_0 <- lapply(loc_Z_need,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_other_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_other_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      Q1_current <- left_join(data_temp, library_output)$output
      if (loc_node == length(temp_node_names)) Q1_current <- rep(1, length(Q1_current))
      # if (loc_node == length(temp_node_names)) Q1_current <- rep(0, length(Q1_current))
    }

    # get Q0 current
    {
      data_temp <- df_to_update_0[1:ind_var]
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 1,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 0,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_temp[1:(ind_var)] %>% unique
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
                                       loc_Z_need <- loc_Z[loc_Z > loc_node]  # only integrate previous variables; debug: not including current
                                       temp_list_0 <- lapply(loc_Z_need,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_other_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_other_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      Q0_current <- left_join(data_temp, library_output)$output
      if (loc_node == length(temp_node_names)) Q0_current <- rep(0, length(Q0_current))
    }

    # get H current
    {
      data_temp <- temp_current[1:ind_var]  # take =1 probs
      temp_ind <- ind_var

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
          left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        # these ratios Z probs will be taken as product
        part_Z <- ifelse_vec(length(loc_Z_needed) == 0, rep(1, length(part_A)),
                             lapply(loc_Z_needed, function(k) {
                               left_join(all_observed_0, list_all_predicted_lkd[[k]])$output /
                                 left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
                             }) %>% pmap_dbl(prod))

        H_current <- ifelse(
          # data_temp[last(intervention_variables_loc_needed)] == 1
          apply(data_temp[intervention_variables_loc_needed] == 1, 1, prod) == 1  # this is the A_ind decided by the library tasks
          , 1/part_A*part_Z, 0) %>% as.vector
      }
      # Z nodes
      if (loc_node %in% loc_Z) {
        loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
        loc_RLY_needed <- loc_RLY[loc_RLY < loc_node]

        # these A probs will be taken as product
        part_A <- lapply(loc_A_needed, function(k) {
          left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        # these ratios Z probs will be taken as product
        part_LR <- lapply(loc_RLY_needed, function(k) {
          left_join(all_observed_1, list_all_predicted_lkd[[k]])$output /
            left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        H_current <- ifelse(
          # ZW todo: stochastic intervention
          # data_temp[last(intervention_variables_loc_needed)] == 0
          apply(data_temp[intervention_variables_loc_needed] == 0, 1, prod) == 1
          , 1/part_A*part_LR, 0)
      }
    }

    # get newH current
    newH_current <- H_current * (Q1_current - Q0_current)
    return(newH_current)
  } else {
    return(NULL)
  }
}




#' @export
get_obs_Q <- function(tmle_task, obs_data, list_H,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      list_all_predicted_lkd,
                      lt) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  list_Q <- list()
  for (loc_node in 1:length(list_H)) {
    if(!is.null(list_H[[loc_node]])) {  # set current var as lt, prod and sum all future conditional lkd
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
      all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:(loc_current_var-1) ),
                                           rule_variables = c(last(obs_variable_names), obs_variable_names[loc_current_var], intervention_variables),
                                           rule_values = c(1, lt, intervention_levels_treat))
      all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:(loc_current_var-1) ),
                                           rule_variables = c(last(obs_variable_names), obs_variable_names[loc_current_var], intervention_variables),
                                           rule_values = c(1, lt, intervention_levels_control))

      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- obs_data[, 1:(loc_current_var-1)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- data.frame(unique_input[which_row, ], all_possible_RZLY_0)
                                       temp_all_comb_1 <- data.frame(unique_input[which_row, ], all_possible_RZLY_1)
                                       for (i in 1:length(intervention_variables_loc)) {
                                         temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                         temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_needed <- loc_Z[loc_Z > loc_node]  # only product children variables
                                       # temp_list_0 <- lapply(loc_Z_needed,
                                       temp_list_0 <- lapply(loc_Z_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_RLY_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_RLY_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      list_Q[[loc_node]] <- left_join(obs_data[, 1:(loc_current_var-1)], library_output)$output
      if (loc_node == length(list_H))  list_Q[[loc_node]] <- rep(lt, nrow(obs_data))
      # if (loc_node == length(list_H))  list_Q[[loc_node]] <- rep(0, nrow(obs_data))
    }
  }

  return(list_Q)
}





#' @export
get_obs_H <- function(tmle_task, obs_data, current_likelihood,
                      cf_task_treatment, cf_task_control,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      fold_number = "full"
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

    list_H <- list()  # get a list of corresponding H covariates; ordered by nodes, not variables
  for (temp_ind in loc_RLY) {  # calculate RLY nodes
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_Z_needed <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
    A_ind <- apply(sapply(loc_A_needed, function(k) {
      obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_treat[tmle_task$npsem[[k]]$name]
    }), 1, all)
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)
    }) %>% pmap_dbl(prod)
    if(length(part_Z) == 0) part_Z <- 1
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
  }
  for (temp_ind in loc_Z) {    # calculate Z nodes
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_RLY_needed <- loc_RLY[loc_RLY < temp_ind]
    A_ind <- apply(sapply(loc_A_needed, function(k) {
      obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_control[tmle_task$npsem[[k]]$name]
    }), 1, all)
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)) %>% pmap_dbl(prod)
    part_RLY <- lapply(loc_RLY_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)
    }) %>% pmap_dbl(prod)
    if(length(part_RLY) == 0) part_RLY <- 1
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
  }
  return(list_H)
}



#' @export
get_obs_H_raw <- function(tmle_task, obs_data, current_likelihood,
                      cf_task_treatment, cf_task_control,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      fold_number = "full"
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # get a list of corresponding H covariates; ordered by nodes, not variables
  list_H <- list()
  # calculate RLY nodes
  for (temp_ind in loc_RLY) {
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_Z_needed <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
    # this is the At indicators for H_RLY; now
    A_ind <-
      # obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == intervention_levels_treat[tmle_task$npsem[[last(loc_A_needed)]]$variables]
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_treat[tmle_task$npsem[[k]]$variables]
      }), 1, prod) == 1
    # A_ind <- obs_data[[temp_node_names[last(loc_A_needed)]]]  # using variable names (rather than node names) to inquire obs_data
    # these A probs will be taken as product
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number) ) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)
    }) %>% pmap_dbl(prod)
    if(length(part_Z) == 0) part_Z <- 1

    list_H[[temp_ind]] <- (1/part_A*part_Z) %>% as.vector
  }
  # calculate Z nodes
  for (temp_ind in loc_Z) {
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_RLY_needed <- loc_RLY[loc_RLY < temp_ind]
    A_ind <-
      # obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == intervention_levels_control[tmle_task$npsem[[last(loc_A_needed)]]$variables]
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_control[tmle_task$npsem[[k]]$variables]
      }), 1, prod) == 1
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) ) %>% pmap_dbl(prod)
    part_RLY <- lapply(loc_RLY_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)
    }) %>% pmap_dbl(prod)
    list_H[[temp_ind]] <- (1/part_A*part_RLY) %>% as.vector
  }
  return(list_H)
}






#' @export
get_obs_H_full <- function(tmle_task, obs_data, current_likelihood,
                           cf_task_treatment, cf_task_control,
                           intervention_variables, intervention_levels_treat, intervention_levels_control,
                           fold_number = "full",
                           bound = NULL
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # get a list of corresponding H covariates; ordered by nodes, not variables
  list_H <- list()
  # calculate RLY nodes
  for (temp_ind in loc_RLY) {
    loc_A_needed <- loc_A
    loc_Z_needed <- loc_Z
    # this is the At indicators for H_RLY; now
    A_ind <-
      # obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == intervention_levels_treat[tmle_task$npsem[[last(loc_A_needed)]]$variables]
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_treat[tmle_task$npsem[[k]]$name]
      }), 1, prod) == 1
    # A_ind <- obs_data[[temp_node_names[last(loc_A_needed)]]]  # using variable names (rather than node names) to inquire obs_data
    # these A probs will be taken as product
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number) ) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)
    }) %>% pmap_dbl(prod)
    if(length(part_Z) == 0) part_Z <- 1

    if (!is.null(bound)) part_A[part_A < bound] <- bound
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
  }
  # calculate Z nodes
  for (temp_ind in loc_Z) {
    loc_A_needed <- loc_A
    loc_RLY_needed <- loc_RLY
    A_ind <-
      # obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == intervention_levels_control[tmle_task$npsem[[last(loc_A_needed)]]$variables]
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_control[tmle_task$npsem[[k]]$name]
      }), 1, prod) == 1
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) ) %>% pmap_dbl(prod)
    part_RLY <- lapply(loc_RLY_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)
    }) %>% pmap_dbl(prod)
    if (!is.null(bound)) part_A[part_A < bound] <- bound
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
  }
  return(list_H)
}





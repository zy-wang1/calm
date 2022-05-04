#' Longitudinal Mediation Targets
#'
#' Parameter definition for longitudinal mediation
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
#'   \code{define_param(Param_mediation, observed_likelihood, intervention_list, ..., outcome_node)}
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
Param_mediation <- R6Class(
  classname = "Param_mediation",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control) {
      # outcome_node is used to check self$supports_outcome_censoring; not checked for now
      super$initialize(observed_likelihood, list())
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
      observed_likelihood$get_likelihoods(observed_likelihood$training_task)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", update = F, node = NULL, submodel_type = "EIC") {
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
        full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
        full_variable_names <- colnames(full_data)
        list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
          if (loc_node > 1) {
            # currently only support univariate node for t>0
            current_variable <- full_task$npsem[[loc_node]]$variables
            temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
            temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
            temp_target_node <- intersect(self$update_nodes, full_node_names[loc_node])
            if (length(temp_target_node) == 1) {
              setattr(temp_task, "target_nodes", full_node_names[loc_node])
              temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
            } else {
              setattr(temp_task, "target_nodes", "no_update")
              temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
            }
            data.frame(temp_input, output = temp_output) %>% return
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
          loc_A <- grep("A", temp_node_names)
          loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
          loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))

          obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # note this is compatible if tmle_task is a cf task
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
          list_Q[[length(list_Q)+1]] <- tmle_task$get_tmle_node(length(list_Q))
          list_delta_Q <- lapply(1:length(list_H), function(i) {
            if (is.null(list_Q[[i]]))
              return(NULL)
            else {
              temp_i_plus <- first(which(!sapply(list_Q[(i+1):length(list_Q)], is.null)))  # search for the first non-null loc after i
              return(list_Q[[i+temp_i_plus]] - list_Q[[i]])
            }
          })
          list_EIC <- lapply(1:length(list_H), function(i) {
            if (is.null(list_H)) NULL else
              list_H[[i]]*list_delta_Q[[i]]
          })
          names(list_EIC) <- temp_node_names

          # last column might be needed for some tmle update functions
          list_EIC[[length(list_EIC) + 1]] <- do.call(cbind, list_EIC)
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
        } else {
          # for library tasks; it's only needed in tlik updates, with single node
          if (is.null(node)) stop("Please specify single update node for library tasks")

          tmle_task_backup <- tmle_task
          tmle_task <- self$observed_likelihood$training_task  # let tmle_task be obs task when calculating for library tasks
          loc_node <- which(names(tmle_task$npsem) == node)
          obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
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
                                          list_all_predicted_lkd  # this is decided above by fold_number
          )
          current_Q <- get_current_Q(loc_node, which_Q = 0,
                                     tmle_task, obs_variable_names,
                                     intervention_variables, intervention_levels_treat, intervention_levels_control,
                                     list_all_predicted_lkd  # this is decided above by fold_number
          )
          current_delta_Q <- current_Q_next - current_Q
          current_EIC <- current_H*current_delta_Q

          current_EIC <- list(current_EIC)
          names(current_EIC) <- node

          return(current_EIC)
        }

      }
    },
    estimates = function(tmle_task = NULL, fold_number = "full", update = F) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # todo: extend for stochastic
      cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)
      # ZW todo: to handle long format and wide format

      if (fold_number == "full") {
        list_EIC <- private$.list_EIC
        result <- private$.result
      } else if (fold_number == "validation") {
        list_EIC <- private$.list_EIC_val
        result <- private$.result_val
      }

      # load full_p list first
      full_task <- self$observed_likelihood$training_task
      full_node_names <- names(full_task$npsem)
      full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
      full_variable_names <- colnames(full_data)
      list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- full_task$npsem[[loc_node]]$variables
          temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
          temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
          temp_target_node <- intersect(self$update_nodes, full_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            setattr(temp_task, "target_nodes", full_node_names[loc_node])
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            setattr(temp_task, "target_nodes", "no_update")
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input, output = temp_output) %>% return
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
        # nodes to integrate out in the target identification
        # only support univaraite node for now; assume treatment level is one
        all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                             rule_variables = c(intervention_variables,
                                                                last(obs_variable_names)),
                                             rule_values = c(intervention_levels_treat, 1))
        all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                             rule_variables = c(intervention_variables,
                                                                last(obs_variable_names)),
                                             rule_values = c(intervention_levels_control, 1))

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
                                     pmap_dbl(temp_list, prod) %>% sum %>% return
                                   })
        )
        # substitution estimator
        vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output %>% suppressMessages
        psi <- mean(vec_est)

        list_EIC <- self$clever_covariates(tmle_task, fold_number, submodel_type = "EIC")
        list_EIC[[1]] <- vec_est - psi
        list_EIC <- list_EIC[-which(names(list_EIC) == "IC")]

        vec_D <- list_EIC %>% compact %>% pmap_dbl(sum)
        IC <- vec_D

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
      loc_A <- grep("A", temp_node_names)
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      nodes_to_update <- temp_node_names[if_not_0 & !((1:length(temp_node_names)) %in% loc_A)]
      # nodes_to_update <- nodes_to_update[-length(nodes_to_update)]
      return(nodes_to_update)
    },
    list_EIC = function() {
      return(private$.list_EIC)
    },
    list_EIC_val = function() {
      return(private$.list_EIC_val)
    }
  ),
  private = list(
    .type = "mediation",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .list_EIC = NULL,  # the clever covariates as the EIC
    .list_EIC_val = NULL,
    .result = NULL,
    .result_val = NULL,
    .submodel_type_supported = c("EIC")
  )
)











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
Param_mediation_survival <- R6Class(
  classname = "Param_mediation",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control) {
      # outcome_node is used to check self$supports_outcome_censoring; not checked for now
      super$initialize(observed_likelihood, list())
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
      loc_A <- grep("A", temp_node_names)  # not used here; it can include both A_E and A_C
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

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
          loc_impute <- grep("Y_|A_C_", full_node_names)  # remain alive and uncensored before current variable
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

        loc_impute <- grep("Y_|A_C_", temp_node_names)  # in integral, Y and A_C always 1

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
                                     pmap_dbl(temp_list, prod) %>% sum %>% return
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
    .submodel_type_supported = c("EIC")
  )
)




#' Longitudinal Mediation via projection

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
Param_mediation_projection <- R6Class(
  classname = "Param_mediation_projection",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y", static_likelihood = NULL, n_resampling = NULL) {
      if(inherits(observed_likelihood, "Targeted_Likelihood")){
        fold_number <- observed_likelihood$updater$update_fold
      } else {
        fold_number <- "full"
      }

      temp_node_names <- names(observed_likelihood$training_task$npsem)
      loc_A <- grep("A", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      tau <- last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2]))

      all_nodes <- names(observed_likelihood$training_task$npsem)
      A_nodes <- grep("A", all_nodes, value = T)
      Z_nodes <- grep("Z", all_nodes, value = T)
      RLY_nodes <- grep("(R|L|Y).[1-9]$", all_nodes, value = T)

      private$.static_likelihood <- static_likelihood
      private$.update_nodes <- c(Z_nodes, RLY_nodes)

      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)

      tmle_task <- observed_likelihood$training_task
      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)

      if (!is.null(n_resampling)) {  # use expanded Monte Carlo samples to train the HAL projection
        # temp_input <- generate_Zheng_data(B = n_resampling, tau = 1, if_LY_misspec = T) %>% data.frame()
        # names(temp_input)[grep("L1_1", names(temp_input))] <- "L_1"

        temp_input <- tmle_task$get_tmle_node(temp_node_names[1])
        temp_input <- rbind(
          temp_input[c(
            sample(nrow(temp_input), abs(round(n_resampling)), replace = T)
            # 1:n_resampling
          ), ]
        )
        for (i in 2:length(temp_node_names)) {
          temp_input <- cbind(temp_input, 1) %>% as.data.frame()
          names(temp_input)[ncol(temp_input)] <- temp_node_names[i]
          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:i])
          temp_input[,ncol(temp_input)] <- rbinom(temp_task$nrow, 1, static_likelihood$get_likelihood(temp_task, node = temp_node_names[i]))
        }


        # temp_input <- expand_values(variables = obs_variable_names)  # all possible inputs
        # temp_input <- obs_data
        # temp_input <- rbind(
        #   # obs_data,
        #                     temp_input[c(
        #                            # sample(nrow(temp_input), abs(round(n_resampling)), replace = T)
        #                       1:n_resampling
        #                            ), ]
        # )
        temp_input <- data.frame(temp_input,
                                 id = 1:nrow(temp_input),
                                 t = 0)
        # temp_input <- rbind(temp_input, obs_data)
        temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem)

        # private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(temp_task)[[1]]
        # private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(temp_task)[[1]]
        # Train the gradient
        private$.gradient <- Gradient$new(observed_likelihood,
                                          ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                          cf_likelihood_control = self$cf_likelihood_control,
                                                          intervention_list_treatment = self$intervention_list_treatment,
                                                          intervention_list_control = self$intervention_list_control,
                                                          # cf_task_treatment = self$cf_task_treatment,
                                                          # cf_task_control = self$cf_task_control,
                                                          static_likelihood = self$static_likelihood
                                          ),
                                          projection_task_generator = gradient_generator_middle,
                                          target_nodes =  self$update_nodes)

        private$.gradient$train_projections(temp_task, fold_number = fold_number)
      } else {
        # todo: extend for stochastic
        # private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
        # private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
        # Train the gradient
        private$.gradient <- Gradient$new(observed_likelihood,
                                          ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                          cf_likelihood_control = self$cf_likelihood_control,
                                                          intervention_list_treatment = self$intervention_list_treatment,
                                                          intervention_list_control = self$intervention_list_control,
                                                          # cf_task_treatment = self$cf_task_treatment,
                                                          # cf_task_control = self$cf_task_control,
                                                          static_likelihood = self$static_likelihood
                                          ),
                                          projection_task_generator = gradient_generator_middle,
                                          target_nodes =  self$update_nodes)

        private$.gradient$train_projections(self$observed_likelihood$training_task, fold_number = fold_number)
      }









      setattr(self$observed_likelihood, "target_nodes", self$update_nodes)
      self$observed_likelihood$get_likelihoods(self$observed_likelihood$training_task, fold_number = fold_number)
      for (node in self$update_nodes) {
        temp_long_task <- private$.gradient$expand_task(observed_likelihood$training_task, node)
        self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
        self$observed_likelihood$get_likelihood(observed_likelihood$training_task, node, fold_number)
        # private$.gradient$expand_task(private$.cf_task_treatment, node)
        # private$.gradient$expand_task(private$.cf_task_control, node)
      }

      list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- tmle_task$npsem[[loc_node]]$variables
          temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
          temp_target_node <- intersect(self$update_nodes, temp_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            # for each short task, only the last node (if it is an update_node) needs to be updated
            setattr(temp_task, "target_nodes", temp_target_node)
            for (node in attr(temp_task, "target_nodes")) {
              temp_long_task <- private$.gradient$expand_task(temp_task, node)
              self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
            }
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            # A nodes won't get updated
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input, output = temp_output) %>% return
        }
      })

    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", node = NULL) {
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
        return(self$gradient$compute_component(tmle_task, node, fold_number = fold_number)$EIC)
      })

      names(EICs) <- update_nodes

      EICs[[length(EICs) + 1]] <- do.call(cbind, EICs)
      colnames(EICs[[length(EICs)]]) <- update_nodes
      EICs[[length(EICs)]] <- as.data.frame(EICs[[length(EICs)]])

      names(EICs)[length(EICs)] <- "IC"

      return(EICs)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # clever_covariates happen here (for this param) only, but this is repeated computation
      EIC <- (do.call(cbind, self$clever_covariates(tmle_task, fold_number)$IC))

      #TODO need to montecarlo simulate from likleihood to eval parameter.

      temp_node_names <- names(self$observed_likelihood$training_task$npsem)
      loc_A <- grep("A", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      tau <- last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2]))

      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)
      list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- tmle_task$npsem[[loc_node]]$variables
          temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
          temp_target_node <- intersect(self$update_nodes, temp_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input, output = temp_output) %>% return
        }
      })

      intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
      intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
      intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
      intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
      # nodes to integrate out in the target identification
      # only support univaraite node for now; assume treatment level is one
      all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                           rule_variables = c(intervention_variables,
                                                              last(obs_variable_names)),
                                           rule_values = c(intervention_levels_treat, 1))
      all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                           rule_variables = c(intervention_variables,
                                                              last(obs_variable_names)),
                                           rule_values = c(intervention_levels_control, 1))
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
                                   pmap_dbl(temp_list, prod) %>% sum %>% return
                                 })
      )
      vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output %>% suppressMessages
      psi <- mean(vec_est)

      EIC <- cbind(EIC, vec_est - psi)

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
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
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
    gradient = function(){
      private$.gradient
    },
    static_likelihood = function(){
      private$.static_likelihood
    }
  ),
  private = list(
    .type = "mediation",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .cf_task_treatment = NULL,
    .cf_task_control = NULL,
    .supports_outcome_censoring = FALSE,
    .gradient = NULL,
    .submodel_type_supported = c("EIC"),
    .update_nodes = NULL,
    .static_likelihood = NULL
  )
)












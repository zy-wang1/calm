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
#'   \code{define_param(Param_middle, observed_likelihood, intervention_list, ..., outcome_node)}
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
Param_middle <- R6Class(
  classname = "Param_middle",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
      # observed_likelihood$get_likelihoods(observed_likelihood$training_task)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", update = F) {
      if (fold_number == "full") {  # tmle
        list_newH <- private$.list_newH
      } else if (fold_number == "validation") {  # cvtmle
        list_newH <- private$.list_newH_val
      }
      if (!is.null(list_newH) & update == F) {
        return(list_newH)
      } else {
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
        loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
        if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

        obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
        obs_variable_names <- colnames(obs_data)
        # ZW todo: to handle long format and wide format

        if (fold_number == "full") {
          list_all_predicted_lkd <- self$observed_likelihood$list_all_predicted_lkd # if it's null, it will be recalculated in observed_likelihood
        } else if (fold_number == "validation") {
          list_all_predicted_lkd <- self$observed_likelihood$list_all_predicted_lkd_val
          temp <- self$observed_likelihood$list_all_predicted_lkd  # calculate full version if not already
          rm(temp)
        }

        intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
        intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
        intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
        intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
        names(intervention_levels_treat) <- names(self$intervention_list_treatment)
        names(intervention_levels_control) <- names(self$intervention_list_control)

        list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = self$observed_likelihood,
                            cf_task_treatment, cf_task_control,
                            intervention_variables, intervention_levels_treat, intervention_levels_control,
                            fold_number)
        list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H,
                              intervention_variables, intervention_levels_treat, intervention_levels_control,
                              list_all_predicted_lkd,  # val version decided above for fold_number == "validation"
                              lt = 1)
        list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H,
                              intervention_variables, intervention_levels_treat, intervention_levels_control,
                              list_all_predicted_lkd,
                              lt = 0)

        list_newH <- list()
        for (loc_node in 1:length(list_H)) {
          if(!is.null(list_H[[loc_node]])) {
            list_newH[[loc_node]] <- ( list_H[[loc_node]] * (list_Q_1[[loc_node]] - list_Q_0[[loc_node]]) ) %>% as.matrix
          }
        }
        names(list_newH) <- temp_node_names

        if (fold_number == "full") {
          private$.list_newH <- list_newH
        } else if (fold_number == "validation") {
          private$.list_newH_val <- list_newH
        }

        return(list_newH)
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
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)
      # ZW todo: to handle long format and wide format

      if (fold_number == "full") {
        list_all_predicted_lkd <- self$observed_likelihood$list_all_predicted_lkd
        list_D <- private$.list_D
        result <- private$.result
      } else if (fold_number == "validation") {
        list_all_predicted_lkd <- self$observed_likelihood$list_all_predicted_lkd_val
        temp <- self$observed_likelihood$list_all_predicted_lkd
        rm(temp)
        list_D <- private$.list_D_val
        result <- private$.result_val
      }

      # make sure we force update this after each updating step; this helps speed up updater$check_convergence
      if (!is.null(list_D) & !is.null(result) & update == F) {
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
                                     temp_all_comb_0 <- cbind(unique_L0[which_row, ], all_possible_RZLY_0)
                                     temp_all_comb_1 <- cbind(unique_L0[which_row, ], all_possible_RZLY_1)
                                     # for all non-A, non-0 variables, calculate the variable by rule
                                     # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                     # note that list_all_predicted_lkd is ordered by node
                                     temp_list_0 <- lapply(loc_Z,
                                                           function(each_t) {
                                                             left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                           })
                                     temp_list_1 <- lapply(loc_RLY,
                                                           function(each_t) {
                                                             left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                           })
                                     temp_list <- c(temp_list_0, temp_list_1)
                                     pmap_dbl(temp_list, prod) %>% sum %>% return
                                   })
        )
        # substitution estimator
        vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output
        psi <- mean(vec_est)

        list_newH <- self$clever_covariates(tmle_task, fold_number)

        list_D <- list_D_null <- list_D_fit <- list()
        for (loc_node in 1:length(list_newH)) {
          if(!is.null(list_newH[[loc_node]])) {
            # ZW todo: for discretized variables
            current_ind <- (obs_data[[tmle_task$npsem[[loc_node]]$variables]] == 1)*1
            # temp_p <- self$observed_likelihood$get_likelihoods(tmle_task, temp_node_names[loc_node])
            if (loc_node %in% loc_Z) temp_p <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[loc_node], fold_number) else
              temp_p <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[loc_node], fold_number)
            temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)
            list_D[[loc_node]] <- (current_ind - temp_p) *list_newH[[loc_node]]
            list_D_fit[[loc_node]] <- (1 - temp_p) *list_newH[[loc_node]]  # e.g. L nodes, ind_A all 1's for treat cf task
            list_D_null[[loc_node]] <- rep(0, tmle_task$nrow)
          }
        }
        list_D[[1]] <-
        # list_D_trt[[1]] <- list_D_ctrl[[1]] <-
        vec_est - psi

        names(list_D) <- names(list_D_null) <- names(list_D_fit) <- names(list_newH)

        vec_D <- list_D %>% compact %>% pmap_dbl(sum)
        IC <- vec_D

        result <- list(psi = psi, IC = IC
                       # , full_IC = list_D
                       )

        # these are cached; unless likelihood is updated, or we force it to update, they shouldn't be changed
        if (fold_number == "full") {
          private$.list_D <- list_D
          private$.list_D_null <- list_D_null
          private$.list_D_fit <- list_D_fit
          private$.result <- result
        } else if (fold_number == "validation") {
          private$.list_D_val <- list_D
          private$.list_D_null_val <- list_D_null
          private$.list_D_fit_val <- list_D_fit
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
    list_D = function() {
      return(private$.list_D)
    },
    list_D_null = function() {
      return(private$.list_D_null)
    },
    list_D_fit = function() {
      return(private$.list_D_fit)
    },
    list_D_val = function() {
      return(private$.list_D_val)
    },
    list_D_null_val = function() {
      return(private$.list_D_null_val)
    },
    list_D_fit_val = function() {
      return(private$.list_D_fit_val)
    },
    update_nodes = function() {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      nodes_to_update <- temp_node_names[if_not_0 & !((1:length(temp_node_names)) %in% loc_A)]
      # nodes_to_update <- nodes_to_update[-length(nodes_to_update)]
      return(nodes_to_update)
    },
    list_newH = function() {
      return(private$.list_newH)
    },
    list_newH_val = function() {
      return(private$.list_newH_val)
    }
  ),
  private = list(
    .type = "middle",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .list_newH = NULL,  # the clever covariates
    .list_newH_val = NULL,
    .list_D = NULL,
    .list_D_null = NULL,
    .list_D_fit = NULL,
    .list_D_val = NULL,
    .list_D_null_val = NULL,
    .list_D_fit_val = NULL,
    .result = NULL,
    .result_val = NULL
  )
)


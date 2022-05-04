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
Param_middle_projection <- R6Class(
  classname = "Param_middle_projection",
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
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
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
      vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output
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
    .type = "ATE",
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
Param_middle_projection_survival <- R6Class(
  classname = "Param_middle_projection",
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
      temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
      loc_A_E <- grep("A_E", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      tau <- as.numeric(last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2])))

      A_E_nodes <- grep("A_E", temp_node_names, value = T)
      Z_nodes <- grep("Z", temp_node_names, value = T)
      RLY_nodes <- grep("(R|L|Y).[1-9]$", temp_node_names, value = T)

      private$.static_likelihood <- static_likelihood
      private$.update_nodes <- c(Z_nodes, RLY_nodes)

      private$.supports_outcome_censoring <- T
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)

      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)

      tmle_task <- observed_likelihood$training_task
      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))
      obs_variable_names <- colnames(obs_data)
      loc_A_C_or_Y <- grep("A_C_|Y_", temp_node_names)

      if (!is.null(n_resampling)) {  # use expanded Monte Carlo samples to train the HAL projection
        # temp_input <- generate_Zheng_data(B = n_resampling, tau = 1, if_LY_misspec = T) %>% data.frame()
        # names(temp_input)[grep("L1_1", names(temp_input))] <- "L_1"
        temp_input <- tmle_task$get_tmle_node(temp_node_names[1])
        temp_input <- rbind(
          temp_input[c(
            sample(nrow(temp_input), abs(round(n_resampling)), replace = T)
          ), ]
        )
        for (i in 2:length(temp_node_names)) {
          temp_input <- cbind(temp_input, 1) %>% as.data.frame()
          names(temp_input)[ncol(temp_input)] <- temp_node_names[i]
          delta_vars <- names(sapply(paste0("delta_", temp_node_names[1:i]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
          if (length(delta_vars) > 0) {
            full_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
            colnames(full_input)[(ncol(temp_input) + 1):ncol(full_input)] <- delta_vars
          } else {
            full_input <- temp_input
          }
          temp_task <- tmle3_Task_drop_censored$new(full_input, tmle_task$npsem[c(1:i,
                                                                    sapply(paste0("delta_", temp_node_names[1:i]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
                                                                    )])
          temp_input[,ncol(temp_input)] <- rbinom(temp_task$nrow, 1, static_likelihood$get_likelihood(temp_task, node = temp_node_names[i]))
        }
        loc_A_C_or_Y_needed <- loc_A_C_or_Y[loc_A_C_or_Y < i]
        loc_current_var <- which(colnames(temp_input) == tmle_task$npsem[[i]]$variables)
        if (length(loc_A_C_or_Y_needed) > 0) for (s in loc_A_C_or_Y_needed) {
          loc_var <- which(colnames(temp_input) == tmle_task$npsem[[s]]$variables)
          if_censored <- which(temp_input[, loc_var] == 0)
          if_censored <- if_censored[!is.na(if_censored)]
          temp_input[if_censored, (loc_var+1):loc_current_var] <- NA
        }
        for (name in colnames(temp_input)) {
          if(any(is.na(temp_input[, name]))) {
            temp_input <- cbind(temp_input, !is.na(temp_input[, name]))
            colnames(temp_input)[ncol(temp_input)] <- paste0("delta_", name)
          }
        }
        temp_input <- data.frame(temp_input,
                                 id = 1:nrow(temp_input),
                                 t = 0)
        # temp_input <- rbind(temp_input, obs_data)
        temp_task <- tmle3_Task_drop_censored$new(temp_input, tmle_task$npsem)

        # private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(temp_task)[[1]]
        # private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(temp_task)[[1]]
        # Train the gradient
        private$.gradient <- Gradient_wide$new(observed_likelihood,
                                          ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                          cf_likelihood_control = self$cf_likelihood_control,
                                                          intervention_list_treatment = self$intervention_list_treatment,
                                                          intervention_list_control = self$intervention_list_control,
                                                          # cf_task_treatment = self$cf_task_treatment,
                                                          # cf_task_control = self$cf_task_control,
                                                          static_likelihood = self$static_likelihood
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
                                                          intervention_list_treatment = self$intervention_list_treatment,
                                                          intervention_list_control = self$intervention_list_control,
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
          loc_impute <- grep("Y_|A_C_", temp_node_names)  # remain alive and uncensored before current variable
          loc_impute <- loc_impute[loc_impute < loc_node]
          if (length(loc_impute) == 0) {
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
          } else {
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)],
                                        rule_variables = sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables),
                                        rule_values = rep(1, length(loc_impute))
                                        )  # all possible inputs
          }
          delta_vars <- names(sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
          if (length(delta_vars) > 0) {
            temp_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
            colnames(temp_input)[(ncol(temp_input) - length(delta_vars) + 1):ncol(temp_input)] <- delta_vars
          }

          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[c(1:loc_node,
                                                                    sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
                                                                    )])
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

      EICs_impute <- lapply(update_nodes, function(node){
        temp_vec <- self$gradient$compute_component(tmle_task, node, fold_number = fold_number)$EIC
        if (!is.null(tmle_task$npsem[[node]]$censoring_node)) {
          full_vec <- if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[node]]$censoring_node$variables)
          full_vec[if_observed] <- temp_vec
          full_vec[!if_observed] <- 0
          temp_vec <- full_vec
        }
        return(temp_vec)
      })
      EICs[[length(EICs) + 1]] <- do.call(cbind, EICs_impute)
      EICs[[length(EICs)]]
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
      temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
      loc_A_E <- grep("A_E", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      tau <- last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2]))

      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))
      obs_variable_names <- colnames(obs_data)
      list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- tmle_task$npsem[[loc_node]]$variables
          loc_impute <- grep("Y_|A_C_", temp_node_names)  # remain alive and uncensored before current variable
          loc_impute <- loc_impute[loc_impute < loc_node]
          if (length(loc_impute) == 0) {
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
          } else {
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)],
                                        rule_variables = sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables),
                                        rule_values = rep(1, length(loc_impute))
            )  # all possible inputs
          }
          delta_vars <- names(sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
          if (length(delta_vars) > 0) {
            temp_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
            colnames(temp_input)[(ncol(temp_input) - length(delta_vars) + 1):ncol(temp_input)] <- delta_vars
          }

          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[c(1:loc_node,
                                                                    sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
          )])
          temp_target_node <- intersect(self$update_nodes, temp_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            # A nodes won't get updated
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input[1:which(obs_variable_names == current_variable)], output = temp_output) %>% return
        }
      })

      intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
      intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
      intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
      intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
      # nodes to integrate out in the target identification
      # only support univaraite node for now; assume treatment level is one
      loc_impute <- grep("Y_|A_C_", temp_node_names)  # in integral, Y and A_C always 1
      all_possible_RZLY_1 <- expand_values(variables = obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                  rule_variables = c(intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                  rule_values = c(intervention_levels_treat, rep(1, length(loc_impute)))
      )
      all_possible_RZLY_0 <- expand_values(variables = obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                           rule_variables = c(intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                           rule_values = c(intervention_levels_control, rep(1, length(loc_impute)))
      )
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
      vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output
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
    .type = "middle_survival",
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

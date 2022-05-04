#' Defines a TML Estimator (except for the data)
#'
#' Longitudinal Mediation Targets
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_med <- R6Class(
  classname = "tmle3_Spec_med",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(treatment_level, control_level, ...) {
      super$initialize(
        treatment_level = treatment_level,
        control_level = control_level, ...
      )
    },
    make_tmle_task = function(data, node_list, if_drop_censored = NULL, ...) {
      variable_types <- self$options$variable_types
      if(is.null(if_drop_censored)) tmle_task <- middle_task(data, node_list, variable_types)
      if(if_drop_censored == T) tmle_task <- middle_task_drop_censored(data, node_list, variable_types)
      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # produce trained likelihood when likelihood_def provided

      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- middle_likelihood(tmle_task, learner_list)  # see middle_helper
      }

      return(likelihood)
    },
    make_params = function(tmle_task, likelihood, options = NULL, if_projection = NULL, static_likelihood = NULL, n_resampling = NULL) {
      if (is.null(static_likelihood)) static_likelihood <- likelihood
      if (is.null(options)) options <- list("tc")

      temp_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_names)
      # ZW todo: in future can be dynamic
      A_levels <- tmle_task$npsem[[ temp_names[loc_A[1]] ]]$variable_type$levels

      tmle_params <-lapply(options, function(option) {
        if (option == "tc") {
          treatment_value <- self$options$treatment_level
          control_value <- self$options$control_level
        } else if (option == "tt") {
          treatment_value <- self$options$treatment_level
          control_value <- self$options$treatment_level
        } else if (option == "cc") {
          treatment_value <- self$options$control_level
          control_value <- self$options$control_level
        }  # decide the target of inference
        if (!is.null(A_levels)) {
          treatment_value <- treatment_value
          # factor(treatment_value, levels = A_levels)
          control_value <- control_value
          # factor(control_value, levels = A_levels)
        }
        # list of intervention nodes as LF_static objects
        treatment <- lapply(temp_names[loc_A], function(eachA) {
          define_lf(LF_static, eachA, value = treatment_value)
        })
        control <- lapply(temp_names[loc_A], function(eachA) {
          define_lf(LF_static, eachA, value = control_value)
        })
        names(treatment) <- names(control) <- temp_names[loc_A]
        if (is.null(if_projection)) {
          param <- Param_med$new(likelihood, treatment, control, outcome_node = last(temp_names))
        } else if (if_projection) {
          param <- Param_middle_projection$new(likelihood, treatment, control, outcome_node = last(temp_names), static_likelihood, n_resampling)
        }
        return(param)
      })
      return(tmle_params)
    },
    make_params_survival = function(tmle_task, likelihood, options = NULL, if_projection = NULL, static_likelihood = NULL, n_resampling = NULL) {
      if (is.null(static_likelihood)) static_likelihood <- likelihood
      if (is.null(options)) options <- list("tc")

      temp_names <- names(tmle_task$npsem)
      temp_names <- temp_names[-grep("delta_", temp_names)]  # remove delta nodes for wide format fitting
      loc_A_E <- grep("A_E", temp_names)
      # ZW todo: in future can be dynamic
      A_levels <- tmle_task$npsem[[ temp_names[loc_A_E[1]] ]]$variable_type$levels

      tmle_params <-lapply(options, function(option) {
        if (option == "tc") {
          treatment_value <- self$options$treatment_level
          control_value <- self$options$control_level
        } else if (option == "tt") {
          treatment_value <- self$options$treatment_level
          control_value <- self$options$treatment_level
        } else if (option == "cc") {
          treatment_value <- self$options$control_level
          control_value <- self$options$control_level
        }  # decide the target of inference
        if (!is.null(A_levels)) {
          treatment_value <- treatment_value
          # factor(treatment_value, levels = A_levels)
          control_value <- control_value
          # factor(control_value, levels = A_levels)
        }
        # list of intervention nodes as LF_static objects
        treatment <- lapply(temp_names[loc_A_E], function(eachA) {
          define_lf(LF_static, eachA, value = treatment_value)
        })
        control <- lapply(temp_names[loc_A_E], function(eachA) {
          define_lf(LF_static, eachA, value = control_value)
        })
        names(treatment) <- names(control) <- temp_names[loc_A_E]
        if (is.null(if_projection)) {
          param <- Param_med_survival$new(likelihood, treatment, control, outcome_node = last(temp_names))
        } else if (if_projection) {
          param <- Param_middle_projection_survival$new(likelihood, treatment, control, outcome_node = last(temp_names), static_likelihood, n_resampling)
        }
        return(param)
      })
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' Longitudinal Mediation Targets
#'
#' O=(L0, A1, R1, Z1, L1, (Y1), ..., An, Rn, Zn, Ln, Yn)
#' L0=Baseline covariates
#' A: Treatment (binary or categorical)
#' Z: Mediators
#' R: (time-varying) Covariates before mediator
#' L: (time-varying) Covariates after mediator
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
tmle_med <- function(treatment_level, control_level) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_med$new(treatment_level, control_level)
}

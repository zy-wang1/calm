#' average density, or integral of squared desnity function
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
#'     \item{\code{...}}{Not currently used.
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{...}}{... from \code{\link{LF_base}} ...
#'     }
#' }
#' @export
Param_ave_dens_tilde <- R6Class(
  classname = "Param_ave_dens_tilde",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, P_likelihood, choose_grid = seq(-20, 20, by = 0.01)) {
      super$initialize(observed_likelihood, list())
      private$.choose_grid <- choose_grid
      private$.P_likelihood <- P_likelihood
      # observed_likelihood$get_likelihoods(observed_likelihood$training_task)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", node = NULL, submodel_type = "EIC") {
      if (is.null(tmle_task)) {  # calculate for obs data task if not specified
        tmle_task <- self$observed_likelihood$training_task
      }
      temp_node_names <- names(tmle_task$npsem)

      # prepare pmf
      support_points <- self$choose_grid
      new_data <- data.table(one = 0,
                             outcome = support_points)
      new_task <- tmle3_Task$new(new_data, tmle_task$npsem)

      new_dens <- data.frame(outcome = new_data$outcome,
                             cont_dens = self$observed_likelihood$get_likelihood(new_task, "outcome", fold_number))
      temp <- c(1)
      for (i in 2:nrow(new_data)) {
        if (new_dens$cont_dens[i] != new_dens$cont_dens[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
      }
      new_dens[["group"]] <- temp

      dat_int_tilde <- left_join(new_dens,
                                 data.frame(outcome = support_points,
                                            dens_p = self$P_likelihood$get_likelihood(tmle_task = new_task, fold_number = fold_number, node = "outcome"))
      ) %>% mutate(value = dens_p) %>% mutate(to_sum = cont_dens*value*(diff(self$choose_grid))[1])
      int_tilde <- sum(dat_int_tilde$to_sum)

      p = self$P_likelihood$get_likelihood(tmle_task = tmle_task, fold_number = fold_number, node = "outcome")

      IC <- 2*p - 2*int_tilde

      list_D <- lapply(temp_node_names, function(name) {
        if (length(unique(tmle_task$get_tmle_node(name))) == 1) return(NULL) else {
          return(IC)
        }
      })
      names(list_D) <- temp_node_names

      if (!is.null(node)) return(list_D[node]) else return(list_D)

    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      obs_data <- tmle_task$data %>% dplyr::select(-c(id, t))
      # new_dens <- data.frame(obs_data,
      #                        cont_dens = initial_likelihood$get_likelihood(tmle_task, "outcome")) %>% arrange(outcome)

      # descritize;
      {
        support_points <- self$choose_grid
        new_data <- data.table(one = 0,
                               outcome = support_points)
        new_task <- tmle3_Task$new(new_data, npsem)
        new_dens <- data.frame(outcome = new_data$outcome,
                               cont_dens = self$observed_likelihood$get_likelihood(new_task, "outcome", fold_number))

        temp <- c(1)
        for (i in 2:nrow(new_data)) {
          if (new_dens$cont_dens[i] != new_dens$cont_dens[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
        }
        new_dens[["group"]] <- temp

        psi <- sum((new_dens$cont_dens^2)*(diff(self$choose_grid))[1])
      }

      IC <- self$clever_covariates(tmle_task = tmle_task, fold_number = fold_number, node = "outcome")$outcome
      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("Expectation of p(O) = Int p(o)^2 do")
      return(param_form)
    },
    update_nodes = function() {
      tmle_task <- self$observed_likelihood$training_task
      temp_node_names <- names(tmle_task$npsem)
      if (length(unique(tmle_task$data[[1]])) == 1) nodes_to_update <- temp_node_names[-1] else nodes_to_update <- temp_node_names
      return(nodes_to_update)
    },
    choose_grid = function() {
      return(private$.choose_grid)
    },
    P_likelihood = function() return(private$.P_likelihood)
  ),
  private = list(
    .type = "ave_dens_tilde",
    .submodel_type_supported = c("EIC"),
    .choose_grid = NULL,
    .P_likelihood = NULL
  )
)


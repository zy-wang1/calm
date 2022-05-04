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
Param_ave_dens <- R6Class(
  classname = "Param_ave_dens",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, choose_grid = seq(-20, 20, by = 0.01)) {
      super$initialize(observed_likelihood, list())
      private$.choose_grid <- choose_grid
      # observed_likelihood$get_likelihoods(observed_likelihood$training_task)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", node = NULL, submodel_type = "EIC") {
      if (is.null(tmle_task)) {  # calculate for obs data task if not specified
        tmle_task <- self$observed_likelihood$training_task
      }
      temp_node_names <- names(tmle_task$npsem)

      result <- self$estimates(fold_number = fold_number)
      psi <- result$psi
      IC <- result$IC

      list_D <- lapply(temp_node_names, function(name) {
        if (length(unique(tmle_task$get_tmle_node(name))) == 1) return(NULL) else {
          # if (identical(tmle_task, self$observed_likelihood$training_task)) return(IC) else {
            IC <- self$observed_likelihood$get_likelihood(tmle_task, "outcome", fold_number)*2 - psi*2
            return(IC)
          # }
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

      # descritize;
      {
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

        # new_pmf <- sapply(new_dens$group %>% unique, function(which_group) {
        #   temp <- which(new_dens$group == which_group)
        #   a1 <- first(temp)
        #   a2 <- last(temp)
        #   c(new_dens$outcome[a1],
        #     (new_dens$outcome[a2] - new_dens$outcome[a1]) *
        #       new_dens$cont_dens[a1],
        #     new_dens$cont_dens[a1]
        #     )
        # }) %>% t
        # # temp_loc <- last(which(new_pmf[, 2] !=0))
        # # new_pmf[temp_loc, 2] <- 1-sum(new_pmf[-temp_loc, 2])
        # # new_pmf[, 2] %>% sum
        # # psi <- sum(new_pmf[, 2]^2)
        # psi <- sum(new_pmf[, 2] * new_pmf[, 3])
        psi <- sum((new_dens$cont_dens^2)*(diff(self$choose_grid))[1])
      }

      IC <- self$observed_likelihood$get_likelihood(tmle_task, "outcome", fold_number)*2 - psi*2
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
    }
  ),
  private = list(
    .type = "ave_dens",
    .submodel_type_supported = c("EIC"),
    .choose_grid = NULL
  )
)


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
Param_ave_dens_2 <- R6Class(
  classname = "Param_ave_dens_2",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, targeted_likelihood_1, updater_1, tilde_likelihood, choose_grid = seq(-20, 20, by = 0.01)) {
      super$initialize(observed_likelihood, list())
      private$.choose_grid <- choose_grid
      private$.targeted_likelihood_1 <- targeted_likelihood_1
      private$.updater_1 <- updater_1
      private$.tilde_likelihood <- tilde_likelihood
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
      tilde_dens <- data.frame(outcome = new_data$outcome,
                               cont_dens = self$tilde_likelihood$get_likelihood(new_task, "outcome", fold_number))
      temp <- c(1)
      for (i in 2:nrow(new_data)) {
        if (tilde_dens$cont_dens[i] != tilde_dens$cont_dens[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
      }
      tilde_dens[["group"]] <- temp
      ini_dens <- data.frame(outcome = new_data$outcome,
                             cont_dens = self$observed_likelihood$get_likelihood(new_task, "outcome", fold_number))
      temp <- c(1)
      for (i in 2:nrow(new_data)) {
        if (ini_dens$cont_dens[i] != ini_dens$cont_dens[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
      }
      ini_dens[["group"]] <- temp
      temp <- c(1)
      for (i in 2:length(ini_dens$group)) {
        if (ini_dens$group[i] != ini_dens$group[i-1] | tilde_dens$group[i] != tilde_dens$group[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
      }
      tilde_dens$group <- ini_dens$group <- temp

      # tilde_pmf <- sapply(tilde_dens$group %>% unique, function(which_group) {
      #   temp <- which(tilde_dens$group == which_group)
      #   a1 <- first(temp)
      #   a2 <- last(temp)
      #   c(tilde_dens$outcome[a1],
      #     (tilde_dens$outcome[a2] - tilde_dens$outcome[a1]) *
      #       tilde_dens$cont_dens[a1],
      #     tilde_dens$cont_dens[a1]
      #   )
      # }) %>% t %>% as.data.frame
      # colnames(tilde_pmf) <- c("outcome", "prob", "dens")
      # # est_tilde <- sum(tilde_pmf$prob^2)
      # est_tilde <- sum(tilde_pmf$prob * tilde_pmf$dens)
      est_tilde <- sum((tilde_dens$cont_dens^2)*(diff(self$choose_grid))[1])

      # ini_pmf <- sapply(ini_dens$group %>% unique, function(which_group) {
      #   temp <- which(ini_dens$group == which_group)
      #   a1 <- first(temp)
      #   a2 <- last(temp)
      #   c(ini_dens$outcome[a1],
      #     (ini_dens$outcome[a2] - ini_dens$outcome[a1]) *
      #       ini_dens$cont_dens[a1],
      #     ini_dens$cont_dens[a1]
      #   )
      # }) %>% t %>% as.data.frame
      # colnames(ini_pmf) <- c("outcome", "prob", "dens")
      # # est_1 <- sum(ini_pmf$prob^2)
      # est_1 <- sum(ini_pmf$prob * ini_pmf$dens)
      est_1 <- sum((ini_dens$cont_dens^2)*(diff(self$choose_grid))[1])

      # param_1 <- self$param_1
      # if (is.null(param_1)) param_1 <- Param_ave_dens$new(observed_likelihood = self$observed_likelihood, choose_grid = self$choose_grid)

      # get cnP
      epsilon_1 <- self$updater_1$epsilons %>% unlist %>% sum
      dat_int_tilde <- left_join(tilde_dens,
                                 # tilde_pmf,
                                 data.frame(outcome = support_points,
                                            D1P = self$observed_likelihood$get_likelihood(new_task, "outcome", fold_number)*2 - est_1*2
                                            # param_1$clever_covariates(tmle_task = new_task, fold_number = fold_number, node = "outcome")$outcome
                                            ,
                                            prob_p = self$observed_likelihood$get_likelihood(tmle_task = new_task, fold_number = fold_number, node = "outcome"))
      ) %>% mutate(value = (D1P)^2 / (1 + epsilon_1 * D1P)^2) %>% mutate(to_sum = cont_dens*value*(diff(self$choose_grid))[1])
      cnP <- sum(dat_int_tilde$to_sum)

      # get dP
      # param_tilde <- self$param_tilde  # the param at the updated likelihood
      # if (is.null(param_tilde)) param_tilde <- Param_ave_dens$new(observed_likelihood = self$targeted_likelihood_1, choose_grid = self$choose_grid)
      dat_int_P <- ini_dens %>%
        # ini_pmf %>%
        left_join(
          data.frame(outcome = new_data$outcome,
                     D1P =
                       self$observed_likelihood$get_likelihood(new_task, "outcome", fold_number)*2 - est_1*2
                     # param_1$clever_covariates(tmle_task = new_task, fold_number = fold_number, node = "outcome")$outcome
                     ,
                     D1P1 =
                       self$targeted_likelihood_1$get_likelihood(new_task, "outcome", fold_number)*2 - est_tilde*2
                     # param_tilde$clever_covariates(tmle_task = new_task, fold_number = fold_number, node = "outcome")$outcome
                     ,
                     prob_p_tilde = self$tilde_likelihood$get_likelihood(tmle_task = new_task, fold_number = fold_number, node = "outcome"))
        ) %>% mutate(to_sum = cont_dens*D1P*D1P1*(diff(self$choose_grid))[1])
      dP <- sum(dat_int_P$to_sum) / cnP


      D1P1 <-
        # param_tilde$clever_covariates(tmle_task = tmle_task, fold_number = fold_number, node = "outcome")$outcome
        self$targeted_likelihood_1$get_likelihood(tmle_task, "outcome", fold_number)*2 - est_tilde*2
      D1P <-
        # param_1$clever_covariates(tmle_task = tmle_task, fold_number = fold_number, node = "outcome")$outcome
        self$observed_likelihood$get_likelihood(tmle_task, "outcome", fold_number)*2 - est_1*2
      line1 <- D1P1 * (1 + epsilon_1 * D1P)

      tilde_p <- self$tilde_likelihood$get_likelihood(tmle_task = tmle_task, node = "outcome", fold_number = fold_number)
      p <- self$observed_likelihood$get_likelihood(tmle_task = tmle_task, node = "outcome", fold_number = fold_number)
      dat_int_tilde <- dat_int_tilde %>% mutate(denominator = 1 / (1 + epsilon_1 * D1P)^2) %>% mutate(to_sum_denominator = denominator * cont_dens * (diff(self$choose_grid))[1])
      tilde_Pn_denominator <- sum(dat_int_tilde$to_sum_denominator)
      line2 <- 2 * dP * (tilde_p / (1 + epsilon_1 * D1P)^2 - 2*p*tilde_Pn_denominator)

      dat_int_P <- dat_int_P %>% mutate(value_line3 = prob_p_tilde / (1 + epsilon_1 * D1P)^2) %>% mutate(to_sum_line3 = cont_dens * value_line3 * (diff(self$choose_grid))[1])
      line3 <- -2 * dP * sum(dat_int_P$to_sum_line3)

      line4 <- 4 * dP *
        # param_1$estimates(fold_number = fold_number)$psi
        est_1 * tilde_Pn_denominator

      dat_int_P <- dat_int_P %>% mutate(to_sum_line5 = cont_dens * D1P1 * (diff(self$choose_grid))[1])
      line5 <- 2 * epsilon_1 * p * (D1P1 - 2 * sum(dat_int_P$to_sum_line5))

      dat_int_P <- dat_int_P %>% mutate(value_line6 = cont_dens * (D1P1 - 2 * sum(dat_int_P$to_sum_line5))) %>% mutate(to_sum_line6 = cont_dens * value_line6 * (diff(self$choose_grid))[1])
      line6 <- -2 * epsilon_1 * sum(dat_int_P$to_sum_line6)

      list_D <- lapply(temp_node_names, function(name) {
        if (length(unique(tmle_task$get_tmle_node(name))) == 1) return(NULL) else {
          IC <- line1 + line2 + line3 + line4 + line5 + line6
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

        # new_pmf <- sapply(new_dens$group %>% unique, function(which_group) {
        #   temp <- which(new_dens$group == which_group)
        #   a1 <- first(temp)
        #   a2 <- last(temp)
        #   c(new_dens$outcome[a1],
        #     (new_dens$outcome[a2] - new_dens$outcome[a1]) *
        #       new_dens$cont_dens[a1],
        #     new_dens$cont_dens[a1]
        #   )
        # }) %>% t
        # temp_loc <- last(which(new_pmf[, 2] !=0))
        # new_pmf[temp_loc, 2] <- 1-sum(new_pmf[-temp_loc, 2])
        # new_pmf[, 2] %>% sum
        # psi <- sum(new_pmf[, 2]^2)
        # psi <- sum(new_pmf[, 2] * new_pmf[, 3])
        psi <- sum((new_dens$cont_dens^2)*(diff(self$choose_grid))[1])
      }

      # {
      #   # new_data <- sample_all(tmle_task, initial_likelihood, 10, "one", "outcome")
      #   psi <- self$observed_likelihood$get_likelihood(tmle_task, "outcome", fold_number) %>% mean
      #
      # }

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
    targeted_likelihood_1 = function() return(private$.targeted_likelihood_1),
    updater_1 = function() return(private$.updater_1),
    tilde_likelihood = function() return(private$.tilde_likelihood),
    param_1 = function() return(private$.param_1),
    param_tilde = function() return(private$.param_tilde)
  ),
  private = list(
    .type = "ave_dens_2",
    .submodel_type_supported = c("EIC"),
    .choose_grid = NULL,
    .targeted_likelihood_1 = NULL,
    .updater_1 = NULL,
    .tilde_likelihood = NULL,
    .param_1 = NULL,
    .param_tilde = NULL
  )
)


library(dplyr)
library(purrr)

# B <- 100  # n
# tau <- 4  # last time point
# seed <- 202008  # random seed

generate_Zheng_data <- function(B, tau, seed = NULL, setAM = NULL, if_LY_misspec = F, if_A_misspec = F) {
  if (is.null(seed)) {} else set.seed(seed)

  # time point 0
  t <- 0
  L01 <- rbinom(n = B, size = 1, prob = 0.4)
  L02 <- rbinom(n = B, size = 1, prob = 0.6)
  temp_data <- list()
  temp_data[[1]] <- data.frame(L1 = L01, L2 = L02)
  t <- t + 1

  if(is.null(setAM)) {
    # time point 1~tau
    # t = 1, ..., tau; to update the t-th time point in the t+1 slot
    while(t <= tau) {
      # omit censoring for now
      if (!if_A_misspec) {
        temp_A <- map_dbl(.x = expit(- 0.1 + 1.2*L02 + 0.7*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) - 0.1*ifelse_vec(t>1, temp_data[[t]]$A, 0)),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      } else {
        temp_A <- map_dbl(.x = expit(- 0.1 + 1.2*L02^2 + 0.7*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) - 0.1*ifelse_vec(t>1, temp_data[[t]]$A^2, 0)),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      }

      temp_R <- map_dbl(.x = expit(- 0.8 + temp_A + 0.1*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.3*ifelse_vec(t>1, temp_data[[t]]$R, L02)),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )
      temp_Z <- map_dbl(.x = expit(- 0.5 + 0.8*L02 + 0.8*temp_A + temp_R),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )
      if (!if_LY_misspec) {
        # default, correct data
        temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02 + temp_A + 0.7*temp_Z - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = expit(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A - 0.3*temp_Z
                                     # - 0.2*temp_A*temp_Z
                                     - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0)
                                     ),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      } else {
        # LY misspec with squares
        temp_L <- map_dbl(.x = scale_bounded_relu(- 1 + 0.3*L02 + temp_A^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = scale_bounded_relu(0.2 + 1.5*L02^2 + temp_R + 0.4*temp_L^2 - 0.3*temp_A - 0.3*temp_Z - 0.8*temp_A*temp_Z -
                                          0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      }

      temp_data[[t + 1]] <- data.frame(A = temp_A,
                                       R = temp_R,
                                       Z = temp_Z,
                                       L1 = temp_L,
                                       Y = temp_Y
      )
      t <- t + 1
    }
  } else {
    # for mediation targets
    # time point 1~tau
    # t = 1, ..., tau; to update the t-th time point in the t+1 slot
    while(t <= tau) {
      # omit censoring for now
      temp_A <- rep(setAM[1], B)
      temp_R <- map_dbl(.x = expit(- 0.8 + temp_A + 0.1*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.3*ifelse_vec(t>1, temp_data[[t]]$R, L02)),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )
      # replace with mediator level
      temp_Z <- map_dbl(.x = expit(- 0.5 + 0.8*L02 +
                                     # 0.8*temp_A +
                                     0.8*rep(setAM[2], B) +
                                     temp_R),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )
      if (!if_LY_misspec) {
        # default, correct data
        temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02 + temp_A + 0.7*temp_Z
                                     - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)
                                     ),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = expit(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A - 0.3*temp_Z
                                     # - 0.2*temp_A*temp_Z
                                     - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0)
                                       ),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      } else {
        # LY misspec
        temp_L <- map_dbl(.x = scale_bounded_relu(- 1 + 0.3*L02 + temp_A^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = scale_bounded_relu(0.2 + 1.5*L02^2 + temp_R + 0.4*temp_L^2 - 0.3*temp_A - 0.3*temp_Z - 0.8*temp_A*temp_Z -
                                          0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      }

      temp_data[[t + 1]] <- data.frame(A = temp_A,
                                       R = temp_R,
                                       Z = temp_Z,
                                       L1 = temp_L,
                                       Y = temp_Y
      )
      t <- t + 1
    }
  }

  for (s in 1:(tau + 1)) {
    colnames(temp_data[[s]]) <- paste0(colnames(temp_data[[s]]), "_", s-1)
  }

  return(temp_data)
}







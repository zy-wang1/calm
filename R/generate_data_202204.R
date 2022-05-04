library(dplyr)
library(purrr)

# sample_size <- 100  # n
# tau <- 4  # last time point
# seed <- 202008  # random seed

generate_Zheng_data_survival_202204 <- function(sample_size, tau, seed = NULL, setAM = NULL, if_LY_misspec = F, if_A_misspec = F, event_label = 1,
                                                if_competing_risk = F, if_omit_censoring = F,
                                                A_C_scaling = 1, A_E_scalling = 1, Z_scaling = 1, Y_scaling = 1,
                                                A_E_int = 0, Y_int = 0
                                                ) {
    if (is.null(seed)) {} else set.seed(seed)

    # time point 0
    t <- 0
    L01 <- rbinom(n = sample_size, size = 1, prob = 0.4)
    L02 <- rbinom(n = sample_size, size = 1, prob = 0.6)
    temp_data <- list()
    temp_data[[1]] <- data.frame(L1 = L01, L2 = L02)
    t <- t + 1  # note that t 0 value is associated with the first time point; t is actually the t-1 time point in the paper

    if(is.null(setAM)) {
        # time point 1~tau
        # t = 1, ..., tau; to update the t-th time point in the t+1 slot
        while(t <= tau) {
            if (!if_A_misspec) {
                temp_A_C <- map_dbl(.x = expit(A_C_scaling * (1.5 - 0.8*L02 - 0.4*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.5*ifelse_vec(t>1, temp_data[[t]]$A_E, 0))),
                                    .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                if(if_omit_censoring) temp_A_C <- 1
                temp_A_E <- map_dbl(.x = expit(A_E_int + A_E_scalling * (- 0.1 + 1.2*L02 + 0.7*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) - 0.1*ifelse_vec(t>1, temp_data[[t]]$A_E, 0))),
                                    .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
            } else {  # all quadratic
                temp_A_C <- map_dbl(.x = expit(A_C_scaling*(1.5 - 0.8*L02^2 -
                                                   0.4*ifelse_vec(t > 1, temp_data[[t]]$L1, L01)^2 +
                                                   0.5*ifelse_vec(t>1, temp_data[[t]]$A_E, 0)^2)
                                               ),
                                    .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                if(if_omit_censoring) temp_A_C <- 1
                temp_A_E <- map_dbl(.x = expit(A_E_int + A_E_scalling * (- 0.1 + 1.2*L02^2 + 0.7*ifelse_vec(t > 1, temp_data[[t]]$L1, L01)^2 - 0.1*ifelse_vec(t>1, temp_data[[t]]$A_E^2, 0)^2)),
                                    .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
            }

            temp_R <- map_dbl(.x = expit(- 0.8 + temp_A_E + 0.1*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.3*ifelse_vec(t>1, temp_data[[t]]$R, L02)),
                              .f = ~ rbinom(n = 1, size = 1, prob = .x)
            )
            temp_Z <- map_dbl(.x = expit(Z_scaling * (- 0.5 + 0.8*L02 + 0.8*temp_A_E + temp_R)),
                              .f = ~ rbinom(n = 1, size = 1, prob = .x)
            )

            if (!if_LY_misspec) {
                # default, correct data
                temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02 + temp_A_E + 0.7*temp_Z - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)),
                                  .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                temp_Y <- map_dbl(.x = expit(Y_int + Y_scaling * (0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                              # - 0.2*temp_A_E*temp_Z
                                              - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0))/4 ),
                                  .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                # temp_Y <- 1-temp_Y
                if (if_competing_risk) {
                    temp_Y2 <- map_dbl(.x = expit(-0.5 + 0.1*temp_R + 0.1*temp_L - 0.1*temp_A_E - 0.1*temp_Z
                                                  - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0) ),
                                       .f = ~ rbinom(n = 1, size = 1, prob = .x)
                    )
                }
            } else {
                # LY misspec with squares
                # temp_L <- map_dbl(.x = scale_bounded_relu(- 1 + 0.3*L02 + temp_A_E^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0), upper = 4),
                #                   .f = ~ rbinom(n = 1, size = 1, prob = .x)
                # )
                # temp_Y <- map_dbl(.x = scale_bounded_relu(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                #                                           - 0.2*temp_A_E*temp_Z
                #                                           - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0), upper = 4),
                #                   .f = ~ rbinom(n = 1, size = 1, prob = .x)
                # )
                # temp_Y <- 1-temp_Y
                temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02^2 + temp_A_E^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)^2),
                                  .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                temp_Y <- map_dbl(.x = expit(Y_int + Y_scaling * (0.2 + 1.5*L02^2 + temp_R^2 + 0.2*temp_L^2 - 0.3*temp_A_E^2 - 0.3*temp_Z^2
                                              # - 0.2*temp_A_E*temp_Z^2
                                              - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0)^2)/4 ),
                                  .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                # temp_Y <- 1-temp_Y
                if (if_competing_risk) {
                    temp_Y2 <- map_dbl(.x = expit(-0.5 + 0.1*temp_R^2 + 0.1*temp_L^2 - 0.1*temp_A_E^2 - 0.1*temp_Z^2
                                                  - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0)^2 ),
                                       .f = ~ rbinom(n = 1, size = 1, prob = .x)
                    )
                }

            }
            temp_data[[t + 1]] <- data.frame(A_C = temp_A_C,
                                             A_E = temp_A_E,
                                             R = temp_R,
                                             Z = temp_Z,
                                             L1 = temp_L,
                                             Y = temp_Y
            )
            if (if_competing_risk){
                temp_Y2[which(temp_Y == 1)] <- 0  # if dead, cannot have competing event Y2
                temp_data[[t + 1]]$Y2 <- temp_Y2
            }
            t <- t + 1
        }
    } else {
        # for mediation targets
        # time point 1~tau
        # t = 1, ..., tau; to update the t-th time point in the t+1 slot
        while(t <= tau) {
            # omit censoring for now
            temp_A_E <- rep(setAM[1], sample_size)
            temp_A_C <- 1
            temp_R <- map_dbl(.x = expit(- 0.8 + temp_A_E + 0.1*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.3*ifelse_vec(t>1, temp_data[[t]]$R, L02)),
                              .f = ~ rbinom(n = 1, size = 1, prob = .x)
            )
            temp_Z <- map_dbl(.x = expit(Z_scaling * (- 0.5 + 0.8*L02 + 0.8*rep(setAM[2], sample_size) + temp_R)),
                              .f = ~ rbinom(n = 1, size = 1, prob = .x)
            )
            if (!if_LY_misspec) {
                # default, correct data
                temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02 + temp_A_E + 0.7*temp_Z - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)),
                                  .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                temp_Y <- map_dbl(.x = expit(Y_int + Y_scaling * (0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                              # - 0.2*temp_A_E*temp_Z
                                              - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0))/4 ),
                                  .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                # temp_Y <- 1-temp_Y
                if (if_competing_risk) {
                    temp_Y[temp_data[[t]]$Y2 == 1] <- 0  # remove cause specific death if there has been a competing condition
                    temp_Y2 <- map_dbl(.x = expit(-1.5 + 0.1*temp_R + 0.1*temp_L - 0.1*temp_A_E - 0.1*temp_Z
                                                  - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0) ),
                                       .f = ~ rbinom(n = 1, size = 1, prob = .x)
                    )
                }
            } else {
                # LY misspec
                # temp_L <- map_dbl(.x = scale_bounded_relu(- 1 + 0.3*L02 + temp_A_E^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0), upper = 4),
                #                   .f = ~ rbinom(n = 1, size = 1, prob = .x)
                # )
                # temp_Y <- map_dbl(.x = scale_bounded_relu(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                #                                           - 0.2*temp_A_E*temp_Z
                #                                           - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0), upper = 4),
                #                   .f = ~ rbinom(n = 1, size = 1, prob = .x)
                # )
                # temp_Y <- 1-temp_Y
                temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02^2 + temp_A_E^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)^2),
                                  .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                temp_Y <- map_dbl(.x = expit(Y_int + Y_scaling * (0.2 + 1.5*L02^2 + temp_R^2 + 0.2*temp_L^2 - 0.3*temp_A_E^2 - 0.3*temp_Z^2
                                              # - 0.2*temp_A_E*temp_Z^2
                                              - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0)^2)/4 ),
                                  .f = ~ rbinom(n = 1, size = 1, prob = .x)
                )
                # temp_Y <- 1-temp_Y
                if (if_competing_risk) {
                    temp_Y[temp_data[[t]]$Y2 == 1] <- 0  # remove cause specific death if there has been a competing condition
                    temp_Y2 <- map_dbl(.x = expit(-1.5 + 0.1*temp_R^2 + 0.1*temp_L^2 - 0.1*temp_A_E^2 - 0.1*temp_Z^2
                                                  - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0)^2 ),
                                       .f = ~ rbinom(n = 1, size = 1, prob = .x)
                    )
                }
            }

            temp_data[[t + 1]] <- data.frame(A_C = temp_A_C,
                                             A_E = temp_A_E,
                                             R = temp_R,
                                             Z = temp_Z,
                                             L1 = temp_L,
                                             Y = temp_Y
            )
            if (if_competing_risk){
                temp_Y2[which(temp_Y == 1)] <- 0  # if dead, cannot have competing event Y2
                temp_data[[t + 1]]$Y2 <- temp_Y2
            }
            t <- t + 1
        }
    }

    for (t in 1:tau) {
        temp_if_censored <- temp_data[[t+1]]$A_C == 0
        for (s in t:tau) {
            for (l in 2:ncol(temp_data[[s+1]])) temp_data[[s+1]][[l]][temp_if_censored] <- NA
            temp_data[[s+1]][[1]][temp_if_censored] <- ifelse(s == t, 0, NA)  # censored at the time point t; delete after t
        }
    }
    if (tau > 1) {
        for (t in 1:(tau-1)) {
            temp_if_dead <- temp_data[[t+1]]$Y == event_label & (!is.na(temp_data[[t+1]]$Y))
            if (if_competing_risk) temp_if_compete <- temp_data[[t+1]]$Y2 == event_label & (!is.na(temp_data[[t+1]]$Y2))
            for (s in (t+1):tau) for (l in 1:ncol(temp_data[[s+1]])) {
                if (colnames(temp_data[[s+1]])[l] == "Y"  ) {
                    temp_data[[s+1]][[l]][temp_if_dead] <- 1  # remain dead after death
                } else if (colnames(temp_data[[s+1]])[l] == "Y2") {
                    temp_data[[s+1]][[l]][temp_if_dead] <- 0  # remain no competing condition if dead
                }  else {
                    temp_data[[s+1]][[l]][temp_if_dead] <- NA  # dead at t, delete after t (including censoring node)
                }
                # remain no death if there has been a competing condition
                if (if_competing_risk) {
                    if (colnames(temp_data[[s+1]])[l] == "Y"  ) {
                        temp_data[[s+1]][[l]][temp_if_compete] <- 0  # remain no death after compete
                    } else if (colnames(temp_data[[s+1]])[l] == "Y2") {
                        temp_data[[s+1]][[l]][temp_if_compete] <- 1  # remain compete if already
                    }  else {
                        temp_data[[s+1]][[l]][temp_if_compete] <- NA  # compete at t, delete after t (including censoring node)
                    }
                }

            }
        }
    }

    for (s in 1:(tau + 1)) {
        colnames(temp_data[[s]]) <- paste0(colnames(temp_data[[s]]), "_", s-1)
    }

    return(temp_data)
}


# generate_Zheng_data_survival_202112(sample_size = 100, tau = 2, seed = 123, setAM = NULL, if_LY_misspec = F, if_A_misspec = F, event_label = 1) %>% lapply(function(x) head(x, 10))
# generate_Zheng_data_survival_202112(sample_size = 100, tau = 2, seed = 123, setAM = NULL, if_LY_misspec = F, if_A_misspec = F, event_label = 1)
# generate_Zheng_data_survival_202112(sample_size = 100000, tau = 2, seed = 123, setAM = c(1, 0), if_LY_misspec = F, if_A_misspec = F, event_label = 1)[[2]]$Y_1 %>% table(useNA = "always")
# generate_Zheng_data_survival_202112(sample_size = 100000, tau = 2, seed = 123, setAM = c(1, 0), if_LY_misspec = F, if_A_misspec = F, event_label = 1)[[3]]$Y_2 %>% table(useNA = "always")
# generate_Zheng_data_survival_202112(sample_size = 100000, tau = 3, seed = 123, setAM = c(1, 0), if_LY_misspec = F, if_A_misspec = F, event_label = 1)[[4]]$Y_3 %>% table(useNA = "always")



# generate_Zheng_data_survival_202112(sample_size = 100, tau = 2, seed = 123, setAM = NULL, if_LY_misspec = F, if_A_misspec = F, event_label = 1, if_competing_risk = T) %>% lapply(function(x) head(x, 10))
# generate_Zheng_data_survival_202112(sample_size = 100, tau = 2, seed = 123, setAM = NULL, if_LY_misspec = F, if_A_misspec = F, event_label = 1, if_competing_risk = T)
# generate_Zheng_data_survival_202112(sample_size = 100000, tau = 2, seed = 123, setAM = c(1, 0), if_LY_misspec = F, if_A_misspec = F, event_label = 1, if_competing_risk = T)[[3]]$Y_2 %>% table(useNA = "always")
# generate_Zheng_data_survival_202112(sample_size = 100000, tau = 2, seed = 123, setAM = c(1, 0), if_LY_misspec = F, if_A_misspec = F, event_label = 1, if_competing_risk = T)[[3]]$Y2_2 %>% table(useNA = "always")



# test <- generate_Zheng_data_survival_202204(sample_size = 500, tau = 2, seed = 123, setAM = NULL, if_LY_misspec = F, if_A_misspec = F, event_label = 1,
#                                             A_C_scaling = 1.5, A_E_scalling = 0.3, Z_scaling = 2, Y_scaling = 0.05
#                                             )
# test[[2]]$A_C_1 %>% table
# test[[2]]$Y_1 %>% table
# test[[3]]$A_C_2 %>% table
# test[[3]]$Y_2 %>% table
# test[[2]]$A_E_1 %>% table
# test[[3]]$A_E_2 %>% table


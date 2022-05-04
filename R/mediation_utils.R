# shared functions across versions of mediation branches

#' @export
expand_values <- function(variables, to_drop = NULL, values = NULL, rule_variables = NULL, rule_values = NULL, ...) {
  # variables is the vector of variable names
  # drop is either a vector of which variables to drop, or their indices
  # values are the possible values to expand; if null, generate binary values
  # ... are other value rules

  # input_list <- list(A = 1)

  input_list <- list(...)
  rules_list <- lapply(names(input_list), function(eachName) {
    if (length(grep("_", eachName)) > 0) eachName else variables[grep(eachName, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))]
  })
  if(!is.null(rule_variables)) {
    temp_list <- as.list(rule_values)
    names(temp_list) <- rule_variables
    input_list <- c(input_list, temp_list)
    rules_list <- c(rules_list, rule_variables)
  }

  # drop variables that don't want to generate
  # to_drop <- c("A", "L1_0")
  if (is.null(to_drop)) variables_to_generate <- variables else
    if(is.numeric(to_drop)) variables_to_generate <- variables[-to_drop] else
      if (is.character(to_drop)) {
        ind_to_drop <- map(to_drop, ~which(variables == .x)) %>% compact %>% unlist
        if (length(ind_to_drop) == 0) variables_to_generate <- variables else
          variables_to_generate <- variables[-ind_to_drop]
      } else {
        # consider how to identify error later
        print("invalid variables to drop")
        variables_to_generate <- variables
      }

  all_possible_values <- map(variables_to_generate, function(eachVar) {
    test_rules <- which(map_dbl(rules_list, ~length(grep(eachVar, .x))) != 0)
    if (length(test_rules) == 0) {
      if(is.null(values[[eachVar]])) return(0:1) else return(values[[eachVar]])
    } else {
      if (length(test_rules) > 1) {
        if (all(unlist(input_list[test_rules]) == input_list[test_rules][[1]])) {
          return(input_list[test_rules][1] %>% as.numeric)
        } else {
          stop("Contradicting rules for expand_values. ")
        }
      } else {
        return(input_list[test_rules] %>% as.numeric)
      }
    }
  }) %>% expand.grid() %>% data.frame
  colnames(all_possible_values) <- variables_to_generate

  # to add values
  # to add variable name rules

  return(all_possible_values)
}

# values <- expand_values(variables, A = 1, rule_variables = last(variables), rule_values = 1)
# merge(values[1:5, ], values)

#' @export
which_variable <- function(variables, target_variable, timepoint) {
  grep(paste0(target_variable, "_", timepoint), variables)
}

#' @export
# get the orders of variables after droping by name or time
which_variable_drop <- function(variables, to_drop_variable = NULL, to_drop_time = NULL, to_drop = NULL) {
  if (!is.null(to_drop_variable)) {
    temp_drop_variable <- lapply(to_drop_variable, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))) %>% unlist
  } else
    temp_drop_variable <- NULL
  if (!is.null(to_drop_time)) {
    temp_drop_time <- lapply(to_drop_time, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][2]))) %>% unlist
  }  else temp_drop_time <- NULL
  if (!is.null(to_drop)) {
    temp_drop_both <- lapply(to_drop, function(s) grep(s, variables)) %>% unlist
  } else temp_drop_both <- NULL
  temp_drop <- c(temp_drop_variable, temp_drop_time, temp_drop_both)
  if (!is.null(temp_drop)) (1:length(variables))[-temp_drop] %>% return else 1:length(variables)
}

# which_variable_drop(variables, "L")

#' @export
# take variables by name or time
which_variable_take <- function(variables, to_take_variable = NULL, to_take_time = NULL, logic = "or") {
  if (!is.null(to_take_variable)) {
    temp_take_variable <- lapply(to_take_variable, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))) %>% unlist
  } else
    temp_take_variable <- NULL
  if (!is.null(to_take_time)) {
    temp_take_time <- lapply(to_take_time, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][2]))) %>% unlist
  }  else temp_take_time <- NULL
  # default is to take all selected, by name or by time
  # and is to take the intersection
  if (logic == "or") temp_take <- c(temp_take_variable, temp_take_time) else
    if (logic == "and") temp_take <- intersect(temp_take_variable, temp_take_time)
    if (!is.null(temp_take)) sort(unique(temp_take)) else NULL
}





#' @export
ifelse_vec <- function(condition, out1, out2) {
  if (condition) return(out1) else return(out2)
}

#' @export
logit <- function(x) {
  log(x/(1-x))
}

#' @export
expit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' @export
scale_01 <- function(x) scale(x, center = min(x), scale = max(x) - min(x))

#' @export
scale_bounded_relu <- function(x, move = 1, upper = 0.95, lower = 0.1) {
  x <- x + move
  x[x >= upper] <- upper
  if (upper > 1) x <- x/upper
  x[x <= lower] <- lower
  return(x)
}





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






generate_Zheng_data_survival <- function(sample_size, tau, seed = NULL, setAM = NULL, if_LY_misspec = F, if_A_misspec = F, event_label = 0) {
  if (is.null(seed)) {} else set.seed(seed)

  # time point 0
  t <- 0
  L01 <- rbinom(n = sample_size, size = 1, prob = 0.4)
  L02 <- rbinom(n = sample_size, size = 1, prob = 0.6)
  temp_data <- list()
  temp_data[[1]] <- data.frame(L1 = L01, L2 = L02)
  t <- t + 1

  if(is.null(setAM)) {
    # time point 1~tau
    # t = 1, ..., tau; to update the t-th time point in the t+1 slot
    while(t <= tau) {
      # omit censoring for now
      if (!if_A_misspec) {
        temp_A_C <- map_dbl(.x = expit(1.5 - 0.8*L02 - 0.4*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.5*ifelse_vec(t>1, temp_data[[t]]$A_E, 0)),
                            .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_A_E <- map_dbl(.x = expit(- 0.1 + 1.2*L02 + 0.7*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) - 0.1*ifelse_vec(t>1, temp_data[[t]]$A_E, 0)),
                            .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      } else {
        temp_A_C <- map_dbl(.x = expit(1.5 - 0.8*L02 - 0.4*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.5*ifelse_vec(t>1, temp_data[[t]]$A_E, 0)),
                            .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_A_E <- map_dbl(.x = expit(- 0.1 + 1.2*L02^2 + 0.7*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) - 0.1*ifelse_vec(t>1, temp_data[[t]]$A_E^2, 0)),
                            .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      }

      temp_R <- map_dbl(.x = expit(- 0.8 + temp_A_E + 0.1*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.3*ifelse_vec(t>1, temp_data[[t]]$R, L02)),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )
      temp_Z <- map_dbl(.x = expit(- 0.5 + 0.8*L02 + 0.8*temp_A_E + temp_R),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )

      if (!if_LY_misspec) {
        # default, correct data
        temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02 + temp_A_E + 0.7*temp_Z - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = expit(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                     - 0.2*temp_A_E*temp_Z
                                     - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0) ),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- 1-temp_Y
      } else {
        # LY misspec with squares
        temp_L <- map_dbl(.x = scale_bounded_relu(- 1 + 0.3*L02 + temp_A_E^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = scale_bounded_relu(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                                  - 0.2*temp_A_E*temp_Z
                                                  - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- 1-temp_Y
      }
      temp_data[[t + 1]] <- data.frame(A_C = temp_A_C,
                                       A_E = temp_A_E,
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
      temp_A_E <- rep(setAM[1], sample_size)
      temp_A_C <- 1
      temp_R <- map_dbl(.x = expit(- 0.8 + temp_A_E + 0.1*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.3*ifelse_vec(t>1, temp_data[[t]]$R, L02)),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )
      temp_Z <- map_dbl(.x = expit(- 0.5 + 0.8*L02 + 0.8*rep(setAM[2], sample_size) + temp_R),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )
      if (!if_LY_misspec) {
        # default, correct data
        temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02 + temp_A_E + 0.7*temp_Z - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = expit(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                     - 0.2*temp_A_E*temp_Z
                                     - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0) ),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- 1-temp_Y
      } else {
        # LY misspec
        temp_L <- map_dbl(.x = scale_bounded_relu(- 1 + 0.3*L02 + temp_A_E^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = scale_bounded_relu(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                                  - 0.2*temp_A_E*temp_Z
                                                  - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- 1-temp_Y
      }

      temp_data[[t + 1]] <- data.frame(A_C = temp_A_C,
                                       A_E = temp_A_E,
                                       R = temp_R,
                                       Z = temp_Z,
                                       L1 = temp_L,
                                       Y = temp_Y
      )
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
      for (s in (t+1):tau) for (l in 1:ncol(temp_data[[s+1]])) temp_data[[s+1]][[l]][temp_if_dead] <- NA  # dead at t, delete after t (including censoring node)
    }
  }

  for (s in 1:(tau + 1)) {
    colnames(temp_data[[s]]) <- paste0(colnames(temp_data[[s]]), "_", s-1)
  }

  return(temp_data)
}






generate_Zheng_data_survival_mn <- function(sample_size, tau, seed = NULL, setAM = NULL, if_LY_misspec = F, if_A_misspec = F, event_label = 0) {
  if (is.null(seed)) {} else set.seed(seed)

  # time point 0
  t <- 0
  L01 <- rbinom(n = sample_size, size = 1, prob = 0.4)
  L02 <- rbinom(n = sample_size, size = 1, prob = 0.6)
  temp_data <- list()
  temp_data[[1]] <- data.frame(L1 = L01, L2 = L02)
  t <- t + 1

  if(is.null(setAM)) {
    # time point 1~tau
    # t = 1, ..., tau; to update the t-th time point in the t+1 slot
    while(t <= tau) {
      # omit censoring for now
      if (!if_A_misspec) {
        temp_A_C <- map_dbl(.x = expit(1.5 - 0.8*L02 - 0.4*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.5*ifelse_vec(t>1, temp_data[[t]]$A_E, 0)),
                            .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_A_E <- map_dbl(.x = expit(- 0.1 + 1.2*L02 + 0.7*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) - 0.1*ifelse_vec(t>1, temp_data[[t]]$A_E, 0)),
                            .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      } else {
        temp_A_C <- map_dbl(.x = expit(1.5 - 0.8*L02 - 0.4*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.5*ifelse_vec(t>1, temp_data[[t]]$A_E, 0)),
                            .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_A_E <- map_dbl(.x = expit(- 0.1 + 1.2*L02^2 + 0.7*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) - 0.1*ifelse_vec(t>1, temp_data[[t]]$A_E^2, 0)),
                            .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
      }

      temp_R <- map_dbl(.x = expit(- 0.8 + temp_A_E + 0.1*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.3*ifelse_vec(t>1, temp_data[[t]]$R, L02)),
                        .f = ~ rbinom(n = 3, size = 1, prob = .x) %>% sum
      )
      temp_Z <- map_dbl(.x = expit(- 0.5 + 0.8*L02 + 0.8*temp_A_E + temp_R),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )

      if (!if_LY_misspec) {
        # default, correct data
        temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02 + temp_A_E + 0.7*temp_Z - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = expit(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                     - 0.2*temp_A_E*temp_Z
                                     - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0) ),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- 1-temp_Y
      } else {
        # LY misspec with squares
        temp_L <- map_dbl(.x = scale_bounded_relu(- 1 + 0.3*L02 + temp_A_E^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = scale_bounded_relu(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                                  - 0.2*temp_A_E*temp_Z
                                                  - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- 1-temp_Y
      }
      temp_data[[t + 1]] <- data.frame(A_C = temp_A_C,
                                       A_E = temp_A_E,
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
      temp_A_E <- rep(setAM[1], sample_size)
      temp_A_C <- 1
      temp_R <- map_dbl(.x = expit(- 0.8 + temp_A_E + 0.1*ifelse_vec(t > 1, temp_data[[t]]$L1, L01) + 0.3*ifelse_vec(t>1, temp_data[[t]]$R, L02)),
                        .f = ~ rbinom(n = 3, size = 1, prob = .x) %>% sum
      )
      temp_Z <- map_dbl(.x = expit(- 0.5 + 0.8*L02 + 0.8*rep(setAM[2], sample_size) + temp_R),
                        .f = ~ rbinom(n = 1, size = 1, prob = .x)
      )
      if (!if_LY_misspec) {
        # default, correct data
        temp_L <- map_dbl(.x = expit(- 1 + 0.3*L02 + temp_A_E + 0.7*temp_Z - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0)),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = expit(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                     - 0.2*temp_A_E*temp_Z
                                     - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0) ),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- 1-temp_Y
      } else {
        # LY misspec
        temp_L <- map_dbl(.x = scale_bounded_relu(- 1 + 0.3*L02 + temp_A_E^2 + 0.7*temp_Z^2 - 0.2*ifelse_vec(t>1, temp_data[[t]]$L1, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- map_dbl(.x = scale_bounded_relu(0.2 + 1.5*L02 + temp_R + 0.2*temp_L - 0.3*temp_A_E - 0.3*temp_Z
                                                  - 0.2*temp_A_E*temp_Z
                                                  - 0.1*ifelse_vec(t>1, temp_data[[t]]$R, 0), upper = 4),
                          .f = ~ rbinom(n = 1, size = 1, prob = .x)
        )
        temp_Y <- 1-temp_Y
      }

      temp_data[[t + 1]] <- data.frame(A_C = temp_A_C,
                                       A_E = temp_A_E,
                                       R = temp_R,
                                       Z = temp_Z,
                                       L1 = temp_L,
                                       Y = temp_Y
      )
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
      for (s in (t+1):tau) for (l in 1:ncol(temp_data[[s+1]])) temp_data[[s+1]][[l]][temp_if_dead] <- NA  # dead at t, delete after t (including censoring node)
    }
  }

  for (s in 1:(tau + 1)) {
    colnames(temp_data[[s]]) <- paste0(colnames(temp_data[[s]]), "_", s-1)
  }

  return(temp_data)
}




# independent random uniform, not exactly expand.grid over all possible levels
expand_values_size <- function(variables, to_drop = NULL, values = NULL, rule_variables = NULL, rule_values = NULL, size = 5000, ...) {
  # variables is the vector of variable names
  # drop is either a vector of which variables to drop, or their indices
  # values are the possible values to expand; if null, generate binary values
  # ... are other value rules

  # input_list <- list(A = 1)

  input_list <- list(...)
  rules_list <- lapply(names(input_list), function(eachName) {
    if (length(grep("_", eachName)) > 0) eachName else variables[grep(eachName, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))]
  })
  if(!is.null(rule_variables)) {
    temp_list <- as.list(rule_values)
    names(temp_list) <- rule_variables
    input_list <- c(input_list, temp_list)
    rules_list <- c(rules_list, rule_variables)
  }

  # drop variables that don't want to generate
  # to_drop <- c("A", "L1_0")
  if (is.null(to_drop)) variables_to_generate <- variables else
    if(is.numeric(to_drop)) variables_to_generate <- variables[-to_drop] else
      if (is.character(to_drop)) {
        ind_to_drop <- map(to_drop, ~which(variables == .x)) %>% compact %>% unlist
        if (length(ind_to_drop) == 0) variables_to_generate <- variables else
          variables_to_generate <- variables[-ind_to_drop]
      } else {
        # consider how to identify error later
        print("invalid variables to drop")
        variables_to_generate <- variables
      }

  prep_list <- map(variables_to_generate, function(eachVar) {
    test_rules <- which(map_dbl(rules_list, ~length(grep(eachVar, .x))) != 0)
    if (length(test_rules) == 0) {
      if(is.null(values[[eachVar]])) return(0:1) else return(values[[eachVar]])
    } else {
      if (length(test_rules) > 1) {
        if (all(unlist(input_list[test_rules]) == input_list[test_rules][[1]])) {
          return(input_list[test_rules][1] %>% as.numeric)
        } else {
          stop("Contradicting rules for expand_values. ")
        }
      } else {
        return(input_list[test_rules] %>% as.numeric)
      }
    }
  })

  names(prep_list) <- variables_to_generate

  all_possible_values <- prep_list %>% sapply(function(u) sample(u, size, replace = T)) %>% data.frame
  colnames(all_possible_values) <- variables_to_generate


  # to add values
  # to add variable name rules

  return(list(all_possible_values, prep_list))
}

#' @export
ipw_gen <- function(Y, A, g, W){
  Y*A/g  - Y*(1-A)/g
}

#' @export
generator_ate <-function(tmle_task, lik = NULL, target_param = NULL, node, outcome = T){
  task <- tmle_task$get_regression_task(node)
  A <- task$X$A
  Y <- task$Y
  W <- task$X$W

  g <- lik$get_likelihood(tmle_task, "A")

  IC <- ipw_gen(Y,A,g, W)

  cols <- task$add_columns(data.table(IC = IC))
  task <- task$clone()
  nodes <- task$nodes
  nodes$outcome <- "IC"
  nodes$covariates <- c(nodes$covariates, node)
  task$initialize(
    task$internal_data,
    nodes = nodes,
    folds = task$folds,
    column_names = cols,
    row_index = task$row_index,
    outcome_type = "continuous"
  )
  return(task)
  # task$next_in_chain(column_names = cols, covariates  = c(task$nodes$covariates, task$nodes$outcome), outcome = "IC")
}


#' Class representing a (possibly non-canonical) gradient for some parameter.
#' Currently, this gradient object takes a IPW/unobserved model gradient
#' and numerically projects it onto the tangent space.
#'
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom digest digest
#' @import data.table
#' @import sl3
#' @param likelihood A trained likelihood object from which to compute the gradient from.
#' @param projection_task_generator A function that takes a task, likelihood, target parameter, and target_node
#' and returns a task with the outcome being the evaluated gradient,
#' and covariates being whatever variables the conditional density of the target node likelihood factor depends on.
#' (E.g. covariates + outcome in the regression task)
#' @export
#'
#' @keywords data
#'
#' @return \code{Gradient} object
#'
#' @format \code{\link{R6Class}} object.
#'
#'
#' @export
#
Gradient_wide <- R6Class(
  classname = "Gradient",
  portable = TRUE,
  class = TRUE,
  inherit = Lrnr_base,
  public = list(
    initialize = function(likelihood, ipw_args = NULL, projection_task_generator, target_nodes = "Y"){
      params <- sl3::args_to_list()
      params$target_nodes <- target_nodes
      private$.params <- params
      private$.cache <- new.env()
      private$.learner <- Lrnr_hal9001a$new(max_degree = 2, family = "gaussian")
    },
    train_projections = function(tmle_task, fold_number = private$.fold_number){
      private$.fold_number <- fold_number
      self$train(tmle_task)
    },
    generate_task = function(tmle_task, node, include_outcome = T, fold_number = "full"){
      self$projection_task_generator(tmle_task, self$likelihood, node, include_outcome = include_outcome, self$params$ipw_args, fold_number = fold_number)
    },
    expand_task = function(tmle_task, node, force = F){
      #Computes expanded task where observations are repeated (with fake ids) for all levels of node
      # TODO these expanded tasks should only be targetted for this node.


      if(tmle_task$uuid %in% names(private$.uuid_expanded_history)){


        if(private$.uuid_expanded_history[[tmle_task$uuid]] != node){
          stop("This expanded task does not match its node. You shouldn't be targeting this node for this task ")
        }
        return(tmle_task)
      }

      key <- paste0(tmle_task$uuid, node, sep = "%")

      cached_task <- get0(key, self$cache, inherits = FALSE)
      if(!is.null(cached_task)){
        if(is.null(attr(cached_task, "target_nodes"))) {
          print(node)
          stop("wrong")
        }
        return(cached_task)
      }
      variables <- tmle_task$npsem[[node]]$variables

      if(length(variables) >1) stop("Multivariate nodes not supported")
      data <- tmle_task$data
      data$trueid <- data$id
      time <- tmle_task$npsem[[node]]$time
      levels <- sort(unique(unlist(tmle_task$get_tmle_node(node))))  #data[, variables, with = F])))

      if(!force & length(levels) > 100){
        stop("Too many levels in node.")
      }
      long_data <- rbindlist(lapply(levels, function(level) {
        data <- copy(data)
        set(data , which(data$t == time), variables, level)
        data$levelcopy <- level
        return(data)
      }))

      long_data$id <-  paste(long_data$trueid,long_data$levelcopy, sep = "_")

      # ZW: only keep full entries till the expanded node
      long_data <- na.omit(long_data, cols = 1:which(colnames(long_data) == tmle_task$npsem[[node]]$variables))

      suppressWarnings(long_task <- tmle3_Task$new(long_data, tmle_task$npsem, id = "id", time = "t", force_at_risk = tmle_task$force_at_risk, summary_measure_columns = c(tmle_task$summary_measure_columns, "trueid")))

      setattr(long_task, "target_nodes", node)
      if(is.null(attr(long_task, "target_nodes"))) {
        print(node)
        stop("wrong")
      }
      assign(key, long_task, self$cache)

      private$.uuid_expanded_history[[long_task$uuid]] <- node
      return(long_task)
    },
    compute_component = function(tmle_task, node, fold_number = "full"){
      time <- tmle_task$npsem[[node]]$time  # not used for wide format

      self$assert_trained()
      #Converts squashed basis to R functions of tmle3_tasks

      long_task <- self$expand_task(tmle_task, node)

      IC_task <- self$generate_task(tmle_task, node, include_outcome = F, fold_number = fold_number)

      if (long_task$npsem[[node]]$variable_type$glm_family() %in% c("binomial")) {
        col_index <- which(colnames(IC_task$X) == node )
      } else {
          stop("continuous and categorical time-varying nodes: in development")
        # col_index <- grep(paste0("^", node, ".X"), colnames(IC_task$X))  # for categorical
      }

      long_preds <- NULL

      tryCatch({
        value_step1 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)
        value_step2 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], long_task, fold_number, node = node)
        if(is.null(value_step1)){
          value_step1 <- 0
        }
        if(is.null(value_step2)){
          value_step2 <- 0
        }
        if(value_step1!=value_step2) {
          # if (value_step1 - value_step2 > 1)
            # step 2 original can compute EIC from step 1 expanded; expanded task will be updated in its own turn;
            # expanded tasks will appear later in cache
            stop("Long_task and tmle_task out of sync.")
        }
        long_preds <- self$likelihood$get_likelihood(long_task, node, fold_number = fold_number, drop_id = T, drop_time = T, drop = T  )},
      error = function(e){
        #long_task is probably out of sync with tmle_task
        #Update it to the same level
        value_step <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)

        self$likelihood$sync_task(long_task, fold_number = fold_number, check = F, max_step = value_step)
        value_step1 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)
        value_step2 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], long_task, fold_number, node = node)
        if(is.null(value_step1)){
          value_step1 <- 0
        }
        if(is.null(value_step2)){
          value_step2 <- 0
        }
        if(value_step1!=value_step2) {
          stop("Long_task and tmle_task out of sync.")
        }
        long_preds <<- self$likelihood$get_likelihood(long_task, node, fold_number = fold_number, drop_id = T, drop_time = T, drop = T, check_sync = F  )

      })
      # long_preds should be the conditional probability value at the currently expanded node variable values (the X in (X, Pa(X)))

      data <- IC_task$data

      # TODO

      variables <- node
      #TODO check id order

      data <- merge(long_task$data[,c("id", "trueid", "t")], 
                                     cbind(long_task$get_tmle_node(node, format = T, include_id = T, include_time = T), long_preds), 
                                     by = c("id", "t")
                                     )

      idkey <- data$trueid
      data[, t := NULL]
      data[, id := NULL]
      
      setnames(data, c("id", node, "pred"))

      setkeyv(data, cols = c("id", node))

      #TODO handle updating of expanded task
      #This is done by stacking copies of the cdf

      data <- dcast(data, as.formula(paste0("id ~ ", node)), value.var = "pred")
      id <- data$id
      data[, id := NULL]
      levels <- as.numeric(colnames(data))

      cdf <- as.data.table(t(apply(data, 1, cumsum)))
      setnames(cdf, as.character(levels))
      #print(cdf)

      if(long_task$uuid == tmle_task$uuid){
        #if expanded task is tmle_task then obtain then expand cdf to match
        #This ensures we dont have any recursion errors by expanding an expanded task

        match_index <- match(idkey, id)
        cdf <- cdf[match_index]
      }

      fit_obj <- private$.component_fits[[node]]
      basis_list <- fit_obj$basis_list  # as of 0.4.3 hal9001, I{X >= u} for knot points u (order 0)
      coefs <- fit_obj$coefs

      # focus on binary first; 
      col_index <- which(colnames(IC_task$X) == node )
      {
          # if y is not involved, then centering basis equals 0
        keep <- sapply(basis_list, function(b){
          col_index %in% b$cols
          # any(col_index %in% b$cols)
        }) %>% unlist

        basis_list <- basis_list[keep]
        coefs <- coefs[c(T, keep)]

        #Should already be sorted
        X <- as.matrix(IC_task$X)
        design <- as.data.table(as.matrix(hal9001::make_design_matrix(X, basis_list)))  # this is phi basis value
        
        # diff_map <- sapply(seq_along(basis_list), function(i) {
        #   basis <- basis_list[[i]]
        #   result <- (list(which(levels == basis$cutoffs[which(basis$cols == col_index)])))
        #   # result <- (list(which(levels == basis$cutoffs[which(basis$cols %in% col_index)])))
        # 
        #   return(result)
        # })  # which level x is the cutoff for this indicator basis; for binary variable and I{Y >= u_Y} basis, cutoff is always 1, i.e. the second level in levels
        diff_map <- lapply(seq_along(basis_list), function(x) return(2))

        # I{X >= u_X}I{Y >= u_Y} - P(Y >= u_Y | X); I{X >= u_X} is redundant but okay with clean_design of I{X >= u_X} again
        center_basis <- lapply(seq_along(diff_map), function(i){
          temp_index <- diff_map[[i]]
          diff <- design[[as.integer(i)]] - 1 + cdf[[temp_index - 1]]  # cdf is the value of <= level; we need >= level probability
          set(design, , as.integer(i), diff)
        })
      }

      min_val <- min(IC_task$X[[node]]) - 5
      clean_basis <- function(basis){
        # index = which(basis$cols == col_index)
        index = which(basis$cols %in% col_index)
        basis$cutoffs[index] <- min_val
        return(basis)
      }
      clean_list = lapply(basis_list, clean_basis)  # for involved phi, set cutoff at y as -infty; this is the Pa(X) component in the basis
      clean_design <- as.data.table(as.matrix(hal9001::make_design_matrix(X, clean_list)))

      mid_result <- as.matrix(design * clean_design)  # design is the centered basis (with a redundant Pa(X) component); clean_design is the Pa(X) component
      result =  mid_result %*% coefs[-1]
      out = list(
          # col_index = col_index,Y = IC_task$Y, cdf = cdf,design = design,  mid_result = mid_result, coefs = coefs[-1], 
          EIC = result)
      return(out)

    },
    
    base_train = function(task, pretrain) {
      fit_object <- private$.train(task, pretrain)
      #new_object <- self$clone() # copy parameters, and whatever else
      self$set_train(fit_object, task)
      private$.training_task <- task

      return(self)
    }
  ),
  active = list(
    params = function(){
      private$.params
    },
    likelihood = function(){
      private$.params$likelihood
    },
    target_param = function(){
      private$.params$target_param
    },
    target_nodes = function(){
      private$.params$target_nodes
    },
    projection_task_generator = function(){
      private$.params$projection_task_generator
    },
    learner = function(){
      private$.learner
    },
    basis = function(){
      private$.basis
    },
    training_task = function(){
      private$.training_task
    },
    cache = function(){
      private$.cache
    },
    fold_number = function(){
      private$.fold_number
    },
    component_fits = function(){
      private$.component_fits
    },
    hal_args = function(args_to_add = NULL){
      if(!is.null(args_to_add)){
        for(name in names(args_to_add)){
          arg_value <- args_to_add[[name]]
          private$.learner_args[name] <- arg_value
        }
        args <- private$.learner_args
        #args$learner_class <- Lrnr_hal9001a
        private$.learner <- do.call(Lrnr_hal9001a$new, args)
      }
      return(private$.learner_args)
    }
  ),
  private = list(
    .train_sublearners = function(tmle_task){
      nodes <- c(self$target_nodes)

      projected_fits <- lapply(nodes, function(node){
        task <- self$generate_task(tmle_task, node, fold_number = self$fold_number)
        lrnr <- self$learner$clone()

        return(delayed_learner_train(lrnr, task))
      })

      projected_fits <- bundle_delayed(projected_fits)
      return(projected_fits)
    },
    .train = function(tmle_task, projected_fits){
      #Store hal_fits
      names(projected_fits) <- self$target_nodes

      component_fits <- lapply(projected_fits, `[[`, "fit_object")
      private$.fit_object <- projected_fits

      component_fits <- lapply(component_fits, function(fit){

        basis_list <- fit$basis_list[as.numeric(names(fit$copy_map))]
        coefs <- fit$coefs

        keep <- coefs[-1]!=0

        basis_list <- basis_list[keep]
        coefs_new <- coefs[c(T, keep)]
        if(sum(coefs_new) != sum(coefs)) stop("squash went wrong")
        if(length(coefs_new) != length(basis_list)+1) stop("squash went wrong")
        return(list(basis_list = basis_list, coefs = coefs_new))
      })
      names(component_fits) <- self$target_nodes

      private$.component_fits <- component_fits
      return(private$.fit_object)

    },
    .predict = function(tmle_task) {
      stop("This gradient has nothing to predict.")
    },
    .chain = function(tmle_task) {
      stop("This gradient has nothing to chain")
    },
    .params = NULL,
    .learner = NULL,
    .learner_args = list(max_degree = 3, family = "gaussian"),
    .component_fits = list(),
    .basis = NULL,
    .training_task = NULL,
    .cache = NULL,
    .uuid_expanded_history = list(),
    .fold_number = "full"
  )
)




#' @docType class
#' @importFrom R6 R6Class
#' @importFrom digest digest
#' @import hal9001
#' @import data.table
#' @import sl3
#' @export
Lrnr_hal9001a <- R6::R6Class(
  classname = "Lrnr_hal9001a", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 2,
                          fit_type = "glmnet",
                          n_folds = 10,
                          use_min = TRUE,
                          reduce_basis = NULL,
                          return_lasso = TRUE,
                          return_x_basis = FALSE,
                          basis_list = NULL,
                          cv_select = TRUE,
                          smoothness_orders = 0,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {
      args <- self$params

      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family <- args$family <- outcome_type$glm_family()
      }

      args$X <- as.matrix(task$X)
      args$Y <- outcome_type$format(task$Y)
      args$yolo <- FALSE

      if (task$has_node("weights")) {
        args$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }
      args$standardize <- F
      fit_object <- sl3:::call_with_args(hal9001::fit_hal, args)
      # ZW: try using cv.glmnet directly
      args$lambda <- fit_object$lambda_star
      args$lambda <- args$lambda * 0.15
      args$fit_control <- list(cv_select = F, n_folds = 10, foldid = NULL, use_min = TRUE,
                               lambda.min.ratio = 1e-04, prediction_bounds = "default")
      # args$cv_select <- F  # for hal9001 version update
      if(args$fit_type == "glmnet"){
        fit_object <- sl3:::call_with_args(hal9001::fit_hal, args)
      }

      return(fit_object)
    },
    .predict = function(task = NULL) {
      predictions <- predict(self$fit_object, new_data = as.matrix(task$X))
      return(predictions)
    },
    .required_packages = c("hal9001")
  )
)





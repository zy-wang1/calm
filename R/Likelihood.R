#' Class for Likelihood
#'
#' This object represents an estimate of the relevant factors of the likelihood estimated from data, or based on \emph{a priori} knowledge where appropriate.
#' That is, it represents some subset of $P_n$. This object inherits from \code{\link[sl3]{Lrnr_base}}, and so shares some properties with \code{sl3} learners.
#' Specifically, to fit a likelihood object to data, one calls \code{likelihood$train(tmle3_task)}.
#' Each likelihood factor is represented by an object inheriting from \code{\link{LF_base}}.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base
#' @importFrom assertthat assert_that is.count is.flag
#' @importFrom delayed bundle_delayed
#' @import data.table
#' @family Likelihood objects
#' @export
#'
#' @keywords data
#'
#' @return \code{Likelihood} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template Likelihood_extra
#'
#' @export
Likelihood <- R6Class(
  classname = "Likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Lrnr_base,
  public = list(
    initialize = function(factor_list, cache = NULL, ...) {
      params <- args_to_list()
      if (inherits(factor_list, "LF_base")) {
        factor_list <- list(factor_list)
      }

      factor_names <- sapply(factor_list, `[[`, "name")
      factor_list_unpooled <- list()
      for(i in seq_along(factor_names)){
        names <- c(factor_names[[i]])
        factor <- factor_list[[i]]
        for(name in names){
          factor_list_unpooled[[name]] <- factor
        }
      }



      names(factor_list) <-  sapply(factor_names, function(name) paste(name, collapse = "%"))

      params$factor_list <- factor_list_unpooled
      params$factor_list_pooled <- factor_list

      if (is.null(cache)) {
        cache <- Likelihood_cache$new()
      }
      private$.cache <- cache

      super$initialize(params)
    },
    print = function() {
      lapply(self$factor_list, print)
      invisible(NULL)
    },
    validate_task = function(tmle_task) {
      assert_that(is(tmle_task, "tmle3_Task"))

      factor_list <- self$factor_list
      factor_names <- names(factor_list)
      task_nodes <- names(tmle_task$npsem)
      if (!all(factor_names %in% task_nodes)) {
        stop("factor_list and task$npsem must have matching names")
      }
    },
    get_likelihood = function(tmle_task, node, fold_number = "full",  drop_id = T, drop_time = T, drop = T, to_wide = F, expand = T) {
      # TODO from regression task extract risk_set and handle this if degeneracy option is set.


      # Likelihood factors just compute likelihoods
      # Risk sets are computed here?
      if(length(node) > 1){
        likelihood_factor <- self$factor_list_pooled[[paste0(node, collapse = "%")]]
        warning("You shouldn't be calling node likelihoods in a pooled way as these likelihoods don't get updated in targeting.")

      } else {
        likelihood_factor <- self$factor_list[[node]]

      }
      # first check for cached values for this task
      likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task, fold_number, node = paste0(node, collapse = "%"))
      if(!is.null(likelihood_values) & length(likelihood_values)==0){
        return(likelihood_values)
      }

      if(!expand & !is.null(likelihood_values)) {
        #Only store the full likelihood
        #Regression task should be cached so this is cheap
        keep <- tmle_task$get_regression_task(node, expand = T)$get_data()$at_risk == 1
        likelihood_values <- likelihood_values[keep]

      }
      if(is.null(likelihood_values)) {
        # if not, generate new ones
        #Include id's and time
        args <- formalArgs(likelihood_factor$get_likelihood)
        if(!all(c("expand", "node") %in% args)){
          likelihood_values <- likelihood_factor$get_likelihood(tmle_task, fold_number)
        } else {
          likelihood_values <- likelihood_factor$get_likelihood(tmle_task, fold_number, expand = expand, node = node)

        }

        # if(!is.data.table(likelihood_values)){
        #   likelihood_values <- as.data.table(likelihood_values)
        #   if(!to_wide){
        #     out_name <- paste0(node, collapse = "%")
        #     setnames(likelihood_values, out_name)
        #   }
        #
        # }
        if(expand){
          self$cache$set_values(likelihood_factor, tmle_task, 0, fold_number, likelihood_values, node = paste0(node, collapse = "%"))
        }

      }

      if(length(likelihood_values)==0){
        return(likelihood_values)
      }
     # names_of <- colnames(likelihood_values)
      #keep_cols <- intersect(c("t", "id", grep(node, names_of, value = T)), names_of)


     # likelihood_values <- likelihood_values[,  keep_cols, with = F]

      # if(to_wide & "t" %in% colnames(likelihood_values)  & length(unique(likelihood_values$t))==1){
      #
      #   likelihood_values$t <- NULL
      # }
      # else if(to_wide & "t" %in% colnames(likelihood_values) & "id" %in% colnames(likelihood_values)){
      #   #likelihood_values <- reshape(likelihood_values, idvar = "id", timevar = "t", direction = "wide")
      #   likelihood_values <- dcast(likelihood_values, id ~ t, value.var = setdiff(names(likelihood_values), c("t", "id")))
      #   if(length(node) + 1 == ncol(likelihood_values)){
      #     setnames(likelihood_values, c("id", node))
      #   } else if (length(node)==1){
      #     setnames(likelihood_values, c("id", paste0( node, "_", names(likelihood_values)[-1])))
      #   }
      # }
      #if(drop_id & "id" %in% colnames(likelihood_values)) likelihood_values$id <- NULL
      #if(drop_time & "t" %in% colnames(likelihood_values)) likelihood_values$t <- NULL
      #if(drop & ncol(likelihood_values) == 1) likelihood_values <- likelihood_values[[1]]
      return(likelihood_values)
    },
    get_likelihoods = function(tmle_task, nodes = NULL, fold_number = "full", drop_id = T, drop_time = T, drop = T, expand = T) {
      if (is.null(nodes)) {
        nodes <- self$nodes
      }

      if (length(nodes) > 1) {

        all_likelihoods <- lapply(nodes, function(node) {
          self$get_likelihood(tmle_task, node = node, fold_number = fold_number, drop_id = F, drop_time = F, to_wide = to_wide, expand =  expand, drop=F)
        })
        likelihood_dt <- as.data.table(all_likelihoods)

        # keep <- unlist(lapply(all_likelihoods, function(lik) length(lik) != 0))
        # all_likelihoods <- all_likelihoods[keep]
        # contains_t <- all(unlist(lapply(all_likelihoods, function(lik){
        #   "t" %in% colnames(lik)
        # })))
        # if(contains_t){
        #   likelihood_dt <- all_likelihoods %>% purrr::reduce(merge, c("id", "t"))#as.data.table(all_likelihoods)
        # } else{
        #
        #   all_likelihoods <- lapply(all_likelihoods, function(lik){
        #     if(!(all(c("id", "t") %in% colnames(lik)))){
        #       if("t" %in% colnames(lik)) lik$t <- NULL
        #       return(lik)
        #     }
        #     if(length(unique(lik$t))==1){
        #       if("t" %in% colnames(lik)) lik$t <- NULL
        #       return(lik)
        #     }
        #     if(all(c("id", "t")) %in% colnames(lik)){
        #       return(reshape(lik, idvar = "id", timevar = "t", direction = "wide"))
        #     }
        #     lik
        #   })
        #
        #   likelihood_dt <- NULL
        #   tryCatch(
        #     {
        #       likelihood_dt <- all_likelihoods %>% purrr::reduce(merge, by = "id")
        #     }, error = function(cond) {
        #       # Handle case when one of predictions doesn't include ID.
        #       # This is mainly for backwards compatibility pre-ltmle changes
        #
        #       likelihood_dt <- setDT(unlist(all_likelihoods, recursive = F))[]
        #       # Need to change likelihood_dt object in parent frame
        #       likelihood_dt <<- likelihood_dt[, which(!duplicated(names(likelihood_dt))), with = F]
        #
        #     }
        #   )
        #
        #   #as.data.table(all_likelihoods)
        #   if("t" %in% colnames(likelihood_dt)) likelihood_dt$t <- NULL
        # }
        # TODO Already has names
        #setnames(likelihood_dt, nodes)
        #if(drop_id) likelihood_dt$id <- NULL
        #if(drop_time & "t" %in% colnames(likelihood_dt)) likelihood_dt$t <- NULL
        return(likelihood_dt)
      } else {
        return(self$get_likelihood(tmle_task, nodes[[1]], fold_number, drop_id = drop_id, drop_time =drop_time, drop = drop, expand = expand, to_wide = to_wide ))
      }
    },
    get_possible_counterfactuals = function(nodes = NULL) {

      # get factors for nodes
      factor_list <- self$factor_list
      if (!is.null(nodes)) {
        factor_list <- factor_list[nodes]
      }

      all_levels <- lapply(factor_list, function(likelihood_factor) {
        likelihood_factor$variable_type$levels
      })
      all_levels <- all_levels[!(sapply(all_levels, is.null))]
      level_grid <- expand.grid(all_levels)
      return(level_grid)
    },
    base_train = function(task, pretrain) {
      self$validate_task(task)
      fit_object <- private$.train(task, pretrain)
      new_object <- self$clone() # copy parameters, and whatever else
      new_object$set_train(fit_object, task)
      return(new_object)
    },
    add_factors = function(factor_list) {
      #TODO this will fail if the factor_list contains pooled likelihood factors
      if (inherits(factor_list, "LF_base")) {
        factor_list <- list(factor_list)
      }

      factor_names <- sapply(factor_list, `[[`, "name")

      # train factors if necessary
      factor_list <- lapply(factor_list, train_lf, self$training_task)

      # add factors to list of factors
      private$.params$factor_list[factor_names] <- factor_list
    },
    sample = function(tmle_task = NULL, sample_lib = NULL) {
      # for now assume nodes are in order
      # TODO: order nodes based on dependencies
      stop("This doesn't work")
      if (is.NULL(sample_lib = NULL)) {
        nodes <- names(self$factor_list)
        sample_lib <- rep(list(NULL), length(nodes))
        names(sample_lib) <- nodes
      }

      for (node in names(self$factor_list)) {
        tmle_task <- factor_list$node$sample(tmle_task, sample_lib$node)
      }

      return(tmle_task)
    }
  ),
  active = list(
    factor_list_pooled = function(){
      return(self$params$factor_list_pooled)
    },
    factor_list = function() {
      return(self$params$factor_list)
    },
    nodes = function() {
      return(names(self$factor_list))
    },
    cache = function() {
      return(private$.cache)
    },
    censoring_nodes = function() {
      return(private$.censoring_nodes)
    }
    ,
    list_all_predicted_lkd = function() {
      if (!is.null(private$.list_all_predicted_lkd)) {
        return(private$.list_all_predicted_lkd)
      } else {
        # all not A, not t=0 nodes
        tmle_task <- self$training_task
        temp_node_names <- names(tmle_task$npsem)
        loc_A <- grep("A", temp_node_names)
        loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
        loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
        if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

        # get list of all possible predicted lkds
        obs_data <- tmle_task$data
        obs_variable_names <- colnames(obs_data)
        list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
          if (loc_node > 1) {
            # currently only support univariate node for t>0
            current_variable <- tmle_task$npsem[[loc_node]]$variables
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
            temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
            temp_output <- self$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number = "full")  # corresponding outputs
            data.frame(temp_input, output = temp_output) %>% return
          }
        })
        names(list_all_predicted_lkd) <- temp_node_names
        private$.list_all_predicted_lkd <- list_all_predicted_lkd
        return(list_all_predicted_lkd)
      }
    },
    list_all_predicted_lkd_val = function() {
      if (!is.null(private$.list_all_predicted_lkd_val)) {
        return(private$.list_all_predicted_lkd_val)
      } else {
        # all not A, not t=0 nodes
        tmle_task <- self$training_task
        temp_node_names <- names(tmle_task$npsem)
        loc_A <- grep("A", temp_node_names)
        loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
        loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
        if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

        # get list of all possible predicted lkds
        obs_data <- tmle_task$data
        obs_variable_names <- colnames(obs_data)
        list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
          if (loc_node > 1) {
            # currently only support univariate node for t>0
            current_variable <- tmle_task$npsem[[loc_node]]$variables
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
            temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
            temp_output <- self$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number = "validation")  # corresponding outputs
            data.frame(temp_input, output = temp_output) %>% return
          }
        })
        names(list_all_predicted_lkd) <- temp_node_names
        private$.list_all_predicted_lkd_val <- list_all_predicted_lkd
        return(list_all_predicted_lkd)
      }
    }
  ),
  private = list(
    .train_sublearners = function(tmle_task) {
      factor_fits <- lapply(self$factor_list_pooled, function(factor) factor$delayed_train(tmle_task))
      result <- bundle_delayed(factor_fits)
      return(result)
    },
    .train = function(tmle_task, factor_fits) {
      factor_list <- self$factor_list_pooled
      for (i in seq_along(factor_list)) {
        factor_list[[i]]$train(tmle_task, factor_fits[[i]])
      }
      # TODO: mutating factor list of Lrnr_object instead of returning a fit
      #       which is not what sl3 Lrnrs usually do

      censoring_nodes <- lapply(tmle_task$npsem, function(node) {
        node$censoring_node$name
      })

      names(censoring_nodes) <- names(tmle_task$npsem)
      private$.censoring_nodes <- censoring_nodes
      return("trained")
    },
    .predict = function(tmle_task) {
      stop("predict method doesn't work for Likelihood. See Likelihood$get_likelihoods for analogous method")
    },
    .chain = function(tmle_task) {
      stop("chain method doesn't work for Likelihood. Currently, no analogous functionality")
    },
    .cache = NULL,
    .censoring_nodes = NULL
    ,
    .list_all_predicted_lkd = NULL,
    .list_all_predicted_lkd_val = NULL
  )
)

#' @param ... Passes all arguments to the constructor. See documentation for the
#'  Constructor below.
#'
#' @rdname Likelihood
#'
#' @export
#
make_Likelihood <- Likelihood$new

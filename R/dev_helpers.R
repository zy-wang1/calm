reload_calm <- function(path_folder = "~/Documents/projects/calm/", path_product = "./temp") {
    devtools::build(path_folder, path = path_product)
    install.packages("./temp/calm_1.0.tar.gz", repos = NULL, type = "source")
    detach("package:calm", unload = T)
    library(calm)
}

get_psi <- function(estimates, arms = NULL, times = NULL) {
    if (!is.null(arms)) {
        temp_vec_arm <- lapply(arms, function(arm) grep(arm, names(estimates))) %>% unlist %>% unique %>% sort
    } else temp_vec_arm <- 1:length(estimates)
    if (!is.null(times)) {
        temp_vec_time <- lapply(times, function(u) grep(paste0("^", u, "_"), names(estimates))) %>% unlist %>% unique %>% sort
    } else temp_vec_time <- 1:length(estimates)
    temp_vec <- intersect(temp_vec_arm, temp_vec_time)
    temp_psi <- sapply(estimates[temp_vec], function(u) u$psi)
    return(temp_psi)
}
get_IC <- function(estimates, arms = NULL, times = NULL) {
    if (!is.null(arms)) {
        temp_vec_arm <- lapply(arms, function(arm) grep(arm, names(estimates))) %>% unlist %>% unique %>% sort
    } else temp_vec_arm <- 1:length(estimates)
    if (!is.null(times)) {
        temp_vec_time <- lapply(times, function(u) grep(paste0("^", u, "_"), names(estimates))) %>% unlist %>% unique %>% sort
    } else temp_vec_time <- 1:length(estimates)
    temp_vec <- intersect(temp_vec_arm, temp_vec_time)
    temp_IC <- sapply(estimates[temp_vec], function(u) u$IC)
    return(temp_IC)
}


message_parallel <- function(...){
    system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}
message_temp <- function(...){
    system(sprintf('echo "\n%s" >> ./temp/temp', paste0(..., collapse="")))
}
message_temp_seed <- function(...){
    system(sprintf('echo "\n%s" >> ./temp/temp_seed', paste0(..., collapse="")))
}
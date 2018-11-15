#' Function that takes a list of parameter objects and runs them all, with optional parallelization
#' @import pbapply
#' @importFrom readr write_rds
#' @importFrom parallel makeCluster stopCluster clusterEvalQ clusterExport
#' @export sim_radiation_multi
sim_radiation_multi <- function(parm_list, save_tree = TRUE, progress = TRUE, trait_hist = TRUE, trait_hist_prop = 0.01, ncpus = NULL, save_folder = NULL, save_prefix = "sim_", compress = "gz", return_sim = FALSE) {
  
  run_sim <- function(parms) {
    rep <- parms[[2]]
    sim <- try(sim_radiation(parms[[1]], save_tree = save_tree, progress = FALSE, trait_hist = trait_hist,
                             trait_hist_prop = trait_hist_prop))
    if(!is.null(save_folder)) {
      file_name <- paste0(save_folder, "/", save_prefix, rep, "_", gsub(":", "-", gsub(" ", "_", Sys.time())), ".rds")
      write_rds(sim, file_name, compress = compress)
    }
  }
   if(!is.null(ncpus)) {
     cl <- makeCluster(ncpus)
     clusterEvalQ(cl, library(nichefillr))
     clusterEvalQ(cl, library(readr))
     clusterExport(cl, c("save_folder", "save_tree", "trait_hist", "trait_hist_prop", "save_prefix", "compress", "run_sim"), envir=environment())
   } else {
     cl <- NULL
   }
  
  parm_list <- mapply(function(x, y) list(x, y), parm_list, seq_along(parm_list), SIMPLIFY = FALSE)
  
  sim_list <- pblapply(parm_list, run_sim, cl = cl)
  
  if(!is.null(cl)) {
    stopCluster(cl)
  }
  
  if(return_sim) {
    return(sim_list)
  } else {
    return(NULL)
  }
}
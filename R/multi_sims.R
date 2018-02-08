#' Function that takes a list of parameter objects and runs them all, with optional parallelization
#' @import pbapply
#' @importFrom readr write_rds
#' @importFrom parallel makeCluster stopCluster
sim_radiation_multi <- function(parm_list, save_tree = TRUE, progress = TRUE, trait_hist = TRUE, trait_hist_prop = 0.01, ncpus = NULL, save_folder = NULL, save_prefix = "sim_", compress = "gz") {
  
  run_sim <- function(parms) {
    sim <- try(sim_radiation(parms, save_tree = save_tree, progress = FALSE, trait_hist = trait_hist,
                             trait_hist_prop = trait_hist_prop))
    if(!is.null(save_folder)) {
      file_name <- paste0(save_folder, "/", save_prefix, gsub(":", "-", gsub(" ", "_", Sys.time())), ".rds")
      write_rds(sim, file_name, compress = compress)
    }
  }
   if(!is.null(ncpus)) {
     cl <- makeCluster(ncpus)
     clusterEvalQ(cl, library(nichefillr))
     clusterEvalQ(cl, library(readr))
     clusterExport(cl, c("save_folder", "save_tree", "trait_hist", "trait_hist_prop", "save_prefix", "compress"))
   } else {
     cl <- NULL
   }
  sim_list <- pblapply(parm_list, run_sim, cl = cl)
  
  if(!is.null(cl)) {
    stopCluster(cl)
  }
}
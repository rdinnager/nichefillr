#' Generate Uniform Random Parameters for Simulation
#' 
#' Function that takes a set of minimum and maximum values and generates uniform random values
#' for each parameter
#' 
#' @param reps Number of random sets of paramters to generate
#' @param xx_range Vectors of length two, specifying the minimum and maximum values of
#' parameter \code{xx} in the model
#' 
#' @return A data.frame containing the random parameter values
#' 
#' @import dplyr
#' 
#' @export
generate_random_params <- function(reps = 10000,
                                   potent_vol_range = c(1, 20), total_vol_range = c(1, 1), 
                                   num_peaks_range = c(2, 10),
                                   P_min_range = c(0.9, 0.9),
                                   P_max_range = c(1.1, 1.1),
                                   prop_variance_range = c(0.1, 0.5),
                                   d_range = c(2, 2), P_range = c(1, 1),
                                   e_var_range = c(0.1, 0.4),
                                   gamma_range = c(0.01, 1),
                                   b_rate_range = c(0.001, 0.003),
                                   D_range = c(1, 1),
                                   h_to_sig_ratio_range = c(0.5, 2),
                                   a_prop_range = c(0.05, 0.05),
                                   tot_time_range = c(50000, 50000),
                                   V_range = c(0.01, 0.01),
                                   mult_range = c(1, 1)) {
  
  if(length(unique(d_range)) < 2) {
    d_range <- rep(d_range, 2)
  } else {
    d_range <- d_range[1]:d_range[2]
  }
  
  if(length(unique(num_peaks_range)) < 2) {
    num_peaks_range <- rep(num_peaks_range, 2)
  } else {
    num_peaks_range <- num_peaks_range[1]:num_peaks_range[2]
  }
  
  param_df <- data_frame(rep = 1:reps,
                         potent_vol = runif(reps, potent_vol_range[1], potent_vol_range[2]),
                         total_vol = runif(reps, total_vol_range[1], total_vol_range[2]),
                         num_peaks = sample(num_peaks_range, reps, replace = TRUE),
                         P_min = runif(reps, P_min_range[1], P_min_range[2]),
                         P_max = runif(reps, P_max_range[1], P_max_range[2]),
                         prop_variance = runif(reps, prop_variance_range[1], prop_variance_range[2]),
                         d = sample(d_range, reps, replace = TRUE),
                         e_var = runif(reps, e_var_range[1], e_var_range[2]),
                         gamma = runif(reps, gamma_range[1], gamma_range[2]),
                         b_rate = runif(reps, b_rate_range[1], b_rate_range[2]),
                         D = runif(reps, D_range[1], D_range[2]),
                         h_to_sig_ratio = exp(runif(reps, log(h_to_sig_ratio_range[1]),
                                                log(h_to_sig_ratio_range[2]))),
                         P = runif(reps, P_range[1], P_range[2]),
                         a_prop = runif(reps, a_prop_range[1], a_prop_range[2]),
                         tot_time = runif(reps, tot_time_range[1], tot_time_range[2]),
                         V_gi = runif(reps, V_range[1], V_range[2]),
                         mult = runif(reps, mult_range[1], mult_range[2])) %>%
    mutate(dirichlet_param = gen_dirichlet(prop_variance, num_peaks))
    
}

#' Create Lists of Parameters from a Parameter Dataframe
#' 
#' Function to generate parameter lists in list columns using a data.frame with
#' parameter values used to generate random landscapes
#' 
#' @param param_df Dataframe containing parameters to be used to generate random
#' carrying capacity landscapes and for competition and macroevolutionary parts
#' of the simulation
#' 
#' @return A data.frame containing the final parameter lists in a list column. This can be 
#' used by \code{\link{sim_radiation}}
#' 
#' @details Note: this function can take some time as it numerically integrates each landscape
#' to make sure the integral adds up to the correct total volume
#' 
#' @import dplyr
#' 
#' @export
make_params_lists <- function(param_df, trait_hist_prop = 0.3) {
  param_list <- param_df %>%
    rowwise %>%
    do(K_list = generate_landscape(potent_vol = .$potent_vol,
                                   total_vol = .$total_vol,
                                   num_peaks = .$num_peaks,
                                   h_to_sig_ratio = .$h_to_sig_ratio,
                                   P_min_max = c(.$P_min, .$P_max),
                                   dirichlet_param = .$dirichlet_param,
                                   d = .$d,
                                   P = .$P,
                                   a_prop = .$a_prop),
       a_list = list(gamma_i = rep(.$gamma, .$d), D_i = rep(.$D, .$d)),
       macro_list = list(b_rate = .$b_rate, init_traits = rep(0, .$d),
                         e_var = rep(.$e_var, .$d), init_Ns = rep(0.05, 2),
                         init_br = 2, check_extinct = 0.005, tot_time = .$tot_time,
                         V_gi = rep(.$V_gi, .$d),
                         mult = .$mult,
                         m = 2, d = .$d, save_tree = FALSE,
                         save_tree_interval = 100,
                         progress = TRUE,
                         trait_hist = TRUE,
                         trait_hist_prop = trait_hist_prop))
  
  param_df_list <- bind_cols(param_df, param_list)
  
}

#' Run Simulations from a Parameter List Data.frame
#' 
#' Function to take a parameter data.frame with list columns and setup
#' and run the simulation in parallel (using \code{multidplyr})
#' 
#' @param param_df_list A data.frame with list columns containing paramter values to
#' be run in the simulation.
#' @param parallel Logical: Should the simulation be run in parallel?
#' @param cores How many cores to use?
#' @param save_folder Folder to save simulation results to.
#' @param num_tries Number of times to try any simulation if their are errors before
#' giving up.
#' @param prefix Prefix to prepend to filenames.
#' @param compress Compression method to use: "none", "gz" ,"bz", or "xz".
#' 
#' @return A vector of filenames pointing to where each simulation object is saved
#' 
#' @import multidplyr
#' 
#' @export
run_sims <- function(param_df_list, parallel = TRUE, cores = 8, save_folder = "results",
                     num_tries = 10, prefix = "", compress = "gz") {
  
  if(!dir.exists(save_folder)) {
    dir.create(save_folder)
  }
  
  if(parallel) {
    cluster <- create_cluster(cores)
    
    param_df_list <- param_df_list %>%
      mutate(group = gl(cores, ceiling(nrow(param_df_list) / (cores)), nrow(param_df_list)),
             file_name = paste0(save_folder, "/", prefix, "_sim_result_rep_", rep, ".rds"))
    
    param_part <- partition(param_df_list, group, cluster = cluster)
    cluster_library(param_part, "nichefillr")
    cluster_assign_value(param_part, "num_tries", num_tries)
    
    result <- param_part %>%
      group_by(rep) %>%
      do(nichefillr:::run_sim_and_save(parms = list(K_parms = .$K_list[[1]],
                                       a_parms = .$a_list[[1]],
                                       macro_parms = .$macro_list[[1]]),
                          file_name = .$file_name[[1]],
                          num_tries = num_tries))
    
    result_data <- collect(result)
    
    # test <- nichefillr:::run_sim_and_save(parms = list(K_parms = param_df_list$K_list[[1]],
    #                                                    a_parms = param_df_list$a_list[[1]],
    #                                                    macro_parms = param_df_list$macro_list[[1]]),
    #                                       file_name = param_df_list$file_name[[1]],
    #                                       num_tries = 10)
    # 
    # test <- nichefillr:::run_sim_and_save(parms = example_parms,
    #                                       file_name = param_df_list$file_name[[1]])
    
  } else {
    result_data <- param_df_list %>%
      mutate(file_name = paste0(save_folder, "/", prefix, "_sim_result_rep_", rep, ".rds")) %>%
      group_by(rep) %>%
      do(nichefillr:::run_sim_and_save(parms = list(K_parms = .$K_list[[1]],
                                                    a_parms = .$a_list[[1]],
                                                    macro_parms = .$macro_list[[1]]),
                                       file_name = .$file_name[[1]],
                                       num_tries = num_tries, compress = compress))
  }
  result_data
}

#' Helper Function to Run Simulation and Save It
#' 
#' @param parms Simulation Parameters.
#' @param file_name Filename to save simulation results to.
#' @param num_tries Number of times to try each simulation in the case of errors before
#' giving up and moving on.
#' @param compress Compression method to use: "none", "gz" ,"bz", or "xz".
#' 
#' @return A data.frame containing the filename saved to.
#' 
#' @importFrom readr write_rds
run_sim_and_save <- function(parms, file_name, num_tries, compress = "gz") {
  error_test <- TRUE
  try_num <- 1
  while(error_test & try_num <= num_tries) {
    result <- try(sim_radiation(parms), silent = TRUE)
    if(class(result) != "try-error") {
      error_test <- FALSE
    } else {
      try_num <- try_num + 1L
    }
  }
  if(class(result) == "try-error") {
    error_message <- result
    sr <- NA
  } else {
    error_message <- "No errors"
    sr <- sum(result$sim_object$extant)
  }
  write_rds(result, file_name, compress)
  data_frame(file_name = file_name, error = error_message, num_tries = try_num,
             final_SR = sr)
}

#' Run Simulations for Paper #1
#' 
#' Function to run simulations used in Paper #1
#' 
#' @param num_reps Number of simulation reps to run.
#' @param num_cores Number of cores to use in parallel computations.
#' @param folder Name of folder to place output files.
#' @param prefix Prefix of type character to prepend to any filenames.
#' @param ... Other parameters passed onto \code{\link{generate_random_params}}
#' 
#' @return None
#' 
#' @importFrom readr write_csv
#' 
#' @export
run_sims_paper_1 <- function(num_reps, num_cores, folder = "results", prefix, compress = "gz", param_df_list = NULL, ...) {
  
  
  ## generate 10000 sets of random parameters
  if(is.null(param_df_list)) {
    param_df <- generate_random_params(num_reps, ...)
  
  ## Use parameters to generate random landscapes and put into list columns
  ## This could take up to an hour for 10000 reps because each landscape is numerically
  ## integrated to make sure their volumes are standardized
    param_df_list <- make_params_lists(param_df)
  }
  ## run simulations!
  
  results <- run_sims(param_df_list, cores = num_cores, save_folder = folder,
                      prefix = prefix,
                      compress = compress)
  write_csv(results, paste0(folder, "/", prefix, "_results_data.csv"))
}

#' #' Generate Quasi Monte Carlo Random Parameters for Simulation
#' #' 
#' #' Function that takes a set of minimum and maximum values and generates uniform random values
#' #' for each parameter within a hypercube using QuasiMC
#' #' 
#' #' @param reps Number of random sets of paramters to generate
#' #' @param xx_range Vectors of length two, specifying the minimum and maximum values of
#' #' parameter \code{xx} in the model
#' #' 
#' #' @return A data.frame containing the random parameter values
#' #' 
#' #' @import dplyr
#' #' @importFrom randtoolbox sobol
#' #' 
#' #' @export
#' generate_random_params_QMC <- function(reps = 10000,
#'                                    potent_vol_range = c(1, 20), total_vol_range = c(5, 5), 
#'                                    num_peaks_range = c(2L, 10L),
#'                                    P_min_range = c(0.9, 0.9),
#'                                    P_max_range = c(1.1, 1.1),
#'                                    prop_variance_range = c(0.1, 0.5),
#'                                    d_range = c(2L, 2L), P_range = c(1, 1),
#'                                    e_var_range = c(0.1, 0.4),
#'                                    gamma_range = c(0.01, 1),
#'                                    b_rate_range = c(0.001, 0.003),
#'                                    D_range = c(1, 1),
#'                                    h_to_sig_ratio_range = c(0.5, 2),
#'                                    a_prop_range = c(0.05, 0.05),
#'                                    tot_time_range = c(50000, 50000),
#'                                    V_range = c(0.01, 0.01),
#'                                    mult_range = c(1, 1)) {
#'   
#'   h_to_sig_ratio_range <- log(h_to_sig_ratio_range)
#'   
#'   varying <- c(potent_vol = ifelse(potent_vol_range[2] - potent_vol_range[1] > 0, TRUE, FALSE),
#'                total_vol = ifelse(total_vol_range[2] - total_vol_range[1] > 0, TRUE, FALSE),
#'                P_min = ifelse(P_min_range[2] - P_min_range[1] > 0, TRUE, FALSE),
#'                P_max = ifelse(P_max_range[2] - P_max_range[1] > 0, TRUE, FALSE),
#'                prop_variance = ifelse(prop_variance_range[2] - prop_variance_range[1] > 0, TRUE, FALSE),
#'                e_var = ifelse(e_var_range[2] - e_var_range[1] > 0, TRUE, FALSE),
#'                gamma = ifelse(gamma_range[2] - gamma_range[1] > 0, TRUE, FALSE),
#'                b_rate = ifelse(b_rate_range[2] - b_rate_range[1] > 0, TRUE, FALSE),
#'                D = ifelse(D_range[2] - D_range[1] > 0, TRUE, FALSE),
#'                h_to_sig_ratio = ifelse(h_to_sig_ratio_range[2] - h_to_sig_ratio_range[1] > 0, TRUE, FALSE),
#'                a_prop = ifelse(a_prop_range[2] - a_prop_range[1] > 0, TRUE, FALSE),
#'                tot_time = ifelse(tot_time_range[2] - tot_time_range[1] > 0, TRUE, FALSE),
#'                V_gi = ifelse(V_range[2] - V_range[1] > 0, TRUE, FALSE),
#'                mult = ifelse(mult_range[2] - mult_range[1] > 0, TRUE, FALSE),
#'                P = ifelse(P_range[2] - P_range[1] > 0, TRUE, FALSE))
#'   
#'   ranges <- list(potent_vol = potent_vol_range,
#'               total_vol = total_vol_range,
#'               P_min = P_min_range,
#'               P_max = P_max_range,
#'               prop_variance = prop_variance_range,
#'               e_var = e_var_range,
#'               gamma = gamma_range,
#'               b_rate = b_rate_range,
#'               D = D_range,
#'               h_to_sig_ratio = h_to_sig_ratio_range,
#'               a_prop = a_prop_range,
#'               tot_time = tot_time_range,
#'               V_gi = V_range,
#'               mult = mult_range,
#'               P = P_range)
#'   
#'   qmc_dim <- sum(varying)
#'   
#'   rand_pts <- sobol(reps, qmc_dim, scrambling = 3)
#'   colnames(rand_pts) <- names(varying[varying])
#'   
#'   for(i in 1:ncol(rand_pts)) {
#'     rand_pts[ , i] <- (rand_pts[ , i] * (ranges[[colnames(rand_pts)[i]]][2] - ranges[[colnames(rand_pts)[i]]][1])) + 
#'       ranges[[colnames(rand_pts)[i]]][1]
#'   }
#'   
#'   if(length(unique(d_range)) < 2) {
#'     d_range <- rep(d_range, 2)
#'   } else {
#'     d_range <- d_range[1]:d_range[2]
#'   }
#'   
#'   if(length(unique(num_peaks_range)) < 2) {
#'     num_peaks_range <- rep(num_peaks_range, 2)
#'   } else {
#'     num_peaks_range <- num_peaks_range[1]:num_peaks_range[2]
#'   }
#'   
#'   other_dim <- sum(varying == FALSE)
#'   
#'   non_varying_params <- replicate(other_dim, rep(1, reps))
#'   colnames(non_varying_params) <- names(varying[!varying])
#'   
#'   for(i in 1:ncol(non_varying_params)) {
#'     non_varying_params[ , i] <- (non_varying_params[ , i] * (ranges[[colnames(non_varying_params)[i]]][2] - ranges[[colnames(non_varying_params)[i]]][1])) + 
#'       ranges[[colnames(non_varying_params)[i]]][1]
#'   }
#'   
#'   param_df <- data_frame(rep = 1:reps,
#'                          num_peaks = sample(num_peaks_range, reps, replace = TRUE),
#'                          d = sample(d_range, reps, replace = TRUE)) %>%
#'     bind_cols(rand_pts %>% as.data.frame()) %>%
#'     bind_cols(non_varying_params %>% as.data.frame()) %>%
#'     mutate(dirichlet_param = nichefillr:::gen_dirichlet(prop_variance, num_peaks)) %>%
#'     transform(h_to_sig_ratio = exp(h_to_sig_ratio))
#'   
#' }
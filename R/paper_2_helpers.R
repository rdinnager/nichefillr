#' Generate Hopeful Niche Monster parameters
#' 
#' Need to fill this in ... ignore the below, it is wrong.
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
generate_all_params <- function(potent_vol = 1, total_vol = 1, d = 2L, mult = 0.1, gamma_i = 1, 
                                V_gi = 0.01, b_rate = 0.01, total_time = 50000, D = 1.25) {
  
  ## use the equation for the volume of an hyper-cube to get super-gaussian radius equal to desired
  ## total volume
  #total_rad <- (potent_vol / ((pi^(d/2))/(gamma((d/2)+1))))^(1/d)
  
  total_rad <- (potent_vol ^ (1/d)) / 2
  
  ## use the full width at half maximum formula to guestimate sigma for total landscape super-gaussian
  P <- 20
  total_sig <- (2 * total_rad) / (2*((2*log(2))^(1/(2*P)))*sqrt(2))
  
  # integral <- nichefillr:::integrate_fun(nichefillr:::K_func, lower = rep(-total_rad*2, d),
  #                                        upper = rep(total_rad*2, d),
  #                           h0 = 1, sig0i = rep(total_sig, d), P0i = rep(P, d),
  #                           hz = 1, biz = matrix(rep(0, d), nrow = d),
  #                           sigiz = matrix(rep(100, d), nrow = d),
  #                           Piz = matrix(rep(1, d), nrow = d), a = 0)
  
  integr_fun <- function(sig, target, total_rad) {
    abs(target - integrate_fun(K_func, lower = rep(-total_rad, d), 
                               upper = rep(total_rad, d),
                               h0 = 1, sig0i = rep(100, d), P0i = rep(P, d),
                               hz = 1, biz = matrix(rep(0, d), nrow = d),
                               sigiz = matrix(rep(sig, d), nrow = d),
                               Piz = matrix(rep(1, d), nrow = d), a = 0))
  }
  #integr_fun(total_rad*5, 5, total_rad)
  
  ## find sigma for fitness hill that gives volume equal to total_vol
  find_sig <- optim(0.5, integr_fun, lower = 0.001, upper = total_rad*100,
                    method = "L-BFGS-B",
                    target = total_vol, total_rad = total_rad)
  
  # integral2 <- nichefillr:::integrate_fun(nichefillr:::K_func, lower = rep(-total_rad*2, d),
  #                                        upper = rep(total_rad*2, d),
  #                                        h0 = 1, sig0i = rep(total_sig, d), P0i = rep(P, d),
  #                                        hz = 1, biz = matrix(rep(0, d), nrow = d),
  #                                        sigiz = matrix(rep(find_sig$par, d), nrow = d),
  #                                        Piz = matrix(rep(1, d), nrow = d), a = 0)
  
  params <- list(K_parms = list(h0 = 1, hz = 1, biz = matrix(rep(0, d), nrow = d),
                                sigiz = matrix(rep(find_sig$par, d), nrow = d),
                                Piz = matrix(rep(1, d), nrow = d),
                                sig0i = rep(total_sig, d),
                                P0i = rep(P, d),
                                a = 0),
                 a_parms = list(gamma_i = rep(gamma_i, d),
                                D_i = rep(D, d)),
                 macro_parms = list(b_rate = b_rate,
                                    init_traits = rep(0, d),
                                    e_var = rep(total_sig/2, d),
                                    init_Ns = rep(0.05, 2),
                                    init_br = 1,
                                    check_extinct = 0.005,
                                    tot_time = total_time,
                                    V_gi = rep(V_gi, d),
                                    m = 2L,
                                    d = d,
                                    mult = mult,
                                    save_tree = TRUE,
                                    save_tree_interval = 500,
                                    progress = TRUE,
                                    trait_hist = TRUE),
                 other_parms = list(potent_vol = potent_vol, total_vol = total_vol))
  
  params
  
}

#' Generate Hopeful Niche Monster random parameters
#' 
#' Need to fill this in ... ignore the below, it is wrong.
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
#' @import dplyr
#' 
#' @export
generate_random_params_2 <- function(nreps = 10, total_vol_range = c(0.05, 1),
                                     mult_range = c(0.01, 1),
                                     gamma_range = c(0.01, 1),
                                     V_gi_range = c(0.001, 0.1),
                                     d_range = c(1L, 4L),
                                     b_rate_range = c(0.0001, 0.01),
                                     tot_time_range = c(20000, 100000)) {
  
  param_df <- data_frame(rep = 1:nreps, total_vol = runif(nreps, total_vol_range[1], total_vol_range[2]),
                         mult = runif(nreps, mult_range[1], mult_range[2]),
                         gamma_i = runif(nreps, gamma_range[1], gamma_range[2]),
                         V_gi = runif(nreps, V_gi_range[1], V_gi_range[2]),
                         d = sample(d_range[1]:d_range[2], nreps, replace = TRUE),
                         b_rate = runif(nreps, b_rate_range[1], b_rate_range[2]),
                         total_time = runif(nreps, tot_time_range[1], tot_time_range[2]),
                         D = sample(c(0.8, 1.25), nreps, replace = TRUE)) %>%
    rowwise %>%
    do(rep = .$rep, param_list = generate_all_params(total_vol = .$total_vol,
                                        d = .$d,
                                        mult = .$mult,
                                        gamma_i = .$gamma_i,
                                        V_gi = .$V_gi,
                                        b_rate = .$b_rate,
                                        total_time = .$total_time,
                                        D = .$D)) %>%
    transform(rep = unlist(rep))
  
  param_df
  
  
}

#' Run Simulations from a Parameter List Data.frame for Hopeful Niche Monster Paper
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
run_sims_2 <- function(param_df, parallel = TRUE, cores = 8, save_folder = "results",
                     num_tries = 10, prefix = "", compress = "gz") {
  
  if(!dir.exists(save_folder)) {
    dir.create(save_folder)
  }
  
  if(parallel) {
    cluster <- create_cluster(cores)
    
    param_df <- param_df %>%
      mutate(group = gl(cores, ceiling(nrow(param_df) / (cores)), nrow(param_df)),
             file_name = paste0(save_folder, "/", prefix, "_sim_result_rep_", rep, ".rds"))
    
    param_part <- partition(param_df, group, cluster = cluster)
    cluster_library(param_part, "nichefillr")
    cluster_assign_value(param_part, "num_tries", num_tries)
    
    result <- param_part %>%
      group_by(rep) %>%
      do(nichefillr:::run_sim_and_save(parms = .$param_list[[1]],
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
    result_data <- param_df %>%
      mutate(group = gl(cores, ceiling(nrow(param_df) / (cores)), nrow(param_df)),
             file_name = paste0(save_folder, "/", prefix, "_sim_result_rep_", rep, ".rds"))%>%
      group_by(rep) %>%
      do(nichefillr:::run_sim_and_save(parms = .$param_list[[1]],
                                       file_name = .$file_name[[1]],
                                       num_tries = num_tries, compress = compress))
  }
  result_data
}

# test <- sim_radiation(generate_all_params(total_vol = 0.75, mult=0.05, gamma_i = 1, V_gi = 0.001, 
#                                           total_time = 50000, d = 2, b_rate = 0.001))
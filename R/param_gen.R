#' Function to generate example parameter objects.
#' @param name Name of example parameter set. See details for list of available options (not done yet).
#' @return A list containing a set of simulation parameters.
#' @export example_parms
example_parms <- function(name = c("simple_peaks_2", "blackhole_comp")) {
  params <- switch(name, 
                   list(simple_peaks_2 =
                          list(K_parms = list(h0 = 1,
                                              hz = c(3.169580, 2.791599, 3.082317, 1.876961, 2.981982, 3.655125, 2.238074, 3.149038, 1.718011, 2.096480),
                                              biz = matrix(c(9.739022, -0.616485, -1.224424, -2.154302, 0.559290, 1.645315,  5.403136, 0.1929105, -2.329065, -6.092679,
                                                             -5.565801,  1.241130, -3.010253, -1.343778, 3.892568, 3.143613, -6.108992, 3.8414613, -1.944608,  2.926509),
                                                             ncol = 10, byrow = TRUE),
                                              sigiz = matrix(c(1.1738626, 0.9748214, 1.073415, 0.9556417, 1.045984, 0.8519923, 1.3346488, 1.1340722, 0.9071014, 1.2096623,
                                                               0.8475463, 1.1222201, 1.140055, 1.1014532, 1.199568, 0.8232535, 0.9490592, 0.9499395, 0.9043750, 0.8512329),
                                                             ncol = 10, byrow = TRUE),
                                              Piz = matrix(c(1.015310, 0.9854997, 1.0317821, 1.083120, 1.0530023, 0.9852945, 0.9110535, 1.1364469, 1.006924, 1.0311518,
                                                             1.142575, 1.0721543, 0.7942551, 1.052459, 0.9825081, 1.1663972, 0.7700944, 0.9466385, 1.000279, 0.9650881),
                                                           ncol = 10, byrow = TRUE),
                                              sig0i = c(5, 5),
                                              P0i = c(1.5, 1.5),
                                              a = 0.5,
                                              c_var = 0.1),
                               a_parms = list(gamma_i = c(0.9, 0.9),
                                              D_i = c(1, 1),
                                              C = 1),
                               macro_parms = list(b_rate = 0.001,
                                                  init_traits = c(3.866015, 1.141986),
                                                  e_var = c(0.2, 0.2),
                                                  init_Ns = c(0.05, 0.05),
                                                  init_br = 20,
                                                  check_extinct = 0.005,
                                                  tot_time = 50000,
                                                  V_gi = c(0.01, 0.01),
                                                  m = 2,
                                                  d = 2,
                                                  mult = 0.5))))
  if(length(name) == 1) {
    params <- params[[1]]
  }
}


#' Generate Quasi Monte Carlo Random Parameters for Simulation
#' 
#' Function that takes a set of minimum and maximum values and generates uniform random values
#' for each parameter within a hypercube using QuasiMC
#' 
#' @param reps Number of random sets of paramters to generate
#' @param xx_range Vectors of length two, specifying the minimum and maximum values of
#' parameter \code{xx} in the model
#' 
#' @return A data.frame containing the random parameter values
#' 
#' @import dplyr
#' @importFrom randtoolbox sobol
#' 
#' @export
generate_random_params_QMC <- function(reps = 10000,
                                       potent_vol_range = c(1, 20), total_vol_range = c(5, 5), 
                                       num_peaks_range = c(2L, 10L),
                                       P_min_range = c(0.9, 0.9),
                                       P_max_range = c(1.1, 1.1),
                                       prop_variance_range = c(0.1, 0.5),
                                       d_range = c(2L, 2L), P_range = c(1, 1),
                                       e_var_range = c(0.1, 0.4),
                                       gamma_range = c(0.01, 1),
                                       b_rate_range = c(0.001, 0.003),
                                       D_range = c(1, 1),
                                       h_to_sig_ratio_range = c(0.5, 2),
                                       a_prop_range = c(0.05, 0.05),
                                       tot_time_range = c(50000, 50000),
                                       V_range = c(0.01, 0.01),
                                       mult_range = c(1, 1),
                                       C_range = c(1, 1),
                                       c_var_range = c(0.1, 0.1)) {
  
  h_to_sig_ratio_range <- log(h_to_sig_ratio_range)
  
  varying <- c(potent_vol = ifelse(potent_vol_range[2] - potent_vol_range[1] > 0, TRUE, FALSE),
               total_vol = ifelse(total_vol_range[2] - total_vol_range[1] > 0, TRUE, FALSE),
               P_min = ifelse(P_min_range[2] - P_min_range[1] > 0, TRUE, FALSE),
               P_max = ifelse(P_max_range[2] - P_max_range[1] > 0, TRUE, FALSE),
               prop_variance = ifelse(prop_variance_range[2] - prop_variance_range[1] > 0, TRUE, FALSE),
               e_var = ifelse(e_var_range[2] - e_var_range[1] > 0, TRUE, FALSE),
               gamma = ifelse(gamma_range[2] - gamma_range[1] > 0, TRUE, FALSE),
               b_rate = ifelse(b_rate_range[2] - b_rate_range[1] > 0, TRUE, FALSE),
               D = ifelse(D_range[2] - D_range[1] > 0, TRUE, FALSE),
               h_to_sig_ratio = ifelse(h_to_sig_ratio_range[2] - h_to_sig_ratio_range[1] > 0, TRUE, FALSE),
               a_prop = ifelse(a_prop_range[2] - a_prop_range[1] > 0, TRUE, FALSE),
               tot_time = ifelse(tot_time_range[2] - tot_time_range[1] > 0, TRUE, FALSE),
               V_gi = ifelse(V_range[2] - V_range[1] > 0, TRUE, FALSE),
               mult = ifelse(mult_range[2] - mult_range[1] > 0, TRUE, FALSE),
               P = ifelse(P_range[2] - P_range[1] > 0, TRUE, FALSE),
               C = ifelse(C_range[2] - C_range[1] > 0, TRUE, FALSE),
               c_var = ifelse(c_var_range[2] - c_var_range[1] > 0, TRUE, FALSE))
  
  ranges <- list(potent_vol = potent_vol_range,
                 total_vol = total_vol_range,
                 P_min = P_min_range,
                 P_max = P_max_range,
                 prop_variance = prop_variance_range,
                 e_var = e_var_range,
                 gamma = gamma_range,
                 b_rate = b_rate_range,
                 D = D_range,
                 h_to_sig_ratio = h_to_sig_ratio_range,
                 a_prop = a_prop_range,
                 tot_time = tot_time_range,
                 V_gi = V_range,
                 mult = mult_range,
                 P = P_range,
                 C = C_range,
                 c_var = c_var_range)
  
  qmc_dim <- sum(varying)
  
  rand_pts <- sobol(reps, qmc_dim, scrambling = 3)
  colnames(rand_pts) <- names(varying[varying])
  
  for(i in 1:ncol(rand_pts)) {
    rand_pts[ , i] <- (rand_pts[ , i] * (ranges[[colnames(rand_pts)[i]]][2] - ranges[[colnames(rand_pts)[i]]][1])) + 
      ranges[[colnames(rand_pts)[i]]][1]
  }
  
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
  
  other_dim <- sum(varying == FALSE)
  
  non_varying_params <- replicate(other_dim, rep(1, reps))
  colnames(non_varying_params) <- names(varying[!varying])
  
  for(i in 1:ncol(non_varying_params)) {
    non_varying_params[ , i] <- (non_varying_params[ , i] * (ranges[[colnames(non_varying_params)[i]]][2] - ranges[[colnames(non_varying_params)[i]]][1])) + 
      ranges[[colnames(non_varying_params)[i]]][1]
  }
  
  param_df <- data_frame(rep = 1:reps,
                         num_peaks = sample(num_peaks_range, reps, replace = TRUE),
                         d = sample(d_range, reps, replace = TRUE)) %>%
    bind_cols(rand_pts %>% as.data.frame()) %>%
    bind_cols(non_varying_params %>% as.data.frame()) %>%
    mutate(dirichlet_param = nichefillr:::gen_dirichlet(prop_variance, num_peaks)) %>%
    transform(h_to_sig_ratio = exp(h_to_sig_ratio))
  
  param_list <- param_df %>%
    rowwise %>%
    do(params = list(K_parms = c(generate_landscape(potent_vol = .$potent_vol,
                                   total_vol = .$total_vol,
                                   num_peaks = .$num_peaks,
                                   h_to_sig_ratio = .$h_to_sig_ratio,
                                   P_min_max = c(.$P_min, .$P_max),
                                   dirichlet_param = .$dirichlet_param,
                                   d = .$d,
                                   P = .$P,
                                   a_prop = .$a_prop), c_var = .$c_var),
       a_parms = list(gamma_i = rep(.$gamma, .$d), D_i = rep(.$D, .$d), C = .$C),
       macro_parms = list(b_rate = .$b_rate, init_traits = rep(0, .$d),
                         e_var = rep(.$e_var, .$d), init_Ns = rep(0.05, 2),
                         init_br = 2, check_extinct = 0.005, tot_time = .$tot_time,
                         V_gi = rep(.$V_gi, .$d),
                         mult = .$mult,
                         m = 2, d = .$d)))
  
  return(param_list$params)
}
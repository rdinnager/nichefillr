#' Convert ODE Output into a Data.frame
#' 
#' Function to convert the ode objects outputted by \code{deSolve} into a data.frame
#' 
#' @param ode_mat An ode object matrix
#' @param d The number of trait dimensions in the simulation
#' @param m The number of species in the simulation
#' @param thin How much to thin the ode output, returns 1 out of every \code{thin} timesteps
#' 
#' @return A data.frame with simulation output
#' 
#' @import dplyr
make_ode_df <- function(ode_mat, d, m, thin = 100) {
  split_list <- c(0, as.numeric(gl(m, d)), 1:m)
  test_ode_df <- as.data.frame(ode_mat)
  ode_df <- do.call(rbind, lapply(1:max(split_list), function(x) cbind(ode_mat[ , 1], ode_mat[ , split_list == x])))
  ode_df <- as.data.frame(ode_df) %>%
    mutate(Species = paste0("Species_", as.numeric(gl(m, nrow(ode_mat)))))
  colnames(ode_df)[1:(d+2)] <- c("Time", paste0("Niche_Axis_", 1:d), "Population") 
  ## thin times
  ode_df_thinned <- ode_df %>%
    filter(Time %in% seq(1, nrow(ode_mat), by = thin)) %>%
    mutate(SR = m)
  ode_df_thinned
}

#' Extract Data.frame of Simulation Data
#' 
#' Function to extract useful data from a simulation object and return it in a data.frame
#' 
#' @param full_sim_ob A simulation object as produced by \code{\link{sim_radiation}}
#' @param thin_to How many timesteps of data to return. The function will thin the output to match
#' this length.
#' 
#' @return A data.frame with simulation data for \code{thin_to} timesteps
#' 
#' @import dplyr
#' 
#' @export
make_sim_df <- function(full_sim_ob, thin_to = 10000, include_trees = TRUE,
                        summarise_trees = FALSE,
                        PD = "PSE") {
  #require(dplyr)
  sim_ob <- full_sim_ob$sim_object
  extant_list <- sim_ob$extant_list[!sapply(sim_ob$extant_list, is.null)]
  full_dat <- sim_ob$full_dat[!sapply(sim_ob$full_dat, is.null)]
  ms <- sapply(extant_list, length)
  dfs <- mapply(make_ode_df, ode_mat = full_dat, d = full_sim_ob$params$macro_parms$d, m = ms, thin = 1, SIMPLIFY = FALSE)
  
  for(i in 2:length(dfs)) {
    dfs[[i]]$Time <- dfs[[i]]$Time + max(dfs[[i - 1]]$Time)
  }
  
  big_df <- bind_rows(dfs)
  
  # if(summarise_trees) {
  #   tip_names <- mapply(function(tree, extant) tree$tip.label[extant],
  #                       )
  # }
  
  ## thin entries
  if(full_sim_ob$params$macro_parms$save_tree) {
    tree_times <- floor(full_sim_ob$saved_trees$tree_times[full_sim_ob$saved_trees$tree_times != 0])
    tree_time_vec <- c(ceiling(seq(1, nrow(big_df), length.out = thin_to)), tree_times)
    big_df <- big_df %>%
      filter(Time %in% tree_time_vec)
  } else {
    big_df <- big_df %>%
      filter(Time %in% ceiling(seq(1, nrow(big_df), length.out = thin_to)))
  }
  
  ## add trees
  if(full_sim_ob$params$macro_parms$save_tree & include_trees) {
    time_vec <- findInterval(tree_time_vec, tree_times, left.open = FALSE) + 1
    tree_df <- data_frame(Time = tree_time_vec, Tree = full_sim_ob$saved_trees$trees[time_vec])
    big_df <- big_df %>%
      left_join(tree_df)
    if(summarise_trees) {
      if(PD == "PSE") {
        
      }
    } else {
      
    }
  }
  big_df
  
}

#' Extract Final Species Richness of Simulation
#' 
#' Function to extract the final species richness in a simulation
#' 
#' @param sim_ob Simulation object from which to extract the final species richness
#' 
#' @return The final species richness value in the simulation
#' 
#' @export
extract_final_SR <- function(sim_ob) {
  sum(sim_ob$sim_object$extant)
}

#' Extract Useful Parameters of Simulation Object
#' 
#' Function to extract the useful parameters used in the simulation
#' 
#' @param sim_ob Simulation object from which to extract the final species richness
#' 
#' @return Single row data.frame of parameter values
#' 
#' @export
extract_params <- function(sim_ob) {
  data_frame(total_rad = sim_ob$params$K_parms$total_rad,
             potent_vol = sim_ob$params$K_parms$potent_vol,
             total_vol = sim_ob$params$K_parms$total_vol,
             num_peaks = sim_ob$params$K_parms$num_peaks,
             h_to_sig_ratio = sim_ob$params$K_parms$h_to_sig_ratio,
             dirichlet_param = sim_ob$params$K_parms$dirichlet_param,
             gamma = sim_ob$params$a_parms$gamma_i[1],
             d = sim_ob$params$macro_parms$d,
             e_var = sim_ob$params$macro_parms$e_var[1],
             b_rate = sim_ob$params$macro_parms$b_rate,
             tot_time = sim_ob$params$macro_parms$tot_time,
             V_gi = sim_ob$params$macro_parms$V_gi[1])
}



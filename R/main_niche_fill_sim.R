#' Main Niche Filling Simulation Function
#' 
#' Function to setup and run the niche filling simulation
#' 
#' @param parms A named list containing three named list elements:
#' \itemize{
#' \item{"K_parms"}{: Named list of parameters relating to the carrying capacity landscape function}
#' \item{"a_parms"}{: Named list of parameters relating to the competition function}
#' \item{"macro_parms"}{: Named list of parameters relating to the macroevolutionary simulation}
#' }
#' 
#' See Parameter List section for details on what parameters are in which lists.
#' 
#' @section Parameter List:
#' There are three named lists in the \code{parms} parameter, relating to the
#' carrying capacity, the competition, and the macroevolutionary model. The parameters found
#' in each one are as follows:
#' @section \code{K_parms}:
#' \itemize{
#' \item{\code{h0}}{: Total maximum height of the carrying capacity landscape}
#' \item{\code{hz}}{: Maximum height of each of u peaks in the landscape; vector of length u}
#' \item{\code{biz}}{: Centre of each of u peaks in the landscape for each of d dimensions; 
#' matrix of dimension d by u}
#' \item{\code{sigiz}}{: Width of each of u peaks for each of d dimensions; matrix of length d by u}
#' \item{\code{Piz}}{: Super-gaussian parameter for each of u peaks for each of d dimensions; 
#' matrix of dimension d by u}
#' \item{\code{sig0i}}{: Total width of landscape in all dimensions; determines how far from zero the carrying capacity
#' drops off to nearly zero; vector of length d}
#' \item{\code{P0i}}{: Total super-gaussian parameter; determines how quickly or gradually the carrying
#' capacity drops off near the landscape borders in each dimension. Higher values give more extreme drop-offs;
#' vector of length d}
#' \item{\code{a}}{: Minimum values of carrying capacity within landscape limits}
#' }
#' @section \code{a_parms}:
#' \itemize{
#' \item{\code{gamma_i}}{: Strength of competition; determines how quickly the competition
#' strength between two species drops off with increasing distance between their niche
#' trait values; vector of length d, where d is the number of niche dimensions}
#' \item{\code{D_i}}{: Super-gaussiain parameter for competition; this parameter controls
#' how the precipitously the competition strength drops off with increasing niche distance.
#' Higher values lead to a more cliff-like competition kernel, with strong competition
#' between highly similar species, then very little competition beyond some threshold}
#' }
#' @section \code{macro_parms}:
#' \itemize{
#' \item{\code{b_rate}}{: Pure birth rate for phylogeny simulation; the rate at which new species
#' form by splitting}
#' \item{\code{init_traits}}{: Initial trait values for ancestral species; vector of length 
#' n_traits, where n_traits is the number of traits being simulated}
#' \item{\code{e_var}}{: Evolutionary rates for all traits; vector of length n_traits}
#' \item{\code{init_Ns}}{: Starting population sizes for two initial species in simulation,
#' after first split with ancestral species; vector of length 2}
#' \item{\code{init_br}}{: Time to ancestral species for initial speciation split}
#' \item{\code{check_extinct}}{: Rate at which to check for species whose population size
#' has dropped below a threshold (currently hard-coded), which are then set as extinct}
#' \item{\code{tot_time}}{: The total amount of time to run the simulation in simulation time
#' steps}
#' \item{\code{V_gi}}{: Mutation rate for niche trait evolution; determines how quickly
#' traits can evolve in the trait evolution simulation; typically this is set very low
#' so that evolution proceeds much slower than population dynamics; vector of length n_trait}
#' #' \item{\code{mult}}{: Multiplier that determines how much traits 'jump' in trait space
#' after a speciation event. This is a multiple of \code{e_var}, and is the standard deviation
#' of a normal deviate that is added to all traits.}
#' \item{\code{save_tree}}{: A logical determining whether the simulation should save
#' intermediate states of the phylogeny during the simulation run}
#' \item{\code{save_tree_interval}}{: How often to save intermediate phylogenies if
#' \code{save_tree} = TRUE; phylogenies are saved every \code{save_tree_interval} time steps}
#' \item{\code{progress}}{: Print progress bar if TRUE.}
#' #' \item{\code{trait_hist}}{: A logical determining whether the simulation should save
#' all trait evolution histories. Set to FALSE if you are only interested in the end state of
#' the simulation in order to save memory.}
#' }
#' 
#' @return A niche_fill_sim object containing the final simulation object,
#' a set of intermediate phylogenies, and the parameters used to run the simulation
#' 
#' @importFrom deSolve ode
#' @import progress
#' 
#' @export
sim_radiation <- function(parms) {
  
  params <- list(K_parms = parms$K_parms, a_parms = parms$a_parms, b_rate = parms$macro_parms$b_rate, 
                 init_traits = parms$macro_parms$init_traits, e_var = parms$macro_parms$e_var, 
                 init_Ns = parms$macro_parms$init_Ns, init_br = parms$macro_parms$init_br, 
                 check_extinct = parms$macro_parms$check_extinct, tot_time = parms$macro_parms$tot_time, 
                 V_gi = parms$macro_parms$V_gi, m = parms$macro_parms$m, d = parms$macro_parms$d, 
                 save_tree = parms$macro_parms$save_tree, save_tree_interval = parms$macro_parms$save_tree_interval, 
                 progress = parms$macro_parms$progress, mult = parms$macro_parms$mult,
                 trait_hist = parms$macro_parms$trait_hist, trait_hist_prop = parms$macro_parms$trait_hist_prop)
  
  tree_list <- list()
  save_tree_times <- vector("numeric", length = 1000000)
  
  tree_ob <- make_tree_ob(b_rate = params$b_rate, params$init_traits, 
                          params$e_var, init_Ns = params$init_Ns, 
                          trait_hist = params$trait_hist)  
  tree_ob <- update_br_lens(tree_ob, params$init_br)
  
  ode_parms <- list(d = params$d, m = params$m, u = length(params$K_parms$hz), 
                    a = params$K_parms$a, h0 = params$K_parms$h0,
                    h_z = params$K_parms$hz,
                    P0_i = params$K_parms$P0i, sigma0_i = params$K_parms$sig0i,
                    P_iz = params$K_parms$Piz, D_i = params$a_parms$D_i,
                    b_iz = params$K_parms$biz, 
                    V_gi = params$V_gi,
                    sigma_iz = params$K_parms$sigiz, gamma_i = params$a_parms$gamma_i)
  
  event_vec <- c(birth = tree_ob$b_rate, check_extinct = params$check_extinct)  
  t <- tree_ob$n_tips_hist[,1][nrow(tree_ob$n_tips_hist)]
  tree_i <- 1L
  i <- 1L
  current_event_vec <- event_vec
  last_print <- 0
  last_tree_save <- 0
  #Ns <- c(0.05, 0.05)
  ## Main Simulation Loop
  
  pr <- progress_bar$new(total = ceiling(params$tot_time), clear = FALSE,
                         format = "Running simulation [:bar] :percent eta: :eta; current species richness: :sr")
  #pr_i <- 0
  print_interval <- params$tot_time / 100
  
  while (t < params$tot_time & sum(tree_ob$extant) > 0) {
    i <- i + 1L
    # generate next event ----------------------------
    #current_event_vec[1] <- event_vec[1]*sum(tree_ob$extant)
    next_event <- rexp(1, sum(current_event_vec))
    next_event_type <- sample(names(current_event_vec), 1, prob = current_event_vec/sum(current_event_vec)) 
    
    ## stuff to do no matter what the event is
    ## setup adaptive dynamics -----------------------------------------------
    
    if(floor(next_event) > 1) {
      ode_parms$m <- sum(tree_ob$extant)
      ode_init <- c(as.vector(tree_ob$traits[ , tree_ob$extant]), tree_ob$Ns[tree_ob$extant])
      #test_fun <- adapt_landscape_comp_dyn_de(1, ode_init, ode_parms)
      #test_fun
      
      ## run ODE
      #print(tree_ob$Ns)
      
      next_ode <- ode(ode_init, 1:floor(next_event), trait_pop_sim_de, ode_parms)
      tree_ob <- update_br_lens(tree_ob, next_event, 0)
      
      #print(tree_ob$Ns)
      
      tree_ob$traits[ , tree_ob$extant] <- matrix(next_ode[nrow(next_ode), -c(1, (ncol(next_ode) - ode_parms$m + 1):ncol(next_ode))], nrow = ode_parms$d)
      tree_ob$Ns[tree_ob$extant] <- next_ode[nrow(next_ode), c((ncol(next_ode) - ode_parms$m + 1):ncol(next_ode))]
      
      #print(tree_ob$Ns)
      
      if(params$trait_hist) {
        if(params$trait_hist_prop < 1) {
          n_sample_rows <- ceiling(nrow(next_ode) * params$trait_hist_prop)
          sample_rows <- seq(1, nrow(next_ode), by = floor(nrow(next_ode) / n_sample_rows))
          next_ode <- next_ode[sample_rows, ]
        }
        tree_ob$full_dat$add(next_ode)
      }
      tree_ob$extant_list$add(which(tree_ob$extant))
      
      #print(tree_ob)
      
      ## birth events --------------------------------------
      if(next_event_type == "birth"){
        #tree_ob <- update_tree_ob_remove_zero_radius(tree_ob, radius_cutoff)
        if(any(tree_ob$Ns < 0.0001 & tree_ob$extant)) {
          extinct <- which(tree_ob$Ns < 0.0001 & tree_ob$extant)
          tree_ob <- update_tree_ob_extinction(tree_ob, extinct)
        }
        split_tip <- sample(tree_ob$tips[tree_ob$extant], 1)
        tree_ob <- update_tree_ob_speciate(tree_ob, split_tip, mult = params$mult)
      }
      
      if(next_event_type == "check_extinct") {
        if(any(tree_ob$Ns < 0.001 & tree_ob$extant)) {
          extinct <- which(tree_ob$Ns < 0.001 & tree_ob$extant)
          tree_ob <- update_tree_ob_extinction(tree_ob, extinct)
        }
      }
      
      t <- t + next_event
      tree_ob$n_tips_hist <- rbind(tree_ob$n_tips_hist, c(t, sum(tree_ob$extant)))
      
      if((t - last_print) > print_interval) {
        pr$tick(t - last_print, tokens = list(sr = sum(tree_ob$extant)))
        last_print <- t
      }
      
      if (params$save_tree) {
        if((t - last_tree_save) > params$save_tree_interval) {
          tree_list[[tree_i]] <- tree_ob$phylo
          save_tree_times[tree_i] <- t
          tree_i <- tree_i + 1L
          last_tree_save <- t
        }
      }
    }
  }
  pr$tick(t - last_print, tokens = list(sr = sum(tree_ob$extant)))
  
  tree_ob$full_dat <- tree_ob$full_dat$as.list()
  tree_ob$extant_list <- tree_ob$extant_list$as.list()
  list(sim_object = tree_ob, saved_trees = list(tree_times = save_tree_times, trees = tree_list), params = parms)
}

#' Continue Niche Filling Simulation Function
#' 
#' Function to continue running the niche filling simulation on a simulation that has already been run,
#' picking up where it left off.
#' 
#' @param sim_ob A simulation object created by a previous run of \code{\link{sim_radiation}}
#' @param time_steps Number of additional time steps to run the simulation.
#' @param arrest_speciation Logical. TRUE if the speciation part of the simulation should be 
#' arrested. This is usually done to allow the simulation to come to some equilibrium state.
#' 
#' @return A niche_fill_sim object containing the final simulation object,
#' a set of intermediate phylogenies, and the parameters used to run the simulation
#' 
#' @importFrom deSolve ode
#' @import progress
#' 
#' @export
continue_sim_radiation <- function(sim_ob, time_steps = 10000, arrest_speciation = FALSE) {
  
  parms <- sim_ob$params
  
  params <- list(K_parms = parms$K_parms, a_parms = parms$a_parms, b_rate = parms$macro_parms$b_rate, 
                 init_traits = parms$macro_parms$init_traits, e_var = parms$macro_parms$e_var, 
                 init_Ns = parms$macro_parms$init_Ns, init_br = parms$macro_parms$init_br, 
                 check_extinct = parms$macro_parms$check_extinct, tot_time = parms$macro_parms$tot_time, 
                 V_gi = parms$macro_parms$V_gi, m = parms$macro_parms$m, d = parms$macro_parms$d, 
                 save_tree = parms$macro_parms$save_tree, save_tree_interval = parms$macro_parms$save_tree_interval, 
                 progress = parms$macro_parms$progress, mult = parms$macro_parms$mult,
                 trait_hist = parms$macro_parms$trait_hist, trait_hist_prop = parms$macro_parms$trait_hist_prop)
  
  
  tree_list <- sim_ob$saved_trees$trees
  save_tree_times <- sim_ob$saved_trees$tree_times
  
  tree_ob <- sim_ob$sim_object
  
  if(params$trait_hist){
    full_dat2 <- expandingList()
  }
  extant_list2 <- expandingList()
  
  ode_parms <- list(d = params$d, m = params$m, u = length(params$K_parms$hz), 
                    a = params$K_parms$a, h0 = params$K_parms$h0,
                    h_z = params$K_parms$hz,
                    P0_i = params$K_parms$P0i, sigma0_i = params$K_parms$sig0i,
                    P_iz = params$K_parms$Piz, D_i = params$a_parms$D_i,
                    b_iz = params$K_parms$biz, 
                    V_gi = params$V_gi,
                    sigma_iz = params$K_parms$sigiz, gamma_i = params$a_parms$gamma_i)
  
  event_vec <- c(birth = tree_ob$b_rate, check_extinct = params$check_extinct)  
  t <- tree_ob$n_tips_hist[,1][nrow(tree_ob$n_tips_hist)]
  tree_i <- length(tree_list) + 1L
  i <- 1L
  current_event_vec <- event_vec
  last_print <- params$tot_time
  last_tree_save <- 0
  #Ns <- c(0.05, 0.05)
  ## Main Simulation Loop
  
  params$tot_time <- params$tot_time + time_steps
  
  pr <- progress_bar$new(total = ceiling(time_steps), clear = FALSE,
                         format = "Running simulation [:bar] :percent eta: :eta; current species richness: :sr")
  #pr_i <- 0
  print_interval <- params$tot_time / 100
  
  while (t < params$tot_time & sum(tree_ob$extant) > 0) {
    i <- i + 1L
    # generate next event ----------------------------
    #current_event_vec[1] <- event_vec[1]*sum(tree_ob$extant)
    next_event <- rexp(1, sum(current_event_vec))
    next_event_type <- sample(names(current_event_vec), 1, prob = current_event_vec/sum(current_event_vec)) 
    
    ## stuff to do no matter what the event is
    ## setup adaptive dynamics -----------------------------------------------
    
    if(floor(next_event) > 1) {
      ode_parms$m <- sum(tree_ob$extant)
      ode_init <- c(as.vector(tree_ob$traits[ , tree_ob$extant]), tree_ob$Ns[tree_ob$extant])
      #test_fun <- adapt_landscape_comp_dyn_de(1, ode_init, ode_parms)
      #test_fun
      
      ## run ODE
      #print(tree_ob$Ns)
      
      next_ode <- ode(ode_init, 1:floor(next_event), trait_pop_sim_de, ode_parms)
      tree_ob <- update_br_lens(tree_ob, next_event, 0)
      
      #print(tree_ob$Ns)
      
      tree_ob$traits[ , tree_ob$extant] <- matrix(next_ode[nrow(next_ode), -c(1, (ncol(next_ode) - ode_parms$m + 1):ncol(next_ode))], nrow = ode_parms$d)
      tree_ob$Ns[tree_ob$extant] <- next_ode[nrow(next_ode), c((ncol(next_ode) - ode_parms$m + 1):ncol(next_ode))]
      
      #print(tree_ob$Ns)
      
      if(params$trait_hist) {
        if(params$trait_hist_prop < 1) {
          n_sample_rows <- ceiling(nrow(next_ode) * params$trait_hist_prop)
          sample_rows <- seq(1, nrow(next_ode), by = floor(nrow(next_ode) / n_sample_rows))
          next_ode <- next_ode[sample_rows, ]
        }
        full_dat2$add(next_ode)
      }
      extant_list2$add(which(tree_ob$extant))
      
      #print(tree_ob)
      
      ## birth events --------------------------------------
      if(next_event_type == "birth" & !arrest_speciation){
        #tree_ob <- update_tree_ob_remove_zero_radius(tree_ob, radius_cutoff)
        if(any(tree_ob$Ns < 0.0001 & tree_ob$extant)) {
          extinct <- which(tree_ob$Ns < 0.0001 & tree_ob$extant)
          tree_ob <- update_tree_ob_extinction(tree_ob, extinct)
        }
        split_tip <- sample(tree_ob$tips[tree_ob$extant], 1)
        tree_ob <- update_tree_ob_speciate(tree_ob, split_tip, mult = params$mult)
      }
      
      if(next_event_type == "check_extinct") {
        if(any(tree_ob$Ns < 0.001 & tree_ob$extant)) {
          extinct <- which(tree_ob$Ns < 0.001 & tree_ob$extant)
          tree_ob <- update_tree_ob_extinction(tree_ob, extinct)
        }
      }
      
      t <- t + next_event
      tree_ob$n_tips_hist <- rbind(tree_ob$n_tips_hist, c(t, sum(tree_ob$extant)))
      
      if((t - last_print) > print_interval) {
        pr$tick(t - last_print, tokens = list(sr = sum(tree_ob$extant)))
        last_print <- t
      }
      
      if (params$save_tree) {
        if((t - last_tree_save) > params$save_tree_interval) {
          tree_list[[tree_i]] <- tree_ob$phylo
          save_tree_times[tree_i] <- t
          tree_i <- tree_i + 1L
          last_tree_save <- t
        }
      }
    }
  }
  pr$tick(t - last_print, tokens = list(sr = sum(tree_ob$extant)))
  
  full_dat2 <- full_dat2$as.list()
  extant_list2 <- extant_list2$as.list()
  
  tree_ob$full_dat <- c(tree_ob$full_dat, full_dat2)
  tree_ob$extant_list <- c(tree_ob$extant_list, extant_list2)
  parms$macro_parms$tot_time <- params$tot_time
  list(sim_object = tree_ob, saved_trees = list(tree_times = save_tree_times, trees = tree_list), params = parms)
}
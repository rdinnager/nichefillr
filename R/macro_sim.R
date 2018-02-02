#' Initialize Tree Object for Macroevolution Simulation
#' 
#'  Initializes an object used by the macroevolutionary simulation to store data
#'  about species, their relationships, and their population sizes and trait values 
#'  
#' @param b_rate Pure birth rate for new species, here just to get first branching;
#' @param init_traits Vector of initial trait values of length n_trait, where n_trait is
#' the number of traits to be simulated
#' @param e_var Vector of gaussian variances determining rate of evolution along branches; of 
#' length n_trait; This is the Brownian motion rate parameter
#' 
#' @return Initialized object with macroevolutionary simulation setup
make_tree_ob <- function(b_rate = 0.001, init_traits, e_var, init_Ns, trait_hist = TRUE){
  tree_ob <- list()
  start_tree <- list()
  start_tree$edge <- matrix(integer(), nrow = 0, ncol = 2)
  start_tree$edge.length <- numeric()
  start_tree$tip.label <- character()
  start_tree$Nnode <- integer()
  class(start_tree) <- "phylo"
  tips <- NA
  n_tips <- 0
  
  first_birth <- rexp(1, b_rate)
  
  start_tree$edge <- rbind(start_tree$edge, matrix(c(3, 3, 1, 2), nrow = 2, ncol = 2))
  start_tree$edge.length <- c(start_tree$edge.length, 0, 0)
  start_tree$tip.label <- c(start_tree$tip.label, paste0("t_", c(1, 2)))
  start_tree$Nnode <- 1
  start_tree$root.edge <- first_birth
  
  tips <- c(1, 2)
  n_tips <- 2
  nodes <- c(3)
  n_nodes <- 1
  
  tree_ob$phylo <- start_tree
  tree_ob$tips <- tips
  tree_ob$n_tips <- n_tips
  tree_ob$nodes <- nodes
  tree_ob$n_nodes <- n_nodes
  tree_ob$b_rate <- b_rate
  
  ## add trait data
  n_traits <- length(init_traits)
  tree_ob$n_traits <- n_traits
  tree_ob$traits <- matrix(init_traits, ncol = 2, nrow = length(init_traits), byrow = FALSE)
  #tree_ob$traits <- tree_ob$traits + rnorm(n_traits, mean = 0,
  #                                         sd = rep(e_var, n_tips))*first_birth
  
  tree_ob$node_traits <- cbind(3, matrix(init_traits, nrow = 1))
  
  ## add extant data
  tree_ob$extant <- c(TRUE, TRUE)
  
  ## save history of species richness
  tree_ob$n_tips_hist <- cbind(0, n_tips)
  tree_ob$e_var <- e_var
  tree_ob$Ns <- init_Ns
  
  if(trait_hist){
    tree_ob$full_dat <- expandingList()
  }
  tree_ob$extant_list <- expandingList()
  
  tree_ob
}

#' Function to Update Branch-lengths (and Traits) in Macroevolutionary Simulation
#' 
#' Updates branch-length by adding to time to every branch; also updates trait values using
#' Brownian motion simulation of the same time period. This is used at the beginning of
#' the simulation to setup the first split of two species, and to grow branches within
#' the main loop (with no trait evolution).
#' 
#' @param tree_ob The simulation object to update
#' @param add_br_len How much branch-length should be added (in units of simulation time)
#' @param e_var Vector of evolutionary rates of length n_trait, where n_trait is the number
#' of traits being simulated
#' 
#' @return An updated simulation object
update_br_lens <- function(tree_ob, add_br_len, e_var = NULL) {
  
  if(is.null(e_var)) {
    e_var <- tree_ob$e_var
  }
  
  ## find tips that are extant
  extant_tips <- 1:tree_ob$n_tips %in% which(tree_ob$extant)
  extant_tips <- tree_ob$phylo$edge[ , 2] %in% which(extant_tips)
  tree_ob$phylo$edge.length[extant_tips] <- tree_ob$phylo$edge.length[extant_tips] + add_br_len
  
  tree_ob$traits[ , tree_ob$extant] <- tree_ob$traits[ , tree_ob$extant] + rnorm(sum(tree_ob$extant)*tree_ob$n_traits, mean = 0,
                                                                                 sd = rep(e_var, sum(tree_ob$extant)))*add_br_len
  
  tree_ob
}

#' Update Simulation Object with Speciation Event 
#' 
#' Function to update simulation object by adding add_br_len to the tips, 
#' then splitting a tip into two species
#' 
#' @param tree_ob The simulation object to update
#' @param split_tip Integer expressing which tip number to split into two sister species
#' @param mult Multiplier which determines the proportional increase in trait evolution
#' rate at a speciation event. This causes species to diverge more rapidly than usual if
#' greater than one. Generally used to make sure a speciation is accompanied by a 'jump'
#' in trait values, which helps them not immediately go extinct from competition
#' 
#' @return An updated simulation object
update_tree_ob_speciate <- function(tree_ob, split_tip, mult = 1) {
  
  #tree_ob <- update_br_lens(tree_ob, add_br_len, e_var)
  
  start_tree <- tree_ob$phylo
  tips <- tree_ob$tips
  n_tips <- tree_ob$n_tips
  nodes <- tree_ob$nodes
  n_nodes <- tree_ob$n_nodes
  
  ## chose tip to split
  curr_tip <- split_tip
  next_tip <- max(tips) + 1L
  
  start_tree$edge[start_tree$edge > n_tips] <- start_tree$edge[start_tree$edge > n_tips] + 1L
  tree_ob$node_traits[ , 1] <- tree_ob$node_traits[ , 1] + 1L
  next_node <- max(start_tree$edge[ , 1]) + 1L
  
  start_tree$edge[start_tree$edge == curr_tip] <- next_node
  start_tree$edge <- rbind(start_tree$edge, matrix(c(next_node, next_node, curr_tip, next_tip), ncol = 2))
  start_tree$edge.length <- c(start_tree$edge.length, 0, 0)
  #start_tree$edge[ , 1] <- start_tree$edge[ , 1] - min(start_tree$edge[ , 1]) + n_tips + 1
  nodes <- unique(start_tree$edge[ , 1])
  n_nodes <- length(nodes)
  start_tree$tip.label <- c(start_tree$tip.label, paste0("t_", next_tip))
  tree_ob$extant <- c(tree_ob$extant, TRUE)
  start_tree$Nnode <- n_nodes
  
  ##update traits
  
  tree_ob$traits <- cbind(tree_ob$traits, matrix(NA, ncol = 1, nrow = tree_ob$n_traits))
  tree_ob$traits[ , next_tip] <- tree_ob$traits[ , curr_tip] + 
    rnorm(tree_ob$n_traits, mean = 0, sd = tree_ob$e_var*mult)
  
  tree_ob$node_traits <- rbind(tree_ob$node_traits, c(next_node, tree_ob$traits[ , curr_tip]))
  tree_ob$Ns <- c(tree_ob$Ns, 0.005)
  
  tips <- c(tips, next_tip)
  n_tips <- length(tips)
  
  tree_ob$phylo <- start_tree
  tree_ob$tips <- tips
  tree_ob$n_tips <- n_tips
  tree_ob$nodes <- nodes
  tree_ob$n_nodes <- n_nodes
  #tree_ob$b_rate <- b_rate
  # tree_ob$Ns[curr_tip] <- tree_ob$Ns[curr_tip]*0.75
  
  tree_ob
}

#' Update Simulation Object with Extinction Event 
#' 
#' Function to update simulation object by removing a species from the simulation
#' 
#' @param tree_ob The simulation object to update
#' @param extinct_tip Integer expressing which tip number to go extinct (stop updating)
#' 
#' @return An updated simulation object
update_tree_ob_extinction <- function(tree_ob, extinct_tip) {
  # tree_ob <- update_br_lens(tree_ob, add_br_len, e_var)
  tree_ob$extant[extinct_tip] <- FALSE
  tree_ob$Ns[extinct_tip] <- 0
  tree_ob
}
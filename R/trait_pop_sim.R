#' Calculate Niche Filling Trait and Population Simulation Derivative
#' 
#' Calculates the derivative of all traits in the niche-filling simulation
#' For use with \code{deSolve} functions to run the trait evolution simulation. User probably 
#' won't have reason to use this directly 
#' 
#' @param t Time step of the trait evolution simulation; this is for compatability with 
#' \code{deSolve} and is ignored
#' @param y Current value of the traits and populations
#' @param parms Named list containing all necessary parameters
#' 
#' @return List of derivatives compatible with \code{deSolve}
trait_pop_sim_de <- function(t, y, parms) {
  list(pkg.env$adapt_landscape_comp_dyn_cmp(d = parms$d, m = parms$m, u = parms$u, 
                                    a = parms$a, h0 = parms$h0,
                                    h_z = parms$h_z,
                                    P0_i = parms$P0_i, sigma0_i = parms$sigma0_i,
                                    P_iz = parms$P_iz, D_i = parms$D_i,
                                    b_iz = parms$b_iz, state = y, 
                                    V_gi = parms$V_gi,
                                    sigma_iz = parms$sigma_iz, gamma_i = parms$gamma_i,
                                    c_r = parms$c_r, C = parms$C))
}

trait_pop_sim_de_ellip <- function(t, y, parms) {
  list(pkg.env$adapt_landscape_comp_dyn_cmp_ellip(d = parms$d, m = parms$m, u = parms$u, 
                                            a = parms$a, h0 = parms$h0,
                                            h_z = parms$h_z,
                                            P0 = parms$P0, sigma0_i = parms$sigma0_i,
                                            P_z = parms$P_z, D0 = parms$D0,
                                            b_iz = parms$b_iz, state = y, 
                                            V_gi = parms$V_gi,
                                            sigma_iz = parms$sigma_iz, gamma_i = parms$gamma_i,
                                            c_r = parms$c_r, C = parms$C))
}

temp_fitness_interface <- function(...) {
  pkg.env$adapt_landscape_fitness_cmp(...)
}

#' Carrying Capacity Landscape Function (R version)
#' 
#' Carrying capacity function written in R, mostly for visualization or testing purposes.
#' Returns the value of the carrying capacity landscape at coordinate \code{xi}.
#' 
#' @param xi X values to evaluate under carrying capacity function; a vector of length d, 
#' where d is the number of dimensions of the landscape
#' @param h0 Total maximum height of the carrying capacity landscape
#' @param sig0i Total width of landscape in all dimensions; determines how far from zero the carrying capacity
#' drops off to nearly zero; vector of length d
#' @param P0i Total super-gaussian parameter; determines how quickly or gradually the carrying
#' capacity drops off near the landscape borders in each dimension. Higher values give more extreme drop-offs;
#' vector of length d
#' @param hz Maximum height of each of u peaks in the landscape; vector of length u
#' @param biz Centre of each of u peaks in the landscape for each of d dimensions; 
#' matrix of dimension d by u
#' @param sigiz Width of each of u peaks for each of d dimensions; matrix of length d by u
#' @param Piz Super-gaussian parameter for each of u peaks for each of d dimensions; 
#' matrix of dimension d by u
#' @param a Minimum values of carrying capacity within landscape limits
#' 
#' @return Value of carrying capacity at coordinates \code{xi} 
#' 
#' @export
K_func <- function(xi, h0 = 1, sig0i = c(5, 5), P0i = c(1.5, 1.5),
                   hz = c(1, 1), biz = matrix(c(-1, -1, 1, 1), nrow = 2, ncol = 2),
                   sigiz = matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2),
                   Piz = matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2), a = 0.1,
                   use_ellip = TRUE) {
  
  
  if(use_ellip) {
    P0i <- P0i[1]
    if(length(dim(Piz)) > 1) {
      Piz <- Piz[1, , drop = TRUE]
    }
    term_1 <- h0*exp(-sum((xi^2/(2*sig0i^2)))^P0i)
    term_2 <- sum(hz*exp(-colSums(((xi - biz)^2/(2*sigiz^2)))^Piz)) + a
  } else {
    term_1 <- h0*exp(-sum((xi^2/(2*sig0i^2))^P0i))
    term_2 <- sum(hz*exp(-colSums(((xi - biz)^2/(2*sigiz^2))^Piz))) + a
  }

  term_1*term_2  
}

#' @export
fit_func <- function(x, K_parms, a_parms, spec_traits, spec_Ns, res_x, d = 2) {
  
  spec_traits <- cbind(spec_traits, res_x)
  spec_Ns <- c(spec_Ns, K_func(res_x, h0 = K_parms$h0, 
                               sig0i = K_parms$sig0i, 
                               P0i = K_parms$P0i,
                               hz = K_parms$hz, 
                               biz = K_parms$biz,
                               sigiz = K_parms$sigiz,
                               Piz = K_parms$Piz, 
                               a = K_parms$a))
  
  xmat <- apply(spec_traits, 2, function(y) (((y - x)^2) / (2 * a_parms$gamma_i))^a_parms$D_i[1]) 
  comp <- exp(-apply(xmat, 2, sum)) * spec_Ns
  
  fitness <- 1 - ((sum(comp)) / K_func(x, h0 = K_parms$h0, 
                         sig0i = K_parms$sig0i, 
                         P0i = K_parms$P0i,
                         hz = K_parms$hz, 
                         biz = K_parms$biz,
                         sigiz = K_parms$sigiz,
                         Piz = K_parms$Piz, 
                         a = K_parms$a))
  
  fitness
}

#' Generate Random Carrying Capacity Landscape
#' 
#' Function to generate a random carrying capacity landscape as a mix of super-gaussians
#' 
#' @param u Number of carrying capacity peaks to generate
#' @param h_mean The mean height of the carrying capacity peaks
#' @param h_var The variance in height of the peaks
#' @param bi_mean The mean peak centre value for each of d dimensions; vector of length d
#' @param bi_var The variance in peak centre values for each of d dimensions; vector of length d
#' @param h0 The maximum height of the carrying capacity on the landscape
#' @param sig0i The width of the carrying capacity landscape, e.g. how far it extends before dropping off
#' in each of d dimensions; vector of length d
#' @param P0i Total super-gaussian parameter; determines how quickly or gradually the carrying
#' capacity drops off near the landscape borders in each dimension. Higher values give more extreme drop-offs;
#' vector of length d
#' @param sigi_mean The mean width of peaks for each of d dimensions; vector of length d
#' @param sigi_var The variance in peak widths for each of d dimensions; vector of length d
#' @param Pi_mean The mean super-gaussian parameter of peaks for each of d dimensions; determines
#' the sharpness of peak drop-off; vector of length d
#' @param a Minimum values of carrying capacity within landscape limits 
#' @param hz_min Minimum height allowed for any peak
#' 
#' @return List of parameter values that can be passed on to other functions such as 
#' \code{\link{K_func}}
#' 
#' @export
generate_landscape_old <- function(u = 2, h_mean = 1, h_var = 0.5, bi_mean = c(0, 0), bi_var = c(3, 3), 
                               h0 = 1, sig0i = c(5, 5), P0i = c(1.5, 1.5), sigi_mean = c(1, 1), 
                               sigi_var = c(1, 1), Pi_mean = c(1, 1), Pi_var = c(0.25, 0.25), 
                               a = 0.5, Piz_min = 0.25, hz_min = 0.25) {
  
  hz <- rnorm(u, h_mean, h_var)
  hz[hz <= hz_min] <- hz_min
  biz <- matrix(rnorm(length(bi_mean)*u, bi_mean, bi_var), nrow = length(bi_mean), ncol = u)
  sigiz <- matrix(rlnorm(length(sigi_mean)*u, sigi_mean, sigi_var), nrow = length(sigi_mean), ncol = u)
  Piz <- matrix(rnorm(length(Pi_mean)*u, Pi_mean, Pi_var), nrow = length(Pi_mean), ncol = u)
  Piz[Piz <= Piz_min] <- Piz_min
  list(h0 = h0, hz = hz, biz = biz, sigiz = sigiz, Piz = Piz, sig0i = sig0i, P0i = P0i, a = a)
  
}

diffeqr_selection_interface <- function(y, parms2, t) {
  scalar_len <- 8
  peaks <- as.integer(parms2[3])
  dims <- as.integer(parms2[1])
  peaks_len <- peaks + peaks + dims*peaks + dims*peaks
  dims_len <- dims + dims + dims
  specs <- as.integer(parms2[2])
  specs_len <- specs
  parms <- list()
  parms$d <- as.integer(parms2[1])
  parms$m <- as.integer(parms2[2])
  parms$u <- as.integer(parms2[3])
  parms$a <- parms2[4]
  parms$h0 <- parms2[5]
  parms$P0 <- parms2[6]
  parms$D0 <- parms2[7] 
  parms$C <- parms2[8]
  spot <- 9
  spot2 <- spot + peaks - 1
  parms$h_z <- parms2[spot:spot2]
  spot <- spot2 + 1
  spot2 <- spot2 + peaks
  parms$P_z <- parms2[spot:spot2] 
  spot <- spot2 + 1
  spot2 <- spot2 + peaks*dims
  parms$b_iz <- matrix(parms2[spot:spot2], nrow = dims, ncol = peaks)
  spot <- spot2 + 1
  spot2 <- spot2 + peaks*dims
  parms$sigma_iz <- matrix(parms2[spot:spot2], nrow = dims, ncol = peaks)
  spot <- spot2 + 1
  spot2 <- spot2 + dims
  parms$sigma0_i <- parms2[spot:spot2]
  spot <- spot2 + 1
  spot2 <- spot2 + dims
  parms$V_gi <- parms2[spot:spot2]
  spot <- spot2 + 1
  spot2 <- spot2 + dims
  parms$gamma_i <- parms2[spot:spot2]
  # spot <- spot2 + 1
  # spot2 <- spot2 + specs
  # parms$c_r <- parms2[spot:spot2]
  #print(parms)
  
  pkg.env$adapt_landscape_comp_dyn_cmp_ellip(d = parms$d, m = parms$m, u = parms$u, 
                                       a = parms$a, h0 = parms$h0,
                                       h_z = parms$h_z,
                                       P0 = parms$P0, sigma0_i = parms$sigma0_i,
                                       P_z = parms$P_z, D0 = parms$D0,
                                       b_iz = parms$b_iz, state = y, 
                                       V_gi = parms$V_gi,
                                       sigma_iz = parms$sigma_iz, gamma_i = parms$gamma_i,
                                       C = parms$C)
  
  #res <- nichefillr:::trait_pop_sim_de_ellip(t, y, parms)
  
}

diffeqr_drift_interface <- function(y, parms2, t) {
  scalar_len <- 8
  peaks <- as.integer(parms2[3])
  dims <- as.integer(parms2[1])
  peaks_len <- peaks + peaks + dims*peaks + dims*peaks
  dims_len <- dims + dims + dims
  specs <- as.integer(parms2[2])
  specs_len <- specs
  parms <- list()
  parms$d <- as.integer(parms2[1])
  parms$m <- as.integer(parms2[2])
  
  spot <- 9
  spot2 <- spot + peaks - 1
  
  spot <- spot2 + 1
  spot2 <- spot2 + peaks
  
  spot <- spot2 + 1
  spot2 <- spot2 + peaks*dims
  
  spot <- spot2 + 1
  spot2 <- spot2 + peaks*dims
  
  spot <- spot2 + 1
  spot2 <- spot2 + dims
  
  spot <- spot2 + 1
  spot2 <- spot2 + dims
  parms$V_gi <- parms2[spot:spot2]
  
  parms$demo <- 0.1
  pkg.env$drift_cmp(y, parms$V_gi, parms$d, parms$m, parms$demo)
  
}

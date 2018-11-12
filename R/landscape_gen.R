#' Generate a Random Carrying Capacity Landscape
#' 
#' Function to generate a random carrying capacity landscape obeying certain constraints. 
#' It distributes peaks uniformly within an N-ball, which is approximated by the overall
#' Super-gaussian multiplier distribution. Peak heights and widths are chosen to
#' satisfy a particular ratio and to add up to the total desired niche 'volume'.
#' 
#' @param potent_vol Total potential volume of the entire landscape.
#' @param total_vol Total approximate volume of the entire landscape.
#' @param num_peaks Number of fitness peaks to place on the landscape.
#' @param h_to_sig_ratio Ratio of the height of each peak to its width (or sigma value).
#' @param P_min_max Vector of length 2 giving the minimum and maximum peak super-gaussian
#' parameters to use.
#' @param dirichlet_param The alpha parameter of the dircihlet distribution for 
#' determining how the volume should be split between the peaks.
#' @param d The number of dimensions of the landscape
#' @param P Super-gaussian parameter for whole landscape multiplier
#' @param a_prop Minimum carrying capacity within the landscape as a proportion of the maximum
#' height
#' 
#' @return List of parameter values that can be used in the carrying capacity function
#' to generate the landscape
#' 
#' @importFrom nimble rdirch
#' 
#' @export
generate_landscape <- function(potent_vol = 1, total_vol = 1, num_peaks = 6, 
                               h_to_sig_ratio = 2, 
                               P_min_max = c(0.9, 1.1), dirichlet_param = 1.5,
                               d = 2, P = 8, a = 0.01) {
  norms <- matrix(rnorm(num_peaks*d), ncol = d)
  #unifs <- runif(num_peaks)
  rexps <- rexp(num_peaks)
  
  ## use the equation for the volume of an N-ball to get radius equal to desired
  ## total volume
  total_rad <- (potent_vol / ((pi^(d/2))/(gamma((d/2)+1))))^(1/d)
  
  ## uniform radius equals unif(0, 1)*radius^(1/3)
  ## uniform within a unit sphere: x1, x2, x3 are three normal deviates (norm(0,1))
  ## (1/sum(xi^2))*c(x1,x2,x3)
  # sphere_points <- ((total_rad*(unifs^(1/3))) * norms) / 
  #   apply(norms^2, 1, function(x) sqrt(sum(x)))
  
  sphere_points <- (total_rad * norms) / 
    sqrt(rexps + apply(norms^2, 1, function(x) sum(x)))
  
  #plot(sphere_points)
  
  res <- list()
  
  vol_split <- rdirch(1, rep(dirichlet_param, num_peaks))
  vols <- total_vol * vol_split
  sigmas <- ((vols / (h_to_sig_ratio*((2 * pi)^((d/2))))))^(1/(d+1))

  hs <- sigmas * h_to_sig_ratio
  ## use the full width at half maximum formula to guestimate sigma for total landscape super-gaussian
  total_sig <- (2 * total_rad) / (2*((2*log(2))^(1/(2*P)))*sqrt(2))
  
  res$h0 <- 1
  res$hz <- hs
  res$biz <- t(sphere_points)
  res$sigiz <- do.call(rbind, replicate(d, sigmas, simplify = FALSE))
  res$Piz <- matrix(runif(d*num_peaks, P_min_max[1], P_min_max[2]), nrow = d)
  res$sig0i <- rep(total_sig, d)
  res$P0i <- rep(P, d)
  res$a <- 0
  
  integral <- integrate_fun(K_func, lower = rep(-total_rad*1, d), upper = rep(total_rad*1, d),
                            h0 = res$h0, sig0i = res$sig0i, P0i = res$P0i,
                            hz = res$hz, biz = res$biz,
                            sigiz = res$sigiz,
                            Piz = res$Piz, a = res$a)
  
  
  miss_vol <- total_vol - integral
  res$a <- miss_vol / potent_vol
  
  #res$h0 <- res$h0 * (total_vol / integral)
  
  res$total_rad <- total_rad
  res$potent_vol <- potent_vol
  res$total_vol <- total_vol
  
  res$num_peaks = num_peaks 
  res$h_to_sig_ratio = h_to_sig_ratio 
  res$P_min_max = P_min_max
  res$dirichlet_param = dirichlet_param
  
  res
  
}

#' Generate a Random Carrying Capacity Landscape
#' 
#' Function to generate a random carrying capacity landscape obeying certain constraints. 
#' It distributes peaks uniformly within an N-ball, which is approximated by the overall
#' Super-gaussian multiplier distribution. Peak heights and widths are chosen to
#' satisfy a particular ratio and to add up to the total desired niche 'volume'.
#' 
#' @param potent_vol Total potential volume of the entire landscape.
#' @param total_vol Total approximate volume of the entire landscape.
#' @param num_peaks Number of fitness peaks to place on the landscape.
#' @param h_to_sig_ratio Ratio of the height of each peak to its width (or sigma value).
#' @param P_min_max Vector of length 2 giving the minimum and maximum peak super-gaussian
#' parameters to use.
#' @param dirichlet_param The alpha parameter of the dircihlet distribution for 
#' determining how the volume should be split between the peaks.
#' @param d The number of dimensions of the landscape
#' @param P Super-gaussian parameter for whole landscape multiplier
#' @param a Minimum carrying capacity within the landscape
#' 
#' @return List of parameter values that can be used in the carrying capacity function
#' to generate the landscape
#' 
#' @importFrom nimble rdirch
#' 
#' @export
generate_landscape_simple <- function(potent_vol = 1, total_vol = 1, num_peaks = 6, 
                               h_to_sig_ratio = 2, 
                               P_min_max = c(0.9, 1.1), dirichlet_param = 1.5,
                               d = 2, P = 8, a = 0.01) {
  norms <- matrix(rnorm(num_peaks*d), ncol = d)
  unifs <- runif(num_peaks)
  exps <- rexp(num_peaks)
  
  ## use the equation for the volume of an N-ball to get radius equal to desired
  ## total volume
  total_rad <- (potent_vol / ((pi^(d/2))/(gamma((d/2)+1))))^(1/d)
  
  ## uniform radius equals unif(0, 1)*radius^(1/3)
  ## uniform within a unit sphere: x1, x2, x3 are three normal deviates (norm(0,1))
  ## (1/sum(xi^2))*c(x1,x2,x3)
  # sphere_points <- ((total_rad*(unifs^(1/3))) * norms) / 
  #   apply(norms^2, 1, function(x) sqrt(sum(x)))
  
  sphere_points <- (total_rad*norms) /
    (exps + apply(norms^2, 1, function(x) sqrt(sum(x))))
  
  #plot(sphere_points)
  
  res <- list()
  
  vol_split <- rdirch(1, rep(dirichlet_param, num_peaks))
  vols <- total_vol * vol_split
  sigmas <- ((vols / (h_to_sig_ratio*((2 * pi)^((d/2))))))^(1/(d+1))
  
  hs <- sigmas * h_to_sig_ratio
  ## use the full width at half maximum formula to guestimate sigma for total landscape super-gaussian
  total_sig <- (2 * total_rad) / (2*((2*log(2))^(1/(2*P)))*sqrt(2))
  
  res$h0 <- 1
  res$hz <- hs
  res$biz <- t(sphere_points)
  res$sigiz <- do.call(rbind, replicate(d, sigmas, simplify = FALSE))
  res$Piz <- matrix(runif(d*num_peaks, P_min_max[1], P_min_max[2]), nrow = d)
  res$sig0i <- rep(total_sig, d)
  res$P0i <- rep(P, d)
  res$a <- 0
  
  integral <- integrate_fun(K_func, lower = rep(-total_rad*1, d), upper = rep(total_rad*1, d),
                            h0 = res$h0, sig0i = res$sig0i, P0i = res$P0i,
                            hz = res$hz, biz = res$biz,
                            sigiz = res$sigiz,
                            Piz = res$Piz, a = res$a)
  
  
  miss_vol <- total_vol - integral
  res$a <- miss_vol / potent_vol
  
  #res$h0 <- res$h0 * (total_vol / integral)
  
  res$total_rad <- total_rad
  res$potent_vol <- potent_vol
  res$total_vol <- total_vol
  
  res$num_peaks = num_peaks 
  res$h_to_sig_ratio = h_to_sig_ratio 
  res$P_min_max = P_min_max
  res$dirichlet_param = dirichlet_param
  
  res
  
}

#' Visualize a Carrying Capacity Landscape using RGL
#' 
#' This function uses the R package \code{rgl} to visualize a carrying capacity lanscape 
#' function. Currently only works for two-dimensional landscapes.
#' 
#' @param K_parms A names list of parameters describing the landscape; currently only
#' excepts two dimensional landscapes
#' 
#' @return None
#' 
#' @import rgl
#' 
#' @export
vis_K_landscape <- function(K_parms, total_rad) {
  
  res <- K_parms
  
  if(dim(res$biz)[1] != 2) {
    stop("error: vis_K_landscape can only visualize two dimensional landscapes")
  }
  
  x1 <- seq(-total_rad*1.5, total_rad*1.5, 0.1) 
  y1 <- x1
  z1 <- expand.grid(x1, y1)
  
  z1_mat <- matrix(apply(z1, 1, function(x) K_func(x, h0 = res$h0, sig0i = res$sig0i, P0i = res$P0i,
                                                   hz = res$hz, biz = res$biz,
                                                   sigiz = res$sigiz,
                                                   Piz = res$Piz, a = res$a)))
  
  open3d()
  surface3d(y1, x1, z1_mat, color = "green", back = "lines")
  axes3d()
}

#' Generate Dirichlet Distribution Parameters for Peaks
#' 
#'  Function to generate the alpha parameters for a dirichlet distribution, which
#'  will provide a particular proportion of the maximum variance for different numbers
#'  of peaks. Uses a simple numerical optimization to find alpha values generating a 
#'  particular variance, as calculated from the formula for the variance in a dirichlet.
#'  
#'  @param prop_variance Vector of proportions of the maximum variance to calculate.
#'  @param num_peaks Vector of number of carrying capacity peaks to use
#'  
#'  @return Data.frame containing the dirichlet parameter in the last column
gen_dirichlet <- function(prop_variance = c(0.1, 0.3, 0.5), num_peaks = c(2, 6, 10)) {
  
  to_min <- function(x, CV, N) {
    abs(x*(((x*N) - x) / (((x*N)^2)*((x*N) + 1))) - CV) ## The dirichlet variance - CV
  }
  ## calculate maximum variances
  max_var <- data_frame(num_peaks = num_peaks,
                        prop_variance = prop_variance,
                        max_var = to_min(0.00000001, 0, num_peaks))
  target_vars <- max_var %>%
    mutate(target_var = prop_variance*max_var) %>%
    group_by(num_peaks, prop_variance) %>%
    do(data_frame(a = optim(1, to_min, method = "L-BFGS-B", lower = 0.00000001, upper = 1000,
             CV = .$target_var[[1]], N = .$num_peaks[[1]])$par)) %>%
    left_join(max_var, .)
  
  target_vars$a
}

#' Numerically Integrate Carrying Capacity Landscape
#' 
#' Integrate a carrying capacity landscape to estimate its 'volume'
#' 
#' @param K_func R function specifying the carrying capacity function.
#' @param lower Lower bound to integrate across; vector of length d, where d is the number of
#' dimensions
#' @param upper Upper bound to integrate across; vector of length d.
#' @param ... Named arguments to be passed to fun.
#' 
#' @return The total volume of the landscape within the lower and upper bounds.
#' 
#' @importFrom cubature adaptIntegrate
integrate_fun <- function(K_func, lower, upper, ...) {
  # lower <- c(-10, -10, -10)
  # upper <- c(10, 10, 10)
  # test <- adaptIntegrate(K_func, lowerLimit = lower, upperLimit = upper,
  #                        h0 = res$h0, sig0i = 100, P0i = res$P0i,
  #                        hz = res$hz[1, drop = FALSE], biz = res$biz[ , 1, drop = FALSE], 
  #                        sigiz = res$sigiz[ , 1, drop = FALSE],
  #                        Piz = matrix(1, nrow = 3, ncol = 1), a = 0)
  # 
  # lower <- c(-10, -10)
  # upper <- c(10, 10)
  # test <- adaptIntegrate(K_func, lowerLimit = lower, upperLimit = upper,
  #                        h0 = res$h0, sig0i = 100, P0i = res$P0i,
  #                        hz = res$hz[1, drop = FALSE], biz = res$biz[ , 1, drop = FALSE], 
  #                        sigiz = res$sigiz[ , 1, drop = FALSE],
  #                        Piz = matrix(1, nrow = 2, ncol = 1), a = 0)
  # 
  # lower <- c(-10, -10, -10)
  # upper <- c(10, 10, 10)
  # test <- adaptIntegrate(K_func, lowerLimit = lower, upperLimit = upper,
  #                        h0 = res$h0, sig0i = res$sig0i, P0i = c(20,20,20),
  #                        hz = res$hz[1, drop = FALSE], biz = matrix(0, nrow=3, ncol=1), 
  #                        sigiz = res$sigiz[ , 1, drop = FALSE],
  #                        Piz = matrix(1, nrow = 3, ncol = 1), a = 0)
  
  test <- adaptIntegrate(K_func, lowerLimit = lower, upperLimit = upper,
                         ...)
  
  test$integral
  
  # (xi, h0 = 1, sig0i = c(5, 5), P0i = c(1.5, 1.5),
  #          hz = c(1, 1), biz = matrix(c(-1, -1, 1, 1), nrow = 2, ncol = 2),
  #          sigiz = matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2),
  #          Piz = matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2), a = 0.1)
}
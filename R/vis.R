#' Function to extract trait history from simulation object as a data_frame
#' @param sim_ob nichefillr_sim object containing the results of a simulation.
#' @import dplyr
#' @import tidyr
#' @export make_df_from_sim
make_df_from_sim <- function(sim_ob) {
  
  full_dat <- sim_ob$sim_object$full_dat$as.list()
  extant_list <- sim_ob$sim_object$extant_list$as.list()
  #tree_ob$tree_list <- tree_ob$tree_list$as.list()
  
  ms <- sapply(extant_list[!sapply(extant_list, is.null)], length)
  d <- length(sim_ob$sim_params$K_parms$sig0i)
  
  # trait_mat <- sim_ob$sim_object$full_dat[[1]]
  trait_mat_as_df <- function(trait_mat, d, m) {
    if(is.null(dim(trait_mat))) {
      trait_mat <- rbind(NULL, trait_mat)
    }
    cols <- c("Time", paste0("spectrait_", paste(gl(m, d), rep(1:d, m), sep = "_")), paste0("pop_", 1:m))
    colnames(trait_mat) <- cols
    trait_mat <- as.data.frame(trait_mat)
    traits <- trait_mat %>%
      select(Time, starts_with("spectrait_")) %>%
      gather(spectrait, val, -Time) %>%
      separate(spectrait, c("junk", "spec", "trait"), "_") %>%
      select(-junk) %>%
      spread("trait", "val") %>%
      transform(spec = paste0("Species_", spec))
    colnames(traits) <- c("Time", "Species", paste0("Niche_Axis_", 1:d))
    pops <- trait_mat %>%
      select(Time, starts_with("pop_")) %>%
      gather(Species, Population, -Time) %>%
      transform(Species = gsub("pop_", "Species_", Species))
    
    new_df <- traits %>%
      left_join(pops, by = c("Time", "Species"))
  }
  
  full_df <- mapply(trait_mat_as_df, full_dat[!sapply(full_dat, is.null)], 
                    d, ms, SIMPLIFY = FALSE) 
  
  for(i in 2:length(full_df)) {
    full_df[[i]]$Time <- full_df[[i]]$Time + max(full_df[[i - 1]]$Time)
  }
  
  full_df <- bind_rows(full_df)
  full_df
}


#' Function to animate the results of a simulation. Currently only supports simulations
#' done in two dimensions. Will eventually support higher dimensions using t-sne dimension
#' reduction.
#' @param sim_ob nichefillr_sim object containing the results of a simulation.
#' @param file_name (optional) Path to to save the animation to. If NULL, and
#' view = TRUE, animation will be output to the rstudio viewer.
#' @param expand_factor Factor by which to expand x and y axes beyond the minimum and
#' maximum found in the simulated species (allows more of the fitness landscape context to ve
#' shown)
#' @param contour_res Resolution for the fitness contours, in terms of how many points per
#' axis to calculate fitness for.
#' @import gganimate 
#' @import tidyr 
#' @import ggplot2
#' @export sim_animation
sim_animation <- function(sim_ob, file_name = NULL, view = FALSE, expand_factor = 0.1, contour_res = 25, height = 300, width = 400, res = 72) {
  message("Extracting trait history data from simulation (this might take awhile).... ")
  plot_df <- make_df_from_sim(sim_ob)
  
  x_lims <- c(min(plot_df$Niche_Axis_1), max(plot_df$Niche_Axis_1))
  y_lims <- c(min(plot_df$Niche_Axis_2), max(plot_df$Niche_Axis_2))
  
  x_range <- x_lims[2] - x_lims[1]
  y_range <- y_lims[2] - y_lims[1]
  
  x_lims <- x_lims + c(-x_range*expand_factor, x_range*expand_factor) 
  y_lims <- y_lims + c(-y_range*expand_factor, y_range*expand_factor) 
  
  pop_lims <- c(min(plot_df$Population), max(plot_df$Population))
  
  contour_df <- crossing(Niche_Axis_1 = seq(x_lims[1], x_lims[2], length.out = contour_res),
                         Niche_Axis_2 = seq(y_lims[1], y_lims[2], length.out = contour_res)) 
  z <- apply(contour_df %>% as.matrix, 1, 
                     function(x) K_func(x, h0 = sim_ob$sim_params$K_parms$h0, 
                                        sig0i = sim_ob$sim_params$K_parms$sigma0i, 
                                        P0i = sim_ob$sim_params$K_parms$P0i,
                                        hz = sim_ob$sim_params$K_parms$hz, 
                                        biz = sim_ob$sim_params$K_parms$biz,
                                        sigiz = sim_ob$sim_params$K_parms$sigiz,
                                        Piz = sim_ob$sim_params$K_parms$Piz, 
                                        a = sim_ob$sim_params$K_parms$a))
  
  contour_df <- contour_df %>%
    mutate(K = z)
  
  #single_df <- plot_df %>% filter(Time == 1)
  draw_single_plot <- function(single_df, x_lims, y_lims, pop_lims) {
    p <- ggplot(single_df, aes(Niche_Axis_1, Niche_Axis_2)) +
      geom_point(aes(size = Population), shape = 21, fill = NA, stroke = 1.1) +
      ylim(y_lims) +
      xlim(x_lims) +
      scale_size_continuous(limits = pop_lims) +
      ggtitle(paste("Time =", single_df$Time[1])) + 
      geom_contour(aes(z = K), data = contour_df, colour = "grey") +
      theme_minimal()
    print(p)
    return(NULL)
  }
  
  message("Drawing animation frames....")
  
  anim <- image_graph(width = width, height = height, res = 72)
  plot_df %>%
    group_by(Time) %>%
    do(null = draw_single_plot(., y_lims = y_lims, x_lims = x_lims, pop_lims = pop_lims))
  dev.off()
  
  message("\nConverting to animation....")
  animation <- image_animate(anim, fps = 20)
  
  if(view) {
    message("Sending animation to viewer....")
    print(animation)
  }
  
  if(!is.null(file_name)) {
    message("Writing animated gif to disk....")
    image_write(animation, file_name)
  }
  
  return(animation)
  
}

#' Plot a \code{nichefillr_sim} object. If the simulation object contains a trait history,
#' this function will by default plot a 'pollock plot'. If the object contains no trait history it
#' will plot the trait values of the extant species, extinct species (optionally), 
#' and ancestors (optionally). Currently only work for 2 dimensional simulations, eventually
#' will work with higher dimensions through t-sne dimension reduction.
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @export plot.nichfillr_sim
plot.nichfillr_sim <- function(x, fitness_contour = TRUE, contour_res = 100, expand_factor = 0.1, bins = 10) {
  message("Extracting trait history data from simulation (this might take awhile).... ")
  plot_df <- make_df_from_sim(x)
  extant_df <- t(x$sim_object$traits[ , x$sim_object$extant]) %>%
    as.data.frame() %>%
    rename(Niche_Axis_1 = V1, Niche_Axis_2 = V2) %>%
    mutate(Species = paste0("Species_", 1:sum(x$sim_object$extant)),
           Population = x$sim_object$Ns[x$sim_object$extant]) 
  
  x_lims <- c(min(plot_df$Niche_Axis_1), max(plot_df$Niche_Axis_1))
  y_lims <- c(min(plot_df$Niche_Axis_2), max(plot_df$Niche_Axis_2))
  
  x_range <- x_lims[2] - x_lims[1]
  y_range <- y_lims[2] - y_lims[1]
  
  x_lims <- x_lims + c(-x_range*expand_factor, x_range*expand_factor) 
  y_lims <- y_lims + c(-y_range*expand_factor, y_range*expand_factor) 
  
  if(fitness_contour) {
    contour_df <- crossing(Niche_Axis_1 = seq(x_lims[1], x_lims[2], length.out = contour_res),
                           Niche_Axis_2 = seq(y_lims[1], y_lims[2], length.out = contour_res))
    z <- apply(contour_df %>% as.matrix, 1, 
               function(y) K_func(y, h0 = x$sim_params$K_parms$h0, 
                                  sig0i = x$sim_params$K_parms$sigma0i, 
                                  P0i = x$sim_params$K_parms$P0i,
                                  hz = x$sim_params$K_parms$hz, 
                                  biz = x$sim_params$K_parms$biz,
                                  sigiz = x$sim_params$K_parms$sigiz,
                                  Piz = x$sim_params$K_parms$Piz, 
                                  a = x$sim_params$K_parms$a))
    
    contour_df <- contour_df %>%
      mutate(K = z)
  }
  message("Generating plot...")
  pp <- ggplot(plot_df, aes(Niche_Axis_1, Niche_Axis_2))
  if(fitness_contour) {
    pp <- pp + geom_raster(aes(fill = K), data = contour_df) + 
      geom_contour(aes(z = K), data = contour_df, colour = "grey20", bins = bins) +
      scale_fill_scico(palette = "bilbao")
  }
  pp <- pp +
    geom_point(aes(colour = Species, alpha = Time, size = Population)) +
    geom_point(aes(size = Population), data = extant_df, shape = 21, fill = NA) +
    scale_size_area() +
    coord_equal() +
    theme_minimal()
  
  
  pp #theme_void() + theme(legend.position = "none")
  
}

#' @import scico
#' @export plot_K_contour
plot_K_contour <- function(K_parms, x_lims = c(-1, 1), y_lims = c(-1, 1), contour_res = 100, bins = 8, cutoff = 0.01, ...) {
  contour_df <- crossing(Niche_Axis_1 = seq(x_lims[1], x_lims[2], length.out = contour_res),
                         Niche_Axis_2 = seq(y_lims[1], y_lims[2], length.out = contour_res))
  z <- apply(contour_df %>% as.matrix, 1, 
             function(y) K_func(y, h0 = K_parms$h0, 
                                sig0i = K_parms$sig0i, 
                                P0i = K_parms$P0i,
                                hz = K_parms$hz, 
                                biz = K_parms$biz,
                                sigiz = K_parms$sigiz,
                                Piz = K_parms$Piz, 
                                a = K_parms$a,
                                ...))
  
  contour_df <- contour_df %>%
    mutate(K = z) %>%
    dplyr::filter(z > cutoff)
  
  pp <- ggplot(contour_df, aes(Niche_Axis_1, Niche_Axis_2)) +
    geom_raster(aes(fill = K)) +
    geom_contour(aes(z = K), data = contour_df, colour = "grey20", bins = bins) +
    scale_fill_scico(palette = "bilbao") +
    coord_equal() +
    theme_minimal()
  pp
}

#' @import scico
#' @export plot_K_contour3D
plot_K_contour3D <- function(K_parms, x_lims = c(-20, 20), y_lims = c(-20, 20), z_lims = c(-20, 20), contour_res = 100, bins = 8, cutoff = 0.01, slices = 9, facet_cols = 3, ...) {
  
  z_vals <- seq(z_lims[1], z_lims[2], length.out = slices)
  
  contour_df <- crossing(Niche_Axis_1 = seq(x_lims[1], x_lims[2], length.out = contour_res),
                         Niche_Axis_2 = seq(y_lims[1], y_lims[2], length.out = contour_res),
                         Niche_Axis_3 = z_vals)
  z <- apply(contour_df %>% as.matrix, 1, 
             function(y) K_func(y, h0 = K_parms$h0, 
                                sig0i = K_parms$sig0i, 
                                P0i = K_parms$P0i,
                                hz = K_parms$hz, 
                                biz = K_parms$biz,
                                sigiz = K_parms$sigiz,
                                Piz = K_parms$Piz, 
                                a = K_parms$a#,
                                #...
                                ))
  
  contour_df <- contour_df %>%
    mutate(K = z) %>%
    dplyr::filter(z > cutoff)
  
  pp <- ggplot(contour_df, aes(Niche_Axis_1, Niche_Axis_2)) +
    geom_raster(aes(fill = K)) +
    geom_contour(aes(z = K), data = contour_df, colour = "grey20", bins = bins) +
    facet_wrap(~Niche_Axis_3, ncol = facet_cols) + 
    scale_fill_scico(palette = "bilbao") +
    coord_equal() +
    theme_minimal()
  pp
}
#' Function to extract trait history from simulation object as a data_frame
#' @param sim_ob nichefillr_sim object containing the results of a simulation.
#' @import dplyr
#' @importFrom tidyr spread, gather, separate
make_df_from_sim <- function(sim_ob) {
  ms <- sapply(sim_ob$sim_object$extant_list[!sapply(sim_ob$sim_object$extant_list, is.null)], length)
  d <- length(sim_ob$params$K_parms$sig0i)
  
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
  
  full_df <- mapply(trait_mat_as_df, sim_ob$sim_object$full_dat[!sapply(sim_ob$sim_object$full_dat, is.null)], 
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
#' @importFrom magick image_graph
sim_animation <- function(sim_ob, file_name = NULL, expand_factor = 0.1, contour_res = 25, height = 300, width = 400, res = 72) {
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
                     function(x) K_func(x, h0 = sim_ob$params$K_parms$h0, 
                                        sig0i = sim_ob$params$K_parms$sigma0i, 
                                        P0i = sim_ob$params$K_parms$P0i,
                                        hz = sim_ob$params$K_parms$hz, 
                                        biz = sim_ob$params$K_parms$biz,
                                        sigiz = sim_ob$params$K_parms$sigiz,
                                        Piz = sim_ob$params$K_parms$Piz, 
                                        a = sim_ob$params$K_parms$a))
  
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
  
  anim <- image_graph(width = width, height = height, res = 72)
  plot_df %>%
    group_by(Time) %>%
    do(draw_single_plot(., y_lims = y_lims, x_lims = x_lims, pop_lims = pop_lims))
  dev.off()
  
  animation <- image_animate(anim, fps = 20)
  
  if(view) {
    print(animation)
  }
  
  if(!is.null(file_name)) {
    image_write(animation, file_name)
  }
  
  return(animation)
  
}
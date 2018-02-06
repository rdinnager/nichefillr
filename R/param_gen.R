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

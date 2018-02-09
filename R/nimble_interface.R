#' @useDynLib nichefillr
#' @importFrom Rcpp sourceCpp
adapt_landscape_comp_dyn_cmp <- function(d, m, u, a, h0, h_z, P0_i, sigma0_i, P_iz, D_i, b_iz, 
                                         state, V_gi, sigma_iz, gamma_i, c_r, C) {
  ans <- .Call("CALL_rcFun_1", d, m, u, a, h0, h_z, P0_i, 
               sigma0_i, P_iz, D_i, b_iz, state, V_gi, sigma_iz, gamma_i, c_r, C)
  ans <- ans[[18]]
  ans
}
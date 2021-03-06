#ifndef __P_1_RCFUN_1
#define __P_1_RCFUN_1
#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <nimble/NimArr.h>
#include <Rinternals.h>
#include <nimble/accessorClasses.h>
#include <nimble/nimDists.h>
#include <nimble/nimOptim.h>

NimArr<1, double>  rcFun_1 ( int ARG1_d_, int ARG2_m_, int ARG3_u_, double ARG4_a_, double ARG5_h0_, NimArr<1, double> & ARG6_h_z_, NimArr<1, double> & ARG7_P0_i_, NimArr<1, double> & ARG8_sigma0_i_, NimArr<2, double> & ARG9_P_iz_, NimArr<1, double> & ARG10_D_i_, NimArr<2, double> & ARG11_b_iz_, NimArr<1, double> & ARG12_state_, NimArr<1, double> & ARG13_V_gi_, NimArr<2, double> & ARG14_sigma_iz_, NimArr<1, double> & ARG15_gamma_i_, NimArr<1, double> & ARG16_c_r_, double ARG17_C_ );

extern "C" SEXP  CALL_rcFun_1 ( SEXP S_ARG1_d_, SEXP S_ARG2_m_, SEXP S_ARG3_u_, SEXP S_ARG4_a_, SEXP S_ARG5_h0_, SEXP S_ARG6_h_z_, SEXP S_ARG7_P0_i_, SEXP S_ARG8_sigma0_i_, SEXP S_ARG9_P_iz_, SEXP S_ARG10_D_i_, SEXP S_ARG11_b_iz_, SEXP S_ARG12_state_, SEXP S_ARG13_V_gi_, SEXP S_ARG14_sigma_iz_, SEXP S_ARG15_gamma_i_, SEXP S_ARG16_c_r_, SEXP S_ARG17_C_ );
#endif

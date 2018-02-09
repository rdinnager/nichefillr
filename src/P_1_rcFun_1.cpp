#ifndef __P_1_RCFUN_1_CPP
#define __P_1_RCFUN_1_CPP
#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <nimble/EigenTypedefs.h>
#include <Rmath.h>
#include <math.h>
#include <nimble/Utils.h>
#include <nimble/accessorClasses.h>
#include <iostream>
#include <nimble/RcppUtils.h>
#include "P_1_rcFun_1.h"

NimArr<1, double>  rcFun_1 ( int ARG1_d_, int ARG2_m_, int ARG3_u_, double ARG4_a_, double ARG5_h0_, NimArr<1, double> & ARG6_h_z_, NimArr<1, double> & ARG7_P0_i_, NimArr<1, double> & ARG8_sigma0_i_, NimArr<2, double> & ARG9_P_iz_, NimArr<1, double> & ARG10_D_i_, NimArr<2, double> & ARG11_b_iz_, NimArr<1, double> & ARG12_state_, NimArr<1, double> & ARG13_V_gi_, NimArr<2, double> & ARG14_sigma_iz_, NimArr<1, double> & ARG15_gamma_i_, NimArr<1, double> & ARG16_c_r_, double ARG17_C_ )  {
NimArr<1, double> new_state;
double Interm_1;
int Interm_2;
NimArr<1, double> N;
NimArr<2, double> new_N_part;
NimArr<2, double> Tr;
NimArr<2, double> beta_1_ri;
NimArr<1, double> beta_1_r;
NimArr<3, double> beta_2_rzi;
NimArr<2, double> beta_2_rz;
NimArr<3, double> beta_3_rsi;
NimArr<2, double> beta_3_rs;
NimArr<3, double> beta_4_rsi;
NimArr<3, double> beta_5_rzi;
NimArr<3, double> beta_full_ris;
NimArr<3, double> K_num_1_riz;
NimArr<2, double> K_num_1_ri;
NimArr<2, double> a_sum_rz;
NimArr<1, double> a_sum_r;
NimArr<2, double> K_num_2_ri;
NimArr<3, double> alpha_num_1_rsi;
double r;
double i;
double p1;
int Interm_3;
double z;
double p2;
int Interm_4;
int Interm_5;
double Interm_6;
int Interm_7;
double s;
double p3;
int Interm_8;
int Interm_9;
double Interm_10;
int Interm_11;
double Interm_12;
Map<MatrixXd> Eig_N(0,0,0);
EigenMapStrd Eig_ARG12_state_Interm_13(0,0,0, EigStrDyn(0, 0));
Map<MatrixXd> Eig_beta_1_ri(0,0,0);
EigenMapStrd Eig_beta_2_rziInterm_14(0,0,0, EigStrDyn(0, 0));
Map<MatrixXd> Eig_a_sum_rz(0,0,0);
EigenMapStrd Eig_K_num_1_rizInterm_15(0,0,0, EigStrDyn(0, 0));
EigenMapStrd Eig_beta_3_rsiInterm_16(0,0,0, EigStrDyn(0, 0));
EigenMapStrd Eig_beta_full_risInterm_17(0,0,0, EigStrDyn(0, 0));
Map<MatrixXd> Eig_new_N_part(0,0,0);
new_state.initialize(0, true, true, true, static_cast<int>(ARG12_state_.size()));
Interm_1 = (ARG1_d_ * ARG2_m_ + 1);
Interm_2 = (ARG1_d_ * ARG2_m_ + ARG2_m_);
N.setSize((Interm_2 - Interm_1) + 1, 0, 0);
new (&Eig_N) Map< MatrixXd >(N.getPtr(),(Interm_2 - Interm_1) + 1,1);
new (&Eig_ARG12_state_Interm_13) EigenMapStrd(ARG12_state_.getPtr() + static_cast<int>(ARG12_state_.getOffset() + static_cast<int>(0)),ARG12_state_.dim()[0],1,EigStrDyn(0, ARG12_state_.strides()[0]));
Eig_N = (Eig_ARG12_state_Interm_13).block(Interm_1 - 1, 0, (Interm_2 - Interm_1) + 1, 1);
new_N_part.initialize(0, true, true, true, ARG2_m_, ARG2_m_);
Tr.initialize(0, true, true, true, ARG2_m_, ARG1_d_);
beta_1_ri.initialize(0, true, true, true, ARG2_m_, ARG1_d_);
beta_1_r.initialize(0, true, true, true, ARG2_m_);
beta_2_rzi.initialize(0, true, true, true, ARG2_m_, ARG3_u_, ARG1_d_);
beta_2_rz.initialize(0, true, true, true, ARG2_m_, ARG3_u_);
beta_3_rsi.initialize(0, true, true, true, ARG2_m_, ARG2_m_, ARG1_d_);
beta_3_rs.initialize(0, true, true, true, ARG2_m_, ARG2_m_);
beta_4_rsi.initialize(0, true, true, true, ARG2_m_, ARG2_m_, ARG1_d_);
beta_5_rzi.initialize(0, true, true, true, ARG2_m_, ARG3_u_, ARG1_d_);
beta_full_ris.initialize(0, true, true, true, ARG2_m_, ARG1_d_, ARG2_m_);
K_num_1_riz.initialize(0, true, true, true, ARG2_m_, ARG1_d_, ARG3_u_);
K_num_1_ri.initialize(0, true, true, true, ARG2_m_, ARG1_d_);
a_sum_rz.initialize(0, true, true, true, ARG2_m_, ARG3_u_);
a_sum_r.initialize(0, true, true, true, ARG2_m_);
K_num_2_ri.initialize(0, true, true, true, ARG2_m_, ARG1_d_);
alpha_num_1_rsi.initialize(0, true, true, true, ARG2_m_, ARG2_m_, ARG1_d_);
for(r=1; r<= static_cast<int>(ARG2_m_); ++r) {
 for(i=1; i<= static_cast<int>(ARG1_d_); ++i) {
  p1 = ARG12_state_[(ARG1_d_ * (r - 1) + i) - 1];
  Tr((r) - 1, (i) - 1) = p1;
  beta_1_ri((r) - 1, (i) - 1) = pow( static_cast<double>(((pow( static_cast<double>(p1),2)) / static_cast<double>((2 * (pow( static_cast<double>(ARG8_sigma0_i_[(i) - 1]),2)))))),ARG7_P0_i_[(i) - 1]);
 }
 Interm_3 = beta_1_ri.dim()[1];
 new (&Eig_beta_1_ri) Map< MatrixXd >(beta_1_ri.getPtr(),beta_1_ri.dim()[0],beta_1_ri.dim()[1]);
 beta_1_r[(r) - 1] = (((Eig_beta_1_ri).block(r - 1, 0, 1, beta_1_ri.dim()[1])).transpose()).sum();
}
for(r=1; r<= static_cast<int>(ARG2_m_); ++r) {
 for(z=1; z<= static_cast<int>(ARG3_u_); ++z) {
  for(i=1; i<= static_cast<int>(ARG1_d_); ++i) {
   p2 = Tr((r) - 1, (i) - 1) - ARG11_b_iz_((i) - 1, (z) - 1);
   beta_2_rzi((r) - 1, (z) - 1, (i) - 1) = pow( static_cast<double>(((pow( static_cast<double>((p2)),2)) / static_cast<double>((2 * pow( static_cast<double>((ARG14_sigma_iz_((i) - 1, (z) - 1))),2))))),ARG9_P_iz_((i) - 1, (z) - 1));
   beta_5_rzi((r) - 1, (z) - 1, (i) - 1) = p2;
  }
  Interm_4 = beta_2_rzi.dim()[2];
  new (&Eig_beta_2_rziInterm_14) EigenMapStrd(beta_2_rzi.getPtr() + static_cast<int>(static_cast<int>((r - 1) * beta_2_rzi.strides()[0] + (z - 1) * beta_2_rzi.strides()[1])),beta_2_rzi.dim()[2],1,EigStrDyn(0, beta_2_rzi.strides()[2]));
  beta_2_rz((r) - 1, (z) - 1) = (Eig_beta_2_rziInterm_14).sum();
  a_sum_rz((r) - 1, (z) - 1) = exp(-(beta_2_rz((r) - 1, (z) - 1))) * ARG6_h_z_[(z) - 1];
 }
 Interm_5 = a_sum_rz.dim()[1];
 new (&Eig_a_sum_rz) Map< MatrixXd >(a_sum_rz.getPtr(),a_sum_rz.dim()[0],a_sum_rz.dim()[1]);
 Interm_6 = (((Eig_a_sum_rz).block(r - 1, 0, 1, a_sum_rz.dim()[1])).transpose()).sum();
 a_sum_r[(r) - 1] = (Interm_6 + ARG4_a_) + ARG16_c_r_[(r) - 1];
}
for(r=1; r<= static_cast<int>(ARG2_m_); ++r) {
 for(i=1; i<= static_cast<int>(ARG1_d_); ++i) {
  for(z=1; z<= static_cast<int>(ARG3_u_); ++z) {
   K_num_1_riz((r) - 1, (i) - 1, (z) - 1) = -(((((exp(-(beta_2_rz((r) - 1, (z) - 1))) * beta_5_rzi((r) - 1, (z) - 1, (i) - 1)) * ARG6_h_z_[(z) - 1]) * ARG9_P_iz_((i) - 1, (z) - 1)) * (pow( static_cast<double>((pow( static_cast<double>(beta_5_rzi((r) - 1, (z) - 1, (i) - 1)),2) / static_cast<double>((2 * pow( static_cast<double>(ARG14_sigma_iz_((i) - 1, (z) - 1)),2))))),(ARG9_P_iz_((i) - 1, (z) - 1) - 1))))) / static_cast<double>((pow( static_cast<double>(ARG14_sigma_iz_((i) - 1, (z) - 1)),2)));
  }
  Interm_7 = K_num_1_riz.dim()[2];
  new (&Eig_K_num_1_rizInterm_15) EigenMapStrd(K_num_1_riz.getPtr() + static_cast<int>(static_cast<int>((r - 1) * K_num_1_riz.strides()[0] + (i - 1) * K_num_1_riz.strides()[1])),K_num_1_riz.dim()[2],1,EigStrDyn(0, K_num_1_riz.strides()[2]));
  K_num_1_ri((r) - 1, (i) - 1) = (Eig_K_num_1_rizInterm_15).sum();
  K_num_2_ri((r) - 1, (i) - 1) = (((((ARG5_h0_ * exp(-(beta_1_r[(r) - 1]))) * ARG7_P0_i_[(i) - 1]) * Tr((r) - 1, (i) - 1)) * (pow( static_cast<double>(((pow( static_cast<double>(Tr((r) - 1, (i) - 1)),2)) / static_cast<double>((2 * pow( static_cast<double>(ARG8_sigma0_i_[(i) - 1]),2))))),(ARG7_P0_i_[(i) - 1] - 1)))) * a_sum_r[(r) - 1]) / static_cast<double>((pow( static_cast<double>(ARG8_sigma0_i_[(i) - 1]),2)));
 }
}
for(r=1; r<= static_cast<int>(ARG2_m_); ++r) {
 for(s=1; s<= static_cast<int>(ARG2_m_); ++s) {
  for(i=1; i<= static_cast<int>(ARG1_d_); ++i) {
   p3 = Tr((s) - 1, (i) - 1) - Tr((r) - 1, (i) - 1);
   beta_3_rsi((r) - 1, (s) - 1, (i) - 1) = pow( static_cast<double>(((pow( static_cast<double>(p3),2)) / static_cast<double>((2 * (pow( static_cast<double>(ARG15_gamma_i_[(i) - 1]),2)))))),ARG10_D_i_[(i) - 1]);
   beta_4_rsi((r) - 1, (s) - 1, (i) - 1) = p3;
  }
  Interm_8 = beta_3_rsi.dim()[2];
  new (&Eig_beta_3_rsiInterm_16) EigenMapStrd(beta_3_rsi.getPtr() + static_cast<int>(static_cast<int>((r - 1) * beta_3_rsi.strides()[0] + (s - 1) * beta_3_rsi.strides()[1])),beta_3_rsi.dim()[2],1,EigStrDyn(0, beta_3_rsi.strides()[2]));
  beta_3_rs((r) - 1, (s) - 1) = (Eig_beta_3_rsiInterm_16).sum();
  if(r != s) {
   new_N_part((r) - 1, (s) - 1) = (N[(s) - 1] * ARG17_C_) * exp(-(beta_3_rs((r) - 1, (s) - 1)));
  } else {
   new_N_part((r) - 1, (s) - 1) = N[(s) - 1] * exp(-(beta_3_rs((r) - 1, (s) - 1)));
  }
 }
}
for(r=1; r<= static_cast<int>(ARG2_m_); ++r) {
 for(s=1; s<= static_cast<int>(ARG2_m_); ++s) {
  for(i=1; i<= static_cast<int>(ARG1_d_); ++i) {
   if(r != s) {
    alpha_num_1_rsi((r) - 1, (s) - 1, (i) - 1) = ((((beta_4_rsi((r) - 1, (s) - 1, (i) - 1) * ARG17_C_) * ARG10_D_i_[(i) - 1]) * exp(beta_1_r[(r) - 1] - beta_3_rs((r) - 1, (s) - 1))) * (pow( static_cast<double>(((pow( static_cast<double>(beta_4_rsi((r) - 1, (s) - 1, (i) - 1)),2)) / static_cast<double>((2 * pow( static_cast<double>(ARG15_gamma_i_[(i) - 1]),2))))),(ARG10_D_i_[(i) - 1] - 1)))) / static_cast<double>((ARG5_h0_ * a_sum_r[(r) - 1]));
    beta_full_ris((r) - 1, (i) - 1, (s) - 1) = N[(s) - 1] * ((((ARG17_C_ * exp(2 * beta_1_r[(r) - 1] - beta_3_rs((r) - 1, (s) - 1))) * ((ARG5_h0_ * (exp(-(beta_1_r[(r) - 1])))) * (K_num_1_ri((r) - 1, (i) - 1)) - K_num_2_ri((r) - 1, (i) - 1))) / static_cast<double>((pow( static_cast<double>((ARG5_h0_ * a_sum_r[(r) - 1])),2)))) - alpha_num_1_rsi((r) - 1, (s) - 1, (i) - 1));
   } else {
    alpha_num_1_rsi((r) - 1, (s) - 1, (i) - 1) = 0;
    beta_full_ris((r) - 1, (i) - 1, (s) - 1) = N[(s) - 1] * (((exp(2 * beta_1_r[(r) - 1] - beta_3_rs((r) - 1, (s) - 1)) * ((ARG5_h0_ * (exp(-(beta_1_r[(r) - 1])))) * (K_num_1_ri((r) - 1, (i) - 1)) - K_num_2_ri((r) - 1, (i) - 1))) / static_cast<double>((pow( static_cast<double>((ARG5_h0_ * a_sum_r[(r) - 1])),2)))) - alpha_num_1_rsi((r) - 1, (s) - 1, (i) - 1));
   }
  }
 }
}
for(r=1; r<= static_cast<int>(ARG2_m_); ++r) {
 for(i=1; i<= static_cast<int>(ARG1_d_); ++i) {
  Interm_9 = beta_full_ris.dim()[2];
  new (&Eig_beta_full_risInterm_17) EigenMapStrd(beta_full_ris.getPtr() + static_cast<int>(static_cast<int>((r - 1) * beta_full_ris.strides()[0] + (i - 1) * beta_full_ris.strides()[1])),beta_full_ris.dim()[2],1,EigStrDyn(0, beta_full_ris.strides()[2]));
  Interm_10 = (Eig_beta_full_risInterm_17).sum();
  new_state[(ARG1_d_ * (r - 1) + i) - 1] = (N[(r) - 1] * ARG13_V_gi_[(i) - 1]) * Interm_10;
 }
 Interm_11 = new_N_part.dim()[1];
 new (&Eig_new_N_part) Map< MatrixXd >(new_N_part.getPtr(),new_N_part.dim()[0],new_N_part.dim()[1]);
 Interm_12 = (((Eig_new_N_part).block(r - 1, 0, 1, new_N_part.dim()[1])).transpose()).sum();
 new_state[((ARG2_m_ * ARG1_d_ + r)) - 1] = N[(r) - 1] * (1 - ((exp(beta_1_r[(r) - 1]) * Interm_12) / static_cast<double>((ARG5_h0_ * a_sum_r[(r) - 1]))));
}
return(new_state);
}

SEXP  CALL_rcFun_1 ( SEXP S_ARG1_d_, SEXP S_ARG2_m_, SEXP S_ARG3_u_, SEXP S_ARG4_a_, SEXP S_ARG5_h0_, SEXP S_ARG6_h_z_, SEXP S_ARG7_P0_i_, SEXP S_ARG8_sigma0_i_, SEXP S_ARG9_P_iz_, SEXP S_ARG10_D_i_, SEXP S_ARG11_b_iz_, SEXP S_ARG12_state_, SEXP S_ARG13_V_gi_, SEXP S_ARG14_sigma_iz_, SEXP S_ARG15_gamma_i_, SEXP S_ARG16_c_r_, SEXP S_ARG17_C_ )  {
int ARG1_d_;
int ARG2_m_;
int ARG3_u_;
double ARG4_a_;
double ARG5_h0_;
NimArr<1, double> ARG6_h_z_;
NimArr<1, double> ARG7_P0_i_;
NimArr<1, double> ARG8_sigma0_i_;
NimArr<2, double> ARG9_P_iz_;
NimArr<1, double> ARG10_D_i_;
NimArr<2, double> ARG11_b_iz_;
NimArr<1, double> ARG12_state_;
NimArr<1, double> ARG13_V_gi_;
NimArr<2, double> ARG14_sigma_iz_;
NimArr<1, double> ARG15_gamma_i_;
NimArr<1, double> ARG16_c_r_;
double ARG17_C_;
SEXP S_returnValue_1234;
NimArr<1, double> LHSvar_1234;
SEXP S_returnValue_LIST_1234;
ARG1_d_ = SEXP_2_int(S_ARG1_d_);
ARG2_m_ = SEXP_2_int(S_ARG2_m_);
ARG3_u_ = SEXP_2_int(S_ARG3_u_);
ARG4_a_ = SEXP_2_double(S_ARG4_a_);
ARG5_h0_ = SEXP_2_double(S_ARG5_h0_);
SEXP_2_NimArr<1>(S_ARG6_h_z_, ARG6_h_z_);
SEXP_2_NimArr<1>(S_ARG7_P0_i_, ARG7_P0_i_);
SEXP_2_NimArr<1>(S_ARG8_sigma0_i_, ARG8_sigma0_i_);
SEXP_2_NimArr<2>(S_ARG9_P_iz_, ARG9_P_iz_);
SEXP_2_NimArr<1>(S_ARG10_D_i_, ARG10_D_i_);
SEXP_2_NimArr<2>(S_ARG11_b_iz_, ARG11_b_iz_);
SEXP_2_NimArr<1>(S_ARG12_state_, ARG12_state_);
SEXP_2_NimArr<1>(S_ARG13_V_gi_, ARG13_V_gi_);
SEXP_2_NimArr<2>(S_ARG14_sigma_iz_, ARG14_sigma_iz_);
SEXP_2_NimArr<1>(S_ARG15_gamma_i_, ARG15_gamma_i_);
SEXP_2_NimArr<1>(S_ARG16_c_r_, ARG16_c_r_);
ARG17_C_ = SEXP_2_double(S_ARG17_C_);
GetRNGstate();
LHSvar_1234 = rcFun_1(ARG1_d_, ARG2_m_, ARG3_u_, ARG4_a_, ARG5_h0_, ARG6_h_z_, ARG7_P0_i_, ARG8_sigma0_i_, ARG9_P_iz_, ARG10_D_i_, ARG11_b_iz_, ARG12_state_, ARG13_V_gi_, ARG14_sigma_iz_, ARG15_gamma_i_, ARG16_c_r_, ARG17_C_);
PutRNGstate();
PROTECT(S_returnValue_LIST_1234 = Rf_allocVector(VECSXP, 18));
PROTECT(S_returnValue_1234 = NimArr_2_SEXP<1>(LHSvar_1234));
PROTECT(S_ARG1_d_ = int_2_SEXP(ARG1_d_));
PROTECT(S_ARG2_m_ = int_2_SEXP(ARG2_m_));
PROTECT(S_ARG3_u_ = int_2_SEXP(ARG3_u_));
PROTECT(S_ARG4_a_ = double_2_SEXP(ARG4_a_));
PROTECT(S_ARG5_h0_ = double_2_SEXP(ARG5_h0_));
PROTECT(S_ARG6_h_z_ = NimArr_2_SEXP<1>(ARG6_h_z_));
PROTECT(S_ARG7_P0_i_ = NimArr_2_SEXP<1>(ARG7_P0_i_));
PROTECT(S_ARG8_sigma0_i_ = NimArr_2_SEXP<1>(ARG8_sigma0_i_));
PROTECT(S_ARG9_P_iz_ = NimArr_2_SEXP<2>(ARG9_P_iz_));
PROTECT(S_ARG10_D_i_ = NimArr_2_SEXP<1>(ARG10_D_i_));
PROTECT(S_ARG11_b_iz_ = NimArr_2_SEXP<2>(ARG11_b_iz_));
PROTECT(S_ARG12_state_ = NimArr_2_SEXP<1>(ARG12_state_));
PROTECT(S_ARG13_V_gi_ = NimArr_2_SEXP<1>(ARG13_V_gi_));
PROTECT(S_ARG14_sigma_iz_ = NimArr_2_SEXP<2>(ARG14_sigma_iz_));
PROTECT(S_ARG15_gamma_i_ = NimArr_2_SEXP<1>(ARG15_gamma_i_));
PROTECT(S_ARG16_c_r_ = NimArr_2_SEXP<1>(ARG16_c_r_));
PROTECT(S_ARG17_C_ = double_2_SEXP(ARG17_C_));
SET_VECTOR_ELT(S_returnValue_LIST_1234, 0, S_ARG1_d_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 1, S_ARG2_m_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 2, S_ARG3_u_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 3, S_ARG4_a_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 4, S_ARG5_h0_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 5, S_ARG6_h_z_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 6, S_ARG7_P0_i_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 7, S_ARG8_sigma0_i_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 8, S_ARG9_P_iz_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 9, S_ARG10_D_i_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 10, S_ARG11_b_iz_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 11, S_ARG12_state_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 12, S_ARG13_V_gi_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 13, S_ARG14_sigma_iz_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 14, S_ARG15_gamma_i_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 15, S_ARG16_c_r_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 16, S_ARG17_C_);
SET_VECTOR_ELT(S_returnValue_LIST_1234, 17, S_returnValue_1234);
UNPROTECT(19);
return(S_returnValue_LIST_1234);
}
#endif

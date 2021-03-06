## This code was used to generate the C++ files for the simulation using the
## nimble package

library(nimble)

adapt_landscape_comp_dyn <- nimbleFunction(
  run = function(d = integer(0), m = integer(0), u = integer(0), a = double(0),
                 h0 = double(0),
                 h_z = double(1),
                 P0_i = double(1),
                 sigma0_i = double(1),
                 P_iz = double(2),
                 D_i = double(1),
                 b_iz = double(2),
                 state = double(1), V_gi = double(1),
                 sigma_iz = double(2), gamma_i = double(1),
                 c_r = double(1), C = double(0)) {
    
    new_state <- numeric(length = length(state), value = 0, init = TRUE)
    N <- state[(d*m+1):(d*m+m)]
    new_N_part <- matrix(value = 0, nrow = m, ncol = m)
    #new_N <- N
    Tr <- matrix(value = 0, nrow = m, ncol = d, init = TRUE, type = "double")
    #new_Tr <- matrix(value = 0, nrow = m, ncol = d, init = TRUE, type = "double")
    beta_1_ri <- matrix(value = 0, nrow = m, ncol = d, init = TRUE, type = "double")
    beta_1_r <- numeric(length = m, value = 0, init = TRUE)
    beta_2_rzi <- array(value = 0, dim = c(m, u, d), init = TRUE, type = "double")
    beta_2_rz <- matrix(value = 0, nrow = m, ncol = u, init = TRUE, type = "double")
    beta_3_rsi <- array(value = 0, dim = c(m, m, d), init = TRUE, type = "double")
    beta_3_rs <- matrix(value = 0, nrow = m, ncol = m, init = TRUE, type = "double")
    beta_4_rsi <- array(value = 0, dim = c(m, m, d), init = TRUE, type = "double")
    beta_5_rzi <- array(value = 0, dim = c(m, u, d), init = TRUE, type = "double")
    #beta_r <- numeric(length = m, value = 0, init = TRUE)
    beta_full_ris <- array(value = 0, dim = c(m, d, m), init = TRUE, type = "double")
    K_num_1_riz <- array(value = 0, dim = c(m, d, u), init = TRUE, type = "double")
    #a_sum_riz <- array(value = 0, dim = c(m, d, u), init = TRUE, type = "double")
    K_num_1_ri <- matrix(value = 0, nrow = m, ncol = d, init = TRUE, type = "double")
    a_sum_rz <- matrix(value = 0, nrow = m, ncol = u, init = TRUE, type = "double")
    a_sum_r <- numeric(length = m, value = 0, init = TRUE)
    K_num_2_ri <- matrix(value = 0, nrow = m, ncol = d, init = TRUE, type = "double")
    alpha_num_1_rsi <- array(value = 0, dim = c(m, m, d), init = TRUE, type = "double")
    #print("done_1")
    ## Do r by i values
    for (r in 1:m){
      for(i in 1:d){
        ## trait states are arranged to vary by dimension first, then species
        p1 <- state[d * (r - 1) + i] ## get the right trait value from the state vector
        Tr[r, i] <- p1 ## put it into a matrix for later use
        beta_1_ri[r, i] <- ((p1^2) / (2*(sigma0_i[i]^2)))^P0_i[i]
        #b4_ri[r, i] <- p1 / (sigma[i]^2)
      }
      beta_1_r[r] <- sum(beta_1_ri[r, ]) ## sum across d
    }
    #print("done_1")
    # Do r by z by i values
    for (r in 1:m) {
      for (z in 1:u) {
        for (i in 1:d) {
          p2 <- Tr[r, i] - b_iz[i, z]
          beta_2_rzi[r, z, i] <- (((p2)^2) / (2*(sigma_iz[i, z])^2))^P_iz[i, z]
          beta_5_rzi[r, z, i] <- p2
          #K_num_1_rzi[r, z, i] <- exp(-beta_2_rzi[r, z, i])*beta_5_rzi[r, z, i]*h_z[z]*P_iz[i, z]*((beta_5_rzi[r, z, i]^2/sigma_iz[i, z]^2)^(P_iz[i, z] - 1))
        }
        beta_2_rz[r, z] <- sum(beta_2_rzi[r, z, ]) ## sum across d
        #print("done_2")
        a_sum_rz[r, z] <- exp(-beta_2_rz[r, z])*h_z[z]
      }
      #beta_r[r] <- sum(exp(-beta_2_rz[r, ])*h_z[r])
      a_sum_r[r] <- sum(a_sum_rz[r, ]) + a + c_r[r]
    }
    #print("done_2")
    # Do r by i by z values
    
    for (r in 1:m) {
      for (i in 1:d) {
        for (z in 1:u) {
          K_num_1_riz[r, i, z] <- -(exp(-beta_2_rz[r, z])*beta_5_rzi[r, z, i]*h_z[z]*P_iz[i, z]*((beta_5_rzi[r, z, i]^2/(2*sigma_iz[i, z]^2))^(P_iz[i, z] - 1)))/(sigma_iz[i, z]^2)
          #a_sum_riz[r, i, z] <- exp(beta_2_rzi[r, z, i])*h_z[z]
        }
        K_num_1_ri[r, i] <- sum(K_num_1_riz[r, i, ])
        K_num_2_ri[r, i] <- (h0*exp(-beta_1_r[r])*P0_i[i]*Tr[r, i]*(((Tr[r, i]^2)/(2*sigma0_i[i]^2))^(P0_i[i] - 1))*a_sum_r[r]) / (sigma0_i[i]^2)
      }
    }
    
    # Do r by s by i values
    for (r in 1:m) {
      for (s in 1:m) {
        for (i in 1:d) {
          p3 <- Tr[s, i] - Tr[r, i]
          beta_3_rsi[r, s, i] <- ((p3^2) / (2*(gamma_i[i]^2)))^D_i[i]
          beta_4_rsi[r, s, i] <- p3 
        }
        beta_3_rs[r, s] <- sum(beta_3_rsi[r, s, ]) ## sum across d
        if(r != s) {
          new_N_part[r, s] <- N[s]*C*exp(-beta_3_rs[r, s]) ## for population dynamics
        } else {
          new_N_part[r, s] <- N[s]*exp(-beta_3_rs[r, s]) ## for population dynamics
        }
      }
    }
    
    ## add them together
    for (r in 1:m) {
      for(s in 1:m) {
        for(i in 1:d){
          #alpha_num_1_rsi[r, s, i] <- (beta_4_rsi[r, s, i]*D_i[i]*exp(beta_1_r[r] - beta_3_rs[r, s])*(((beta_4_rsi[r, s, i]^2)/(2*gamma_i[i]^2))^(D_i[i] - 1))) / (h0*a_sum_r[r])
          if(r != s) {
            alpha_num_1_rsi[r, s, i] <- (beta_4_rsi[r, s, i]*C*D_i[i]*exp(beta_1_r[r] - beta_3_rs[r, s])*(((beta_4_rsi[r, s, i]^2)/(2*gamma_i[i]^2))^(D_i[i] - 1))) / (h0*a_sum_r[r])
            beta_full_ris[r, i, s] <- N[s]*(((C*exp(2*beta_1_r[r] - beta_3_rs[r, s])*(h0*(exp(-beta_1_r[r]))*(K_num_1_ri[r, i]) - K_num_2_ri[r, i])) / ((h0*a_sum_r[r])^2)) - alpha_num_1_rsi[r, s, i])
          } else {
            alpha_num_1_rsi[r, s, i] <- 0
            beta_full_ris[r, i, s] <- N[s]*(((exp(2*beta_1_r[r] - beta_3_rs[r, s])*(h0*(exp(-beta_1_r[r]))*(K_num_1_ri[r, i]) - K_num_2_ri[r, i])) / ((h0*a_sum_r[r])^2)) - alpha_num_1_rsi[r, s, i])
          }
        }
      }
    }
    #print("done_2")
    # #print("done_3")
    
    # #print("done_4")
    ## last step, sum across s to get r by i values
    for(r in 1:m) {
      for(i in 1:d) {
        #new_Tr[r, i] <- N[r]*V_gi*sum(b[r, i, ])
        new_state[d * (r - 1) + i] <- N[r]*V_gi[i]*sum(beta_full_ris[r, i, ])
        #print(N[r]*V_gi[d]*sum(b[r, i, ]))
      }
      ## calculate new Ns
      #new_N[r] <- N[r]*(1 - ((beta_1_r[r]*sum(new_N_part[r, ])) / (h0*a_sum_r[r])))
      new_state[(m * d + r)] <- N[r]*(1 - ((exp(beta_1_r[r])*sum(new_N_part[r, ])) / (h0*a_sum_r[r]))) ## population dynamics
      #print(N[r]*exp(b2_r[r])*(1 - sum(new_N_part[r, ])))
    }
    # #print(N)
    return(new_state)
    returnType(double(1))
  }
)

cplus1 <- compileNimble(adapt_landscape_comp_dyn, dirName = '.')
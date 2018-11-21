## Test diffeqr

library(diffeqr)
diffeq_setup()

data("example_parms")
parms <- list()

parms$d <- example_parms$macro_parms$d
parms$m <- example_parms$macro_parms$m 
parms$u <- length(example_parms$K_parms$hz)
parms$a <-  example_parms$K_parms$a
parms$h0 <- example_parms$K_parms$h0
parms$h_z <- example_parms$K_parms$hz
parms$P0 <- example_parms$K_parms$P0i[1]
parms$sigma0_i <- example_parms$K_parms$sig0i
parms$P_z <- example_parms$K_parms$Piz[1, , drop = TRUE]
parms$D0 <- example_parms$a_parms$D_i[1]
parms$b_iz <- example_parms$K_parms$biz 
parms$V_gi <- example_parms$macro_parms$V_gi
parms$sigma_iz <- example_parms$K_parms$sigiz
parms$gamma_i <- example_parms$a_parms$gamma_i
parms$c_r <- rep(example_parms$K_parms$c_var, parms$m)
parms$C <- example_parms$a_parms$C

scalar_len <- 8
peaks <- parms$u
dims <- parms$d
peaks_len <- peaks + peaks + dims*peaks + dims*peaks
dims_len <- dims + dims + dims
specs <- parms$m
specs_len <- specs
parms2 <- numeric(length = scalar_len + peaks_len + dims_len + specs_len)
parms2[1] <- example_parms$macro_parms$d
parms2[2] <- example_parms$macro_parms$m
parms2[3] <- length(example_parms$K_parms$hz)
parms2[4] <- example_parms$K_parms$a
parms2[5] <- example_parms$K_parms$h0
parms2[6] <- example_parms$K_parms$P0i[1]
parms2[7] <- example_parms$a_parms$D_i[1]
parms2[8] <- example_parms$a_parms$C

spot <- 9
spot2 <- spot + peaks - 1
parms2[spot:spot2] <- example_parms$K_parms$hz
spot <- spot2 + 1
spot2 <- spot2 + peaks
parms2[spot:spot2] <- example_parms$K_parms$Piz[1, , drop = TRUE]
spot <- spot2 + 1
spot2 <- spot2 + peaks*dims
parms2[spot:spot2] <- as.vector(example_parms$K_parms$biz)
spot <- spot2 + 1
spot2 <- spot2 + peaks*dims
parms2[spot:spot2] <- as.vector(example_parms$K_parms$sigiz)
spot <- spot2 + 1
spot2 <- spot2 + dims
parms2[spot:spot2] <- example_parms$K_parms$sig0i
spot <- spot2 + 1
spot2 <- spot2 + dims
parms2[spot:spot2] <- example_parms$macro_parms$V_gi
spot <- spot2 + 1
spot2 <- spot2 + dims
parms2[spot:spot2] <- example_parms$a_parms$gamma_i
spot <- spot2 + 1
spot2 <- spot2 + specs
parms2[spot:spot2] <- rep(example_parms$K_parms$c_var, parms$m)

state = c(-5, 0, 5, 5, 0.1, 0.1) 

abstol = 1e-8
reltol = 1e-8
saveat = 0:10000/100

tspan <- list(0.0,10000.0)

flipper <- function(y, parms2, t) {
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
  spot <- spot2 + 1
  spot2 <- spot2 + specs
  parms$c_r <- parms2[spot:spot2]
  #print(parms)
  res <- nichefillr:::trait_pop_sim_de_ellip(t, y, parms)
  res[[1]]
}

sol = diffeqr::ode.solve(flipper,
                         state,
                         tspan,
                         p = parms2,
                         abstol = abstol,
                         reltol = reltol)
udf = as.data.frame(sol$u)
plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V5, type = 'scatter3d', mode = 'lines')
plotly::plot_ly(udf, x = ~V3, y = ~V4, z = ~V6, type = 'scatter3d', mode = 'lines')

library(nimble)
drifter <- nimbleFunction(
  run = function(state = double(1), V_gi = double(1), d = integer(0), m = integer(0), demo = double(0)) {
    drift_state <- numeric(length = length(state), value = 0, init = TRUE)
    N <- state[(d*m+1):(d*m+m)]
    drift_N <- demo*N
    for (r in 1:m){
      for(i in 1:d){
        drift_state[d * (r - 1) + i] <- 2*V_gi[i]
      }
      drift_state[(m * d + r)] <- demo*N[r]
    }
    return(drift_state)
    returnType(double(1))
  }
  
  
  # #print(N)

)

drifter_cmp <- compileNimble(drifter, showCompilerOutput = FALSE)

drifter_int <- function(y, parms, t) {
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
  drifter_cmp(y, parms$V_gi, parms$d, parms$m, parms$demo)
  
}

drifter_int(state, parms2, 0)
tspan <- list(0.0,5000.0)
sol2 = diffeqr::sde.solve(flipper,
                         drifter_int,
                         state,
                         tspan,
                         p = parms2)

udf = as.data.frame(sol2$u)
plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V5, type = 'scatter3d', mode = 'lines')
plotly::plot_ly(udf, x = ~V3, y = ~V4, z = ~V6, type = 'scatter3d', mode = 'lines')

plot(udf$V1, udf$V2, type = "l", col = "red", xlim = c(-7, 7), ylim = c(-7, 7))
points(udf$V3, udf$V4, type = "l", col = "blue")


f <- function(u,p,t) {
  du1 = p[1]*(u[2]-u[1])
  du2 = u[1]*(p[2]-u[3]) - u[2]
  du3 = u[1]*u[2] - p[3]*u[3]
  return(c(du1,du2,du3))
}
g <- function(u,p,t) {
  return(c(0.3*u[1],0.3*u[2],0.3*u[3]))
}
u0 = c(1.0,0.0,0.0)
tspan <- list(0.0,100.0)
p = c(10.0,28.0,8/3)
sol = diffeqr::sde.solve(f,g,u0,tspan,p=p,saveat=0.005)
udf = as.data.frame(sol$u)
plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')


######## test new SDE version of nichefillr #########
library(diffeqr)
diffeq_setup()
library(nichefillr)
data("example_parms")

example_parms$macro_parms$demo <- 0.2
example_parms$macro_parms$tot_time <- 50000

test <- sim_radiation(example_parms, trait_hist_prop = 0.25)
plot(test, type = "trace")

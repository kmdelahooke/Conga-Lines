##NUMERIC SIMULATIONS OF CONGA PROBABILITY

#depends on 'conga_finding_functions.R'
# source("./conga_finding_functions.R")

#uses h5 point pattern (window with holes)
library(parallel)

h5 <- point_patterns[[15]]

#----------------------------------------------------------------------------------------
# (1) Assuming homogeneous Poisson distribution
#----------------------------------------------------------------------------------------

#find parameters (mm)
#lambda <- intensity(h5)


csr_monte_carlo <- function(win, r, theta, lambda){
  #list of CSR distributed point patterns
  list <- rpoispp(lambda = lambda, win = win, nsim = 10000)
  
  #apply conga finding function to each of the simulated point patterns (in parallel)
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {library(spatstat)})
  
  dat <- parLapply(cl = cl, X = list , fun = congas, r = r, theta = theta) #iterate over r and theta c+p within cluster
  nconga <- lapply(dat, n_conga)
  print(table(unlist(nconga))) # add to spreadsheet
  
  stopCluster(cl)
}

csr_monte_carlo(win = h5$window, r = 10, theta = 30, lambda = 0.00014) #1
csr_monte_carlo(win = h5$window, r = 12, theta = 30, lambda = 0.00014) #2
csr_monte_carlo(win = h5$window, r = 15, theta = 30, lambda = 0.00014) #3
csr_monte_carlo(win = h5$window, r = 20, theta = 30, lambda = 0.00014) #4 

csr_monte_carlo(win = h5$window, r = 15, theta = 20, lambda = 0.00014) #5
csr_monte_carlo(win = h5$window, r = 15, theta = 25, lambda = 0.00014) #6
csr_monte_carlo(win = h5$window, r = 15, theta = 30, lambda = 0.00014) #7
csr_monte_carlo(win = h5$window, r = 15, theta = 45, lambda = 0.00014) #8

csr_monte_carlo(win = h5$window, r = 15, theta = 30, lambda = 0.00001) #9 
csr_monte_carlo(win = h5$window, r = 15, theta = 30, lambda = 0.00005) #10
csr_monte_carlo(win = h5$window, r = 15, theta = 30, lambda = 0.0001) #11
csr_monte_carlo(win = h5$window, r = 15, theta = 30, lambda = 0.00015) #12

#----------------------------------------------------------------------------------------
# (2) Assuming Thomas clustered distribution
#----------------------------------------------------------------------------------------

#find parameters (mm)
#t <- kppm(h5 ~ 1, clusters = "Thomas", method = "mincon", statistic = "pcf")
#t$clustpar[1] #kappa
#t$clustpar[2] #scale
#t$mu #mu

tc_monte_carlo <- function(win, r, theta, kappa, mu, sigma){
  #list of TC distributed point patterns
  list <- rThomas(kappa = kappa, scale = sigma, mu = mu, win = win, nsim = 10000)
  
  #apply conga finding function to each of the simulated point patterns (in parallel)
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {library(spatstat)})
  
  dat <- parLapply(cl = cl, X = list , fun = congas, r = r, theta = theta) #iterate over r and theta c+p within cluster
  nconga <- lapply(dat, n_conga)
  print(table(unlist(nconga))) # add to spreadsheet
  
  stopCluster(cl)
}

tc_monte_carlo(win = h5$window, r = 10, theta = 30, kappa = 55e-06 , mu = 3, sigma = 45) #1
tc_monte_carlo(win = h5$window, r = 12, theta = 30, kappa = 55e-06 , mu = 3, sigma = 45) #2
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 3, sigma = 45) #3
tc_monte_carlo(win = h5$window, r = 20, theta = 30, kappa = 55e-06 , mu = 3, sigma = 45) #4

tc_monte_carlo(win = h5$window, r = 15, theta = 20, kappa = 55e-06 , mu = 3, sigma = 45) #5
tc_monte_carlo(win = h5$window, r = 15, theta = 25, kappa = 55e-06 , mu = 3, sigma = 45) #6
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 3, sigma = 45) #7
tc_monte_carlo(win = h5$window, r = 15, theta = 45, kappa = 55e-06 , mu = 3, sigma = 45) #8

tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 10e-06 , mu = 3, sigma = 45) #9
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 75e-06 , mu = 3, sigma = 45) #10
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 150e-06 , mu = 3, sigma = 45) #11
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 200e-06 , mu = 3, sigma = 45) #12

tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 2, sigma = 45) #13
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 3, sigma = 45) #14 =18
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 4, sigma = 45) #15
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 6, sigma = 45) #16

tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 3, sigma = 20) #17
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 3, sigma = 45) #18
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 3, sigma = 75) #19
tc_monte_carlo(win = h5$window, r = 15, theta = 30, kappa = 55e-06 , mu = 3, sigma = 120) #20

#----------------------------------------------------------------------------------------
# (3) Probability distribution
#----------------------------------------------------------------------------------------

#Read in results
df <- read.csv("D:/monte_carlo_csr_out.csv", header = T)

#poisson distributed variable
pois <- function(x, lambda) {
  lambda^x * exp(-lambda)/factorial(x)
}

#get expected number of congas from the monte-carlo simulations
get_ex <- function(x){
  x <- x[!is.na(x)]
  n <- seq(0, length(x)-1, 1)
  xn <- x*n
  return(sum(xn)/10000)
}

#plot example probability distribution
plot(sim4/10000 ~ nconga, data = df, xlab = "no. congas", ylab = "Probability")
x <- seq(0, 8, 0.1)
y <- pois(x, get_ex(df$sim4))
lines(y ~ x, col = "skyblue")



#----------------------------------------------------------------------------------------
# (4) Comparison of analytical vs numeric methods
#----------------------------------------------------------------------------------------

## A. CSR distribution (df as above)
##-----------------------------------

# Get numeric E(x) for each simulation
expected_congas <- lapply(df[,-1], get_ex)


# Get analytic E(x)
expected_poisson <- function(theta, r, lambda, A){
  (theta * pi^2 * r^4 *lambda^3 * A)/360
}

A <- area(h5$window)
l <- 0.00014 # control parameter

# Plot
par(mfrow = c(1,3))

#lambda
lambda <- c(0.00001, 0.00005, 0.0001, 0.00015)
ex_l <- as.numeric(unlist(expected_congas[9:12])) #numeric

l2 <- seq(0.00001, 0.00015, 1/A)
lp <- expected_poisson(30, 15, l2, A) #analytic

plot(ex_l ~ lambda, xlab = "intensity", ylab = "E(x)" )
lines(lp ~  l2, col = "red")

#theta
theta <- c(20, 25, 30, 45)
ex_t <- as.numeric(unlist(c(expected_congas[5:8])))

t2 <- seq(20, 45, 0.5)
tp <- expected_poisson(t2, 15, l, A)

plot(ex_t ~ theta, xlab = "Theta", ylab = "E(x)" )
lines(tp ~ t2, col = "orange")

#r 
r <- c(10, 12, 15, 20)
ex_r <- as.numeric(unlist(c(expected_congas[1:4])))

r2 <- seq(10, 20, 0.5)
rp <- expected_poisson(30, r2, l, A)

plot(ex_r ~ r, xlab = "r (mm)", ylab = "E(x)" )
lines(rp ~ r2, col = "skyblue" )


## B. TC distribution
#----------------------------

df <- read.csv("D:/monte_carlo_tc_out.csv", header = T)

# Get numeric E(x) for each simulation
expected_congas <- lapply(df[,-1], get_ex)

# Get analytic E(x)
expected_clustered <- function(r, theta, A, kappa, mu, sigma){
  ex <- theta * A * kappa^3 * mu^3/360 *(pi * r^2 + (1-exp(-(r^2)/(4*sigma^2)))/kappa)^2
  return(ex)
}

#control parameters
kp <- 55e-06
m <- 3
sc <- 45

par(mfrow = c(2,3))

#kappa
kappa <- c(10, 75, 150, 200)*10^(-6)
ex_k <- as.numeric(unlist(c(expected_congas[9:12])))

k2 <- seq(10, 200, 0.5)*10^(-6)
kpp <- expected_clustered(15, 30, A, k2, m, sc)

plot(ex_k ~ kappa, xlab = "Kappa", ylab = "E(x)" )
lines(kpp ~ k2, col = "darkblue")

#mu
mu <- c(2, 3, 4, 6)
ex_m <- as.numeric(unlist(c(expected_congas[13:16])))

m2 <- seq(2, 6, 0.2)
mp <- expected_clustered(15, 30, A, kp, m2, sc)

plot(ex_m ~ mu, xlab = "Mu", ylab = "E(x)" )
lines(mp ~ m2, col = "green")

#sigma
sigma <- c(20, 45, 75, 120)
ex_s <- as.numeric(unlist(c(expected_congas[17:20])))

s2 <- seq(20, 120, 0.5)
sp <- expected_clustered(15, 30, A, kp, m, s2)

plot(ex_s ~ sigma, xlab = "Sigma", ylab = "E(x)" )
lines(sp ~ s2, col = "mediumvioletred")


#theta
theta <- c(20, 25, 30, 45)
ex_t <- as.numeric(unlist(c(expected_congas[5:8])))

t2 <- seq(20, 45, 0.5)
tp <- expected_clustered(15, t2, A, kp, m, sc)

plot(ex_t ~ theta, xlab = "Theta", ylab = "E(x)" )
lines(tp ~ t2, col = "orange")


#r
r <- c(10, 12, 15, 20)
ex_r <- as.numeric(unlist(c(expected_congas[1:4])))

r2 <- seq(10, 20, 0.5)
rp <- expected_clustered(r2, 30, A, kp, m, sc)

plot(ex_r ~ r, xlab = "r (mm)", ylab = "E(x)" )
lines(rp ~ r2, col = "skyblue" )


##Finding probability of congas given a random distribution
library(spatstat)
library(plyr)
library(parallel)

##Numerically
#-----------------
#initialise cluster
cl <- makeCluster(detectCores())
clusterEvalQ(cl, {library(spatstat)})
             
#conga finding functions
congas <- function (pp, r, theta){
  coords <- c(1, 1, 1)
  
  #neighbourhood of r around every point in point pattern
  nd <- applynbd(pp,  R = r, function(Y, ...){return(Y)}, exclude = TRUE)
  
  #function to calculate bearing between 2 points
  angle <- function(cur, Y){
    r <- atan2(y = (Y$y - cur$y), x = (Y$x - cur$x)) * 180 / pi - 90
    r[r < 0] <- r + 360
    return(r)
  }
  # plot base map
  #plot(pp, cols = "darkslategray4")
  
  #loop though each neighbourhood that contains points
  for(i in 1:length(nd)) {
    if(nd[[i]]$n > 0){
      
      #define 1st and 2nd points and their indices
      p2 <- nd[[i]]
      ind2 <- nncross(p2, pp, k = 1)$which
      p1 <- pp[i] #i = ind1
      
      #neighbourhood of each element
      for(j in 1:p2$n){
        pp2 <- pp[c(-i, -ind2[j])]
        p3 <- nncross(p2[j], pp2, k = seq(1, 5, 1))
        
        #find p3 within r mm of p2
        for(l in 1:5){
          if(p3[, l] < r) {
            #bearings
            b2 <- angle(pp2[p3[, l + 5]], p2[j])
            b1 <- angle(p2[j], p1)
            
            #check if aligned, allowing tolerance theta
            if(abs(b2 - b1) < theta) {
              #plot(pp2[p3[, l + 5]], col = "red4", add = TRUE) #plot p3s
              #plot(p2[j], col = "red4", add = TRUE) #plot p2s associated with those p3s
              
              #get  x-coords of each of 3 points as a point id.
              cds <- c(pp2[p3[, l + 5]]$x, pp[ind2[j]]$x, pp[i]$x)
              coords <- rbind(coords, cds)
              
            }
          }
        }
      }
    }
  }
  return(coords)
}

n_conga <- function(coords){
  ifelse(is.null(nrow(coords)), return(0), {
    
    coords <- data.frame(coords[-1, ])
    
    #remove rows of x-coords which contain the same x-coords as other rows
    for(i in seq_len(nrow(coords))){
      coords[i,] <- sort(coords[i,])
    }
    uni <- unique(coords)
    for(i in 1:nrow(uni)) {
      if(any(uni[i, ] %in% as.vector(data.matrix(uni[-i,])))){
        uni <- uni[-i,]
      }
    }
    return(nrow(uni)) #number of rows = number of congas
  })
}

#get intensity of the point pattern 
#lambda <- intensity(point_patterns[[5]])
lambda <- 0.0001
cat(paste("lambda :",lambda), sep = "\n")

#list of poisson distributed point patterns
list <- rpoispp(lambda = lambda, win = windows[[5]], nsim = 10000)

#Find number of congas
#dat <- lapply(list, FUN = congas, r = 15, theta = 25)

cat("r: 15, theta: 15", sep = "\n")
dat <- parLapply(cl = cl, X = list , fun = congas, r = 15, theta = 15)
nconga <- lapply(dat, n_conga)
table(unlist(nconga))

cat("r: 15, theta: 20", sep = "\n")
dat <- parLapply(cl = cl, X = list , fun = congas, r = 15, theta = 20)
nconga <- lapply(dat, n_conga)
table(unlist(nconga))

cat("r: 15, theta: 45", sep = "\n")
dat <- parLapply(cl = cl, X = list , fun = congas, r = 15, theta = 45)
nconga <- lapply(dat, n_conga)
table(unlist(nconga))

cat("r: 5, theta: 30", sep = "\n")
dat <- parLapply(cl = cl, X = list , fun = congas, r = 5, theta = 30)
nconga <- lapply(dat, n_conga)
table(unlist(nconga))

cat("r: 12, theta: 30", sep = "\n")
dat <- parLapply(cl = cl, X = list , fun = congas, r = 12, theta = 30)
nconga <- lapply(dat, n_conga)
table(unlist(nconga))

cat("r: 20, theta: 30", sep = "\n")
dat <- parLapply(cl = cl, X = list , fun = congas, r = 20, theta = 30)
nconga <- lapply(dat, n_conga)
table(unlist(nconga))

stopCluster(cl)


##Read in results
df <- read.csv("E:/random_dist_congas.csv", header = T)
str(df)

#probability distribution

pois <- function(x, lambda) {
  lambda^x * exp(-lambda)/factorial(x)
}

get_ex <- function(x){
  x <- x[!is.na(x)]
  n <- seq(0, length(x)-1, 1)
  xn <- x*n
  return(sum(xn)/10000)
}

plot(sim4/10000 ~ nconga, data = df, xlab = "no. congas", ylab = "Probability")
x <- seq(0, 8, 0.1)
y <- pois(x, get_ex(df$sim4))
lines(y ~ x, col = "skyblue")


##Get E(x) for each simulation
expected_congas <- lapply(df, get_ex)

#Relationships between E(x) parameter
par(mfrow = c(1,1))
#intensity/lambda
lambda <- c(0.00005, 0.0001, 0.00015, 0.0002)
ex_l <- as.numeric(unlist(expected_congas[6:9]))
plot(ex_l ~ lambda, xlab = "intensity", ylab = "E(x)" )

#theta
theta <- c(15, 20, 30, 45)
ex_t <- as.numeric(unlist(c(expected_congas[10:11], expected_congas[7], expected_congas[12])))
plot(ex_t ~ theta, xlab = "Theta", ylab = "E(x)" )


# r (mm)
r <- c(5, 12, 15, 20)
ex_r <- as.numeric(unlist(c(expected_congas[13:14], expected_congas[7], expected_congas[15])))
plot(ex_r ~ r, xlab = "r (mm)", ylab = "E(x)" )

#Analytically
#-------------

#get n from intensity
A <- area(windows[[5]])
n <- round(lambda*A)

pois_prob <- function(theta, r, n, A){
  (theta/180)*pi^2*r^4*(n*(n-1)/A)*(n/(2*A))
}
#lambda
n2 <- seq(237, 949, 1)
l2 <- seq(0.00005, 0.0002, 2.103787e-07)
lp <- pois_prob(30,15, n2, A)

lines(lp ~  l2, col = "red")

#theta
t2 <- seq(15, 45, 0.5)
tp <- pois_prob(t2, 15, 474, A)

lines(tp ~t2, col = "orange")

#r
r2<-seq(5,20,0.5)
rp <- pois_prob(30, r2, 474, A)
lines(rp ~r2, col = "skyblue" )


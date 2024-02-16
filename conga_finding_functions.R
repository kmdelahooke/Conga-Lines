##CONGA FINDING FUNCTIONS##

#includes:
# 1. congas
# 2. n_congas
# 3. compare_congas
# 4: plot_expected
# 5: plot_probability

#Libraries
library(spatstat)
library(plyr)
library(parallel)
library(gmp)

#-----------------------------------------------------------------------------------------------
# (1) Find 3-point alignments (r, theta) within a point pattern
#-----------------------------------------------------------------------------------------------
# pp: spatstat point pattern
# r: max. spacing between points (in units of point pattern)
# theta: max. bearing between points (in degrees)

#if un-commented maps showing position of the 3-point alignments will be plotted

# returns a list of coordinates of points that are aligned

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


#-----------------------------------------------------------------------------------------------
# (2) Count the number of 3-point alignments found
#-----------------------------------------------------------------------------------------------
# coords: list of coordinate outputted by "congas" function
# returns the number of 3-point alignments on the surface

n_conga <- function(coords){
  ifelse(is.null(nrow(coords)), return(0), {
    
    coords <- coords[-1, ]
    
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


#------------------------------------------------------------------------------------------------
# (3) Find the number of 3-point alignments on the surface for a range of r (5-50mm), alongside 
# the probability that they are formed by chance
#------------------------------------------------------------------------------------------------

# pp: spatstat point pattern
# theta: max. bearing between points in degrees

# fits CSR and TC processes to the point pattern, finds the best fit model, and uses that as a background distribution

# returns a data frame containing the following parameters:
## rs: value of r
## values: number of 3-point alignments found 
## ex_c: expected number of 3-point alignments assuming TC distribution
## ex_p: expected number of 3-point alignments assuming CSR
## p_c: probability of values, given ex_c
## p_p: probability of values, given ex_p
## pd_c: Diggle's goodness of fit test score for TC distribution model fit
## pd_p: Diggle's goodness of fit test score for CSR model fit

compare_congas <- function(pp, theta){
  
  #range of r (5 - 50mm)
  rs <- seq (5, 50, 1)
  
  #find congas along range of r
  cat("Finding congas.....", sep = "\n")
  
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {library(spatstat)})
  ac <- parLapply(cl = cl, X = rs, fun = congas, pp = pp, theta = theta)
  stopCluster(cl)
  values <- unlist(lapply(ac, n_conga))
  cat("Done.", sep = "\n")
  
  #fit thomas cluster model to point pattern
  library(spatstat)
  t <- kppm(pp ~ 1, clusters = "Thomas", method = "mincon", statistic = "pcf")
  
  #fit CSR model to point pattern
  p <- ppm(pp ~ 1)
  
  #find goodness of fit
  pd_t <- dclf.test(t)$p.value
  pd_r <- dclf.test(p)$p.value
  
  cat(paste("Thomas Cluster Model pd =", pd_t), sep = "\n")
  cat(paste("CSR pd =", pd_r), sep = "\n")
  
  #get p_d values for df
  pd_c <- rep(pd_t, length(rs))
  pd_p <- rep(pd_r, length(rs))
  
  #extract parameters
  A <- area(pp$window)
  kappa <- unname(t$clustpar[1])
  mu <- unname(t$mu)
  sigma <- unname(t$clustpar[2])
  int <- intensity(pp)
  
  
  #Expected congas - clustered distribution
  expected_clustered <- function(r, theta, A, kappa, mu, sigma){
    ex <- theta * A * kappa^3 * mu^3/360 *(pi * r^2 + (1-exp(-(r^2)/(4*sigma^2)))/kappa)^2
    return(ex)
  }
  
  ex_c <- unlist(lapply(rs, expected_clustered, theta = theta, A = A, kappa = kappa, mu = mu, sigma = sigma))
  
  #Expected congas - random distribution
  expected_poisson <-function(r, theta, lambda, A){
    (theta * pi^2 * r^4 *lambda^3 * A)/360
  }
  
  ex_p <- unlist(lapply(rs, expected_poisson, theta = theta, lambda = int, A = A))
  
  
  #calculate probability and add column to dataframe
  pois <- function(x, lambda) {
    lambda^x * exp(-lambda)/factorial(x)
  }
  
  p_c <- pois(x = values, lambda = ex_c)
  p_p <- pois(x = values, lambda = ex_p)
  
  
  #create dataframe
  df <- data.frame(cbind(rs, values, ex_c, ex_p, p_c, p_p, pd_c, pd_p))
  
  
  return(df)
}

#--------------------------------------------------------------------------------------------------------------------
# (4) Plot expected versus actual number of conga lines for range of r
#--------------------------------------------------------------------------------------------------------------------
# df: dataframe output from 'compare_congas'

plot_expected <- function(df){
  
  #work out which distribution to assume
  diff <- df$pd_c[1] - df$pd_p[1]
  
  ifelse(diff>0, ex <- df$ex_c, ex <-df$ex_p)
  ifelse(diff>0, print("assuming clustered distribution", sep = "\n"), print("assuming random distribution", sep = "\n"))
  
  plot(values ~ rs, type = "n", xlab = "r (mm)", ylab = "Number of congas", data = df)
  lines(values ~ rs, col = "skyblue", data = df)
  lines(ex ~ df$rs, col = "red")
  #lines(ex_p ~ rs, col = "orange", data = df)
  #legend("bottomright", legend = c("Actual", "Expected"), lty = 1, col = c("skyblue", "red"))
}


#--------------------------------------------------------------------------------------------------------------------
# (5) Plot probability of actual congas found for range of r
#--------------------------------------------------------------------------------------------------------------------
# df: dataframe output from 'compare_congas'

plot_probability <- function(df){
  
  #work out which distribution to assume
  diff <- df$pd_c[1] - df$pd_p[1]
  ifelse(diff>0, df$ex <- df$ex_c, df$ex <-df$ex_p)
  ifelse(diff>0, df$p <- df$p_c, df$p <-df$p_p)
  ifelse(diff>0, print("assuming clustered distribution", sep = "\n"), print("assuming random distribution", sep = "\n"))
  
  #split based on if more or less congas than expected
  more <- subset(df, values > ex)
  less <- subset(df, ex > values)
  
  #calculate values between you get more or less congas than expected
  sig_m <- subset(more, p < 0.05)
  sig_l <- subset(less, p < 0.05)
  
  #print ranges where significant difference
  cat("Assuming a clustered distribution:", sep = "\n")
  cat("----------------------------------", sep = "\n")
  cat(paste("More congas than expected (p < 0.05):", range(sig_m$rs)[1], "mm -", range(sig_m$rs)[2], "mm"), sep = "\n")
  cat(paste("Fewer congas than expected (p < 0.05):", range(sig_l$rs)[1], "mm -", range(sig_l$rs)[2], "mm"), sep = "\n")
  
  #plot
  plot(log(p) ~ rs, type = "n", xlab = "r (mm)", ylab = "log probability of actual number of congas", data = df)
  lines(log(p) ~ rs, col = "snow3", data = df)
  points(log(p) ~ rs, col = "red", pch = 20, data = more)
  points(log(p) ~ rs, col = "skyblue", pch = 20, data = less)
  abline(h = log(0.05), col = "snow4", lty = 2)
  text(47, -2.8, labels = "p = 0.05", col = "snow4")
  #legend("topright", legend = c("fewer congas", "more congas"), lty = 1, col = c("skyblue", "red"))
  
}


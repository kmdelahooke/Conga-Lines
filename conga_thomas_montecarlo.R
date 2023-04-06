#THOMAS CLUSTER MONTE-CARLO 
library(tictoc)

#1. Fit thomas cluster model to all erect fronds on surface
#----------------------------------------------------------
source("C:/Users/kmd47/OneDrive - University of Cambridge/PhD/dex_repo/dex/dex_outline.R")
source("C:/Users/kmd47/OneDrive - University of Cambridge/PhD/dex_repo/dex/dex.R")
H5 <- "C:/Users/kmd47/OneDrive - University of Cambridge/PhD/H5_mapping/H5_complete/H5_complete_dex.svg"

win.test <- dex_outline(H5, retrodeform = FALSE)
h5 <- dex(H5, retrodeform = FALSE)
h5 <- subset(h5, taxon=='primocandelabrum'|taxon=='taxona'|taxon=='taxonb'|taxon=='frond'|taxon=='charniodiscus'|taxon=='bradgatia'|taxon=='parviscopa'| taxon == 'feffaced')
h5 <- ppp(h5$x, h5$y, window = win.test)

t <- kppm(h5 ~ 1, clusters = "Thomas", method = "mincon")
t
#2. Function to detect conga lines of 3 fronds
#---------------------------------------------
#(1) pp = point pattern or fitted model
#(2) r = distance between fronds (mm)
#(3) theta = leeway in bearing allowed

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
  plot(pp, pch = 20, cols = "darkslategray4")
  
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
              plot(pp2[p3[, l + 5]], col = "red4", cex = 6, pch = "*", add = TRUE) #plot p3s
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

#3. Count the number of conga lines detected using function "congas"
#------------------------------------------------------------------
#(1) coords: returned dataframe of x coordinates from function "congas"

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

#4. Check that it can correctly locate and count congas on h5
#------------------------------------------------------------
h5_coords <- congas(h5, r = 15, theta = 30)
n_conga(h5_coords)

#5. Function to be used in monte carlo simulation 
#------------------------------------------------
#(1) pp: point pattern or fitted thomas cluster model

fun <- function(pp) {
  coords <- congas(pp, r = 15, theta = 25)
  nconga <<- c(nconga, n_conga(coords)) #global assignment
}

#6. Run monte carlo simulation and count number of 3 frond congas
#----------------------------------------------------------------
#library(parallel)
#detectCores()
#tic()
list <- rThomas(kappa = t$clustpar[1], scale = t$clustpar[2], mu = t$mu, win = win.test, nsim = 10)
nconga <- 0
#fun <-mclapply(list, FUN = fun, mc.cores = 12)
fun <- lapply(list, FUN = fun)
nconga <- nconga[-1]
table(nconga)
toc()

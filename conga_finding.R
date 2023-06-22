##FINDING CONGA LINES ON MULTIPLE SURFACES

#depends on 'conga_finding_functions.R'
# source("./conga_finding_functions.R")

#workflow given folder of .csvs of surface markups

#---------------------------------------------------------------------------------------
# (1) Import surface data
#---------------------------------------------------------------------------------------

filelist <- list.files(path = "D:/surface_data", pattern = ".*.csv")

get_df <- function(filename){
  x <- (paste0("D:/surface_data/", filename))
  df <- read.table(x, sep = ",", header = T)
  return(df)
}

surface_data <- lapply(filelist, get_df)


## CREATE WINDOWS

#create_windows <- function(j){
#  plot(j$y ~j$x)
#  win <- clickpoly(add = T)
#  return(win)
#}

#windows <- lapply(surface_data, create_windows)
#coords <- lapply(windows, function(x){data.frame(x$bdry[[1]])})
#windowpath <- unlist(lapply(filelist, function(x){paste("D:/windows/",x)}))
#mapply(write.csv, coords, windowpath)

##IMPORT WINDOWS

window_files <- list.files(path = "D:/windows", pattern = ".*.csv")

get_df2 <- function(filename){
  x <- (paste0("D:/windows/", filename))
  df <- read.table(x, sep = ",", header = F)
  df <- df[,2:3] # remove labels
  df <- df[-1,] # remove header
  return(df)
}

coords <- lapply(window_files, get_df2)

make_windows <- function(x){
  l <- list(x= as.numeric(x[,1]), y = as.numeric(x[,2]))
  win <- owin(poly = l)
  return(win)
}

windows <- lapply(coords, make_windows)


#Add in H5 windows with holes

source("D:/dex_outline.R")
h5 <- "D:/H5_complete_dex.svg"

h5h_win <- dex_outline(h5, retrodeform = FALSE)
h5hr_win <- dex_outline(h5, retrodeform = TRUE)

windows <- append(windows, list(h5h_win, h5hr_win))
surface_data <- append(surface_data, surface_data[10:11])


#-----------------------------------------------------------------------------------------------
# (3) Create point patterns
#-----------------------------------------------------------------------------------------------

point_patterns <- vector("list", length(surface_data))

for(i in 1:length(surface_data)){
  point_patterns[[i]] <- ppp(surface_data[[i]]$x, surface_data[[i]]$y, window = windows[[i]])
}

#-----------------------------------------------------------------------------------------------
# (4) Find conga lines
#-----------------------------------------------------------------------------------------------

cf <- lapply(point_patterns, compare_congas, theta = 30)

#Export
filenames <- gsub(".csv", "", filelist)
filenames <- c(filenames, "h5h", "h5hr")

for(i in 1:length(cf)){
  write.csv(cf[[i]], paste0("D:/congas_found/", filenames[[i]], "_cf.csv"))
}

#--------------------------------------------------------------------------------------------------------------------
# (5) Plotting results
#--------------------------------------------------------------------------------------------------------------------

#Read data back in if required
filelist <- list.files(path = "D:/congas_found/", pattern = ".*.csv")

get_df <- function(filename){
  x <- (paste0("D:/congas_found/", filename))
  df <- read.table(x, sep = ",", header = T)
  return(df)
}

cf <- lapply(filelist, get_df)


#view goodness of fit
for(i in 1:16){
  print(cf[[i]]$pd_p[1])
}

#plot results
par(mfrow = c(4,4))
lapply(cf, plot_expected)
lapply(cf, plot_probability)



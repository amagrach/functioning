#test simmons algorithm

#First I use Song example to keep it simple----
# I also aim to load as less packages as possile to avoid conflicts

#Song approach
source('Data_and_code_Song_et_al_2017/toolbox.R') #load the toolbox

web <- read.csv("Data_and_code_Song_et_al_2017/network.csv") #load network.csv
web <- as.matrix(web)
head(web)
print(NODF <- nestedness_NODF(web)) # this calculates the raw value of NODF
print(max_NODF_S <- max_nest(web)) # this calculates the maximum value of NODF for that network
print(combined_NODF_S <- comb_nest(web,NODF,max_NODF_S)) # this calculates the combined NODF statistic as described in the manuscript
#works as expected.

#Now Simmons approach. 
#NODF <- nestedness_NODF(web) # this calculates the raw value of NODF #as in song
devtools::install_github("CHoeppke/maxnodf")
library(maxnodf)
max_NODF <- maxnodf(web=web, quality=2) # this calculates the maximum value of NODF for that network, based on Simmons correction of Song
#Error in maxnodf(web = web, quality = 2) : 
 #Invalid network. Ensure all row and column totals are >= 1.
colSums(web)
web2 <- web[,-which(colSums(web) < 0.1)]
#start again
nestedness_NODF(web)
(NODF <- nestedness_NODF(web2)) # Fuck, this value changes!!
max_NODF <- maxnodf(web=web2, quality=2) # this calculates the maximum value of NODF for that network, based on Simmons correction of Song
#interestin max_nest(web2) has the same exact value!!
max_NODF_S == max_NODF$max_nodf
max_nest(web2) == max_NODF$max_nodf
#combine
combined_NODF <- comb_nest(web2,NODF,max_NODF$max_nodf) # this calculates the combined NODF statistic as described in the manuscript
#compatre
combined_NODF_S == combined_NODF #FALSE
comb_nest(web2,NODF,max_nest(web2)) == combined_NODF #TRUE

#test issue links > ncol + nrows
sum(web2)
ncol(web2) + nrow(web2)

#for failing in network analysis.R
sum(web4)
ncol(web4) + nrow(web4) #small difference!

#opening up the black box

#debugged using web4, which should fail.
web = web4
max_NODF_tuneada <- function (web, quality = 0) {
  NodesA <- -1
  NodesB <- -1
  Edges <- -1
  if (is.matrix(web)) {
    if (all(is.numeric(web))) {
      web[web > 0] <- 1
      if (any(web < 0)) {
        stop("Invalid network. Ensure all elements of web >= 0.")
      }
      mt_0 <- computeMT0(web) #not in package :(
      mt_t <- computeMTt(web)
      if (any(mt_0 < 1)) {
        stop("Invalid network. Ensure all row and column totals are >= 1.")
      }
      if (any(mt_t < 1)) {
        stop("Invalid network. Ensure all row and column totals are >= 1.")
      }
      NodesA <- nrow(web)
      NodesB <- ncol(web)
      Edges <- sum(web)
    }
    else {
      stop("Parameter 'web' is expected to be a numeric matrix or a numeric vector.")
    }
  }
  else if (is.vector(web)) {
    if (length(web) == 3) {
      if (all(round(web) == web)) {
        NodesA <- web[[1]]
        NodesB <- web[[2]]
        Edges <- web[[3]]
      }
      else {
        stop("The vector 'web' is expected to have three integer members")
      }
    }
    else {
      stop("The vector 'web' is expected to have three integer members")
    }
  }
  else {
    stop("Parameter 'web' is expected to either be a matrix or a vector containing the matrix dimensions and number of links.")
  }
  #NACHO comments this out to trick the function. NOT WORKING because it calls Rccp later on.
  #if (Edges <= NodesA + NodesB) {
   # stop("Number of links needs to satisfy 'Links > nrow(web) + ncol(web).")
  #}
  #if (Edges > NodesA * NodesB) {
   # stop("Number of links needs to satisfy 'Links <= nrow(web) * ncol(web).")
  #}
  if (!quality %in% 0:2) {
    stop("Please chose a valid quality parameter. Options: \n\tquality = 0 -> Use a very fast greedy algorithm.\n\tquality = 1 -> Improved result using hillclimbing in combination with greedy.\n\tquality = 2 -> Use a simulated annealing algorithm. Best results but requires the most computation time.")
  }
  if (quality == 0) {
    mtx <- greedy_solve2(NodesA, NodesB, Edges)
    cat("\n")
  }
  else if (quality == 1) {
    mtx <- greedy_solve2(NodesA, NodesB, Edges)
    cat("\n")
    mtx <- full_hill_climb_cpp(mtx)
  }
  else if (quality == 2) {
    mtx <- greedy_solve2(NodesA, NodesB, Edges) #loaded via copy paste on github. BUt nodf_cpp not available :(
    cat("\n")
    mtx <- sim_anneal_opt_cpp(mtx)
  }
  return(list(max_nodf = nodf_cpp(mtx), max_nodf_mtx = mtx))
}


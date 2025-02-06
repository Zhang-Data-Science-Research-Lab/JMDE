## Simulate 3D structures and corresponding contact maps
## Contact: Yuping Zhang (yuping.zhang@uconn.edu) and Zhengqing Ouyang ( ouyang@schoolph.umass.edu ) 

## Create the directories in the code first if they do not exist. 
## n.track = 1 for MDE
## n.track > 1 for JMDE 

library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(plot3D)

# Structure
set.seed(21)
# obtain the skeleton first
GetStructure <- function(n){
  x <- cos((1:n) / 3)
  y <- sin((1:n) / 3)
  z <- 1:n / 20
  return(cbind(x, y, z))
}

stru <- GetStructure(100) # 100 loci
distmat <- as.matrix(dist(stru))   ### The lower triangle of the distance matrix stored by columns in a vector, say do. If n is the number of observations, i.e., n <- attr(do, "Size"), then for i < j ≤ n, the dissimilarity between (row) i and j is do[n*(i-1) - i*(i-1)/2 + j-i]. The length of the vector is n*(n-1)/2, i.e., of order n^2.

distsq = distmat^2 ### squared distances


Simu <- function(beta0, beta, distsq){
  n <- dim(distsq)[1]
  contactmat <- matrix(0, n, n)
  for (i in 1:(n - 1)){
    for (j in (i+1):n){
      lambda = exp(beta0 - 0.5* beta *log(distsq[i, j]))
      contactmat[i, j] <- rpois(1, lambda)  ### lambda mean parameter
      contactmat[j, i] <- contactmat[i, j]
    }
  }
  diag(contactmat) <- 0
  print(sum(contactmat != 0) / (n ^ 2 - n))
  return(contactmat)
}
# a = 0.95 ===> 30%
# a = 2 ===> 50%
# a = 4 ===> 70%
# a = 9 ==> 90%

dir = "./Simulations/Helix/"

if (! file.exists(dir)){
  dir.create(file.path(getwd(), dir))
}

pdf(paste0(dir, "trueHelix.pdf"))
scatter3D(as.matrix(stru[, 1]), as.matrix(stru[, 2]), as.matrix(stru[, 3]), type = "l", lwd=3,font=2, font.tex =2, main = "True Structure for Helix Simulations", colkey = F)
dev.off()


png(paste0(dir, "trueHelix.png"))
scatter3D(as.matrix(stru[, 1]), as.matrix(stru[, 2]), as.matrix(stru[, 3]), type = "b", lwd=3, font=2, font.tex =2, main = "True Structure for Helix Simulations", colkey = F)
dev.off()


# Simulate 100 contactmaps

n.track = 2
betas = c(1.5, 1.6)
beta0s = c(2, 2)

sink(paste0(dir, "parameters.txt"))
cat("betas \t", betas, "\n")
cat("beta0s \t", beta0s, "\n")
sink()

for (m in 1:n.track){
  beta0 = beta0s[m]
  beta = betas[m]
  sink(file=paste0(dir, "track", m, "_data/", "coverage.txt" ))
  for (i in 1:100){
    set.seed(i)
    simulated_contactmap = Simu(beta0, beta, distsq)
    filename <- paste0(dir, "track", m, "_data/", i, ".txt", sep = "")
    write.table(simulated_contactmap, file = filename, row.names = F, col.names = F)
  }
  sink()
}

# Save structure
write.table(stru, file = paste0(dir, "Helix_Structure.txt"), row.names = F, col.names = F)


###############
# Random Walk #
###############

GaussianRW <- function(n.times, dimension = 3){
  increments <- matrix(rnorm(n.times * dimension), n.times, dimension)
  coordinates <- apply(increments, 2, cumsum)
  return(coordinates)
}

set.seed(1)
stru <- GaussianRW(100)
distmat <- as.matrix(dist(stru))


distsq = distmat^2 ### squared distances

dir = "./Simulations/Randomwalk/"


# a = 5.5 ===> 30%
# a = 13 ==> 50%
# a = 25 ==> 70%
# a = 60 ==> 90%

pdf(paste0(dir, "trueRandomWalk.pdf"))
scatter3D(as.matrix(stru[, 1]), as.matrix(stru[, 2]), as.matrix(stru[, 3]), type = "l", lwd=3,font=2, font.tex =2, main = "True Structure for Random-Walk Simulations", colkey = F)
dev.off()


png(paste0(dir, "trueRandomWalk.png"))
scatter3D(as.matrix(stru[, 1]), as.matrix(stru[, 2]), as.matrix(stru[, 3]), type = "b", lwd=3, font=2, font.tex =2, main = "True Structure for Random-Walk Simulations", colkey = F)
dev.off()

# Simulate 100 contactmaps

n.track = 2

betas = c(1.5, 1.6)
beta0s = c(4, 4)


sink(paste0(dir, "parameters.txt"))
cat("betas \t", betas, "\n")
cat("beta0s \t", beta0s, "\n")
sink()

for (m in 1:n.track){
  beta0 = beta0s[m]
  beta = betas[m]
  sink(file=paste0(dir, "track", m, "_data/", "coverage.txt" ))
  for (i in 1:100){
    set.seed(i)
    simulated_contactmap = Simu(beta0, beta, distsq)
    filename <- paste0(dir, "track", m, "_data/", i, ".txt", sep = "")
    write.table(simulated_contactmap, file = filename, row.names = F, col.names = F)
  }
  sink()
}

# Save structure
write.table(stru, file = paste0(dir, "RandomWalk_Structure.txt"), row.names = F, col.names = F)


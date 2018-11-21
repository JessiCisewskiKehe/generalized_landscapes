### Pickup Sticks Simulator (STIX)
### By Jessi Cisewski-Kehe on November 20, 2018

### Install/Load R packages
# install.packages("TDA")
# install.packages("imager")
library(TDA)
library(imager)


### Read in functions
source("functions.R") #generalized landscape functions
source("stix_simulator.R") #Pickup sticks simulator


### Generate a sample of STIX data
n <- 100
thickness <- 6
xrange <- 1
yrange <- 1
image.name <- "sample"
set.seed(1234)
# An image is generated and saved to your working directory
getPickUpSticks(n, thickness, image.name, xrange, yrange)


### Read-in image
im1 <- load.image(paste0(image.name, ".png"))
x1 <- im1[,,1,1]
rm(im1); gc()  #remove unnecessary variables


### Apply a mild amount of smoothing to image (for sharp contrasts)
# (this takes a few seconds to run)
gridX1 <- expand.grid(1:dim(x1)[1], 1:dim(x1)[2]) #define grid
cX1 <- as.vector(x1) #vectorize the image
fit1 <- loess(cX1 ~ gridX1[,1] + gridX1[,2], degree = 2, span = .001) #local polynomial smoothing
y1 <- matrix(fit1$fitted, nrow = dim(x1)[1], ncol = dim(x1)[2], byrow = FALSE) #turn back into matrix
rm(gridX1, cX1, fit1); gc()  #remove unnecessary variables


### Get persistence diagram
# (this takes a few seconds to run)
Diag1 <- gridDiag(FUNvalues = y1, sublevel = FALSE,location = FALSE, printProgress = TRUE, maxdimension = 2)$diagram 	


### Plot persistence diagram
plot(Diag1[,2], Diag1[,3], col = Diag1[,1]+1, pch = Diag1[,1]+1, xlab = "Birth", ylab = "Death")
abline(a = 0, b = 1, lwd = 1)
legend("bottomright", legend = c(expression(H[0]), expression(H[1])), pch = 1:2, col = 1:2)


### Compute/plot landscape function (TDA package)
which_hom <- 1
tseq <- seq(min(Diag1[Diag1[,1]==which_hom,2:3]), max(Diag1[Diag1[,1]==which_hom,2:3]), length.out = 1000)
KK <- 10 #number of function orders
land <- landscape(Diag1, dimension = which_hom, tseq, KK = 1:KK)
plot(tseq, land[,1], "l", xlab = "t", ylab = "Landscape function")
out <- sapply(2:KK, function(ii) lines(tseq, land[,ii], lty = ii, col = ii))


### Compute/plot silhouette function (TDA package)
which_hom <- 1
tseq <- seq(min(Diag1[Diag1[,1]==which_hom,2:3]), max(Diag1[Diag1[,1]==which_hom,2:3]), length.out = 1000)
silh <- silhouette(Diag1, p = 1, dimension = which_hom, tseq)
plot(tseq, silh, "l", xlab = "t", ylab = "Silhouette function")


### Compute/plot generalized landscape function (from functions.R)
which_hom <- 1
tseq <- seq(min(Diag1[Diag1[,1]==which_hom,2:3]), max(Diag1[Diag1[,1]==which_hom,2:3]), length.out = 1000)
KK <- 10 #number of function orders
kern_fn = kern_land  #kern_epa (Epanechnikov), kern_land (triangle kernel), kern_tri (tricubic kernel)
h = 0.10 #select a bandwidth
rot_diagram <- RotateDiagram(Diag1, which.dim = which_hom)      
glandscape <- KernelMatrix(rot_diagram, tseq=tseq, kern_function = kern_fn, h) 

plot(tseq, glandscape[,1], "l", xlab = "t", ylab = "Generalized landscape function")
out <- sapply(2:10, function(ii) lines(tseq, glandscape[,ii], lty = ii, col = ii))














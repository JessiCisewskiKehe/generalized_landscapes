### Pickup Sticks Simulator (STIX)
### By Jessi Cisewski-Kehe on March 8, 2018
### Please send questions to jessica.cisewski@yale.edu

### Run the function below
getPickUpSticks <- function(n, thickness, image.name, xrange, yrange){
# n = number of sticks 
# thickness = degrees of freedom of X^2 for thickness of sticks'
# image.name = name to save the file as
# xrange = constant > 0 for specifying the range for points to be uniformly sampled in the horizontal direction 
# yrange = constant > 0 for specifying the range for points to be uniformly sampled in the horizontal direction 
# Note:  the xrange and yrange is arbitrary; the images are set to appear square
# Note:  the edges of the image are trimmed

if (!is.numeric(n) || length(n) != 1 || n <= 0) {
        stop("n should be a positive integer")
}

if (!is.numeric(thickness) || length(thickness) != 1 || thickness <= 0) {
        stop("thickness should be a positive integer")
}

if (!is.numeric(xrange) || length(xrange) != 1 || xrange <= 0) {
        stop("xrange should be a positive integer")
}

if (!is.numeric(yrange) || length(yrange) != 1 || yrange <= 0) {
        stop("yrange should be a positive integer")
}

# Define range for random sampling of points
gridx = c(-xrange,xrange)
gridy = c(-yrange,yrange)

# Randomly sample the points on grid
getPoints1 <- cbind(runif(n, gridx[1], gridx[2]), runif(n, gridy[1], gridy[2]))
getPoints2 <- cbind(runif(n, gridx[1], gridx[2]), runif(n, gridy[1], gridy[2]))

# Sample the thickness of each line from a Chi-square distribution with "thickness" degrees of freedom
lwd0 <- rchisq(n, thickness)

# Grayscale for color of lines
col0 <- gray.colors(20, start = .5, end = 1)[sample(1:20,n,replace = TRUE)]

png(paste0(image.name,".png"), width = 600, height = 600)
par(mar = c(0,0,0,0), oma = c(0,0,0,0), bg = "black")
plot(getPoints1[,1], getPoints1[,2], xlab = "", ylab = "", pch = 19, cex = .05, xlim = gridx*.5, ylim = gridy*.5, xaxt = "n", yaxt = "n", axes = FALSE)
points(getPoints2[,1], getPoints2[,2], pch = 19, cex = .05)
segments(getPoints1[,1], getPoints1[,2], getPoints2[,1], getPoints2[,2], lwd = lwd0, col = col0)
dev.off()
}




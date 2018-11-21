### Included are functions and sub-functions for computing generalized landscapes.
library(TDA)
library(pracma)
library(parallel)

###------------------------------------------------------
RotateDiagram <- function(Diag1, which.dim){
	# Used for rotating persistence diagrams
	# Diag1 = persistence diagram object
	# which.dim = which homology dimension to rotate
	d <- Diag1[Diag1[,1]==which.dim,2:3]
	### Rotated diagram
	rot1 <- cbind((d[,1] + d[,2])/2, (d[,2] - d[,1])/2)
	return(rot1)
	}
	
###------------------------------------------------------
RotateDiagramRips <- function(Diag1, which.dim){
	# Used for rotating persistence diagrams
	# Diag1 = persistence diagram object
	# which.dim = which homology dimension to rotate
	if(which.dim == 0){
		Diag1 <- Diag1[-1, ] # Remove the point at infinity for dim = 0
	}
	d <- Diag1[Diag1[,1]==which.dim,2:3]
	### Rotated diagram
	if(length(d)>2){
		rot1 <- cbind((d[,1] + d[,2])/2, (d[,2] - d[,1])/2)
		}else{rot1 <- cbind((d[1] + d[2])/2, (d[2] - d[1])/2)}
	return(rot1)
	}

###------------------------------------------------------
kern_tri <- function(x0, y0, x, h){
	###  kernel function (Tricube)
	# (x0, y0) are the observations
	# x is the location to compute the kernel
	# h is the kernel bandwidth
	u <- (x-x0)/h
	out <- ifelse(abs(u)<=1, y0*(1 - abs(u)^3)^3, 0)
	return(out)
}

###------------------------------------------------------
kern_epa <- function(x0, y0, x, h){
	###  kernel function (Epanechnikov)
	# (x0, y0) are the observations
	# x is the location to compute the kernel
	# h is the kernel bandwidth
	u <- (x-x0)/h
	out <- ifelse(abs(u)<=1, y0*(1 - abs(u)^2), 0)
	return(out)
}

###------------------------------------------------------
kern_land <- function(x0, y0, x, h){
	###  kernel function (triangle/landscape)
	# (x0, y0) are the observations
	# x is the location to compute the kernel
	# h is the kernel bandwidth
	u <- (x-x0)/h
	out <- ifelse(abs(u)<=1, y0*(1-sign(u)*u), 0)
	return(out)
}

###------------------------------------------------------
kern_land_original <- function(x0, y0, x, h){
	###  kernel function (usual landscape)
	# (x0, y0) are the observations
	# x is the location to compute the kernel
	out <- ifelse(x-x0 < y0 - x, x-x0, y0-x)
	out[out<0] <- 0
	return(out)
}


###------------------------------------------------------
### Get the kernels for each point
KernelMatrix <- function(rot, tseq, kern_function, h){
	#rot = rotated diagram
	kernel_matrix <- matrix(0, nrow = length(tseq), ncol = nrow(rot))
		for(ii in 1:nrow(rot)){
			kernel_matrix[,ii] <- kern_function(rot[ii, 1], rot[ii, 2], tseq, h)
	}
	lambda <- sapply(seq(along = tseq), FUN = function(i) {
        sort(kernel_matrix[i, ], decreasing = TRUE)
    })
	return(t(lambda))
}


###------------------------------------------------------
### Get the kernels for each point
KernelMatrixGland <- function(rot, tseq, kern_function, h){
	#rot = rotated diagram
	if(length(rot)>2){
		kernel_matrix <- matrix(0, nrow = length(tseq), ncol = nrow(rot))
			for(ii in 1:nrow(rot)){
				kernel_matrix[,ii] <- kern_function(rot[ii, 1], rot[ii, 2], tseq, h)
		}
		lambda <- sapply(seq(along = tseq), FUN = function(i) {
        	sort(kernel_matrix[i, ], decreasing = TRUE)
    	})
    	}else{
			lambda <- matrix(kern_function(rot[1], rot[2], tseq, h), ncol = 1)
    }
	return(t(lambda))
}

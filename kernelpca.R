# AA2, kernel project
# Cristina Aguiera, Jesus Antonanzas
# GCED, November 2019

# ---------------------------------- KPCA class -------------------------------- #

# Create kernel matrix from a data frame. For now, kernels must be available in "kernlab", 
# but other kernels can be implemented.
kmat <- function(data, kernel = "vanilladot", params = list()) {
  n <- dim(data)[1]
  kernel <- base::do.call(kernel, params)
  return(kernelMatrix(kernel, data))
}

# Given a kernel matrix compute the "ncomp" principal components in the feature space
# with associated eigenvalue > "thres". "ncomp" = -1 means compute all of them.
p.comp.vect <- function(kmatrix, thres, ncomp) {
  n = dim(kmatrix)[1]
  
  # recall thory: eigenvals and eigenvects are of K/n.
  eig <- eigen(kmatrix/n, symmetric = T)
  
  # Pick valid components (>"thres"). Threshold can be user defined.
  nonzero <- which(eig$values > thres)
  if (ncomp != -1) {
    if (ncomp < length(nonzero)) {
      nonzero <- nonzero[1:ncomp]
    } else if (ncomp > length(nonzero)) {
      stop("'ncomp' bigger than number of valid PC's. Please pick a lower 'ncomp' value")
    }
  }
  evalues <- eig$values[nonzero]
  
  # Normalize eigenvectors of the feature space.
  evectors <- t(t(eig$vectors[,nonzero])/sqrt(eig$values[nonzero]))
  colnames(evectors) <- paste("PC", nonzero, sep = "")
  return(list("eigenvalues" = evalues, "princomp" = evectors))
}

# Project a point 'x' into the FS PCs as f(x) = t(\alpha)*g(x), 
# where g(x) = t(k(x, x1), ..., k(x, x_n)). Can operate on set of points and project
# if g(x) not specified (centered) by the user, it is computed given the kernel specified and its parameters.
# Again, compatible with kernlab kernels for now. 'tpoints' are the points used to compute the PC's from. 
# Project them onto "ncomp" <= "number of training points" feature principal components.
project <- function(x = NULL, kpca = NULL, gx = NULL, ncomp = NULL) {
  tpoints <- kpca$train.data
  kernel <- kpca$kernel
  par <- kpca$par
  if (is.null(kpca)) { stop("please specify 'kpca'") }
  
  # If gx null we have to compute it.
  else if (is.null(gx)) {
    kern <- base::do.call(kernel, par)
    n <- dim(tpoints)[1]
    m <- dim(x)[1]
    
    # If we are projecting just one point, 'm' will be NULL.
    if (is.null(m)) { 
      gx <- matrix(nrow = n, ncol = 1)
      for (i in 1:n) {
        gx[i,1] = kern(x, tpoints[i,])
      }
    } 
    else {
      gx <- matrix(nrow = m, ncol = n)
      for (i in 1:m) {
        for (j in 1:n) {
          gx[i,j] = kern(x[i,], tpoints[j,])
        }
      }
    }
    km <- kmat(kpca$train.data, kernel = kernel, params = par)
    n <- dim(km)[1]
    
    # Center the points in the FS according to theory.
    gx <- t(t(gx - rowSums(gx)/n) - rowSums(km)/n) + sum(km)/(n^2)
  }
  if (is.null(ncomp)) { ncomp = length(kpca$eigenvalues) }
  projection <- (gx) %*% (kpca$princomp)[,1:ncomp]
  colnames(projection) <- paste("PC", 1:dim(projection)[2], sep = "")
  return(projection)
}

# Fully integrated KPCA method. The kernel matrix can be already specified. 
kernpca <- function(data, kmatrix = NULL, kernel = "vanilladot", par = list(), thres = 1e-4, ncomp = -1) {
  if (is.null(kmatrix)) {
    kmatrix <- kmat(data, kernel = kernel, params = par)
  }
  else if ((tail(eigen(kmatrix)$values,1) < 0)) { stop("kernel matrix non psd") }
  n <- dim(kmatrix)[1]
  
  # Center training points in the FS.
  k <- t(t(kmatrix - colSums(kmatrix)/n) - rowSums(kmatrix)/n) + sum(kmatrix)/(n^2)
  k.princomp <- p.comp.vect(k, thres, ncomp)
  
  # Add data so methods are more comfortable to use.
  k.princomp$train.data <- data
  k.princomp$kernel <- kernel
  k.princomp$par <- par
  
  # Project training points onto "ncomp" or all kernel principal components.
  projection <- project(kpca = k.princomp, gx = k)
  return(list("princomp" = k.princomp$princomp, "eigenvalues" = k.princomp$eigenvalues, 
              "projected.training.data" = projection, "train.data" = k.princomp$train.data,
              "kernel" = k.princomp$kernel, "par" = k.princomp$par))
}

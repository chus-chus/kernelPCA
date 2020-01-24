# AA2, kernel project
# Cristina Aguiera, Jesus Antonanzas
# GCED, November 2019
# Experimenting with KPCA

library(datasets)
library(kernlab)
library(mlbench)
library(BKPC) # for time execution comparison

source('kernelpca.R')

set.seed(888)

# --------------------------------- How to make it work ---------------------------------  #

# Data must be a matrix.
iris <- data.matrix(iris)

# Train KPCA on "train" data. Compatible with kernel functions in library
# 'kernlab': "rbfdot", "polydot", "tanhdot", "vanilladot", "laplacedot" and more. See the documentation
# for info on the parameters: https://rdrr.io/cran/kernlab/man/dots.html

# par is a list() object containg parameters:

# "rbfdot" -> "sigma" (inverse kernel width)
# "laplacedot" -> "sigma"
# "polydot" -> "degree", "scale", "offset"
# "tanhdot" -> "scale", "offset"

train <- iris[1:20,]
testing <- kernpca(train[,-5], kernel = "rbfdot", par = list(sigma = 0.1))

# Project new points.
newdata <- iris[51:61,]
project(x = newdata[,-5], kpca = testing)

# We can project into less feature principal components.
project(x = newdata[,-5], kpca = testing, ncomp = 2)

# And choose different kernels.
  kernpca(train[,-5], kernel = "polydot", par = list(offset = 0, scale = 2, degree = 2))
  project(x = newdata[,-5], kpca = testing, ncomp = 2)
  
  # Linear kernel -> standard PCA
  kernpca(train[,-5], kernel = "vanilladot", par = list())

# --------------------------------- Execution time ---------------------------------  #

# Ours
start_time <- Sys.time()
testing <- kernpca(iris, kernel = "rbfdot", par = list(sigma = 0.1))
end_time <- Sys.time()

end_time - start_time

# kernlab::kpca
start_time <- Sys.time()
testing <- kernlab::kpca(iris, kernel = "rbfdot", par = list(sigma = 0.1))
end_time <- Sys.time()

end_time - start_time

# BKPC::kPCA
start_time <- Sys.time()
test <- BKPC::kPCA(kmat(iris, "rbfdot", list(sigma = 0.1)))
end_time <- Sys.time()

end_time - start_time

# --------------------------------- Visualizing ---------------------------------  #

# Plots into the first two feature principal components
plotkpca <-function (kernelfu, kerneln, color = 1) {
    xpercent <- kernelfu$eigenvalues[1]/sum(kernelfu$eigenvalues)*100
    ypercent <- kernelfu$eigenvalues[2]/sum(kernelfu$eigenvalues)*100
    
    plot(kernelfu$projected.training.data, col=color,
         main=paste(paste("Kernel PCA (", kerneln, ")", format(xpercent+ypercent,digits=3)), "%"),
         xlab=paste("1st PC -", format(xpercent,digits=3), "%"),
         ylab=paste("2nd PC -", format(ypercent,digits=3), "%"))
}
  
# We also want to plot the found principal components in the input space, to see how well
# the shape of the data is explained. "pca.obj" must be an object from "prcomp".
plotpca <- function(pca.obj, data, color) {
    vects <- pca.obj$rotation
    plot(data, col = color, ylab = "y", xlab = "x")
    for (i in 1:2) {
      sl <- (vects[2,i]/vects[1,i])
      f <- function(x) x*sl
      curve(f, type="l", col="black", add=TRUE)
    }
}

# ---- Simple function
N <- 250
f <- function(x) x
x <- runif(N, min = -2, max = 2)                     
t <- f(x) + rnorm(N, sd=0.5)
data <- data.frame(x, t)

plot(data)
curve(f, type="l", col="blue", add=TRUE)

# PCA
pca <- prcomp(data)
plotpca(pca, data, 1)

# kPCA
kpcairis <- kernpca(iris[,1:4], kernel = "rbfdot", par = list(sigma=0.005))
plotkpca(kpcairis, "rbfdot, sigma = 0.005", color = iris[,5])  

# Notice how the first feature principal component explains almost all the structure.
plot(kpcairis$eigenvalues, type = "b")

# ---- Banana dataset
N <- 500
f <- function(x) (x^2)  
x <- runif(N, min = -2, max = 2)                     
t <- f(x) + rnorm(N, sd=0.2)      
bana <- data.frame(x, t)

plot(bana)
curve(f, type="l", col="green", add=TRUE)

# PCA
pca <- prcomp(bana)
plotpca(pca, bana, 1)

# kPCA
kpcabana <- kernpca(data.matrix(bana), kernel = "polydot", par = list(degree = 2, scale = 1, offset = 0))
plotkpca(kpcabana, "polydot, deg = 2, scale = 1, off = 0")
plot(kpcabana$eigenvalues, type = "b")

# Observe that the variability is almost completely explained by the fisrt two principal components.
# The fact that the plot is almost the same as the original means that the structure of the 
# data is almost completely explained in the feature space.

# ---- Circular dataset
N <- 500
c1 <- function(x) sqrt(1-x^2)  
c2 <- function(x) -sqrt(1-x^2) 
x <- runif(N, min = -1, max = 1)                   
circle <- matrix(nrow = 500, ncol = 2)
t1 <- c1(x[1:(N/2)]) + rnorm(N/2, sd=0.1)      
t2 <- c2(x[((N/2)+1):N]) + rnorm(N/2, sd=0.1)   
circle[,1] <- x
circle[1:250,2] <- t1
circle[251:500,2] <- t2

plot(circle)

# PCA
pca <- prcomp(circle)
plotpca(pca, circle, 1)

# kPCA
circle <- data.matrix(circle)
kpcacirc <- kernpca(circle, kernel = "rbfdot", par = list(sigma = 0.01))
plotkpca(kpcacirc, "rbfdot, sigma = 0.01", color = 1)  
plot(kpcacirc$eigenvalues, type = "b")


# ---- Half moon

N <- 200

make.sinusoidals <- function(m,noise=0.2) 
{
  x <- c(1:2*m)
  y <- c(1:2*m)
  
  for (i in 1:m) {
    x[i] <- (i/m) * pi
    y[i] <- sin(x[i]) + rnorm(1,0,noise)
  }
  
  for (j in 1:m) {
    x[m+j] <- (j/m + 1/2) * pi
    y[m+j] <- cos(x[m+j]) + rnorm(1,0,noise)
  }
  
  target <- as.factor(c(rep(+1,m),rep(-1,m)))
  
  return(data.frame(x,y,target))
}

dataset <- make.sinusoidals (N)

plot(dataset$x,dataset$y,col=dataset$target)

# PCA
pca <- prcomp(dataset[,1:2])
plotpca(pca, dataset[,1:2], col = dataset$target)

# kPCA
kpca_sinus <- data.matrix(dataset)
kpca_s <- kernpca(kpca_sinus, kernel = "rbfdot", par = list(sigma=2))
plotkpca(kpca_s, "rbfdot, sigma = 2", color = dataset$target)  
plot(kpca_s$eigenvalues, type = "b", ylab = "eigenvalues")


# --------------------------------- Advanced Visualization: Isolines ---------------------------------  #

# Let's plot isolines corresponding to some principal components of the feature space.
# These isolines are going to be for 2d data.

library(plotly)

# Returns a matrix containing the point representation of the grid containing some data.
# Granularity is 50x50
create_grid <- function(original_data, extra_margin = 0.1) {
  # extract grid limits
  minx <- min(original_data[,1])
  maxx <- max(original_data[,1])
  miny <- min(original_data[,2])
  maxy <- max(original_data[,2])
  
  # generate points
  epsx <- extra_margin*(maxx - minx)
  epsy <- extra_margin*(maxy - miny)
  x <- seq(minx-epsx, maxx+epsx, length.out = 50)
  y <- seq(miny-epsy, maxy+epsy, length.out = 50)
  
  # bind test points
  grid <- matrix(ncol = (length(y)+1), nrow = length(x))
  grid[,1] <- x
  for (i in 2:(length(y)+1)) {
      grid[,i] <- rep((y[i-1]), length(x))
  }
  return(grid)
}

# Plots isolines for some data onto a given feature principal component and
# overlays the training data.
plot_isolines <- function(grid, original_data, projections) {
  plot_ly(
      x = grid[,1], 
      y = grid[1,2:ncol(grid)], 
      z = t(projections), 
      type = "contour",
      colorscale = "YlOrRd") %>% 
  add_trace(x = original_data[,1], 
            y = original_data[,2], 
            type = "scatter",
            color = I("black"),
            mode = "markers") %>%
  layout(showlegend = FALSE, showscale = FALSE)
}

# Returns a 3d array: "nrow(grid)" arrays of projections onto "ncomp" principal component features.
project_grid <- function(grid, kpca_obj, numcomp = NULL) {
  if (is.null(numcomp)) {numcomp = length(kpca_obj$eigenvalues)}
  p_grid <- array(dim = c(nrow(grid), numcomp, nrow(grid)))
  for (i in 1:(ncol(grid)-1)) {
      p_grid[,,i] <- project(x = cbind(grid[,1], grid[,i]), kpca = kpca_obj, ncomp = numcomp)
  }
  return(p_grid)
}


# ---- Simple dataset
N <- 250
f <- function(x) x  
x <- runif(N, min = -2, max = 2)                    
t <- f(x) + rnorm(N, sd=0.5)    
data <- data.frame(x, t)
plot(data)

kpca_test <- kernpca(data.matrix(data), kernel = "polydot", par = list(offset = 0, degree = 2))

grid <- create_grid(data)

# Can take a minute.
proj <- project_grid(grid, kpca = kpca_test, numcomp = 2)

subplot(plot_isolines(grid, data, proj[,1,]), plot_isolines(grid, data, proj[,2,]))


# ---- Banana dataset
set.seed(888)
N <- 250
f <- function(x) (x^2) 
x <- runif(N, min = -1, max = 1)    
t <- f(x) + rnorm(N, sd=0.2)
bana <- data.frame(x, t)

plot(bana)

kpcabana <- kernpca(data.matrix(bana), kernel = "tanhdot", par = list(scale = 1, offset = 0))

grid_bana <- create_grid(bana, extra_margin = 0.1)

p_grid_bana <- project_grid(grid_bana, kpca = kpcabana, numcomp = 3)

subplot(plot_isolines(grid_bana, bana, p_grid_bana[,1,]), 
        plot_isolines(grid_bana, bana, p_grid_bana[,2,]),
        plot_isolines(grid_bana, bana, p_grid_bana[,3,]),
        shareY = TRUE, titleX = TRUE)


# ---- Circular dataset
N <- 250
c1 <- function(x) sqrt(1-x^2) 
c2 <- function(x) -sqrt(1-x^2) 
x <- runif(N, min = -1, max = 1)   
circle1 <- matrix(nrow = 2*N, ncol = 2)
t1 <- c1(x[1:(N/2)]) + rnorm(N/2, sd=0.1)  
t2 <- c2(x[((N/2)+1):N]) + rnorm(N/2, sd=0.1)
circle1[,1] <- x
circle1[1:250,2] <- t1
circle1[251:500,2] <- t2

plot(circle1)

kpcacirc <- kernpca(data.matrix(circle1), kernel = "rbfdot", par = list(sigma = 0.1))

grid_circ <- create_grid(circle1, extra_margin = 0.2)

proj <- project_grid(grid_circ, kpca = kpcacirc, numcomp = 4)

subplot(plot_isolines(grid_circ, circle1, proj[,1,]), plot_isolines(grid_circ, circle1, proj[,2,]),
        plot_isolines(grid_circ, circle1, proj[,3,]), plot_isolines(grid_circ, circle1, proj[,4,]),
        shareY = TRUE, shareX = TRUE, nrows = 2)


# ---- Clusters

# RBF Gaussia 
set.seed(777)
clusters <- mlbench.2dnormals(n = 500, cl = 3, r = 4, sd = 0.8)

clust_kpca <- kernpca(clusters$x, kernel = "rbfdot", par = list(sigma=0.3))

grid_clus <- create_grid(clusters$x, extra_margin = 0.15)

proj_clus <- project_grid(grid_clus, kpca = clust_kpca, numcomp = 8)

subplot(plot_isolines(grid_clus, clusters$x, proj_clus[,1,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,2,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,3,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,4,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,5,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,6,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,7,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,8,]),
        nrows = 2, shareX = T, shareY = T)

# Hyperbolic tangent kernel
clusters <- mlbench.2dnormals(n = 250, cl = 3, sd = 0.25)

clust_kpca <- kernpca(clusters$x, kernel = "tanhdot", par = list(scale=2, offset = 1))

grid_clus <- create_grid(clusters$x, extra_margin = 0.15)

proj_clus <- project_grid(grid_clus, kpca = clust_kpca)

subplot(plot_isolines(grid_clus, clusters$x, proj_clus[,1,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,2,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,3,]))

# see that vanilladot (linear kernel) is the same as linear pca
clust_kpca <- kernpca(clusters$x, kernel = "vanilladot", par = list())

grid_clus <- create_grid(clusters$x, extra_margin = 0.15)

proj_clus <- project_grid(grid_clus, kpca = clust_kpca)

subplot(plot_isolines(grid_clus, clusters$x, proj_clus[,1,]),
        plot_isolines(grid_clus, clusters$x, proj_clus[,2,]),
        shareX = T, shareY = T)

pcaclust <- prcomp(clusters$x)
plotpca(pcaclust, data = clusters$x, color = 1)

# ---- Beans & Circle
set.seed(888)
cass <- mlbench.cassini(n = 500, relsize = c(2,2,1))

plot(cass)

cass_kpca <- kernpca(cass$x, kernel = "rbfdot", par = list(sigma = 0.2))

grid_cass <- create_grid(cass$x, extra_margin = 0.15)

proj_cass <- project_grid(grid_cass, kpca = cass_kpca, numcomp = 8)

subplot(plot_isolines(grid_cass, cass$x, proj_cass[,1,]),
        plot_isolines(grid_cass, cass$x, proj_cass[,2,]),
        plot_isolines(grid_cass, cass$x, proj_cass[,3,]),
        plot_isolines(grid_cass, cass$x, proj_cass[,4,]),
        plot_isolines(grid_cass, cass$x, proj_cass[,5,]),
        plot_isolines(grid_cass, cass$x, proj_cass[,6,]),
        plot_isolines(grid_cass, cass$x, proj_cass[,7,]),
        plot_isolines(grid_cass, cass$x, proj_cass[,8,]),
        nrows = 2, shareX = T, shareY = T)

plot(cass_kpca$eigenvalues, type = "b", ylab = "eigenvalues")

# see that vanilladot (linear kernel) is the same as linear pca
cass_kpca <- kernpca(cass$x, kernel = "vanilladot", par = list())

proj_cass <- project_grid(grid_cass, kpca = cass_kpca)

subplot(plot_isolines(grid_cass, cass$x, proj_cass[,1,]),
        plot_isolines(grid_cass, cass$x, proj_cass[,2,]),
        shareX = T, shareY = T)

pcacass <- prcomp(cass$x)
plotpca(pcacass, data = cass$x, color = cass$classes)

# --------------------------------- Predicting ---------------------------------  #

# Now let's do some performance benchmarking.

measures <- function(real, pred) {
  t <- table(truth=real, predicted = pred)
  print(t)
  err <- (1-sum(diag(t))/sum(t))
  return(err)
}

# ------------ Pixel Classification ------------- #

pix_tr <- read.csv("px_train.csv", skip = 1)
pix_te <- read.csv("px_test.csv", skip = 1)

# Non relevant variables: they are constant.
pix_tr <- pix_tr[,-4]
pix_te <- pix_te[,-4]

# Mix observations
set.seed(888)
pix_tr <- pix_tr[sample.int(nrow(pix_tr)),]
pix_te <- pix_te[sample.int(nrow(pix_te)),]

# kPCA (Gaussian RBF)
kpca_pix_rbf <- kernpca(data.matrix(pix_tr[,-1]), kernel = "rbfdot", par = list(sigma=0.01))
plot(kpca_pix_rbf$eigenvalues)

plot(kpca_pix_rbf$projected.training.data[,1:2], col = pix_tr[,1])

# kPCA (linear kernel -> PCA)
kpca_pix_lin <- kernpca(data.matrix(pix_tr[,-1]), kernel = "vanilladot", par = list())
plot(kpca_pix_lin$eigenvalues)

plot(kpca_pix_lin$projected.training.data[,1:2], col = pix_tr[,1])

# kPCA (Laplace kernel)
kpca_pix_lap <- kernpca(data.matrix(pix_tr[,-1]), kernel = "laplacedot", par = list(sigma = 0.01))
plot(kpca_pix_lap$eigenvalues)

plot(kpca_pix_lap$projected.training.data[,1:2], col = pix_tr[,1])


# ------- Train a linear predictor on all kernel principal components (Laplace kernel)
lap_svm <- ksvm(pix_tr[,1]~., data = kpca_pix_lap$projected.training.data, kernel = "vanilladot", par = list())

pred_lap <- predict(object = lap_svm, newdata = kpca_pix_lap$projected.training.data)
measures(pix_tr[,1], pred_lap) # training error

# project test points into the feature pcs
test_proj_lap <- project(data.matrix(pix_te[,-1]), kpca_pix_lap)

# and predict
pred_lap <- predict(object = lap_svm, newdata = test_proj_lap)
measures(pix_te[,1], pred_lap) # test error

# ------- We now try the same for linear pca
lin_svm <- ksvm(pix_tr[,1]~., data = kpca_pix_lin$projected.training.data, kernel = "vanilladot", par = list())

pred_lin <- predict(object = lin_svm, newdata = kpca_pix_lin$projected.training.data)
measures(pix_tr[,1], pred_lin) # training error

# project test points into the feature pcs
test_proj_lin <- project(data.matrix(pix_te[,-1]), kpca_pix_lin)

# and predict
pred_lin <- predict(object = lin_svm, newdata = test_proj_lin)
measures(pix_te[,1], pred_lin)

# ------- with Gaussian RBF
rbf_svm <- ksvm(pix_tr[,1]~., data = kpca_pix_rbf$projected.training.data, kernel = "vanilladot", par = list())

pred_rbf <- predict(object = rbf_svm, newdata = kpca_pix_rbf$projected.training.data)
measures(pix_tr[,1], pred_rbf) # training error.

# project test points into the feature pcs
test_proj_rbf <- project(data.matrix(pix_te[,-1]), kpca_pix_rbf)

# and predict.
pred_rbf <- predict(object = rbf_svm, newdata = test_proj_rbf)
measures(pix_te[,1], pred_rbf) # test error.

# ------- and with nonlinear SVM directly on the raw data
nonlin_svm <- ksvm(pix_tr[,1]~., data = pix_tr[,-1], kernel = "laplacedot", par = list(sigma = 0.01))

# training error.
pred_svm <- predict(object = nonlin_svm, newdata = pix_tr[,-1])
measures(pix_tr[,1], pred_svm)

# test error.
pred_svm <- predict(object = nonlin_svm, newdata = pix_te[,-1])
measures(pix_te[,1], pred_svm)


# let's visualize how many components of the "Laplace kernel" option beat linear pca in test:
dims <- seq(5, 100, 5)
lap_err <- array(dim = 20)
for (i in 1:dim(lap_err)) {
  lap_svm <- ksvm(pix_tr[,1]~., data = kpca_pix_lap$projected.training.data[,1:dims[i]], kernel = "vanilladot", par = list())
  test_proj_lap <- project(data.matrix(pix_te[,-1]), kpca_pix_lap, ncomp = dims[i])
  pred_lap <- predict(object = lap_svm, newdata = test_proj_lap)
  lap_err[i] <- measures(pix_te[,1], pred_lap)
}

dims_less <- seq(1, 14)
lap_err_less <- array(dim = 14)
lin_err_less <- array(dim = 14)
for (i in 1:dim(lap_err_less)) {
  # a bug in ksvm
  if (i == 1) {
    lap_svm <- ksvm(pix_tr[,1]~PC1, data = kpca_pix_lap$projected.training.data[,1:2], kernel = "vanilladot", par = list())
    test_proj_lap <- project(data.matrix(pix_te[,-1]), kpca_pix_lap, ncomp = dims_less[i])
    pred_lap <- predict(object = lap_svm, newdata = test_proj_lap)
    lap_err_less[i] <- measures(pix_te[,1], pred_lap)
  
    lin_svm <- ksvm(pix_tr[,1]~PC1, data = kpca_pix_lin$projected.training.data[,1:2], kernel = "vanilladot", par = list())
    test_proj_lin <- project(data.matrix(pix_te[,-1]), kpca_pix_lin, ncomp = dims_less[i])
    pred_lin <- predict(object = lin_svm, newdata = test_proj_lin)
    lin_err_less[i] <- measures(pix_te[,1], pred_lin)
  } else {
    lap_svm <- ksvm(pix_tr[,1]~., data = kpca_pix_lap$projected.training.data[,1:dims_less[i]], kernel = "vanilladot", par = list())
    test_proj_lap <- project(data.matrix(pix_te[,-1]), kpca_pix_lap, ncomp = dims_less[i])
    pred_lap <- predict(object = lap_svm, newdata = test_proj_lap)
    lap_err_less[i] <- measures(pix_te[,1], pred_lap)
    
    lin_svm <- ksvm(pix_tr[,1]~., data = kpca_pix_lin$projected.training.data[,1:dims_less[i]], kernel = "vanilladot", par = list())
    test_proj_lin <- project(data.matrix(pix_te[,-1]), kpca_pix_lin, ncomp = dims_less[i])
    pred_lin <- predict(object = lin_svm, newdata = test_proj_lin)
    lin_err_less[i] <- measures(pix_te[,1], pred_lin)
  }
}

# and plot it
plot_ly(
x = dims_less, 
y = lap_err_less, 
type = "scatter",
mode = "lines") %>%
add_trace(
  x = dims_less, 
  y = lin_err_less, 
  type = "scatter",
  mode = "lines")

plot_ly(x = c(dims_less[1:5], dims), y = c(lap_err_less[1:5], lap_err), type = "scatter", mode = "lines") %>%
add_trace(x = dims_less, y = lin_err_less, type = "scatter", mode = "lines")
    
# Summary of arguments, for details, see man/optspace.Rd
#
# x: an \eqn{(n\times m)} matrix whose missing entries should be flagged as NA.
#
# ropt: pre-defined rank, positive integer (default: 3); or logical (FALSE to guess the rank)
#
# niter: maximum number of iterations allowed
#
# tol: Stopping criterion for reconstruction in Frobenius norm.
#
# verbose: a logical value; \code{TRUE} to show progress, \code{FALSE} otherwise.
`optspace`  <-
  function(x, ropt = 3, niter = 5, tol = 1e-5, verbose = FALSE)
{

  ## Preprocessing : x : partially revealed matrix
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  idxna <- is.na(x)

  if (!is.matrix(x)) {
    stop("* optspace : input 'x' should be a matrix")
  }
  if (!any(idxna)) {
    x # no NAs to be filled in and can be returned immediately
  }  
  if (any(is.infinite(x))) {
    stop("* optspace : infinite values are not allowed in 'x'")
  }
  
  m_e <- array(0, c(nrow(x), ncol(x)))
  m_e[!idxna] <- x[!idxna]
  
  ## Preprocessing : size information
  n <- nrow(x)
  m <- ncol(x)
  
  ## Preprocessing : other sparse-related concepts
  nnz_e <- sum(!idxna)
  E <- array(0, c(nrow(x), ncol(x)))
  E[!idxna] <- 1
  eps <- nnz_e / sqrt(m * n)
  
  ## Preprocessing : ropt  : implied rank
  if (ropt) {
    r <- round(ropt)
    if (!is.numeric(ropt) || (!is.numeric(r)) || (r < 1) || (r > m) || (r > n)) {
      stop("* optspace: value of argument 'ropt' should be integer
            in [1, min(nrow(x), ncol(x))]")
    }
  } else {
    r <- min(max(round(.guess_rank(m_e, nnz_e)), 2), m - 1)
    if (verbose) {
      message(paste0("* optspace: Guessing an implicit rank: Estimated rank 'ropt': ", r))
    }
  }

  ## Preprocessing : niter : maximum number of iterations
  if ((is.infinite(niter)) || (niter <= 1) || (!is.numeric(niter))) {
    stop("* optspace: invalid number provided for argument 'niter'")
  }
  niter <- round(niter)
  rho <-  eps * n
  
  ## Main Computation
  rescal_param <- sqrt(nnz_e * r / (norm(m_e, 'f')^2))
  m_e <- m_e * rescal_param

  # 1. SVD
  if (verbose) {
    message("* optspace: Step 2: SVD ...")
  }

  svdEt <- svd(m_e)
  X0 <- matrix(svdEt$u[, seq_len(r)], ncol=r)
  X0 <- matrix(X0[, rev(seq_len(ncol(X0)))], ncol=r)
  S0 <- diag(rev(svdEt$d[seq_len(r)]))
  Y0 <- matrix(svdEt$v[, seq_len(r)], ncol=r)
  Y0 <- matrix(Y0[, rev(seq_len(ncol(Y0)))], ncol=r)

  # 3. Initial Guess
  if (verbose) {
    message("* optspace: Step 3: Initial Guess ...")
  }
  X0 <- X0 * sqrt(n)
  Y0 <- Y0 * sqrt(m)
  S0 <- S0 / eps

  # 4. Gradient Descent
  if (verbose) {
    message("* optspace: Step 4: Gradient Descent ...")
  }
  X <- X0
  Y <- Y0
  S <- .aux_getoptS(X, Y, m_e, E)
  # initialize
  dist <- array(0, c(1, (niter + 1)))
  dist[1] <- norm((m_e - (X %*% S %*% t(Y))) * E, 'f') / sqrt(nnz_e)
  # Resolution/regularization for gradient calculation
  m0 <- 10000
  
  for (i in seq_len(niter)) {
    # compute the gradient
    tmpgrad <- .aux_gradF_t(X, Y, S, m_e, E, m0, rho)
    W <- tmpgrad$W
    Z <- tmpgrad$Z
    # line search for the optimum jump length
    t <- .aux_getoptT(X, W, Y, Z, S, m_e, E, m0, rho)
    X <- X + t * W
    Y <- Y + t * Z
    S <- .aux_getoptS(X, Y, m_e, E)
    # compute the distortion
    dist[i + 1] <- norm(((m_e - X %*% S %*% t(Y)) * E), 'f') / sqrt(nnz_e)    
    if (dist[i + 1] < tol) {
      dist <- dist[seq_len(i + 1)]
      break
    }
  }
  S <- S / rescal_param  
  # Return Results
  out <- list()

  # re-order Optspace may change order during iters
  index_order <- order(diag(S), decreasing = TRUE)
  X <- matrix(X[, index_order], ncol=length(index_order))
  Y <- matrix(Y[, index_order], ncol=length(index_order))
  S <- matrix(S[index_order, index_order], ncol=length(index_order))
  out$X <- X
  out$S <- S
  out$Y <- Y
  out$dist <- dist
  if (verbose) {
    message('* optspace: estimation finished.')
  }

  # -------------------------------------------

  # This part is not in the Python / Gemelli implementation
  # but has been added in R to provide more direct access
  # to the imputed matrix.

  # Reconstruct the matrix
  M <- X %*% S %*% t(Y)

  # Centering is common operation supporting output visualization
  # Center cols to 0
  M <- as.matrix(scale(M, center = TRUE, scale = FALSE))
  # Center rows to 0
  M <- as.matrix(t(scale(t(M), center = TRUE, scale = FALSE)))

  # Add imputed matrix to the output
  out$M <- M

  # -------------------------------------------

  out
}


# Estimate Matrix Rank 
#
# x: numeric matrix for which the rank is to be estimated
#
# nnz: estimated number of non-zero entries to derive noise threshold
#
# keywords internal
.guess_rank <- function(x, nnz)
{
  maxiter <- 10000
  n <- nrow(x)
  m <- ncol(x)
  epsilon <- nnz / sqrt(m * n)
  svdX <- svd(x)
  S0 <- svdX$d
  
  nsval0 <- length(S0)
  S1 <- S0[seq_len(nsval0 - 1)] - S0[seq(2, nsval0)]  
  nsval1 <- length(S1)
  if (nsval1 > 10) {
    S1_ <- S1 / mean(S1[seq((nsval1 - 10), nsval1)])
  } else {
    S1_ <- S1 / mean(S1[seq_len(nsval1)])
  }
  r1 <- 0
  lam <- 0.05
  
  itcounter <- 0
  while (r1 <= 0) {
    itcounter <- itcounter + 1
    cost <- array(0, c(1, length(S1_)))
    for (idx in seq_len(length(S1_))) {
      cost[idx] <- lam * max(S1_[seq(idx, length(S1_))]) + idx
    }
    v2 <- min(cost)
    i2 <- which(cost == v2)
    if (length(i2) == 1) {
      r1 <- i2 - 1
    } else {
      r1 <- max(i2) - 1
    }
    lam <- lam + 0.05
    if (itcounter > maxiter) {
      break
    }
  }
  
  if (itcounter <= maxiter) {
    cost2 <- array(0, c(1, (length(S0) - 1)))
    for (idx in seq_len(length(S0) - 1)) {
      cost2[idx] <- (S0[idx + 1] + sqrt(idx * epsilon) * S0[1] / epsilon) / S0[idx]
    }
    v2 <- min(cost2)
    i2 <- which(cost2 == v2)
    if (length(i2) == 1) {
      r2 <- i2
    } else {
      r2 <- max(i2)
    }
    
    if (r1 > r2) {
      r <- r1
    } else {
      r <- r2
    }
  } else {
    r <- min(nrow(x), ncol(x))
  }
  r
}


# Auxiliary Gradient Contribution Function
#
# Compute distortion. Computes nonlinear transformation of the squared
# row norms of a matrix, with thresholding based on scaled values.
# Typically used as part of a gradient computation or optimization routine,
# where it selectively activates rows based on their scaled magnitude.
#
# x: numeric matrix; the function computes the squared norm of each row.
#
# m0 positive scalar that acts as a scaling parameter in the nonlinear regularization term
#
# r: numeric scalar; another scaling factor applied alongside m0 in the normalization of row norms.
#' Its effect is similar to m0, the two are often coupled in optimization problems.
#
# keywords internal
.aux_G <- function(x, m0, r)
{
  z <- rowSums(x^2) / (2 * m0 * r)
  y <- exp((z - 1)^2) - 1
  idxfind <- (z < 1)
  y[idxfind] <- 0
  out <- sum(y)
  out
}

# Total Loss Function with Regularization
#
# Composite loss function combining a weighted squared error term and 
# nonlinear regularization penalties based on the row norms of the input matrices. 
# This function is commonly used in matrix factorization or low-rank modeling 
# frameworks where latent matrices are learned under structural constraints.
#
# x numeric matrix, typically representing latent factors or components.
#
# y numeric matrix, typically representing latent factors or components.
#
# s numeric matrix used for linear transformation between \code{x} and \code{y}.
#
# m_e numeric matrix of size representing observed or target values to approximate.
#
# e numeric matrix of the same dimension as \code{m_e}, used as a weight or mask matrix. Values typically range from 0 to 1.
#
# m0 positive scalar that acts as a scaling parameter in the nonlinear regularization term
#
# rho positive scalar controlling the strength of the regularization terms
#
# keywords internal
.aux_F_t <- function(x, y, s, m_e, e, m0, rho)
{
  n <- nrow(x)
  r <- ncol(x)
  out1 <- sum((((x %*% s %*% t(y)) - m_e) * e)^2) / 2
  out2 <- rho * .aux_G(y, m0, r)
  out3 <- rho * .aux_G(x, m0, r)
  out  <- out1 + out2 + out3
  out
}



# Gradient of Auxiliary Regularization Term
#
# Computes the gradient of the nonlinear regularization function \code{.aux_G} 
# with respect to its matrix input \code{x}. 
#
# x numeric matrix of size \code{n x r}. Each row is treated as a vector
# whose squared norm determines the activation of the gradient.
#
# m0 positive scalar that acts as a scaling parameter in the nonlinear regularization term
#
# r: numeric scalar; another scaling factor applied alongside m0 in the normalization of row norms.
# Its effect is similar to m0, the two are often coupled in optimization problems.
#
# keywords internal
.aux_Gp <- function(x, m0, r)
{
  z <- rowSums(x^2) / (2 * m0 * r)
  z <- 2 * exp((z - 1)^2) / (z - 1)
  idxfind <- (z < 0)
  z[idxfind] <- 0  
  out <- x * matrix(z, nrow = nrow(x), ncol = ncol(x), byrow = FALSE) / (m0 * r)
}


# Gradient of Composite Loss Function
#
# Computes the gradient of the total loss. The result includes both the data 
# reconstruction gradient and the nonlinear regularization gradient.
#
# x numeric matrix of size \code{n x r}, representing one set of latent 
# variables or factor matrix.
#
# y numeric matrix of size \code{m x r}, representing the counterpart latent 
# variable matrix.
#
# s numeric matrix used in the bilinear transformation between \code{x} and \code{y}.
#
# m_e A numeric matrix representing the observed data or target matrix 
#
# e numeric matrix used as a weighting or masking matrix.
#
# m0 positive scalar for gradient regularization 
#
# rho weight applied to the regularization gradient terms.
#
# @keywords internal
.aux_gradF_t <- function(x, y, s, m_e, e, m0, rho)
{
  n <- nrow(x)
  r <- ncol(x)
  m <- nrow(y)
  if (ncol(y) != r) {
    stop("dimension error from the internal function .aux_gradF_t")
  }
  
  XS  <- x %*% s
  YS  <- y %*% t(s)
  XSY <- XS %*% t(y)
  
  Qx <- t(x) %*% ((m_e - XSY) * e) %*% YS / n
  Qy <- t(y) %*% t((m_e - XSY) * e) %*% XS / m
  
  W <- ((XSY - m_e) * e) %*% YS  + (x %*% Qx) + rho * .aux_Gp(x, m0, r)
  Z <- t((XSY - m_e) * e) %*% XS + (y %*% Qy) + rho * .aux_Gp(y, m0, r)
  
  resgrad <- list()
  resgrad$W <- W
  resgrad$Z <- Z
  resgrad
  
}


# Solve for Optimal Transformation Matrix S
#
# Computes the optimal transformation matrix that minimizes the squared 
# error, optionally weighted by importance matrix \code{e}.
# This function forms and solves a linear least-squares problem.
#
# x A numeric matrix representing one side of the bilinear transformation.
#
# y A numeric matrix representing the other side of the transformation.
#
# m_e A numeric matrix representing the observed data or target matrix to approximate.
#
# e numeric matrix used as a weighting or masking matrix.
#
# # @keywords internal
.aux_getoptS <- function(x, y, m_e, e)
{
  n <- nrow(x)
  r <- ncol(x)  
  C <- t(x) %*% (m_e) %*% y
  C <- matrix(as.vector(C))  
  nnrow <- ncol(x) * ncol(y)
  A <- matrix(NA, nrow = nnrow, ncol = (r^2))
  
  for (i in seq_len(r)) {
    for (j in seq_len(r)) {
      ind <- (j - 1) * r + i
      tmp <- t(x) %*% (outer(x[, i], y[, j]) * e) %*% y      
      A[, ind] <- as.vector(tmp)
    }
  }
  
  S <- solve(A, C)
  out <- matrix(S, nrow = r)
  out
}



# Backtracking Line Search for Step Size Optimization
#
# Optimal line search.
# Determines an optimal step size \code{t} for descending along a direction 
# defined by perturbations \code{w} and \code{z} for the matrices \code{x} and \code{y}. 
#
# x numeric matrix representing current values of one latent factor.
#
# w numeric matrix representing the search direction (gradient or descent vector) for \code{x}.
#
# y numeric matrix representing current values of the second latent factor.
#
# z numeric matrix representing the search direction for \code{y}.
#
# s numeric matrix used in the bilinear product in \code{.aux_F_t}.
#
# m_e numeric the observed or target matrix.
#
# e numeric mask or weighting matrix.
#
# m0 positive scalar regularization parameter.
#
# rho positive scalar weighting the regularization term in the total loss function.
#
# keywords internal
.aux_getoptT <- function(x, w, y, z, s, m_e, e, m0, rho)
{
  norm2WZ <- norm(w, 'f')^2 + norm(z, 'f')^2
  f <- array(0, c(1, 21))
  f[1] <- .aux_F_t(x, y, s, m_e, e, m0, rho)
  t <- -1e-1
  for (i in seq_len(21)) {
    f[i + 1] <- .aux_F_t(x + t * w, y + t * z, s, m_e, e, m0, rho)
    if ((f[i + 1] - f[1]) <= 0.5 * t * norm2WZ) {
      out <- t
      break
    }
    t <- t / 2
  }
  out <- t
  t
}


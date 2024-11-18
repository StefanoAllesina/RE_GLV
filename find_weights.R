source("check_copositivity.R") # edit as needed 

# this takes a matrix A # and the Nash equilibrium p
# and transforms it into a matrix B with 
# * Nash equilibrium 1/n for all components
# * It is defined by three coefficients (a,b,c) as well as either k = {1,-1,0}
abck_form <- function(A, p){
  n <- nrow(A)
  B1 <- A %*% diag(p)
  B2 <- B1 - outer(rep(1, n), diag(B1))
  k <- rowSums(B2)[1]
  if (k != 0) B2 <- B2 / abs(k)
  # return the matrix as well as the values a, b, c, k
  return(list(B = B2, a = B2[1,2], b = B2[1,3], c = B2[3,1], k = rowSums(B2)[1]))
}

# Form matrix H(Btilde) given a matrix in abck form
# and weights summing to unity
get_HBtilde <- function(w, B, n){
  # form B tilde
  Btilde <- (outer(rep(1, n), w) -diag(n)) %*% B %*% diag(1/w)
  # and its symmetric part
  HBtilde <- (Btilde + t(Btilde)) / 2
  return(HBtilde)
}

# Try to find weights numerically
find_weights <- function(abck, 
                         npoints = 100,
                         method = "Nelder-Mead",
                         maxit = 5000){
  to_minimize <- function(pars, B, n){
    w <- exp(pars) # make sure w are nonnegative
    # normalize to keep them in the simplex
    w <- w / sum(w)
    HBtilde <- get_HBtilde(w, B, n)
    # compute x^T HBtilde x for each x in X
    quadratic_form <- rowSums((X %*% HBtilde) * X)
    return(-sum(quadratic_form[quadratic_form < 0]))
  }
  B <- abck$B
  if (abck$k != 1) {
    print("The equilibrium cannot be globally stable, because it is not locally stable")
    return(NULL)
  }
  n <- nrow(B)
  # take values in the simplex
  X <- matrix(rexp(npoints * n), npoints, n)
  X <- X / rowSums(X)
  # start with parameters wi = 1/n for all i
  # choose weights such that they are in the simplex
  # try to maximize the minimum of the quadratic form
  tmp <- optim(par = rep(0, n), 
               fn = to_minimize, 
               method = method,
               control = list(maxit = maxit),
               B = B, n = n)
  # check whether we have found a solution
  if (tmp$value == 0){
    print("Found a solution")
    w <- exp(tmp$par)
    w <- w / sum(w)
    print(w)
    return(w)
  }
  print("Could not find a solution")
  return(NULL)
}

# Example from Taylor and Yonker 1978
TY_A <- matrix(c(
  2, 1, 5, 
  5, 1, 0, 
  1, 4, 3),
  3, 3, byrow = TRUE)

TY_p <- c(15, 11, 9) / 35

abck <- abck_form(TY_A, TY_p)
w <- find_weights(abck, npoints = 1000)
if (!is.null(w)){
  # verify the solution
  print("Verifying the solution")
  HBtilde <- get_HBtilde(w, abck$B, nrow(abck$B))
  print(check_copositivity(HBtilde))
}
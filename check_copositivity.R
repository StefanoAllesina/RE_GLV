# These functions test whether a given matrix is copositive
# They follow the approach detailed in
#
# J.-B. Hiriart-Urruty & A. Seeger
# A Variational Approach to Copositive Matrices
# SIAM Review 52(4) 593-629, 2010

# In particular, Proposition 2.2 provides 
# necessary and sufficient conditions 
# for 3x3 matrices

test_3x3 <- function(A){
  # 1 by 1
  if (any(diag(A) < 0)) return(FALSE)
  # 2 by 2
  baraij <- c(
    A[1,2] + sqrt(A[1,1] * A[2,2]),
    A[1,3] + sqrt(A[1,1] * A[3,3]),
    A[2,3] + sqrt(A[2,2] * A[3,3])
  )
  if (any(baraij < 0)) return(FALSE)
  # now the last condition
  return(
    sqrt(A[1,1] * A[2,2] * A[3,3]) + 
      A[1,2] * sqrt(A[3,3]) + A[1,3] * sqrt(A[2,2]) + A[2,3] * sqrt(A[1,1]) + 
      sqrt(2 * prod(baraij)) >= 0
  )
}

test_3x3_det0 <- function(A){
  # 1 by 1
  if (any(diag(A) < 0)) return(FALSE)
  # 2 by 2
  baraij <- c(
    A[1,2] + sqrt(A[1,1] * A[2,2]),
    A[1,3] + sqrt(A[1,1] * A[3,3]),
    A[2,3] + sqrt(A[2,2] * A[3,3])
  )
  if (any(baraij < 0)) return(FALSE)
  return(TRUE)
}

# For larger matrices, it uses the approach
# by Cottle–Habetler–Lemke, in Theorem 3.3
test_larger <- function(A){
  # a matrix with copositive principal 
  # submatrices is copositive if its determinant
  # is nonnegative
  if (det(A) >= 0) return(TRUE)
  # or the adjugate matrix contains a negative entry
  B <- det(A) * solve(A)
  if (any(B < 0)) return(TRUE)
  return(FALSE)
}

# This function tests whether x^T A x >= 0 for every x>0
# where the determinant of A is zero
# The matrix has size at least n = 3 (because this is the first interesting case)
check_copositivity <- function(A){
  # the tests are meant for symmetric matrices---take the symmetric part
  A <- (A + t(A)) / 2
  n <- nrow(A)
  # it assumes that the determinant is zero
  if (n == 3) return(test_3x3_det0(A))
  # Start by checking that all 3x3 submatrices are copositive
  all_3 <- combn(1:n, 3)
  test_3 <- apply(all_3, 2, function(x) test_3x3(A[x,x]))
  if (any(test_3 == FALSE)) return(FALSE)
  # now move to larger sizes
  if(n > 4) for (k in 4:(n-1)){
    all_k <- combn(1:n, k)
    test_k <- apply(all_k, 2, function(x) test_larger(A[x,x]))  
    if (any(test_k == FALSE)) return(FALSE)
  }
  return(TRUE)
}

# To show an interesting example, consider Horn's matrix
# a 5x5 matrix that is copositive but cannot be written
# as the sum of a PSD matrix and a matrix with nonnegative entries
test_horn <- function(){
  A <- matrix(c(
    1,-1,1,1,-1,
    -1,1,-1,1,1,
    1,-1,1,-1,1,
    1,1,-1,1,-1,
    -1,1,1,-1,1), 5, 5, byrow = TRUE)
  return(list(A, check_copositivity(A)))
}



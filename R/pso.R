#' Particle Swarm Optimizer
#'
#' This particle swarm optimizer takes in a function and outputs
#' the optimal parameter values. The function must be formatted in such
#' a way that it takes a single vector input for all parameters (i.e.
#' VEC = c(1,2,3) for x=1, y=2, z=3).
#'
#' @param fn Function to be minimized
#' @param S Number of particles in swarm
#' @param lb Lower bounds of inputs
#' @param ub Upper bounds of inputs
#' @param w Weight coefficient
#' @param max_iter Maximum number of iterations
#' @return List containing the optimal parameters, optimal function value, and number of iterations
#' @export
pso <- function(fn, S, lb, ub, w=.1, max_iter = 1000){

  d = length(lb) # No. dimensions
  X = matrix(0, nrow = d, ncol = S)
  V = matrix(0, nrow = d, ncol = S)
  ranges = ub - lb
  for (i in 1:d) {
    # Initilaize particle position matrix, x
    X[i,] = runif(S, lb[i], ub[i])
    # Initialize each particle's velocity
    V[i,] = runif(S, -ranges[i], ranges[i])
  }

  # Initialize matrix to track each particle's best known point
  P = X

  # Calculate swarm's best known position
  f_x = apply(X, MARGIN = 2, FUN = fn)
  g = X[,which.min(f_x)]

  iter = 0
  # Add convergence criteria
  while(iter < max_iter){

    # Update velocities
    # V = w*V + r_p*(P-X) - r_g * -t(t(X)-g)
    for(i in 1:S){
      for(j in 1:d){
        # Choose random r values
        r_p = runif(1)
        r_g = runif(1)
        V[j,i] = w*V[j,i] + r_p*(P[j,i]-X[j,i]) - r_g*(g[j] - X[j,i])
      }
    }

    # Update particles' position
    X = X + V

    # Update swarm's best position
    f_x = apply(X, MARGIN = 2, FUN = fn)
    g = X[,which.min(f_x)]

    # Count iterations
    iter = iter + 1
  }

  # Add counts vector
  return(list(par = g, value = fn(g), iterations = iter))
}


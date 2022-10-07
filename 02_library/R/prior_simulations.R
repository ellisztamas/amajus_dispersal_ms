#' PDF of the generalised normal distribution.
#'
#' @param x Vector of deviates (from zero).
#' @param scale: Scale parameter.
#' @param shape: Shape parameter
#'
#' @return Vector of probabilities
dgennorm <- function(x, scale, shape, b){
  (shape/(2 * scale * gamma(1/shape) )) * exp(-(x/scale)^shape)
}

#' Samples from the generalised normal distribution
#' 
#' Draw samples from a vector deviates in proportion to the generalised normal
#' distribution
#' 
#' Since there is no built-in function to draw numbers from the generalised 
#' normal distribution (GND), and this is non-trivial to implement, this instead 
#' takes a vector of input distances, and draws samples with replacement from
#' these distances in proportion to their probability under the GND. If values
#' in the input are evenly spaced and n is large enough, this will approximate
#' drawing from a continuous number line (but neither assumption is mandatory).
#' 
#' @param x Vector of deviates (from zero)
#' @param scale Scale parameter
#' @param shape Shape parameter
#' @param n Number of samples to draw
#' 
#' @return Vector of size n of samples (with replacement) from x.
#' 
#' @author Tom Ellis
#' @seealso dgennorm
rgennorm <- function(x, scale, shape, n){
  p <- dgennorm(x, scale, shape)
  sample(x, size = n, prob = p, replace = TRUE)
}

#' Summarise prior simulation
#' 
#' Helper function to summarise the output of a vector.
#' 
#' @param draws Vector of prior draws.
#' @return Vector giving the mean, median, proportion of values > 500, the
#' proportion of values > 1000 and maximum value in `draws`.
summarise_draws <- function(draws){
  c(
    mean(draws),
    median(draws),
    mean(draws > 500),
    mean(draws > 1000),
    max(draws)
  )
}

#' Prior draws from the generalised normal
#' 
#' Simulate values from vectors of input parameters to the generalised normal
#' distribution (GND) and return a summary of each parameter set.
#' 
#' @param x Vector of deviates (positive values representing distances travelled)
#' @param scale Vector of input values for the scale parameter of the GND. 
#' Should be the same length as `shape`.
#' @param shape Vector of input values for the shape parameter of the GND.
#' Should be the same length as `scale`.
#' @param n Number of values to draw for each set of input values.
prior_sim <- function(x, scale, shape, n){
  if( length(scale) != length(shape) ){
    stop("scale and lenhgth are not the same length.")
  }
  # Matrix of simulated dispersal distance
  # Each column is for a separate simulation pair of parameter values
  draws <- mapply(rgennorm, scale=scale, shape=shape, n=n, MoreArgs = list(x=x))
  # Summarise each simulation
  sim_table <- apply(draws, 2, summarise_draws)
  sim_table <- as.data.frame(t(sim_table))
  names(sim_table) <- c("mean", "median", "geq500", "geq1000", "max")
  
  return(
    list(draws = draws, summary = sim_table)
  )
}
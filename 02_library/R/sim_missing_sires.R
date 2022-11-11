#' Estimate the number of unsampled fathers
#' 
#' This uses information from a table of mating events with identifiable fathers,
#' the sizes of each paternal family and the overall maternal family size to
#' predict how many families are present within groups of offspring whose
#' father could not be identified.
#' 
#' This uses a Poisson GLM to predict (log) number of paternal families as a
#' function of maternal family size. This model must go through zero when the
#' total maternal array size is zero, otherwise estimates are inflated by many 
#' apparent mating events with one offspring but poor support. The proportion of
#' mating events, p, from unsampled fathers is then the pedicted number divided
#' by total array size.
#' 
#' In a second step, we use p to draw a binomially distributed number of paternal
#' families from unsampled fathers. This is repeated for each iteration from 
#' the mating events file, and a single value is drawn for each.
#' 
#' @param mating_events A table listing mother, father (with unsampled fathers
#' given as NA), posterior probability for the mating event, weighted-mean 
#' number of offspring in the family, total maternal array size and iteration in
#' the MCMC. This is given by the output of `03_analysis/04_mating_events/`.
#' @return For each iteration in the input file, a data.frame giving simulated
#' values for the number of paternal families with an unsampled father, and the
#' regression coefficient
#' @author Tom Ellis
sim_missing_sires <- function(mating_events){
  
  # Labels for iterations
  ix <- unique(mating_events$iter)
  # Vectors to hold the results
  out <- data.frame(
    n_mating_events   = rep(NA, length(ix)),
    sim_missing_sires = rep(NA, length(ix)),
    sim_coef          = rep(NA, length(ix))
  )
  
  pb <- txtProgressBar(0, length(ix), style = 3)
  # Loop over results for each iteration in the input file
  for(i in 1:length(ix)){
    setTxtProgressBar(pb, i)
    # Table of mating events for which the father COULD be identified.
    families_with_sires <- mating_events %>% 
      # Keep only mating events with an unsampled father, and with offspring number
      # distinguishable from zero, because this causes NAs.
      # (offspring number is weighted by the probability of having occurred, so
      # this may be less than 1.)
      filter(!is.na(father), iter == ix[i], offspring > 0) %>% 
      # Get the weight mean number of mating events for each mother
      group_by(mother, array_size) %>% 
      summarise(
        n = sum(prob),
        n_plus_1 = n + 1,
        .groups = "drop_last"
      ) %>% 
      ungroup()
    
    # Poisson GLM for the number of mating events as a function of the number of 
    # progeny in the maternal family.
    # Notice the -1, to remove the intercept, meaning the curve will go through 
    # 0 on the link function scale, and 1 on the response scale
    mod <- glm(n_plus_1 ~ -1 + array_size, data = families_with_sires, family = poisson)
    
    # A table of groups of progeny in each maternal family for whom no father was
    # found.
    # Then add the expected number of full sibships expected within each group from
    # the GLM above.
    missing <- mating_events %>%
      # Keep only mating events with an unsampled father, and with offspring number
      # distinguishable from zero, because this causes NAs.
      # (offspring number is weighted by the probability of having occurred, so
      # this may be less than 1.)
      filter( is.na(father), offspring > 0, iter == ix[i]) %>%
      # select(mother, prob, offspring) %>% 
      # mutate(array_size = offspring) %>% 
      # Get predicted number of families from the GLM
      # The -1 corrects for the +1 correction in the prediction model
      mutate(
        offspring = round(offspring), # We need integers to draw binomial draws from
        pred = predict(object = mod, newdata = ., type = "response") -1,
        n_obs = 1,
        p = pred / array_size
      )
    
    # Simulate a number of missing fathers and send it to `sim_missing_sires`
    out$sim_missing_sires[i] <- sum(
      rbinom(
        n = missing$n_obs, size = round(missing$offspring), prob = missing$p
      )
    )
    # Save the regression coefficient
    out$sim_coef[i] <- coef(mod)
    # Save the number of mating events with identifiable fathers
    out$n_mating_events[i] <- sum( families_with_sires$n )
  }
  close(pb)
  
  # Return a dataframe of results
  out
}

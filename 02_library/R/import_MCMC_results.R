#' Import MCMC results from independent chains
#' 
#' Import and concatenate MCMC results from independent chains from within a 
#' single folder. This is is usually the output of the Python function 
#' `amajusdispersal.mcmc.run_MCMC()`.
#' 
#' @param dir Directory where MCMC results can be found. Files should be
#' labelled as chains 1 to 4 with the suffix `.out`. The directory should 
#' contain only the results of a single analysis.
#' @return A tibble containing the full output of each chain (no trimming for
#' burnin iterations is done).
import_mcmc_results <- function(dir){
  # cat(paste("Importing MCMC results from", dir, "\n"))
  list(
    read_tsv(paste0(dir, "/output/chain1.out"), col_types = 'innnnnnnnf') %>% 
      mutate(chain = "1"),
    read_tsv(paste0(dir, "/output/chain2.out"), col_types = 'innnnnnnnf') %>% 
      mutate(chain = "2"),
    read_tsv(paste0(dir, "/output/chain3.out"), col_types = 'innnnnnnnf') %>% 
      mutate(chain = "3"),
    read_tsv(paste0(dir, "/output/chain4.out"), col_types = 'innnnnnnnf') %>% 
      mutate(chain = "4")
  ) %>% 
    do.call(what="rbind") %>% 
    mutate(scenario = as.factor(basename(dir)))
}

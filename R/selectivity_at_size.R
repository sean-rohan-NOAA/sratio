#' Selectivity at size
#' 
#' Selectivity-at-size for dome-shaped or asymptotic selectivity functions.
#' 
#' @param size Vector of fork sizes
#' @param selectivity_opts List containing parameters for selectivity functions. Current support for double normal (pattern 24 in SS3) or logistic (pattern 1 in SS3).
#' @export

selectivity_at_size <- function(size, 
                           selectivity_opts) {
  
  if(selectivity_opts$type == "doublenormal") {
    #Double normal selectivity
    beta1 <- selectivity_opts$begin_top
    beta2 <- inv_logit(selectivity_opts$logit_top)
    beta3 <- selectivity_opts$ln_sd1
    beta4 <- selectivity_opts$ln_sd2
    beta5 <- inv_logit(selectivity_opts$logit_start)
    beta6 <- inv_logit(selectivity_opts$logit_end)
    sel <- doubleNorm24.fn(x = size, a=beta1, b=beta2, c=beta3, d=beta4, e=beta5, f=beta6)
  } else if(selectivity_opts$type == "asymptotic") {
    # Logistic selectivity
    sel <- (1 + exp(-log(19) * (size - selectivity_opts$begin_top)/selectivity_opts$ln_sd1))^-1
  }
  
  return(sel)
  
}
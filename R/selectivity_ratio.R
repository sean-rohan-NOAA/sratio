#' Calculate selectivity catch comparison rate and selectivity ratio
#' 
#' @param size Optional numeric vector of sizes/ages
#' @param count1 Count for gear #1
#' @param count2 Count for gear #2
#' @param effort1 Effort for gear #1
#' @param effort2 Effort for gear #2
#' @param sampling_factor1 Sampling factor for gear #1 (multiplied by count to estimate total count)
#' @param sampling_factor2 Sampling factor for gear #2 (multiplied by count to estimate total count)
#' @export
#' @return A data.frame containing catch proportions for gear 1 (r1) and gear 2(2), catch comparison rate (p12), selectivity ratio between gear 1 and gear 2 (s12)

selectivity_ratio <- function(size = NULL,
                              count1, 
                              count2, 
                              effort1 = 1, 
                              effort2 = 1, 
                              sampling_factor1 = 1, 
                              sampling_factor2 = 1) {
  stopifnot(
    "selectivity_ratio: count1 and count2 are not the same length." = 
      length(count1) == length(count2)
  )
  stopifnot(
    "selectivity_ratio: effort1 must be 1L or the same length as count1." = 
              length(effort1) == 1 | length(effort1) == length(count1)
    )
  stopifnot(
    "selectivity_ratio: effort2 must be 1L or the same length as count2." = 
              length(effort2) == 1 | length(effort2) == length(count2)
    )
  stopifnot(
    "selectivity_ratio: sampling_factor1 must be 1L or the same length as count1." = 
              length(sampling_factor1) == 1 | length(sampling_factor1) == length(count1)
    )
  stopifnot(
    "selectivity_ratio: sampling_factor2 must be 1L or the same length as count2." = 
              length(sampling_factor2) == 1 | length(sampling_factor2) == length(count2)
    )
  
  expanded_count1 <- count1*sampling_factor1/effort1
  
  expanded_count2 <- count2*sampling_factor2/effort2
  
  r1 <- expanded_count1/sum(expanded_count1)
  
  r2 <- expanded_count2/sum(expanded_count2)
  
  p12 <- r1/(r1+r2)
  
  s12 <- p12/(1-p12)
  
  output <- data.frame(count1 = count1,
                       count2 = count2,
                       sampling_factor1 = sampling_factor1,
                       sampling_factor2 = sampling_factor2,
                       effort1 = effort1,
                       effort2 = effort2,
                       r1 = r1, 
                       r2 = r2, 
                       p12 = p12, 
                       s12 = s12)
    
  if(!is.null(size)) {
    
    output$size <- size
    
  }

  
  return(output)
  
}

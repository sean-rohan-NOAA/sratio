#' Calculate selectivity catch comparison rate and selectivity ratio
#' 
#' @param count1 Count for gear #1
#' @param count2 Count for gear #2
#' @param effort1 Effort for gear #1
#' @param effort2 Effort for gear #2
#' @export

selectivity_ratio <- function(count1, count2, effort1 = NULL, effort2 = NULL) {
  
  if(!is.null(effort1) & !is.null(effort2)) {
    count1 <- count1/effort1
    count2 <- count2/effort2
  } else {
    relative_selectivity <- TRUE
  }
  
  
  r1 <- count1/sum(count1)
  r2 <- count2/sum(count2)
  
  
  p12 <- r1/(r1+r2)
  
  s21 <- 1/p12-1
  
  return(data.frame(r1 = r1, r2 = r2, s21 = s21, p12 = p12))
  
}

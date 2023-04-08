#' Transform [0,1] bounded variable for beta regression
#' 
#' Transform 0 and 1 values to meet conditions for beta regression.
#' 
#' @param x Numeric vector
#' @param method "constant" adds/subtracts 1e-7 from 0/1; "sv" uses the transformation from Smithson and Verkulien (2006)
#' @export
#' @references Smithson, M. & Verkuilen, J. A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychol. Methods 11, 54–71 (2006). DOI: 10.1037/1082-989X.11.1.54


scale_for_betareg <- function(x, method = "constant") {
  
  # Add/subtract a constant.
  if(method == "constant") {
    x[x == 1] <- x[x == 1] - 1e-7
    x[x == 0] <- x[x == 0] + 1e-7
  }

  # Smithson and Verkulien (2006) - http://dx.doi.org/10.1037/1082-989X.11.1.54
  if(method == "sv") {
    x <- x * ((length(x) - 1) + 0.5) / length(x)
  }

  return(x)
}

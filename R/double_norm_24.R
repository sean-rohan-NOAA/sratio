#' Double normal selectivity function
#' 
#' Selectivity pattern 24 in SS3. Copied from https://github.com/r4ss/r4ss/blob/414f60b88f404a92eb7eab323c6cc580798565ea/inst/Shiny/Selectivities/server.R on October 26, 2021.
#'
#' @param x Vector of lengths
#' @param a Parameter a
#' @param b Parameter b
#' @param c Parameter c
#' @param d Parameter d
#' @param e Parameter e
#' @param f Parameter f
#' @export

double_norm_24 <- function(x, a, b, c, d, e, f) {
  #UPDATED: - input e and f on 0 to 1 scal and transfrom to logit scale
  #         - changed bin width in peak2 calculation
  #         - updated index of sel when j2 < length(x)
  #	  - renamed input parameters, cannot have same names as the logitstic function
  #         - function not handling f < -1000 correctly
  if(e == 0) {#Avoid errors on the bounds
    e <- 1-0.999955 #an input that results in approx -10
  } 
  if(e == 1) {
    e <- 0.999955  #an input that results in approx 10
  }
  e <- log(e/(1-e)) #transform input to logit
  
  if(f == 0) {#Avoid errors on the bounds
    f <- 1-0.999955 #an input that results in approx -10
  } 
  if(f == 1) {
    f <- 0.999955  #an input that results in approx 10
  }
  f <- log(f/(1-f)) #transform input to logit
  
  sel <- rep(NA, length(x))
  startbin <- 1
  peak <- a
  upselex <- exp(c)
  downselex <- exp(d)
  final <- f
  if (e < -1000) {
    j1 <- -1001 - round(e)
    sel[1:j1] <- 1e-06
  }
  if (e >= -1000) {
    j1 <- startbin - 1
    if (e > -999) {
      point1 <- 1/(1 + exp(-e))
      t1min <- exp(-(x[startbin] - peak)^2/upselex)
    }
  }
  if (f < -1000)
    j2 <- -1000 - round(f)
  if (f >= -1000)
    j2 <- length(x)
  bin_width <- x[2] - x[1]
  peak2 <- peak + bin_width + (0.99 * x[j2] - peak - bin_width)/(1 +
                                                                   exp(-b))
  if (f > -999) {
    point2 <- 1/(1 + exp(-final))
    t2min <- exp(-(x[j2] - peak2)^2/downselex)
  }
  t1 <- x - peak
  t2 <- x - peak2
  join1 <- 1/(1 + exp(-(20/(1 + abs(t1))) * t1))
  join2 <- 1/(1 + exp(-(20/(1 + abs(t2))) * t2))
  if (e > -999)
    asc <- point1 + (1 - point1) * (exp(-t1^2/upselex) -
                                      t1min)/(1 - t1min)
  if (e <= -999)
    asc <- exp(-t1^2/upselex)
  if (f > -999)
    dsc <- 1 + (point2 - 1) * (exp(-t2^2/downselex) -
                                 1)/(t2min - 1)
  if (f <= -999)
    dsc <- exp(-(t2)^2/downselex)
  idx.seq <- (j1 + 1):j2
  sel[idx.seq] <- asc[idx.seq] * (1 - join1[idx.seq]) + join1[idx.seq] * (1 -
                                                                            join2[idx.seq] + dsc[idx.seq] * join2[idx.seq])
  if (startbin > 1 && e >= -1000) {
    sel[1:startbin] <- (x[1:startbin]/x[startbin])^2 *
      sel[startbin]
  }
  if (j2 < length(x))
    sel[(j2 + 1):length(x)] <- sel[j2]
  return(sel)
}
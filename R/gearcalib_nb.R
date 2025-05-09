#' Non-parametric relative selectivity log-Gaussian Cox process model
#'
#' Estimates the relative selectivity ratio following Thygesen et al. (2019), can use a ZIP or negative binomial instead of a Poisson.
#'
#' @param d A list containing input data with elements:
#'   \itemize{
#'     \item \code{N}: Matrix of counts (hauls x size bins).
#'     \item \code{SweptArea}: Numeric vector of swept areas per haul.
#'     \item \code{group}: Factor indicating group membership.
#'     \item \code{Gear}: Factor indicating gear type.
#'     \item \code{Lvec} (optional): Numeric vector of size bin indices (defaults to \code{1:ncol(d$N)} if missing).
#'   }
#' @param fit0 Logical; if \code{TRUE}, fits a reduced model with no gear size structure and performs a likelihood ratio test.
#' @param logsd Optional numeric value to fix the standard deviation (log-scale) of the random walk increments.
#' @param phi Optional numeric value to fix the AR(1) correlation parameter.
#' @param logsdnug Optional numeric value to fix the standard deviation (log-scale) of the nugget effect.
#' @param logsdres Optional numeric value to fix the standard deviation (log-scale) of the residuals.
#' @param logsdGearRW Optional numeric value to fix the standard deviation (log-scale) of the gear random walk.
#' @param logalpha Optional numeric value to fix the exponent of the swept area effect (default is 0).
#' @param logtheta Optional numeric value to fix the overdispersion parameter for negative binomial and ZIP models.
#' @param logitpi Optional numeric value to fix the logit of the zero-inflation probability for the ZIP model.
#' @param model Character; model to use. One of \code{"poisson"}, \code{"nb"}, or \code{"zip"}. Note that zip and nb are not stable.
#'
#' @return An object of class \code{"gearcalibFit"} containing:
#'   \item{d}{Input data list.}
#'   \item{Pvalue}{P-value from likelihood ratio test (if \code{fit0 = TRUE}), otherwise \code{NA}.}
#'   \item{rep}{TMB \code{sdreport} object containing estimated parameters and standard deviations.}
#'   \item{opt}{Optimization results from \code{nlminb}.}
#'   \item{obj}{TMB objective function.}
#'   \item{est}{Estimated gear effect sizes.}
#'   \item{sd}{Standard deviations of the gear effect estimates.}
#'
#' @details
#' Model and code adapted from Thygesen et al. (2019; https://github.com/Uffe-H-Thygesen/Intercalibration).
#' @references Thygesen, U.H., Kristensen, K., Jansen, T., Beyer, J.E., 2019. Intercalibration of survey methods using paired fishing operations and log-Gaussian Cox processes. ICES J. Mar. Sci. 76, 1189–1199. https://doi.org/10.1093/icesjms/fsy191
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb pchisq
#' @export
gearcalib_fit <- 
  function(d,
           fit0=FALSE,
           logsd=NA,
           phi=NA,
           logsdnug=NA,
           logsdres=NA,
           logsdGearRW=NA,
           logalpha=0,
           logtheta=NA,
           logitpi=NA,
           model = "poisson") {
    
    nsize <- ncol(d$N)
    ngroup <- nlevels(d$group)
    nhaul <- nrow(d$N)
    ngear <- nlevels(d$Gear)
    
    if(!("Lvec" %in% names(d))) {
      d$Lvec <- 1:ncol(d$N)
    }
    
    ## Set default RW order
    rw_order <- c(1,1)
    
    data <- list(
      N=d$N,
      SweptArea=d$SweptArea,
      group=d$group,
      Gear=d$Gear,
      huge=10,
      tiny=0.01,
      rw_order=rw_order
    )
    
    data$model_type <- 
      switch(
        model,
        "poisson" = 1,
        "nb" = 2,
        "zip" = 3
      )
    
    
    parameters <- 
          list(
            logspectrum=matrix(0,ngroup,nsize),
            nugget=matrix(0,nhaul,nsize),
            residual=matrix(0,nhaul,nsize),
            loggear=numeric(nsize),
            logsd=-1,
            phi=0.9,
            logsdnug=-1,
            logsdres=-1,
            logsdGearRW=-1,
            logalpha = 0,
            logitpi = 0,
            logtheta = -1
          )
    
    random <- c("logspectrum","residual","loggear","nugget")
  
    map <- list()
    
    ## If any parameters have been specified in the call, do not estimate those
    ## parameters but fix them to the specified value
    setparameter <- function(name) {
      var <- get(name)
      if(!is.na(var))
      {
        map[[name]] <<- factor(NA)
        parameters[[name]] <<- var
      }
    }
    parameternames <- c("logsd","phi","logsdnug","logsdres","logsdGearRW","logalpha","logitpi","logtheta")
    sapply(parameternames, setparameter)
    
    # Map parameters
    map <- 
      switch(
          model,
          "poisson" = 
            c(map, 
              list(logtheta = factor(NA)),
              list(logitpi = factor(NA))
              ),
          "nb" = 
            c(map, 
              list(logitpi = factor(NA))
            ),
          "zip" = 
            c(map, 
              list(logtheta = factor(NA))
            )
        )
    
    file.copy(from = system.file("executables", "gearcalib_switch.cpp", package = "sratio"),
              to = paste0(getwd(), "/gearcalib_switch.cpp"))
    compile(paste0(getwd(), "/gearcalib_switch.cpp"))
    dyn.load(paste0(getwd(), "/", dynlib("gearcalib_switch")))
    
    obj <- 
      TMB::MakeADFun(
      data = data,
      parameters = parameters,
      map = map, 
      random = random,
      DLL="gearcalib_switch"
    )
    
    
    lower <- 0*obj$par-Inf
    upper <- 0*obj$par+Inf
    if(any("phi" == names(obj$par)))
    {
      lower["phi"] <- 0
      upper["phi"] <- 0.99
    }
    
    
    obj$env$tracepar <- TRUE
    
    system.time( opt <- nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper) )
    Pvalue <- NA
    
    if(fit0)
    {
      parameters0 <- parameters
      parameters0$logsdGearRW <- -10
      
      obj0 <- 
        TMB::MakeADFun(
        data = data,
        parameters = parameters0,
        map = c(map,list(logsdGearRW=factor(NA))),
        random = random,
        DLL="gearcalib_switch"
      )
      
      obj0$env$tracepar <- TRUE
      
      system.time( opt0 <- nlminb(obj0$par,obj0$fn,obj0$gr,lower=lower,upper=upper) )
      
      Pvalue <- pchisq(2*(opt0$objective-opt$objective),df=1)
      
      print(paste("Chisq-test of no size structure in gear effect:",Pvalue))
      
      rm(obj0)
    }
    
    rep <- sdreport(obj)
    repsum <- summary(rep,"random")
    s <- repsum[grep("gear",rownames(repsum)),]
    est <- 2*s[,1]
    sd <- 2*s[,2]
    
    ret <- list(d=d,
                Pvalue=Pvalue,rep=rep,opt=opt,obj=obj,est=est,sd=sd)
    
    class(ret)<-"gearcalibFit" 
    
    return(ret)
    
  }


#' Bootstrap Estimates of Gear Efficiency Ratio
#'
#' Performs bootstrap resampling to estimate the uncertainty of relative gear efficiency across size bins.
#'
#' @param d A list containing input data with elements:
#'   \itemize{
#'     \item \code{N}: Matrix of counts (hauls x size bins).
#'     \item \code{SweptArea}: Numeric vector of swept areas per haul.
#'     \item \code{group}: Factor indicating group membership.
#'     \item \code{Gear}: Factor indicating gear type (must have exactly two levels).
#'     \item \code{Lvec}: Numeric vector of size bin labels.
#'   }
#' @param quantiles Numeric vector of quantiles to compute for bootstrap estimates (default is \code{c(0.025, 0.16, 0.5, 0.84, 0.975)}).
#' @param nboot Integer; number of bootstrap replicates to perform (default is 1000).
#'
#' @return A list with components:
#'   \item{RawEstimate}{Vector of raw relative gear efficiency estimates across size bins.}
#'   \item{BootEstimate}{Matrix of bootstrap estimates (rows = replicates, columns = size bins).}
#'   \item{BootQuantiles}{Matrix of bootstrap quantiles (rows = quantiles, columns = size bins).}
#'   \item{Density1}{Vector of catch density estimates for gear 1.}
#'   \item{Density2}{Vector of catch density estimates for gear 2.}
#'   \item{GearNames}{Character vector of the two gear level names.}
#'
#' @details
#' The function estimates relative gear efficiency as the ratio of catch densities per swept area between two gear types.
#' Bootstrap resampling is done at the group level (e.g., strata), maintaining group structure while resampling.
#' A small constant is added to avoid division by zero.
#'
#' @importFrom stats quantile
#' @export
gearcalib_boot <- function(d,quantiles = c(0.025,0.16,0.5,0.84,0.975),nboot=1000)
{
  GearNames <- levels(d$Gear)
  
  ## Raw data: Total catches
  I1 <- d$Gear==GearNames[1]
  I2 <- d$Gear==GearNames[2]
  
  totcatch1 <- apply(d$N[I1,],2,sum)
  totcatch2 <- apply(d$N[I2,],2,sum)
  
  Density1 <- totcatch1 / sum(d$SweptArea[I1])
  Density2 <- totcatch2 / sum(d$SweptArea[I2])
  
  tiny <- 1e-6
  
  RawEstimate <- (Density2 + tiny) / (Density1 + tiny)
  
  BootEstimate <- array(NA,c(nboot,length(RawEstimate)))
  
  for(i in 1:nboot)
  {
    ## Select random groups with replacement
    GroupBoot <- sample(1:length(levels(d$group)),length(levels(d$group)),replace=TRUE)
    
    Haul1 <- as.numeric(sapply(GroupBoot,function(g) which( (as.numeric(d$group) == g) & I1)))
    Haul2 <- as.numeric(sapply(GroupBoot,function(g) which( (as.numeric(d$group) == g) & I2)))
    
    Nb1 <- d$N[Haul1,]
    Nb2 <- d$N[Haul2,]
    
    Db1 <- apply(Nb1,2,sum) / sum(d$SweptArea[Haul1])
    Db2 <- apply(Nb2,2,sum) / sum(d$SweptArea[Haul2])
    
    BootEstimate[i,] <- (Db2 + tiny) / (Db1 + tiny) 
  }
  
  BootQuantiles <- apply(BootEstimate,2,function(x)quantile(x,quantiles))
  
  colnames(BootQuantiles) <- d$Lvec
  
  return(list(RawEstimate=RawEstimate,
              BootEstimate=BootEstimate,
              BootQuantiles=BootQuantiles,
              Density1=Density1,
              Density2=Density2,
              GearNames=GearNames))
  
}

#' Plot log-Gaussian Cox process model and bootstrap results
#'
#' Produces plots of relative selectivity ratios, sample sizes, and bootstrap density estimates by gear.
#'
#' @param fit An object of class \code{"gearcalibFit"} as returned by \code{\link{gearcalib_fit}}.
#' @param select Character vector specifying which plots to return. Options are \code{"relsel"} (relative selectivity) and \code{"density"} (CPUE by gear). Default is \code{c("relsel", "density")}.
#' @param boot Optional list of bootstrap results from \code{\link{gearcalib_boot}}. Required for plotting observed ratios and densities.
#' @param xlab Character; label for the x-axis (default is \code{"Length group"}).
#' @param add_bootquantiles Logical; if \code{TRUE}, overlay bootstrap quantiles on the relative selectivity plot. Default is \code{FALSE}.
#' @param log_transform_cpue Logical; if \code{TRUE}, log-transform CPUE values (base 10) in the density plot. Default is \code{FALSE}.
#'
#' @return A named list of \code{ggplot} objects. Contains:
#'   \item{p_fit}{Relative selectivity plot (if \code{"relsel"} is selected).}
#'   \item{p_cpue}{CPUE by gear plot (if \code{"density"} is selected and \code{boot} is provided).}
#'
#' @details
#' The relative selectivity plot shows model-based estimates (with ±2 SD ribbons) and optionally raw ratios and bootstrap quantiles.
#' The CPUE plot shows catch densities per swept area for both gears, with optional log-scale.
#'
#' @import ggplot2
#' @importFrom scales squish_infinite
#' @export
gearcalib_plot <- function(fit, select = c("relsel", "density"), boot = NULL, xlab = "Length group", add_bootquantiles = FALSE, log_transform_cpue = FALSE) {
  
  Lmin <- min(fit$d$Lvec)
  Lmax <- max(fit$d$Lvec)
  
  output <- list()
  
  if("relsel" %in% select) {
    
    p_fit <- 
      ggplot() +
      geom_ribbon(
        mapping = aes(
          x = fit$d$Lvec, 
          ymin = exp(fit$est + fit$sd * -2), 
          ymax = exp(fit$est + fit$sd *2)
        ),
        alpha = 0.7,
        fill = "grey70"
      ) +
      geom_hline(yintercept = 1, linetype = 2, linewidth = rel(1.1)) +
      geom_path(
        mapping = aes(x = fit$d$Lvec, y = exp(fit$est)), 
        linewidth = rel(1.1)
      )
    
    if(!is.null(boot)) {
      p_fit <- 
        p_fit +
        geom_point(
          mapping = aes(x = fit$d$Lvec, y = boot$RawEstimate)
        )
      
      if(add_bootquantiles) {
        
        boot_quantiles <-
          boot$BootQuantiles |>
          as.table() |>
          as.data.frame()
        
        names(boot_quantiles) <- c("quantile", "l", "rel_s")
        boot_quantiles$l <- as.numeric(as.character(boot_quantiles$l))
        
        p_fit <- 
          p_fit +
          geom_path(data = boot_quantiles,
                    mapping = aes(x = l,
                                  y = rel_s,
                                  group = quantile),
                    linetype = 3
          )
        
      }
    }
    
    
    p_fit <- 
      p_fit +
      scale_x_continuous(
        name = xlab, 
        limits = c(Lmin, Lmax)
      ) +
      scale_y_continuous(
        name = "Relative selectivity",
        # name = paste(levels(fit$d$Gear)[2]," vs. ", levels(fit$d$Gear)[1]),
        limits = c(0, 2),
        oob = scales::squish_infinite
      ) +
      theme_bw()
    
    output <- c(output, list(p_fit = p_fit))
    
  }
  
  
  if(!is.null(boot) && "density" %in% select) {
    
    if(log_transform_cpue) {
      p_cpue <- 
        ggplot() +
        geom_point(
          mapping =
            aes(
              x = fit$d$Lvec,
              y = log10(boot$Density1+1),
              shape = boot$GearNames[1])
        ) +
        geom_path(
          mapping =
            aes(
              x = fit$d$Lvec,
              y = log10(boot$Density1+1),
              linetype = boot$GearNames[1])
        ) +
        geom_point(
          mapping =
            aes(
              x = fit$d$Lvec,
              y = log10(boot$Density2+1),
              shape = boot$GearNames[2])
        ) +
        geom_path(
          mapping =
            aes(
              x = fit$d$Lvec,
              y = log10(boot$Density2+1),
              linetype = boot$GearNames[2])
        ) +
        scale_shape(solid = FALSE) +
        scale_x_continuous(name = xlab, limits = c(Lmin, Lmax)) +
        scale_y_continuous(name = expression('CPUE ('*log[10]*'(N/A+1))')) +
        theme_bw() +
        theme(legend.title = element_blank())
    } else {
      p_cpue <- 
        ggplot() +
        geom_point(
          mapping = 
            aes(
              x = fit$d$Lvec,
              y = boot$Density1,
              shape = boot$GearNames[1])
        ) +
        geom_path(
          mapping = 
            aes(
              x = fit$d$Lvec,
              y = boot$Density1,
              linetype = boot$GearNames[1])
        ) +
        geom_point(
          mapping = 
            aes(
              x = fit$d$Lvec,
              y = boot$Density2,
              shape = boot$GearNames[2])
        ) +
        geom_path(
          mapping = 
            aes(
              x = fit$d$Lvec,
              y = boot$Density2,
              linetype = boot$GearNames[2])
        ) +
        scale_shape(solid = FALSE) +
        scale_x_continuous(name = xlab, limits = c(Lmin, Lmax)) +
        scale_y_log10(name = expression('CPUE (#'%.%km^-2*')')) +
        theme_bw() +
        theme(legend.title = element_blank())
    }
    
    
    output <- c(output, list(p_cpue = p_cpue))
    
  }
  
  return(output)
  
}

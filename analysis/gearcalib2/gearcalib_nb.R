##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Fit relative catch effiency by length group between to gears
##' @param d A list with four elements: N (matrix of integers), SweptArea(numeric vector), group(factor vector), and Gear(factor vector). 
##' @param fit0 If TRUE a Chisq-test of no size structure in gear effect is performed.
##' @param linearEffort If TRUE, catch is assumed proportional to swept area, otherwise proportional to a power function of swept area.
##' @param model "poisson" = Poisson, "nb" = Negative binomial, "zip" = zero-inflated Poisson
##' @return A list
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
    
    # check.gearcalib.data(d)
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
    
    obj <- MakeADFun(
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
      
      obj0 <- MakeADFun(
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Empirical estimate with bootstrap confidence regions
##' @param d gearcalib data
##' @param quantiles report these quantiles
##' @param nboot number of bootstrap replicates
##' @return a list
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


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Empirical estimate with bootstrap confidence regions
##' @param d gearcalib data
##' @param quantiles report these quantiles
##' @param Nboot number of bootstrap replicates
##' @return a list
boot <- function(d,quantiles = c(0.025,0.16,0.5,0.84,0.975),Nboot=1000)
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
  
  BootEstimate <- array(NA,c(Nboot,length(RawEstimate)))
  
  for(i in 1:Nboot)
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plots a gear calibration fit
##' @param fit an onject of class 'gearcalibFit' as fitted by gearcalib()
##' @param select a vector specifying which plots are wanted. Can be "relsel" for "Relative selectivity" or "density".
##' @param boot (optional) list with bootstrap estimates as produced by boot()
##' @param Lvec (optional) vector with labels for each length group
##' @return nothing
plot.gearcalibFit <- function(fit,select=c("relsel","density"),boot=NULL, Lvec=NULL)
{
  
  if(is.null(Lvec)) Lvec <- 1:(ncol(fit$d$N))
  Lmin <- min(Lvec)
  Lmax <- max(Lvec)
  
  with(fit,{
    
    if("relsel" %in% select){
      plot(range(Lvec),c(0,max(exp(est + sd %o% c(0,-2,2)))),type="n",
           ylim=c(0,2),
           xlim=range(Lvec),
           xlab="Length group",
           ylab=paste(levels(fit$d$Gear)[2]," vs. ",levels(fit$d$Gear)[1]),
           main="Relative selectivity")
      
      lines(range(Lvec),rep(1,2),col="darkgrey",lwd=3,lty="dashed")
      
      polygon(c(Lvec,rev(Lvec)),c(exp(est-2*sd),rev(exp(est+2*sd))),
              col="grey",border=NA)
      lines(Lvec,exp(est),lwd=3)
      if(!is.null(boot)){
        points(Lvec,boot$RawEstimate)
        apply(boot$BootQuantiles,1,function(x)lines(Lvec,x))
      }
      grid()
    }
    
    
    if(!is.null(boot) && "density" %in% select){
      with(boot,{
        plot(Lvec,log10(1+Density1),
             ylim=log10(1+range(c(Density1,Density2))),
             xlim=c(Lmin,Lmax),
             type="l",lty="dashed",
             xlab="Length [cm]",ylab="Density (log10(N/A+1))")
        points(Lvec,log10(1+Density1),pch="o")
        lines(Lvec,log10(1+Density2))
        points(Lvec,log10(1+Density2),pch="+")
        legend("topright",legend=GearNames,lty=c("dashed","solid"),pch=c("o","+"))
        
        grid()
      })
    }
  })
}

plot_gearcalibFit <- function(fit, select = c("relsel", "density"), boot = NULL, xlab = "Length group", add_bootquantiles = FALSE) {
  
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
        linewidth = rel(1.5)
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
        name = paste(levels(fit$d$Gear)[2]," vs. ", levels(fit$d$Gear)[1])
      ) +
      theme_bw()
    
    output <- c(output, list(p_fit = p_fit))
    
  }

  
  if(!is.null(boot) && "density" %in% select) {
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
    scale_y_continuous(name = "CPUE (log10(N/A+1))") +
    theme_bw() +
    theme(legend.title = element_blank())
    
    output <- c(output, list(p_cpue = p_cpue))
    
  }
  
  return(output)
  
}

test <- plot_gearcalibFit(
  fit = out3,
  select = c("relsel", "density"),
  boot=out_boot,
  xlab = "Length",
  add_bootquantiles = TRUE
)

test$p_fit

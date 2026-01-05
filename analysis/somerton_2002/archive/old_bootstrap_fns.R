boot_pois_nb <- 
  function(model_list, boot_samples_list) {
    
    model_fits <- vector(mode = "list", length = length(model_list))
    names(model_fits) <- names(model_list)
    
    bootstrap_results <- data.frame()
    
    for(ii in seq_along(model_list)) {
      
      cat(names(model_list)[ii], "\n")
      
      mod <- model_list[[ii]]
      
      fits_list <- lapply(boot_samples_list, function(samp) {
        boot_fit <- update(mod, data = samp)
        
        coef <- 
          fixef(boot_fit) |> 
          as.list() |>
          do.call(what = rbind) |>
          as.data.frame()
        
        type <- rownames(coef)
        rownames(coef) <- NULL
        coef$type <- type
        coef$model_name <- names(model_list)[ii]
        
        coef
        
      })
      
      fits <- do.call(rbind, fits_list)
      
      value_cols <- names(fits)[!(names(fits) %in% c("type", "model_name"))]
      
      ci <- 
        fits %>%
        dplyr::group_by(
          dplyr::across(
            dplyr::all_of(c("type", "model_name")))
        ) %>%
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(value_cols),
            list(
              median = ~median(.x, na.rm = TRUE),
              lci = ~quantile(.x, probs = 0.025),
              uci = ~quantile(.x, probs = 0.975)
            ),
            .names = "{.col}_{.fn}"
          ),
          .groups = "drop"
        )
      
      model_fits[[ii]] <- fits
      
    }
    
    output <- 
      list(
        ci = bootstrap_results,
        fits = model_fits
      )
    
  }


boot_betareg <- 
  function(model_list, boot_samples_list) {
    
    model_fits <- vector(mode = "list", length = length(model_list))
    names(model_fits) <- names(model_list)
    
    bootstrap_results <- data.frame()
    
    for(ii in seq_along(model_list)) {
      
      cat(names(model_list)[ii], "\n")
      
      mod <- model_list[[ii]]
      
      fits <- data.frame()
      
      for(jj in 1:length(boot_samples_list)) {
        
        if(jj%%100 == 0) {
          cat(jj, "\n")
        }
        
        boot_fit <- update(mod, data = boot_samples_list[[jj]])
        
        coef <- 
          fixef(boot_fit) |> 
          as.list() |>
          do.call(what = rbind) |>
          as.data.frame()
        
        type <- rownames(coef)
        rownames(coef) <- NULL
        coef$type <- type
        coef$model_name <- names(model_list)[ii]
        
        fits <- dplyr::bind_rows(
          fits,
          coef
        )
        
      }
      
      value_cols <- names(fits)[!(names(fits) %in% c("type", "model_name"))]
      
      ci <- 
        fits %>%
        dplyr::group_by(
          dplyr::across(
            dplyr::all_of(c("type", "model_name")))
        ) %>%
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(value_cols),
            list(
              median = ~median(.x, na.rm = TRUE),
              lci = ~quantile(.x, probs = 0.025),
              uci = ~quantile(.x, probs = 0.975)
            ),
            .names = "{.col}_{.fn}"
          ),
          .groups = "drop"
        )
      
      model_fits[[ii]] <- fits
      
    }
    
    output <- 
      list(
        ci = bootstrap_results,
        fits = model_fits
      )
    
  }


boot_ols <- 
  function(model_list, boot_samples_list) {
    
    model_fits <- vector(mode = "list", length = length(model_list))
    names(model_fits) <- names(model_list)
    
    bootstrap_results <- data.frame()
    
    for(ii in seq_along(model_list)) {
      
      mod <- model_list[[ii]]
      
      fits <- data.frame()
      
      for(jj in 1:length(boot_samples_list)) {
        
        boot_fit <- update(mod, data = boot_samples_list[[jj]])
        
        boot_fpc <- miller_bias_correct(mod = boot_fit)
        
        boot_fpc$model_name <- names(model_list)[ii]
        
        fits <- dplyr::bind_rows(
          fits,
          boot_fpc
        )
        
      }
      
      model_fits[[ii]] <- fits
      
      bootstrap_median <-
        data.frame(
          model_name = names(model_list)[ii],
          method = "median",
          
          ratio     = quantile(fits$ratio, 0.5),
          ratio_lci        = quantile(fits$ratio, 0.025),
          ratio_uci        = quantile(fits$ratio, 0.975)
        )
      
      bootstrap_mean <-
        data.frame(
          model_name = names(model_list)[ii],
          method = "mean",
          
          ratio  = quantile(fits$ratio_bc, 0.5),
          ratio_lci     = quantile(fits$ratio_bc, 0.025),
          ratio_uci     = quantile(fits$ratio_bc, 0.975)
        )
      
      bootstrap_results <-
        dplyr::bind_rows(
          bootstrap_median,
          bootstrap_mean,
          bootstrap_results
        )
      
    }
    
    output <- 
      list(
        ci = bootstrap_results,
        fits = model_fits
      )
    
  }
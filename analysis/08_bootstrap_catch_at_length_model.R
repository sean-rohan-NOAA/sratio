library(sratio)

n_cores <- 4

treatments <- factor(c(15,30))

bootstrap_sample_path <- list.files(here::here("output"), 
                                    recursive = TRUE, 
                                    pattern = "bootstrap_samples_", 
                                    full.names = TRUE)

cal_rmse_df <- read.csv(file = here::here("output", "cal_model_rmse.csv"))

for(ii in 1:length(bootstrap_sample_path)) {
  
  boot_dat <- readRDS(bootstrap_sample_path[ii])
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_sample_path[ii])))
  
  best_family <- cal_rmse_df$name[cal_rmse_df$LOWEST_RMSE & cal_rmse_df$SPECIES_CODE == sp_code]
  
  gam_family <- switch(best_family,
                        tw = tw(link = "log"),
                        nb = nb(link = "log"),
                        poisson = poisson(link = "log"))
  
  gam_knots <- cal_rmse_df$gam_knots[cal_rmse_df$LOWEST_RMSE & cal_rmse_df$SPECIES_CODE == sp_code]
  
  
  lengths <- seq(min(unlist(lapply(boot_dat, FUN = function(x) {min(x$SIZE_BIN)}))), 
                 max(unlist(lapply(boot_dat, FUN = function(x) {max(x$SIZE_BIN)}))),
                 by = 1)
  
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  bootstrap_output <- foreach::foreach(iter = 1:length(boot_dat), .packages = "mgcv") %dopar% {
    
    boot_df <- boot_dat[[iter]]
    boot_df$dummy_var <- 1
    
    model <- mgcv::gam(formula = FREQ_EXPANDED ~ s(SIZE_BIN, bs = "tp", k = gam_knots, by = TREATMENT) + s(MATCHUP, bs = "re", by = dummy_var) + offset(I(log(AREA_SWEPT_KM2))),
                       family = gam_family,
                       data = boot_df)
    
    fit_df <- expand.grid(SIZE_BIN = lengths,
                         MATCHUP = boot_dat[[iter]]$MATCHUP[1],
                         TREATMENT = treatments,
                         AREA_SWEPT_KM2 = 1,
                         dummy_var = 0,
                         draw = iter) # random effects off
    
    fit_df$fit <- predict(model, newdata = fit_df, type = "response")
    
    return(fit_df)
    
    
  }
  
  doParallel::stopImplicitCluster()
  
  bootstrap_df <- do.call("rbind", bootstrap_output)
  
  saveRDS(bootstrap_df, file = here::here("output", sp_code,  paste0("cal_model_bootstrap_results_", sp_code, ".rds")))
  
}



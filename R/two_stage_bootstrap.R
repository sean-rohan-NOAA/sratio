#' Two-stage bootstrap from paired samples
#' 
#' Two-stage bootstrap resampling (by block/pair then individual within each pair), where samples are drawn from paired treatments.
#' 
#' @param count1 Numeric vector of catch-at-size/age for gear #1
#' @param count2 Numeric vector of catch-at-size/age for gear #2
#' @param size1  Numeric vector of sizes or ages for gear #1
#' @param size2  Numeric vector of sizes or ages for gear #2
#' @param block1 Sample block for gear #1
#' @param block2 Sample block for gear #2
#' @param n_draws Number of samples to draw.
#' @param seed random number generator seed
#' @param treatment_name1 Name of the treatment for gear #1.
#' @param treatment_name2 Name of the treatment for gear #2.
#' @export

two_stage_bootstrap <- function(count1,
                                count2,
                                size1,
                                size2,
                                block1,
                                block2,
                                add_columns1,
                                add_columns2,
                                n_draws,
                                seed = NULL,
                                treatment_name1 = 1,
                                treatment_name2 = 2) {
  
  boot_samples <- vector(mode = "list", length = n_draws)
  
  unique_blocks <- unique(c(block1, block2))
  n_block <- length(unique_blocks)
  
  set.seed(seed)
  sel_block <- sample(x = unique_blocks,
                       size = n_draws*n_block,
                       replace = TRUE)
  
  trt1 <- data.frame(size = size1,
                     block = block1,
                     count = count1,
                     treatment = treatment_name1) |>
    dplyr::group_by(size, block, treatment) |>
    dplyr::mutate(count = sum(count)) |>
    as.data.frame()
  
  trt2 <- data.frame(size = size2,
                     block = block2,
                     count = count2,
                     treatment = treatment_name2) |>
    dplyr::group_by(size, block, treatment) |>
    dplyr::mutate(count = sum(count)) |>
    as.data.frame()
  
  set.seed(seed)
  for(ii in 1:length(boot_samples)) {
    
    group_sample <- data.frame(
      size = vector(), 
      treatment = vector(),
      count = vector(),
      block = vector(),
      new_block = vector()
    )
    
    for(jj in 1:n_block) {
      
      set_block <- sel_block[jj + (ii-1)*n_block]
      trt1_index <- trt1$block == set_block
      trt2_index <- trt2$block == set_block
      
      block_sample <- rbind(
        data.frame(size = 
                     sample(
                       x = trt1$size[trt1_index], 
                       size = sum(trt1$count[trt1_index]),
                       replace = TRUE,
                       prob = trt1$count[trt1_index]
                     ),
                   treatment = treatment_name1),
        data.frame(size = 
                     sample(x = trt2$size[trt2_index], 
                            size = sum(trt2$count[trt2_index]),
                            replace = TRUE,
                            prob = trt2$count[trt2_index]
                     ),
                   treatment = treatment_name2)
      ) |>
        dplyr::group_by(treatment, size) |>
        dplyr::summarise(count = n(), .groups = "keep") |>
        as.data.frame()
      
      block_sample$original_block <- set_block
      block_sample$new_block <- paste0(set_block, "_", jj)
      
      group_sample <- rbind(group_sample, block_sample)
      
    }
    
    boot_samples[[ii]] <- group_sample
    
  }
  
  return(boot_samples)
  
}

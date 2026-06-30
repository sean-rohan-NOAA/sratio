catch_data <- readRDS(
  here::here("analysis", "shelf_slope_gear_selectivity", "data", "catch_data.rds")
)

dplyr::group_by(catch_data,
                COMMON_NAME, GEAR) |>
  dplyr::summarise(N_HAULS = n()) |>
  dplyr::ungroup() |>
  tidyr::pivot_wider(
    names_from = "GEAR", values_from = "N_HAULS"
  )



aa <- rpois(30, 100)
bb <- rpois(30, 100)
bb_frac <- round(runif(n = 30, 0.1, 1), 1)
bb_samp <- round(bb*bb_frac)

test1 <- 
  data.frame(
    aa = aa,
    bb = bb,
    bb_frac = bb_frac,
    bb_samp = bb_samp,
    total = aa + bb_samp,
    p = bb_samp/(bb_samp+aa)
  )

m1 <- glmmTMB::glmmTMB(p~1 + offset(log(bb_frac)),
                 weights = total,
                 family = binomial(), 
                 data = test1)

predict(m1, type = "response", newdata = data.frame(aa = 1, bb = 1, bb_frac = 1, total = 2))


test2 <- 
  data.frame(
    aa = rpois(100, 100),
    bb = rpois(100, 50),
    aa_effort = 1,
    bb_effort = 0.5,
    bb_frac = round(runif(n = 100, 0.1, 1), 1)
  ) |>
  dplyr::mutate(
    bb_samp = round(bb*bb_frac),
    total = aa + bb_samp,
    p = bb_samp/total
  )

m2 <- glmmTMB::glmmTMB(p~1 + offset(log(bb_effort*bb_frac)),
                       weights = total,
                       family = binomial(), 
                       data = test2)

predict(m2, type = "response", newdata = data.frame(aa = 1, bb = 1, bb_frac = 1, total = 2, bb_effort = 1))
predict(m2, type = "response", newdata = data.frame(aa = 1, bb = 1, bb_frac = 1, total = 2, bb_effort = 0.5))
predict(m2, type = "response", newdata = data.frame(aa = 1, bb = 1, bb_frac = 1, total = 2, bb_effort = 2))


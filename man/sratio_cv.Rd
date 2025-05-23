% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sratio_cv.R
\name{sratio_cv}
\alias{sratio_cv}
\title{Selectivity ratio cross-validation}
\usage{
sratio_cv(
  model_type = "binomial",
  size,
  count1,
  count2,
  effort1,
  effort2,
  sampling_factor1 = 1,
  sampling_factor2 = 1,
  obs_weight_control = list(method = "none", max_count = Inf, residual_type = NA,
    normalize_weight = FALSE),
  block,
  k = NULL,
  scale_method = "sv",
  sratio_type = "absolute",
  n_cores = 1
)
}
\arguments{
\item{model_type}{Model type. Options are "binomial" or "beta"}

\item{size}{Numeric vector of sizes or ages}

\item{count1}{Numeric vector of catch-at-size/age for gear #1}

\item{count2}{Numeric vector of catch-at-size/age for gear #2}

\item{effort1}{Numeric vector of effort for for gear #1}

\item{effort2}{Numeric vector of effort for for gear #2}

\item{sampling_factor1}{Numeric vector of sampling factors for count from gear #1 (i.e., the inverse of the proportion of the catch that was sampled).}

\item{sampling_factor2}{Numeric vector of sampling factors for count from gear #2 (i.e., the inverse of the proportion of the catch that was sampled).}

\item{obs_weight_control}{A list indicating which method to use to weight observations. See help documentation `sratio_fit_gamm()` for information about options (?sratio_fit_gamm) .}

\item{block}{Sample block (i.e. paired sample)}

\item{k}{k to use for GAMs. Automatically set to the minimum of 8 or 3 less than the number of unique values in size.}

\item{scale_method}{Method to use for scaling the catch comparison rate for beta regression. See ?scale_for_betareg}

\item{sratio_type}{Which selectivity ratio calculation should be used? Absolute ("absolute") or relative ("relative")?}

\item{n_cores}{Number of cores to use for parallel processing.}
}
\description{
Block-level cross validation for selectivity ratio binomial and beta binomial generalized additive models.
}

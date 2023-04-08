#' Make GAM prediction data.frame
#' 
#' @param model GAM model
#' @param lengths lengths to predict
#' @param type type of prediction to return
#' @param model_name name/type of model used to generate the prediction
#' @export

make_prediction_df <- function(model, lengths = 10:120, type = "link", model_name) {
  
  out <- expand.grid(LEN_MIDPOINT = lengths, MATCHUP = 1, dummy_var = 0)
  out$fit <- predict(model, newdata = out, type = type)
  out$se.fit <- predict(model, newdata = out, type = type, se.fit = TRUE)$se.fit
  out$model <- model_name
  
  if(type == "link") {
    
    out <- out |> 
      dplyr::mutate(p = inv_logit(fit),
                    p_minus1_se = inv_logit(fit-se.fit),
                    p_plus1_se = inv_logit(fit+se.fit),
                    sratio = 1/inv_logit(fit)-1,
                    sratio_minus1_se = 1/inv_logit(fit-se.fit)-1,
                    sratio_plus1_se = 1/inv_logit(fit+se.fit)-1)
  }
  
  return(out)
}

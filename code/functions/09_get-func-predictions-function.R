# Helper function to generate model-averaged predictions
get_ma_predictions_func <- function(model_avg, raw_data, age_mean, age_sd, order_q_label) {
  
  # Create prediction grid across Age range, for each Treatment
  # Use modal/reference values for other covariates
  pred_grid <- expand.grid(
    Age_raw = seq(min(raw_data$Age.raw, na.rm=T), max(raw_data$Age.raw, na.rm=T), length.out = 50),
    Treatment = levels(raw_data$Treatment),
    Season = levels(raw_data$Season),        # reference level
    Site = levels(raw_data$Site)[1],            # reference level
    crop.stage = levels(raw_data$crop.stage)[1],# reference level
    PC1.climate = 0                             # standardised mean = 0
  )
  
  # Re-standardise Age the same way arm::rescale does (divides by 2*SD)
  pred_grid$Age <- (pred_grid$Age_raw - age_mean) / (2 * age_sd)
  
  # Manually construct interaction terms to ensure MuMIn captures them
  pred_grid$`Age:TreatmentControl` <- pred_grid$Age * (pred_grid$Treatment == "Control")
  pred_grid$`Age:SeasonLate` <- pred_grid$Age * (pred_grid$Season == "Late")
  
  # Predict on log scale then back-transform
  pred_grid$log_qFD <- predict(model_avg, newdata = pred_grid, full = FALSE)
  pred_grid$qFD_pred <- exp(pred_grid$log_qFD)
  pred_grid$Order.q <- order_q_label
  
  return(pred_grid)
}

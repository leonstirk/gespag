#' Return models and bounds for input into runGeneAlgorithm
#' @param gdf data.frame
#' @param mf data.frame
#' @param time_var character
#' @return list
#' @export
makeOrganismModels <- function(gdf, response, control_vars = c()) {

  model <- as.formula(paste0(c(response,
                               paste0(c(
                                 'cell_id',
                                 control_vars
                                 ), collapse = '+')), collapse = ' ~ '))

  theoretical_max <- summary(stats::lm(model, data = gdf))$r.squared

  if(length(control_vars) == 0) { baseline <- 0 }
  else {
    model <- as.formula(paste0(c(response,
                                 paste0(c(
                                   control_vars
                                 ), collapse = '+')), collapse = ' ~ '))
    baseline <- summary(stats::lm(model, data = gdf))$r.squared
  }

  print(paste(c('Theoretical max R2: ', theoretical_max), collapse = ''))
  print(paste(c('Baseline R2: ', baseline), collapse = ''))

  model <- as.formula(paste0(c(response,
                               paste0(c(
                                 'spag_class',
                                 control_vars
                               ), collapse = '+')), collapse = ' ~ '))

  return(list('theoretical_max' = theoretical_max, 'baseline' = baseline, 'model' = model))
}

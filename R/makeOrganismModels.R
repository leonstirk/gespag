#' Return models and bounds for input into runGeneAlgorithm
#' @param gdf data.frame
#' @param mf data.frame
#' @param time_var character
#' @return list
#' @export
makeOrganismModels <- function(gdf, mf, time_var = 'sale_quarter') {
  model <- mf %>% modman::modelSpec('resp', c('cell_id', 'property_controls', time_var))
  theoretical_max <- summary(stats::lm(model, data = gdf))$r.squared
  model <- mf %>% modman::modelSpec('resp', c('property_controls', time_var))
  baseline <- summary(stats::lm(model, data = gdf))$r.squared

  print(paste(c('Theoretical max R2: ', theoretical_max), collapse = ''))
  print(paste(c('Baseline R2: ', baseline), collapse = ''))

  model <- mf %>% modman::modelSpec('resp', c('spag_class', 'property_controls', time_var))

  return(list('theoretical_max' = theoretical_max, 'baseline' = baseline, 'model' = model))
}

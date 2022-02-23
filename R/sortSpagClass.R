#' Attach and sort spag class from organism to object to gdf object
#' @param gdf data.frame
#' @param organism list
#' @param mf data.frame
#' @return data.frame
#' @export
sortSpagClass <- function(gdf, organism, organism_models) {

  gdf$spag_class <- plyr::mapvalues(gdf$cell_id, organism[['gene_mapping']]$cell_ids, organism[['gene_mapping']]$genome)

  ## Order the spag_class factor levels
  model <- organism_models$model
  gdf$spag_class <- forcats::fct_relevel(gdf$spag_class, as.character(seq(1:length(levels(gdf$spag_class)))-1))
  a <- lm(model, gdf) %>% summary %>% .$coefficients %>% .[stringr::str_which(rownames(.),'spag_class'),] %>% as.data.frame()
  a <- rbind(a, 'spag_class0' = c(0,0,0,0))
  a <- a[with(a, order(Estimate)),] %>% rownames %>% stringr::str_extract('[0-9]')
  gdf$spag_class <- plyr::mapvalues(gdf$spag_class, a, seq(1:length(a)))
  gdf$spag_class <- forcats::fct_relevel(gdf$spag_class, as.character(seq(1:length(a))))
  rm(a)

  return(gdf)
}

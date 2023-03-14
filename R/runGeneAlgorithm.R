#' Run genetic algorithm
#' @param ttv list
#' @param cell_ids character
#' @param organism_models list
#' @param n_class integer
#' @param n_gen integer
#' @param max_mute_intensity numeric
#' @param mute_decay numeric
#' @param pop integer
#' @param grid sf list
#' @param resolution integer
#' @return list
#' @export
runGeneAlgorithm <- function(ttv, cell_ids, organism_models, n_class = 4, n_gen = 2000, max_mute_intensity = 0.5, mute_decay = 2, pop = 52, grid, resolution) {
  ## problem with double naming the model variable
  model <- organism_models$model
  theoretical_max <- organism_models$theoretical_max
  baseline <- organism_models$baseline

  evolution <- list('train' = c(), 'test' = c(), 'validate' = c())
  gen <- 1
  mdf <- list()

  ## first gen model
  for(i in 1:pop) {
    mdf[[i]] <- firstGenModel(ttv$train, cell_ids, model, n_class)
  }
  mdf <- as.data.frame(data.table::rbindlist(lapply(mdf, as.list)))

  ## n gen models
  while(gen <= n_gen) {
    mdf <- mdf[with(mdf, order(-r2))[1:(pop/2)],]

    ## Record progress
    evolution$train[gen]    <- mdf[1,'r2'] ## r2 of the best performer of the generation
    evolution$test[gen]     <- nGenModel(ttv$test, cell_ids, model, n_class, unlist(strsplit(mdf[1,'genome'], character(0))), gen)[['r2']]
    evolution$validate[gen] <- nGenModel(ttv$validate, cell_ids, model, n_class, unlist(strsplit(mdf[1,'genome'], character(0))), gen)[['r2']]
    ## Print live progress to console
    print(paste('Generation:',
                gen,
                ". Max r2: ",
                mdf[1,'r2'],
                ". % improvement over baseline: ",
                round((mdf[1,'r2'] - baseline)/(theoretical_max - baseline), digits = 4)*100,
                "%",
                sep = ""))
    ## Reset working variables
    gen <- gen + 1
    o_mdf <- list()

    ## Chromosones
    xlen <- ceiling(length(cell_ids)/2)
    x <- lapply(strsplit(mdf[,'genome'], character(0)), function(x) { paste(x[1:xlen], collapse = "") })
    y <- lapply(strsplit(mdf[,'genome'], character(0)), function(x) { paste(x[(xlen+1):length(cell_ids)], collapse = "") })

    offspring_genomes <- strsplit(paste(x,sample(y, length(x), replace = FALSE), sep = ''), character(0)) %>% lapply(as.integer)
    offspring_genomes <- lapply(offspring_genomes, function(x) { mutateGen(x, n_class, max_mute_intensity, n_gen, gen, mute_decay) })

    for(i in 1:length(offspring_genomes)) {
      o_mdf[[i]] <- nGenModel(ttv$train, cell_ids, model, n_class, offspring_genomes[[i]], gen)
    }

    o_mdf <- as.data.frame(data.table::rbindlist(lapply(o_mdf, as.list)))
    mdf <- rbind(mdf, o_mdf)
  }

  ## evolution
  evolution <- do.call(cbind, evolution) %>% as.data.frame()
  evolution$generation <- seq(1:nrow(evolution))

  ## winner
  mdf <- mdf[with(mdf, order(-r2)),]
  winner <- as.numeric(strsplit(mdf[1, 'genome'], character(0))[[1]])
  gene_mapping <- as.data.frame(cbind('genome' = winner, 'cell_ids' = cell_ids))

  return(
    list(
      'grid' = grid,
      'gene_mapping' = gene_mapping,
      'evolution' = evolution,
      'other' = list(
        'resolution' = resolution,
        'theoretical_max' = organism_models$theoretical_max,
        'baseline' = organism_models$baseline,
        'model' = organism_models$model
      )
    )
  )
}

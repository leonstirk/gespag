#' @importFrom magrittr %>%
magrittr::`%>%`

#' @importFrom stats lm
stats::lm

#' Combine x and y location columns into a single identifier column for each grid cell.
#' @param xy matrix
#' @return numeric
cantorPair <- function(xy) {
  a <- xy[,1]
  b <- xy[,2]
  out <- 0.5*(a+b)*(a+b+1)+b
  return(out)
}

#' Convert lon/lat values into grid cell row/col (zero indexed, integer) identifiers using an origin point and cell width/height
#' @param dat data.frame
#' @param origin numeric
#' @param cellsize numeric
#' @return matrix
cellXY <- function(dat, origin, cellsize) {
  apply(dat, 1, function(x) { floor((x - origin)/cellsize) }) %>% t %>% apply(., 2, function(x) {
    x[which(x == max(x))] <- max(x) - 1
    return(x)
    })
}

#' Convert grid cell x/y identifiers and return the original coordinate values of the centre of each grid cell
#' @param xy matrix
#' @param origin numeric
#' @param cellsize numeric
#' @return matrix
cellCoords <- function(xy, origin, cellsize) {
  apply(xy, 1, function(x) { cellsize/2+origin+cellsize*x }) %>% t
}

#' Produce a vector of n classes partitioned by a vector of breaks.
#' Also accepts part = 'equal' for n equal sized classes
#' and part = 'ml' that gives 60%, 20%, 20% partitions intended for test, training, and validation sets.
#' @param df data.frame
#' @param part numeric
#' @param n integer
#' @return integer
bootstrap <- function(df, part, n = 3) {
  if(anyNA(part) | part == 'equal') {
    part <- c(0,seq(1:n)/n)
  } else if (part == 'ml') {
    part <- c(0,0.6,0.8,1)
  }
  return(sample(rep(1:n, diff(floor(nrow(df) * part)))))
}

#' Split one data.frame into list of train/test/validate data.frames
#' @param df data.frame
#' @return list
#' @export
## df %>% TTV()
TTV <- function(df) {
  bs <- bootstrap(df, 'ml') %>% as.factor() %>% forcats::fct_relabel(~c('train', 'test', 'validate'))
  ttv <- levels(as.factor(bs)) %>% lapply(., function(x) { df[bs == x,] })
  names(ttv) <- c('train', 'test', 'validate')
  return(ttv)
}

#' Run the first generation of models
#' @param df data.frame
#' @param cell_ids integer
#' @param model formula
#' @param n_class integer
#' @return list
#' @export
firstGenModel <- function(df, cell_ids, model, n_class) {
  genome <- floor(stats::runif(length(cell_ids), 0, n_class))
  df$spag_class <- plyr::mapvalues(df$cell_id, cell_ids, genome)
  r2 <- summary(stats::lm(model, data = df))$r.squared
  genome <- paste(genome, collapse = "")
  return(list("r2" = r2, "genome" = genome, "gen" = 1))
}

#' Run the nth generation of models
#' @param df data.frame
#' @param cell_ids integer
#' @param model formula
#' @param n_class integer
#' @param genome character
#' @param gen integer
#' @return list
#' @export
nGenModel <- function(df, cell_ids, model, n_class, genome, gen) {
  df$spag_class <- plyr::mapvalues(df$cell_id, cell_ids, genome)
  r2 <- summary(stats::lm(model, data = df))$r.squared
  genome <- paste(genome, collapse = "")
  return(list("r2" = r2, "genome" = genome, "gen" = gen))
}

#' Mutate the offspring
#' @param genome integer
#' @param n_class integer
#' @param mute_intensity numeric
#' @param n_gen numeric
#' @param gen numeric
#' @param mute_decay numeric
#' @return integer
#' @export
mutateGen <- function(genome, n_class, mute_intensity, n_gen = NA, gen = NA, mute_decay = 2) {
  l <- length(genome)
  if(!anyNA(c(n_gen, gen))) {
    n <- ceiling((-(gen/n_gen)^(1/mute_decay)+1)*mute_intensity+0.00001)
  } else {
    n <- round(stats::runif(1,1,ceiling(l*mute_intensity)),digits = 0)
  }
  a <- ceiling(stats::runif(n,0,l))
  b <- floor(stats::runif(n,0,n_class))
  genome[a] <- b
  return(genome)
}

cantorPair <- function(a, b) {
    out <- 0.5*(a+b)*(a+b+1)+b
    return(out)
}

ji <- function(xy, origin, cellsize) {
    t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
}

gridData <- function(data, x, y, res) {

    lon <- data[x]
    lat <- data[y]

    origin <- c(min(lon), min(lat))
    cellsize <- c(diff(range(lon))/res, diff(range(lat))/res)

    xy <- cbind(lon, lat)
    JI <- ji(xy, origin, cellsize)

    X <- JI[,1]
    Y <- JI[,2]
    cell <- paste(X, Y)

    dat <- cbind(data, X, Y, cell)

    f <- c('cell')
    n <- c('X', 'Y')
    dat[f] <- lapply(dat[f], as.factor)
    dat[n] <- lapply(dat[n], function(x) { as.numeric(as.character(x)) })

    dat$cantor <- as.factor(with(dat, cantorPair(X,Y)))

    return(dat)
}

bootstrap <- function(d, n, part) {
    if(anyNA(part)) {
        part <- c(0,seq(1:n)/n)
    }
    return(sample(rep(1:n, diff(floor(nrow(d) * part)))))
}


first_gen_model <- function(ck, data, model, n_sub) {
    genome <- floor(runif(length(ck), 0, n_sub))
    data['genome'] <- plyr::mapvalues(data[,'cantor'], ck, genome)
    r2 <- summary(lm(model, data = data))$r.squared
    genome <- paste(genome, collapse = "")
    return(list("r2" = r2, "genome" = genome, "gen" = 1))
}

n_gen_model <- function(ck, data, model, n_sub, genome, gen) {
    data['genome'] <- plyr::mapvalues(data[,'cantor'], ck, genome)
    r2 <- summary(lm(model, data = data))$r.squared
    genome <- paste(genome, collapse = "")
    return(list("r2" = r2, "genome" = genome, "gen" = gen))
}

mutate <- function(genome, n_sub) {
    l <- length(genome)
    a <- ceiling(runif(3,0,l))
    b <- floor(runif(3,0,n_sub))
    genome[a] <- b
    return(genome)
}

geneObject <- function(dat, res, n_sub, n_gen, vdf) {
    dat_all <- dat
    dat <- gridData(dat, 'lonWGS84', 'latWGS84', 25)[c('X', 'Y', 'cantor', modelVec(c('resp', 'gene_model'), vdf))]
    ck <- levels(dat[,'cantor'])
    bs <- bootstrap(dat, 3, c(0,0.6,0.8,1))
    train <- dat[bs == 1,]
    test <- dat[bs == 2,]
    validate <- dat[bs == 3,]

    model <- modelSpec('resp', c('genome', 'gene_model'), vdf)

    pop <- 52

    train_fit <- c()
    test_fit <- c()
    validation_fit <- c()

    ## first gen model
    gen <- 1
    models <- list()
    for(i in 1:pop) {
        models[[i]] <- first_gen_model(ck, train, model, n_sub)
    }
    m_df <- as.data.frame(data.table::rbindlist(lapply(models, as.list)))

    ## n gen models
    while(gen <= n_gen) {
        m_df <- m_df[with(m_df, order(-r2))[1:(pop/2)],]

        train_fit[gen] <- m_df[1,'r2']
        test_fit[gen] <- n_gen_model(ck, test, model, n_sub, unlist(strsplit(m_df[1,'genome'], character(0))), gen)['r2']
        validation_fit[gen] <- n_gen_model(ck, validate, model, n_sub, unlist(strsplit(m_df[1,'genome'], character(0))), gen)['r2']
        print(paste('Generation:', gen, ". Max r2: ", m_df[1,'r2'], sep = ""))

        gen <- gen + 1

        models <- list()

        xlen <- ceiling(length(ck)/2)
        x <- lapply(strsplit(m_df[,'genome'], character(0)), function(x) { paste(x[1:xlen], collapse = "") })
        y <- lapply(strsplit(m_df[,'genome'], character(0)), function(x) { paste(x[(xlen+1):length(ck)], collapse = "") })
        offspring_genomes <- lapply(strsplit(c(paste(x[1:(pop/4)],y[((pop/4)+1):(pop/2)], sep = ""), paste(y[1:(pop/4)],x[((pop/4)+1):(pop/2)], sep = "")), character(0)), as.numeric)
        offspring_genomes <- lapply(offspring_genomes, function(x) { mutate(x, n_sub) })

        for(i in 1:length(offspring_genomes)) {
            models[[i]] <- n_gen_model(ck, train, model, n_sub, offspring_genomes[[i]], gen)
        }

        o_df <- as.data.frame(data.table::rbindlist(lapply(models, as.list)))
        m_df <- rbind(m_df, o_df)
    }

    ## evolution
    evolution <- as.data.frame(cbind('generation' = seq(1:n_gen), 'train' = unlist(train_fit), 'test' = unlist(test_fit), 'validation' = unlist(validation_fit)))

    model <- modelSpec('resp', c('cantor', 'gene_model'), vdf)
    theoretical_max <- summary(lm(model, data = train))$r.squared

    model <- modelSpec('resp', 'gene_model', vdf)
    baseline <- summary(lm(model, data = train))$r.squared

    ## winner
    winner <- as.numeric(strsplit(m_df[1, 'genome'], character(0))[[1]])
    gene_mapping <- as.data.frame(cbind('genome' = as.numeric(strsplit(m_df[1, 'genome'], character(0))[[1]]), 'cantor_key' = ck))

    ## data
    model <- modelSpec('resp', c('genome', 'gene_model'), vdf)
    dat$genome <- plyr::mapvalues(dat[,'cantor'], ck, winner)
    dat$genome <- plyr::mapvalues(
                            dat[,'genome'],
                            levels(dat$genome)[sort.list(by(dat, dat[,'genome'], function(x) {
                                mean(predict(lm(model, data = dat), x))
                            } ))],
                            seq(1:length(levels(dat[,'genome'])))
                        )

    return(
        list('final_generation' = m_df,
             'evolution' = evolution,
             'plot_lines' = c('baseline' = baseline, 'max' = theoretical_max),
             'winner' = gene_mapping,
             data = cbind(dat_all[which(!names(dat_all) %in% names(dat))], dat)
             )
    )
}


## plots <- lapply(sample(1:length(c_out$mstrataID), 12), function(i) {
##     return(
##         ggmap(map)
##         + geom_tile()
##         + geom_point(data = dat[which(c_out$w > 0 & c_out$mstrata == c_out$mstrataID[i]),], aes(x = lonWGS84, y = latWGS84, color = flood, shape = "1",size = c_out$w[which(c_out$w > 0 & c_out$mstrata == c_out$mstrataID[i])]))
##         + scale_shape(solid = FALSE)
##         + scale_size_area(limits = c(0,max(c_out$w[which(c_out$w > 0 & c_out$mstrata == c_out$mstrataID[i])])))
##         + theme(legend.position="none"))
##     })

plotGrid <- function(data, map, var, stratum = NA, strata = NA, split = 0) {
    alpha <- 0.5
    counts <- as.data.frame(t(matrix(unlist(by(data, data$cantor, function(x) { c(x$X[1], x$Y[1], mean(as.numeric(as.character(x[,var])))) })), nrow = 3)))
    names(counts) <- c('X', 'Y', 'count')

    if(is.na(stratum) & is.na(strata)) {

        plot <- ggmap(map) + geom_tile(data = counts, aes(x = X, y = Y, fill = count), alpha = alpha) + scale_fill_gradientn(colors = c('purple','orange','yellow'), name = 'Spatial submarket') + labs(x = "Longitude (WGS84)", y = "Latitude (WGS84)")

    } else if(is.na(strata) & !is.na(stratum)) { ## One straum

        plot <- ggmap(map) + geom_tile(data = counts, aes(x = X, y = Y, fill = count), alpha = alpha) + scale_fill_gradientn(colors = c('purple','orange','yellow'), name = 'Spatial submarket') + geom_point(data = stratum, aes(lonWGS84, latWGS84, col = flood, size = w, shape = '1')) + labs(x = "Longitude (WGS84)", y = "Latitude (WGS84)") + scale_shape(solid = FALSE)

    } else { ## All strata

        plot <- ggmap(map) + geom_tile(data = counts, aes(x = X, y = Y, fill = count), alpha = alpha) + scale_fill_gradientn(colors = c('purple','orange','yellow'), name = 'Spatial submarket') + geom_point(data = strata, aes(lonWGS84, latWGS84, color = ifelse(stt > split, 1, 0))) + scale_color_gradientn(colors=rainbow(2)) + coord_cartesian(xlim = range(strata$lonWGS84), ylim = range(strata$latWGS84)) + labs(x = "Longitude (WGS84)", y = "Latitude (WGS84)")

        ## plot <- ggplot(counts, aes(X,Y))
        ## + geom_polygon(data = spdf_poly_f, aes( x = long, y = lat, group = group), fill = 'white', color='black', size = 0.2)
        ## + geom_tile(aes(fill = count), alpha = alpha)
        ## + scale_fill_gradientn(colors = c('purple','orange','yellow'))
        ## + geom_point(data = strata, aes(lonWGS84, latWGS84, color = stt))
        ## + coord_cartesian(xlim = range(strata$lonWGS84), ylim = range(strata$latWGS84))

}
return(plot)
}

plotGeneticSimulation <- function(obj) {
    plot(obj$evolution$train ~ obj$evolution$generation,
         type = 'l',
         col = 'red',
         ylim = c(obj$plot_lines['baseline'], obj$plot_lines['max']),
         ylab = 'r2',
         xlab = 'generation')
    with(obj$evolution, lines(test ~ generation, col = 'blue', lty = 2))
    with(obj$evolution, lines(validation ~ generation, col = 'green', lty = 2))
    abline(h = obj$plot_lines['baseline'])
    abline(h = obj$plot_lines['max'])
}

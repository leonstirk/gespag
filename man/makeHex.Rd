% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachCellID.R
\name{makeHex}
\alias{makeHex}
\title{Create hex binnings as an sf object}
\usage{
makeHex(df, lon_colname, lat_colname, resolution, crs)
}
\arguments{
\item{df}{data.frame}

\item{lon_colname}{character}

\item{lat_colname}{character}

\item{resolution}{integer}

\item{crs}{CRS}
}
\value{
SpatialPolygons
}
\description{
Create hex binnings as an sf object
}

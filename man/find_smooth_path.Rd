% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{find_smooth_path}
\alias{find_smooth_path}
\title{Combined Path Processing Pipeline}
\usage{
find_smooth_path(df, na.fill = TRUE, p = 3, n = 13)
}
\arguments{
\item{df}{Input dataframe with raw trajectory data}

\item{na.fill}{Logical passed to smooth_path() for NA handling}

\item{p}{Polynomial order for Savitzky-Golay filter (must be < n)}

\item{n}{Window size (number of points) for Savitzky-Golay filter (must be odd)}
}
\value{
Fully processed dataframe with smoothed coordinates
}
\description{
Wrapper function that sequentially applies path identification and smoothing.
}

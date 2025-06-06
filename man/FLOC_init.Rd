% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLOC.R
\name{FLOC_init}
\alias{FLOC_init}
\title{FLOC initialize}
\usage{
FLOC_init(data, NJ, NK, rhoj, rhok)
}
\arguments{
\item{data}{numerical vector; historical data that is available to estimate
the prechange distribution from.}

\item{NJ}{integer; number of data points in a jump bin}

\item{NK}{integer; number of data points in a kink bin}

\item{rhoj}{numerical; threshold for the jump part of the algorithm.}

\item{rhok}{numerical; threshold for the kink part of the algorithm.}
}
\value{
list with
detect_jump - Boolean; True if a detection was made by the jump test statistic
in this iteration
detect_kink - Boolean; True if a detection was made by the kink test statistic
in this iteration
iteration - integer; number of the iteration
SJ - numerical vector; Sums of the bins for a jump
SK - numerical vector; Sums of the bins for a kink
W - numerical vector; weighted sums of the bins
Tstat - numerical vector; the current test statistic
f - object of class "lm"; estimate of the pre change distribution
NJ - integer; number of data points in a jump bin
NK - integer; number of data points in a kink bin
rhoj - numerical; threshold for the jump part of the algorithm.
rhok - numerical; threshold for the kink part of the algorithm.
}
\description{
Initializes the FLOC, Fast Limited-memory Optimal Change detector.
}

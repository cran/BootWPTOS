\name{WPts}
\alias{WPts}
\title{Compute wavelet packet test statistic for test of stationarity.
}
\description{Computes nondecimated wavelet packet transform of time series.
	Computes b-spectrum (square of nondecimated WP transform) for
	various levels and indices (controlled by \code{levs} and
	\code{indices} arguments). Computes variance (L2 norm) of
	the b-spectra and averages them. Returns the average.
}
\usage{
WPts(x, levs, indices, filter.number = 1, family = "DaubExPhase")
}
\arguments{
  \item{x}{The time series whose statistic you want to compute.
}
  \item{levs}{The b-spectrum levels you want to use.
}
  \item{indices}{The b-spectrum indices you want to use.
}
  \item{filter.number}{The filter number of the underlying wavelet.
}
  \item{family}{The family of the underlying wavelet.
}
}
\details{Description says it all. However, the \code{levs} and
\code{indices} warrant further explanation.

Our code is designed to be used on data sets that are a power
of two, i.e. T = 2^J for some J (note: the test can work on
other values of T but coding is more finickity). Given a series of this
length there are J levels, labelled 0 (coarse) to J-1 fine. Within each
level there are J-j packets indexed 0, 1, ..., J-j-1 for scales
J-1, ..., 0 respectively.

Packet 0 within any scale always corresponds to the father wavelets
at that scale and we don't tend to use this for stationarity testing.
Packet 1 within any scale always corresponds to mother wavelets. We often
use these. Note, at the finest scale J-1 there are only two packets
0 (father wavelet) and 1 (mother wavelet) coefficients.
}
\value{The test statistic value is returned.
}
\references{
Cardinali, A. and Nason, G.P. (2016) Practical Powerful Wavelet Packet Tests
        for Second-Order Stationarity. Applied and Computational Harmonic
        Analysis, 2016.
	\doi{10.1016/j.acha.2016.06.006}
}
\author{
G.P. Nason
}

\seealso{
\code{\link{BootWPTOS}}
}
\examples{
#
# Generate some test data
#
x <- rnorm(512)
#
# Compute the test statistic on mother wavelets and packets from the finest
# scale and the THIRD finest scale
# 
J <- IsPowerOfTwo(length(x))
J
# [1] 9
#
x.ts <- WPts(x, levs=c(8, rep(6, 7)), indices=c(1, 1:7))
x.ts
# [1] 1.792252
}
\keyword{ts}

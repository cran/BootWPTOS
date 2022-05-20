\name{BootWPTOS}
\alias{BootWPTOS}

\title{Compute test of stationarity for time series via bootstrap
	wavelet packet method.
}
\description{Computes b-spectrum (wavelet packet periodogram) for
	predefined set of wavelet packets on time series. Then applies
	bootstrap method to resample new versions of the (assumed for the
	test) stationary series and retests the series. If the value of
	the test statistic is out of line (bigger) than the resampled
	test statistics then the series might well not be stationary.
}
\usage{
BootWPTOS(x, levs, indices, filter.number = 1, family = "DaubExPhase",
	Bsims = 200, lapplyfn = lapply, ret.all = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The time series you wish to test. Length is a power of two
}
  \item{levs}{The levels of the wavelet packets that you want to involve
	in the test.
}
  \item{indices}{The indices of the wavelet packets that you want to involve
	in the test.
}
  \item{filter.number}{The filter number of the wavelet that underlies the
	wavelet packet.
}
  \item{family}{The family of the wavelet that underlies the wavelet
	packet used for the test.
}
  \item{Bsims}{The number of bootstrap simulations.
}
  \item{lapplyfn}{By default this argument is the \code{lapply} function
	which operates sequentially. However, if you have the 
	\pkg{parallel} package you could supply the \code{mclapply}
	function which processes list elements in parallel on a multicore
	machines.
}
  \item{ret.all}{If \code{FALSE} then the results of the test
	are returned in the standard \code{R} \code{htest}
	object, and can be printed and manipulated by those
	standard tools. If \code{TRUE} then a list is returned
	with information about the test.
}
}
\details{Function computes a test statistic for test of stationarity
	on a time series. Then successive bootstrap realizations are
	drawn using the \code{surrogate} function. If the original time
	series WAS stationary then \code{surrogate} causes stationary
	draws with the same spectral characteristic as the data to be
	produced (with Gaussian marginals). Under the null hypothesis
	the time series is assumed stationary and so the distribution of
	all of the test statistics should be the same and the p-value
	of the test statistic be uniformly distributed. If the series
	is nonstationary then the value of the statistic is likely
	to be bigger on the first computed test statistic on the dats
	and much bigger than all the others. We can work out a bootstrap
	p-value by counting how many resampled test statistics are bigger
	than the one computed on the data.
}
\value{Normally a list, cast as a \code{htest} class object with the
following components:
\item{statistic}{The test statistic computed on the data}
\item{p.value}{The bootstrap p.value of the test.}
\item{method}{The name of the method of this test statistic.}
\item{data.name}{The name of the data set tested.}
\item{Bootvals}{The remaining bootstrap generated test statistics.}
}
\references{
Cardinali, A. and Nason, G.P. (2016) Practical Powerful Wavelet Packet Tests
        for Second-Order Stationarity. Applied and Computational Harmonic
        Analysis, 2016
	\doi{10.1016/j.acha.2016.06.006}
}
\author{
G.P. Nason
}

\seealso{\code{\link{WPts}}, \code{\link{WPTOSpickout}}}

\examples{
#
# Generate a stationary time series (e.g. iid standard normals)
#
x <- rnorm(512)
#
# What would be the finest scale?
#
J <- IsPowerOfTwo(length(x))
J
#[1] 9
#
# So, in WaveThresh there are 9 scales indexed 0 to 8.
#
# Let us test x for stationarity
#
# The finest scale wavelets (or packets) are at scale 8
# The next finest scale is 7.
#
# Wavelets themselves are always indexed 1, father wavelets 0.
# We don't tend to use father wavelets for stationary testing.
#
# There are 2^j packets at scale J-j (so 2 at the finest [father and
# mother], 4 at the next finest [father=0, mother=1, packets 2 and 3].
#
# Let's just look at the finest scale wavelet (8,1) and the next finest
# scale wavelet (7,1) and two other wavelet packets  (7,2) and (7,3)
#
x.test <- BootWPTOS(x=x, levs=c(8,7,7,7), indices=c(1,1,2,3), Bsims=30)
#
# Note: Bsims=30 is almost certainly too small (but it is small here because
# on installation R run these examples and I don't want it to take too long.
# 100+ is almost certainly necessary, and probably 500+ useful and 1000+
# to be "sure". If you can load the multicore library then you can
# replace lapplyfn=lapply with lapplyfn=mclapply to get a parallel processing
#
# What are the results of our test?
#
x.test
#
#	WPBootTOS test of stationarity
#
#data:  x
#= 1.8096, p-value = 0.7
#
# So, the p-value is > 0.05 so this test indicates that there is
# no evidence for non-stationarity. Running it for 1000 bootstrap simulations
# gave a p-vale of 0.736.
#
# The next example is nonstationary. However, after the series has been
# generated you should plot it. The second half has a different variance
# to the first half but it is very difficult, usually, to identify the
# different variances on a plot.
#
x2 <- c(rnorm(256), rnorm(256,sd=1.4))
#
# Let's do a test, but involve ALL non-father-wavelet packets from scales
# 8, 7, 6 and 5.
#
x2.test <- BootWPTOS(x=x2, levs=c(8,7,7,7,rep(6,7), rep(5,15)),
	indices=c(1,1,2,3, 1:7, 1:15), Bsims=30)
x2.test
#
#	WPBootTOS test of stationarity
#
#data:  x2
#= 5.4362, p-value < 2.2e-16
#
# So, strong evidence for nonstationarity because p.value < 0.05 (much less
# than!). Again here we've only use 30 bootstrap simulations and this is
# probably too small. Using Bsims=1000 and mclapply (for speed) gave a p-value
# of 0.002, so still assessed to be nonstationary, but we have more confidence
# in the answer.
}
\keyword{ts}

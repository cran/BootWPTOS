\name{WPTOSpickout}
\alias{WPTOSpickout}
\title{Find nonstationarities in a time series
}
\description{The nonstationarities are located by looking for significant
Haar wavelet coefficients of a b-spectrum of a time series (a b-spectrum
is the expectation of a wavelet packet periodgram). The significant
Haar coefficients can locate discontinuities in space and time.
}
\usage{
WPTOSpickout(x, level, index, filter.number = 1, family = "DaubExPhase",
	plot.it = FALSE, verbose = FALSE, lowlev = 3, highlev, nomsize = 0.05)
}

\arguments{
  \item{x}{The time series you wish to analyse.
}
  \item{level}{The level of the b-spectrum you want to examine.
	See help for \code{\link{BootWPTOS}} for more information
	on levels.
}
  \item{index}{The index of the b-spectrum you want to examine.
	See help for \code{\link{BootWPTOS}} for more information
	on indices you can chose.
}
  \item{filter.number}{The filter number of the underlying wavelet
	you wish to examine.
}
  \item{family}{The family of the underlying wavelet you wish to examine.
}
  \item{plot.it}{If \code{TRUE} this plots the time series \code{x} as
	supplied. Then, superimposed in red, is the raw wavelet packet
	periodogram. The Haar wavelet coefficients of the red signal are
	examined for whether they are deemed non-zero or not using
	asymptotic normality results. If \code{FALSE} then no such plot
	is produced.
}
  \item{verbose}{If \code{TRUE} then a single line indicating the
	number of significant coefficients found is printed.
	If \code{FALSE} then the function prints nothing.
}
  \item{lowlev}{Keep away from coarse scales. Typically, Haar wavelet
	coefficients at the coarse scales are contaminated by boundary
	effects. These won't usually cause a problem at scales 3 or higher,
	or maybe 2. Only coefficients at scales \code{lowlev} or finer
	will be selected for testing.
}
  \item{highlev}{Keep away from fine scales. The testing of Haar
	wavelet coeffients depends on utilizing enough data points
	to enable asymptotic normality to kick in. At the finest scale
	only TWO single points are compared and the distribution of each
	point might be far from normality. At coarser scales TWO averages
	are compared and those averages will consist of many points.
	E.g. at the third finest scale each average will consist of four
	points, at the fourth finest scale each average will consist of
	8 points, etc. By default this argument is set to be roughly the
	use a scale one level finer than the halfway number of levels.
	So, if J=10, then highlev is 6. The formula is
	floor(J/2)+1. Note: highlev and lowlev should be specified
	in the WaveThresh scaling (e.g. scale 0 is the coarsest scale) 
}
  \item{nomsize}{The nominal size of the test as a number between
	0 and 1. So, if you want a 5% test then enter 0.05.
}
}
\details{This function computes the nondecimated wavelet packet
	transform of the packet you specified by
	\code{level} and \code{index}. Note: you can only specify
	one number for each of these. Then the b-spectrum (raw
	wavelet packet periodogram) is formed by squaring the
	nondecimated wavelet packet transform. Then the Haar
	wavelet coefficients are obtained for the b-spectrum and
	a multiple hypothesis test is performed on all the Haar
	wavelet coefficients between scales lowlev and highlev.
	The function return information about any significant wavelet
	coefficients.
}
\value{A list of class \code{toswp} containing the following
components:
\item{x}{The time series that was analyzed}
\item{level}{The level of the b-spectrum that we wanted}
\item{index}{The index of the b-spectrum that we wanted}
\item{sigcoefs}{A wd class object containing the significant Haar
	wavelet coefficients, if there are any}
\item{nreject}{The number of significant Haar wavelet coefficients}
\item{ntests}{The total number of hypothesis tests carried out in the
	multiple hypothesis test}
\item{bonsize}{The Bonferroni corrected rate for the multiple hypothesis
	test}
}
\references{
Cardinali, A. and Nason, G.P. (2016) Practical Powerful Wavelet Packet Tests
        for Second-Order Stationarity. Applied and Computational Harmonic
        Analysis, 2016, DOI=10.1016/j.acha.2016.06.006.
}
\author{
G.P. Nason
}


\seealso{
\code{\link{BootWPTOS}}, \code{\link{plot.toswp}}, \code{\link{print.toswp}},
	\code{\link{summary.toswp}}
}
\examples{
#
# Create some simulated data
#
x2 <- c(rnorm(256), rnorm(256, 2))
#
# The following call to BootWPTOS (generic tester) 
#
\dontrun{x2.tos <- BootWPTOS(x2, levs=5, indices=1, Bsims=500)}
\dontrun{x2.tos}
#
#	WPBootTOS test of stationarity
#
#data:  x2
#= 5.2826, p-value < 2.2e-16
#
# So, test indicates that strong evidence for nonstationarity. 
#
# Now let's do the multiple Haar hypothesis test.
#
x2.po <- WPTOSpickout(x=x2, level=7, index=1)
x2.po
#Class 'toswp' : Wavelet Packet Test of Stationarity Object :
#       ~~~~  : List with 7 components with names
#              x level index sigcoefs nreject ntests bonsize 
#
#
#summary(.):
#----------
#Number of individual tests: 56 
#Bonferroni p-value was: 0.0008928571 
#Tests rejected: 2 
#Listing Bonferroni rejects...
#Wavelet Packet  (5,1): HWTlev: 4. Indices: 8  
#Wavelet Packet  (5,1): HWTlev: 5. Indices: 16 
#
# So, this test also shows nonstationarities. For this packet (5,1)
# two significant Haar coefficients were identified. One was at level 4
# position 8 and the other was at scale level 5 position 16.
#
# You can plot them also
#
\dontrun{plot(x2.po)}
#
# You should get a nice plot of the time series with double-headed red
# arrows indicating the location and extent of the nonstationarities.
# For this example, where the spectrum changes dramatically at the halfway
# point - this is where the arrows should be located. Of course, with random
# data you might see other arrows in other locations, but this should be
# unlikely and on repeating the above they should not persist.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

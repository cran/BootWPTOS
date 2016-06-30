\name{BootWPTOS-package}
\alias{BootWPTOS-package}
\docType{package}
\title{Bootstrap Wavelet Packet Tests of Stationarity
}
\description{This package contains two main functions to carry out
	tests of second-order stationarity using wavelet packets.
	One test, \code{\link{BootWPTOS}} carries out the
	bootstrap wavelet packet test of stationarity as described by
	Cardinali and Nason (2016), with algorithm of that same name
	in that paper. The test is carried out respect to a set of
	wavelet packets (one, or more than one).

	The other main function is \code{\link{WPTOSpickout}}.
	Here, the test is carried out using a fixed single wavelet
	packet and inference for test statistics is carried out
	using asymptotic normal approximations as described in
	Cardinali and Nason (2016) but based on ideas in 
	von Sachs and Neumann (2000).
}
\details{
\tabular{ll}{
Package: \tab BootWPTOS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-03-30\cr
License: \tab GPL-2\cr
}
The main functions are documented above. See below for an example of
each.

Both functions require the specification of a set (for \code{\link{BootWPTOS}})
	or a single (for \code{\link{WPTOSpickout}}) wavelet packet.
	This is because the tests use and rely on wavelet packets.

Wavelet packets are indexed by two quantities: scale and index. The scale
is referred to in the functions by the \code{levs} and \code{level} arguments
respectively. Scale can be any scale that you would normally use in
wavethresh. So, for a series of dyadic length, that is T=2^J, the scales
are indexed 0 (coarsest scale) to J-1 (finest scale).

The range of the index argument for a wavelet packet depends on the scale.
Always the scaling function coefficients have index 0 and the regular wavelet
coefficients always have index 1. Note: wavelets are a subset of wavelet
packets. Then, for scale J-j there are 2^j packets.

So, for example, at the finest scale J-1 there are 2^1=2 packets. These
correspond to indices 0 and 1, the father and mother wavelet coefficients
respectively. For the next finest scale J-2 there are 2^2=4 packets. These
are indexed 0, 1, 2 and 3 with 0,1 being the father/mother wavelet
(packet) coefficients at that scale and indices 2 and 3 being the other
two wavelet packets at that scale. Clearly there are many more wavelet
packets at coarser scales. 

In our functions the levs or level contains the scale of any wavelet
packet and the indices or index variable contains the indices.

In \code{\link{BootWPTOS}} you can use any combination of wavelet packets,
but it is important that the entries correspond to each other in the
levs and indices vector. E.g. if you wanted wavelet packet (3,5) and (4,7)
with J=5 then you would use the arguments levs=c(3,4) and indices=c(5,7).
}
\author{
Guy Nason

Maintainer: Who to complain to <g.p.nason@bristol.ac.uk>
}
\references{
Cardinali, A. and Nason, G.P. (2016) Practical Powerful Wavelet Packet Tests
	for Second-Order Stationarity. Applied and Computational Harmonic
	Analysis, 2016, DOI=10.1016/j.acha.2016.06.006. 

Von Sachs, R. and Neumann, M.H. (2000) A Wavelet-Based Test for Stationarity.		Journal of Time Series Analysis, 21, 597-613. 
}
\keyword{ package }
\keyword{ts}
\seealso{
\code{\link{BootWPTOS}},
\code{\link{plot.toswp}},
\code{\link{print.toswp}},
\code{\link{summary.toswp}},
\code{\link{WPTOSpickout}}
}
\examples{
#
# First, we provide an example concerning BootWPTOS
#
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
# So, in WaveThresh there are 9 scaled indexed 0 to 8.
#
# Test x for stationarity
#
# The finest scale wavelets are at 8
# The next finest scale is 7.
#
# Wavelets themselves are always indexed 1, father wavelets 0.
# We don't tend to use father wavelets for stationary testing.
#
# There are 2^j packets at scale J-j (so 2 at the finest [father and
# mother], 4 at the next finest [father=0, mother=1, packets 2 and 3].
#
# Let's just look at the finest scale wavelets and the next finest
# scale wavelets and two other wavelet packets.
#
\dontrun{x.test <- BootWPTOS(x=x, levs=c(8,7,7,7), indices=c(1,1,2,3), Bsims=30)}
#
# Note: Bsims=30 is ALMOST CERTAINLY TOO SMALL (but it is small here because
# on installation R run these examples and I don't want it to take too long.
# 100+ is almost certainly necessary, and probably 500+ useful and 1000+
# to be "sure". If you can load the multicore library then you can
# replace lapplyfn=lapply with lapplyfn=mclapply to get a parallel processing
#
# What are the results of our test?
#
\dontrun{x.test}
#
#       WPBootTOS test of stationarity
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
#
\dontrun{x2.test <- BootWPTOS(x=x2, levs=c(8,7,7,7,rep(6,7), rep(5,15)),
	indices=c(1,1,2,3, 1:7, 1:15), Bsims=500)}
\dontrun{x2.test}
#
#       WPBootTOS test of stationarity
#
#data:  x2
#= 5.4362, p-value < 2.2e-16
#
# So, strong evidence for nonstationarity because p.value < 0.05 (much less
# than!). 
# Using Bsims=1000 and mclapply (for speed) gave a p-value
# of 0.002, so still assessed to be nonstationary, but we have more confidence
# in the answer.
#
# Now we provide an example of using WPTOSpickout
#
#
# Create some simulated data
#
x2 <- c(rnorm(256), rnorm(256, 2))
#
# Note: x2 should be highly nonstationary. The left-hand half of the series
# has variance 1, the right-hand half has variance 2.
#
# First, try out a wavelet packet test of stationarity. This is a check
# on the later test. You should really check both.
#
# We've chosen this packet more or less at random. Its packet at scale 5 and
# index 1 (this happens to be a wavelet)
#
\dontrun{x2.tos <- BootWPTOS(x2, levs=5, indices=1, Bsims=500)}
\dontrun{x2.tos}
#
#       WPBootTOS test of stationarity
#
#data:  x2
#= 5.2826, p-value < 2.2e-16
#
# So, test indicates that strong evidence for nonstationarity.
#
# Now let's do the multiple Haar hypothesis test.
#
\dontrun{x2.po <- WPTOSpickout(x=x2, level=7, index=1)}
\dontrun{x2.po}
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

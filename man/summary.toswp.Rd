\name{summary.toswp}
\alias{summary.toswp}
\title{Summarize a \code{toswp} class object.
}
\description{
This function goes through a \code{toswp} class object and
printing out details of which Haar coefficients were significant.
}
\usage{
\method{summary}{toswp}(object, quiet = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The object you wish to summarize.
}
  \item{quiet}{If \code{TRUE} nothing is printed. However, a list
	of the significant coefficients is returned. This information
	is used, for example, by \code{\link{plot.toswp}}.
}
  \item{\dots}{Other arguments to summary.
}
}
\details{None
}
\value{
A list with the following components:
\item{rejlist}{A list with details on the rejected coefficients.
	Each component of the list is a vector. The first
	element of each vector is the Haar wavelet coefficient
	scale level. The remaining numbers are the indices of
	any significant coefficients at that level.}
\item{nreject}{The number of rejected coefficients}
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
\code{\link{WPTOSpickout}}, \code{\link{plot.toswp}},
	\code{\link{print.toswp}}
}
\examples{
#
# See example of \code{\link{print.toswp}} in the help for
# \code{\link{WPTOSpickout}} which includes a call to this function.
}
\keyword{ts}

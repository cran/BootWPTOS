\name{plot.toswp}
\alias{plot.toswp}
\title{Plot \code{toswp} class object.
}
\description{Plots the time series that was analyzed to produce
	the \code{toswp} class object. Then superimposes the location
	and extent of nonstationarities by means of double-headed red
	arrows. The right-hand axis indicates the scale of the
	significant Haar wavelet coefficients corresponding to the
	nonstationary arrows.
}
\usage{
\method{plot}{toswp}(x, sub = NULL, xlab = "Time", arrow.length = 0.05, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The \code{toswp} class object that you wish to plot.
}
  \item{sub}{A subtitle.
} 
  \item{xlab}{The label for the x-axis.
}
  \item{arrow.length}{
	Length of arrow head.
}
  \item{verbose}{
	If \code{TRUE} prints debugging information from the plotting
	process.
}
  \item{\dots}{Other arguments to the plot function.
}
}
\details{As description says.
}
\value{Nothing of interest.
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
\code{\link{WPTOSpickout}}, \code{\link{print.toswp}},
	\code{\link{summary.toswp}}
}
\examples{
#
# See example in helpfile for \code{\link{WPTOSpickout}}
#
}
\keyword{ts}

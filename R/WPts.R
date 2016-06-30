WPts <-
function(x, levs, indices, filter.number=1, family="DaubExPhase")
{
#
# Computes wavelet packet test statistic on time series x
#
# On wavelet packets indexed by levels in levs and indices in indices.
#

#
# Compute nondecimated wavelet packet transform on series
#
xwpst <- wpst(x, filter.number=filter.number, family=family)

#
# See how many levels and indices we have to compute the statistic on
# and do some other argument checking (ie levels are positive and < J)
#

nlev <- length(levs)
nind <- length(indices)

if (nlev != nind)
	stop("Number of levels and number of indices has to be the same")

J <- nlevelsWT(xwpst)


if (any(levs < 0))
	stop("All levels have to be >= 0")
if (any(levs >= J))
	stop("All levels have to be < log_2(n)")

#
# Now check that indices for each level are correct
#
for(i in 1:nlev)	{

	newlev <- J - levs[i] 
	maxix <- 2^newlev

	if (indices[i] < 0)
		stop("All indices have to be >= 0")
	if (indices[i] >= maxix)
		stop(paste("Index ", indices[i], " is invalid for level ", levs[i], ". Max index is: ", maxix-1))
	}

#
# Now compute and return test statistic
#

the.ts <- 0

for(i in 1:nlev)	{

	xwp <- accessD(xwpst, level=levs[i], index=indices[i])^2

	the.ts <- the.ts + var(xwp)
	}

the.ts <- the.ts/nlev

return(the.ts)

}

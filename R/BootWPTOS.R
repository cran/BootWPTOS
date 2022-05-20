BootWPTOS <-
function(x, levs, indices, filter.number=1, family="DaubExPhase", Bsims=200, lapplyfn=lapply, ret.all=FALSE){

#
# Get the name of the data object, x, to be tested
#

DNAME <- deparse(substitute(x))

#
# Check for any illegal values
#

if (any(is.na(x)) || any(is.nan(x)) || any(is.infinite(x)))
	stop("NA/NaN/Inf found in x")

#
# Compute the wavelet packet test statistic on the actual data
#

TS <- WPts(x=x, levs=levs, indices=indices, filter.number=filter.number, family=family)


#
# Create a function to run the bootstrap
#

bsfn <- function(dummy, x, levs, indices, filter.number, family) {

	# Note: nothing is done with the dummy argument
	#
	# Compute surrogate of time series x
	#
	xs <- as.numeric(surrogate(x=x, ns=1, fft=TRUE, amplitude=TRUE)) 

	if (any(is.na(x)) || any(is.nan(x)) || any(is.infinite(x)))
		stop("NA/NaN/Inf found (in bsfn)")

	#
	# Compute test statistic on surrogate series
	#

	TS <- WPts(x=xs, levs=levs, indices=indices,
		filter.number=filter.number, family=family)
	return(TS)
	}

#
# Create a list where every entry is equal to the single number: Bsims-1
#

dummy.ip <- vector("list", Bsims - 1)

#
# Apply the bootstrap function to every entry in the dummy list
# Note, the bootstrap function takes the object x forms a surrogate
# applies test statistic to the surrogate and returns the test statistic. 
#

ans <- lapplyfn(dummy.ip, bsfn, x=x, levs=levs, indices=indices,
	filter.number=filter.number, family=family)

#
# Convert the answer list to a vector
#

ans <- unlist(ans)

#
# Append the value of the test statistic on the data to all the bootstrap vals
#

TS <- c(TS, ans)

#
# Work out the bootstrap p-value
#

p.value <- sum(TS[1] < TS[-1])/Bsims

#
# For debugging purposes return all the test statistics and the computed
# p-value in a list
#

if (ret.all==TRUE)	{
	l <- list(TS=TS, p.value=p.value)
	return(l)
	}

#
# Otherwise return the information in the form of a standard hypothesis
# test object.
#

	
htest.obj <- list(statistic = TS[1], p.value = p.value,
	method = "WPBootTOS test of stationarity", 
        data.name = DNAME, Bootvals = TS)

class(htest.obj) <- c("BootTOS", "htest")
return(htest.obj)
}

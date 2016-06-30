WPTOSpickout <-
function(x, level, index, filter.number=1, family="DaubExPhase",
	plot.it=FALSE, verbose=FALSE, lowlev=3, highlev, nomsize=0.05)	{

if (missing(highlev))	{
	J <- IsPowerOfTwo(length(x))
	highlev <- floor(J/2)+1
	}
#
# Compute wavelet packet transform of data
#
xwpst <- wpst(x, filter.number=1, family="DaubExPhase")
#
# Extract level and index number relating to packet you're interested in
#
xCoefs <- accessD(xwpst, level=level, index=index)
#
# Form b-spectrum (raw wavelet packet periodogram)
#
Ijk <- xCoefs^2
#
# Plot it if necessary
#
if (plot.it==TRUE)	{
	ts.plot(x)
	lines(Ijk, col=2)
	}

#
# Compute Haar wavelet transform of raw wavelet packet periodogram
#
Ijk.haar <- wd(Ijk, filter.number=1, family="DaubExPhase")
#
# Under null hypothesis each scale levels is N(0, sigma) .
#
# Estimate sigma for each scale.
#
sigma <- rep(0, J-1)
for(j in highlev:lowlev)	{
	sigma[j] <- mad(accessD(Ijk.haar, level=j))
	}
#
# Count how many coefficients we're going to test in total
#
totalc <- 0
for(j in highlev:lowlev)	{
	totalc <- totalc + 2^j
	}

#
# Work out Bonferonni size
#
	
mcsize <- nomsize/totalc

#
# Work out appropriate equivalent Z-value
#

z.mcsize <- abs(qnorm(mcsize/2))

ans.haar <- Ijk.haar

if (lowlev>0)
	ans.haar <- nullevels(ans.haar, 0:(lowlev-1))
if (highlev < J-1)
	ans.haar <- nullevels(ans.haar, (highlev+1):(J-1))

#
# Now do t test for each coefficient
#
survive_count <- 0
for(j in highlev:lowlev)	{
	y <- accessD(Ijk.haar, lev=j)
	z <- y/sigma[j]
	z [ abs(z) < z.mcsize] <- 0
	survive_count <- survive_count + sum(abs(z) > 0)
	ans.haar <- putD(ans.haar, lev=j, v=z)
	}

if (verbose==TRUE)
	cat("Number of Significant Coefficients: ", survive_count, "\n")

l <- list(x=x, level=level, index=index, sigcoefs=ans.haar, nreject=survive_count, ntests=totalc, bonsize=mcsize)
class(l) <- "toswp"
return(l)
}

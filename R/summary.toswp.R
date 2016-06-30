summary.toswp <-
function (object, quiet=FALSE, ...)
{
#
# Identify and return significant Haar wavelet coefficients
#


hwtosop <- object
ntests <- object$ntests
rejpval <- object$bonsize
bonreject <- object$nreject
sigcoefs <- object$sigcoefs
level <- object$level
index <- object$index

if (quiet == FALSE) {
	cat("Number of individual tests:", ntests, "\n")
	cat("Bonferroni p-value was:", rejpval, "\n")
	cat("Tests rejected:", bonreject, "\n")
    }

v <- NULL
count <- 1

if (bonreject != 0) {
	if (quiet == FALSE) 
		cat("Listing Bonferroni rejects...\n")

	sJ <- nlevelsWT(sigcoefs)

	for (k in 0:(sJ - 1)) {

		dk <- accessD(sigcoefs, level=k)

		ix <- which(dk != 0)

		if (length(ix) > 0) {
		  if (quiet == FALSE) {
	  	    cat("Wavelet Packet ", paste("(", level, ",",  index, "):", sep=""), paste("HWTlev: ", k,".", sep=""), "Indices: ")
		    for(i in 1:length(ix))
		      cat(ix[i], " ")
		    cat("\n")
                    }
	    	  v[[count]] <- c(k, ix)
		  count <- count + 1
                  }
		}
	vret <- list(rejlist = v, nreject = bonreject)
	}
else
	vret <- NULL

return(invisible(vret))
}

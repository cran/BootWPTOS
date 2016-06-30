plot.toswp <-
function (x, sub = NULL, xlab = "Time", arrow.length = 0.05, 
    verbose = FALSE, ...) 
{
    object <- x
    x <- x$x
    if (is.null(sub)) 
        sub <- paste("Packet:", paste("(", object$level, ",", object$index,")", sep=""), "Using Bonferonni:", object$nreject, 
            " rejected.")

    #
    # Plot the original time series in gray
    #

    ts.plot(x, xlab = xlab, sub = sub, col = "gray80", ...)

    st <- summary.toswp(object, quiet=TRUE)

    if (is.null(st))
	return(NULL)

    nreject <- st$nreject
    st <- st$rejlist

    stHlevs <- NULL

    for (i in 1:length(st)) {
        stHlevs <- c(stHlevs, st[[i]][1])
    }

    lyn <- min(stHlevs)
    lyx <- max(stHlevs)
    nHlevs <- length(lyn:lyx)
    ry <- range(x)
    mny <- ry[1]
    mxy <- ry[2]
    mainy <- seq(from = mny, to = mxy, length = nHlevs + 1)
    littley <- seq(from = 0, to = (mainy[2] - mainy[1]), length = lyx - 
        lyn + 2)
    if (verbose == TRUE) {
        cat("nHlevs: ", nHlevs, "\n")
        cat("mny, mxy: ", mny, mxy, "\n")
        cat("mainy: ")
        print(mainy)
    }
    abline(h = mainy[1:(length(mainy) - 1)], lty = 2)
    axis(4, at = mainy[1:(length(mainy) - 1)], labels = lyn:lyx)
    J <- IsPowerOfTwo(length(x))
    for (i in 1:length(st)) {
        stH <- st[[i]][1]
        ix <- st[[i]][c(-1)]
        for (j in 1:length(ix)) {
            xl <- 2^(J - stH) * (ix[j] - 1)
            xr <- 2^(J - stH) * (ix[j])
            yy <- mainy[stH - min(stHlevs) + 1]
	#+ littley[stH - lyn + 1]
            arrows(x0 = xl, x1 = xr, y0 = yy, y1 = yy, code = 3, 
                col = 2, length = arrow.length)
            if (verbose == TRUE) {
                cat("stH: ", stH, "\n")
                cat("[xl, xt] ", xl, xr, mainy[stH - min(stHlevs) + 
                  1], "\n")
                scan()
            }
        }
    }
}

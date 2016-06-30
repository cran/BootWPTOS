print.toswp <-
function (x, ...) 
{
    cat("Class 'toswp' : Wavelet Packet Test of Stationarity Object :\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.toswp(x)
}

# analyssi of simulated time series
# identification of inflexion points (storm peaks)
identify.peaks <- function(qsim)
{
	dq<-diff(qsim)
	change.in.sign <- rle(sign(as.vector(dq)))

	peaks <- cumsum(change.in.sign$lengths)
	peaks <- index(qsim)[peaks]
	return(peaks)
}



# identify inflexion point and determine difference between simulated and observed
# flows at these. as may not exactly coincide, need to aggregate and locate nearest value
analyse.peaks <- function(run)
{
	qsim <- aggregate_xts(run$qsim, dt=0.25)
	qobs <- aggregate_xts(run$qobs, dt=0.25)
	
	peaks <- identify.peaks(qsim)
	
	diff <- qsim[peaks]-qobs[peaks]
	
	diff.pc <- 100*(qsim[peaks]-qobs[peaks])/qobs[peaks]

	flood.peaks <- peaks[which(qsim[peaks]>1/1000)]
	
	diff.pc.flood <-  100*(qsim[flood.peaks]-qobs[flood.peaks])/qobs[flood.peaks]

}
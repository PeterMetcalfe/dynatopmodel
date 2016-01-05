
load.source("DynaTopMod.r")
TemporalResponseAnalysis <- function(proj, ntts=c(1, 3, 5, 10, 15), int=NULL, fn=NULL)
{
  qsim.ntt <- NULL #proj$qobs
  run.info <- NULL
  prev<-NULL
  conn <- stdout()
  if(!is.null(fn))
  {
    conn <- fn
   
  }
  if(!is.null(int))  
  {proj$disp.par$graphics.interval<- int}
  
  for(ntt in ntts)
  {    
    start.time <- Sys.time()
    proj$ntt<-ntt
    proj$disp.par$title.main <- paste0("ntt=",ntt)
    run.out <- run.proj(proj)
    run.len <- difftime(Sys.time(), start.time, units="secs")
    eff <- NSE(as.vector(run.out$qobs), as.vector(run.out$qsim))
    sigma.q <- sum(run.out$qsim)*proj$dt
    if(!is.null(prev))
    {
      diff.prev <- mean(abs(run.out$qsim-prev))
      diff.prev.pc <- mean(abs((run.out$qsim-prev)/run.out$qsim))/100
    }
    else
    {
      diff.prev<-NA
      diff.prev.pc <- NA
    }
    # constriuct info for this run
    ri <- c(#"ntt"=ntt, 
    				"dt"=60*proj$dt/ntt, 
    		#		"eff"=round(eff, 2), 
    		#		"run.len"=round(run.len), 
            "wb_pc"=100*run.out$wb/sigma.q, 
    				"sigma.q"=sigma.q, 
    				"ovf"= 1000*run.out$ovf, 
    		#		"deltaqbar"=signif(diff.prev, 3), 
    				"qbqn1_pc"=diff.prev.pc)
 #   ri <- signif(ri,3)
    run.info <- rbind(run.info, ri)
    cat(paste0(ri, "\t"), "\n", file=conn, append=T)
   # write.table(run.info, sep="\t", conn)
    qsim.ntt <-cbind(qsim.ntt, run.out$qsim[,1])
    prev <- run.out$qsim
  } 
  
  # add in the difference between individual runs and the final run
  n.sim <- length(ntts)
  run.info <-data.frame(cbind(run.info, 
                #   "qbqn"= get.qbqn(qsim.ntt),
  								 "qbqn_pc"=get.qbqn(qsim.ntt, pc=T)))
                     
                  #   sapply(1:n.sim, function(i){sum(abs(qsim.ntt[,i]-qsim.ntt[,n.sim]))/nrow(qsim.ntt)}))
  run.info<-signif(run.info, 2)
 	nresp <- ncol(qsim.ntt)
 	run.info$qbqn1_pc <- c(NA, colMeans(abs(qsim.ntt[,1:(nresp-1)]-qsim.ntt[,2:nresp])/qsim.ntt[,2:(nresp)]*100))

 colnames(qsim.ntt)<-ntts 
  write.table(run.info, row.names=F, sep="\t", file=conn)
	

	 if(!is.null(fn))
	 {	
	 	
	 	qsim.fn <- paste0(fn, ".tsv")
	 	write.zoo(qsim.ntt, sep="\t", file=qsim.fn)
	 }
  return(list("proj"=proj, "run.info"=run.info, "qsim"=qsim.ntt))
}


get.qbqn <- function(qsim.ntt, pc=F)
{
  n.sim <- ncol(qsim.ntt)
  qbqn<-sapply(1:n.sim, 
         function(i)
         {
         	if(pc)
         	{
         		return(mean(abs(qsim.ntt[,i]-qsim.ntt[,n.sim])/qsim.ntt[,n.sim])*100)
         	}
           mean(abs(qsim.ntt[,i]-qsim.ntt[,n.sim]))
         }
  )
  return(signif(qbqn, 3))
}


resp.proj <- function(proj, cuts=c(1, 2, 4, 8,  10, 12, 15, 20),  
																	 
														 			rebuild=F,  # use existing discretisations or rebuild from scratch, possibly overwriting parameter
																	 area.thresh=1,...)
{


	# 
#	if(rebuild)
#	{
		proj$disc<-list()
		i.disc <- 0
		for(n in cuts)
		{    
			cuts.n <- list("a"=n)
			disc <- add.disc(proj, cuts=cuts.n, rebuild=rebuild, area.thresh=area.thresh,...) 
			proj$disc[[i.disc+1]] <- disc
			i.disc <- length(proj$disc)    
			
			# build the entire project (exclude reaches once created)
			#	  disc <- disc.from.dir(proj$dir, dem=proj$dem, drn=proj$drn, reaches=proj$reaches,
			#                           cuts=cuts.n, ...)
		}
#	}
	return(proj)
}

sp.resp.analysis <- function(proj, cuts=c(1, 2, 4, 8,  10, 12, 15, 20),  										
														 fn.out=NULL,
														 rebuild=T,  # use existing discretisations or rebuild from scratch, possibly overwriting parameter
														 int=NULL,area.thresh=1,...)
{
	qsim.ntt <- NULL #proj$qobs
	run.info <- NULL
	prev<-NULL
	conn <- stdout()
	if(!is.null(fn.out))
	{
		conn <- fn.out #file(fn.out, open="w")
	#	on.exit(close(conn))
	}	

	runs <- lapply(proj$disc,
		function(disc)
		{        	      	 		
			start.time <- Sys.time()	
			proj$disp.par$title.main <- paste(disc$cuts, collapse=" ")
			run.out <- run.proj(proj, disc=disc)
			run.len <- difftime(Sys.time(), start.time, units="secs")
			eff <- NSE(as.vector(run.out$qobs), as.vector(run.out$qsim))
			sigma.q <- sum(run.out$qsim)
			if(!is.null(prev))
			{
				diff.prev <- mean(abs(run.out$qsim-prev))
			}
			else
			{
				diff.prev<-NA
			}
			#"ntt"=ntt, "dt"=60*proj$dt/ntt, 
			#n.run.info <- c("ncuts"=nrow(run.out$groups)-1, "run.len"=run.len, 
			#								"sigma.q"=sigma.q, "ovf"=run.out$ovf, "deltaqbar"=diff.prev)
			n.run.info <- c("ncuts"=nrow(run.out$groups)-1, 
											"wb_pc"=100*run.out$wb/sigma.q, 
			                "sigma.q"=sigma.q, 
											"ovf"=run.out$ovf)
			
      
#			run.info <- rbind(run.info, n.run.info)
			cat(paste0(n.run.info, collapse="\t"), "\n", file=conn, append=T)
	#		write.table(run.info, sep="\t")
#	write.table(run.info, sep="\t", quote=F, row.names=F, file=conn)
  #		qsim.ntt <-cbind(qsim.ntt, )
#			prev <- run.out$qsim
			return(list(qsim=run.out$qsim[,1],run.info=n.run.info))
		}
	)


	qsim <- lapply(runs, function(run){run$qsim})
	qsim <- do.call(cbind, qsim)
	
	run.info <- data.frame(t(sapply(runs, function(run){run$run.info})))
	run.info <- run.info[order(run.info$ncuts),]					 
	qsim <- qsim[,order(run.info$ncuts)]
	colnames(qsim) <- run.info$ncuts
# add in the difference between each simulation and the final one
#	n.sim <- nrow(run.info)
#	run.info<-cbind(run.info, get.qbqn(qsim.ntt))
browser()
  #"qbqn"=sapply(1:n.sim, function(i){sum(abs(qsim.ntt[,i]-qsim.ntt[,n.sim]))/nrow(qsim.ntt)})
	nresp <- ncol(qsim)
	run.info$qbqn1_pc <- c(NA, colMeans(abs(qsim[,1:(nresp-1)]-qsim[,2:nresp])/qsim[,2:(nresp)]*100))

	run.info$qbfn_pc <- get.qbqn(qsim, pc=T)
	run.info<-signif(run.info,2)
	
	write.table(run.info, sep="\t", quote=F, row.names=F, file=conn)
	if(!is.null(fn.out))
	{	
		
		qsim.fn <- paste0(fn.out, ".tsv")
		write.zoo(qsim, sep="\t", file=qsim.fn)
	}
#	jpeg(dirname(fn), width=2*1024, height=2*768)
#	DisplayResponseResults(proj, qsim.ntt)
#	dev.off()
																#		 sel=c(start(qsim), end(qsim)))	
	
#	names(qsim.ntt)<-ntts
	return(list("proj"=proj, "run.info"=run.info, "qsim"=qsim))
}



# some plots of interest
summarise.response <- function(ri)
{
	layout(matrix(c(1,2),ncol=2))
	par(family="serif")
	par(mar=c(4, 3, 4, 4))
	# how does run time increase with no. time steps
	plot(ri[, c("ntt", "run.len")], type="l", xlab="Number of time steps", ylab="", xaxt="n")
	abline(h=pretty(range(tr[,"run.len"])),lty=2, col="gray")
	axis(side=1, at=ri[,"ntt"], labels=ri[,"ntt"],)
	par(new=T)
	plot(ri[, c("ntt", "eff")], ann=F, xaxt="n", yaxt="n", type="l", lty=2)
	axis(side=4)
	mtext(side=4, text="Efficiency", line=3)
	mtext(side=2, text="Run time (s)", line=2)	
	# grid	
	abline(v=ri[,"ntt"],lty=2, col="gray")
	
	
	plot(ri[,c("dt","ovf")], type="l", ylab="", xlab="Time step (m)", xaxt="n", yaxt="n")
	mtext(side=4, text="Overland flow total (m)", line=2)
	axis(side=4)
	axis(side=1, at=ri[,"dt"], lab=ri[,"dt"])
	abline(v=ri[,"dt"],lty=2, col="gray")
	abline(h=pretty(range(ri[,"ovf"])),lty=2, col="gray")

	
# 	# flow and	
# 	#tr[,"sigma.q"]<-tr[,"sigma.q"]/1000
# 	ovf.pc  <- ri[,"ovf"]
# 	plot(tr[, c("ntt", "sigma.q")], type="l", xlab="Time interval (hr)", ylab="", xaxt="n")
# 	abline(h=pretty(range(tr[,"sigma.q"])),lty=2, col="gray")
# 	#axis(side=1, at=tr[,"ntt"], labels=tr[,"ntt"],)
# 	par(new=T)
# 	plot(ri[, c("ntt....dt", "ovf")], ann=F, xaxt="n", yaxt="n", type="l", lty=2)
# 	axis(side=4)
# 	mtext(side=4, text="Overland flow (% of total)", line=2)
# 	mtext(side=2, text="Total discharge (m)", line=2)	
# 	# grid
# 	
# 	abline(v=ri[,"ntt"],lty=2, col="gray")	


	# how does overland flow increase with size of time step
	
	# how does total flow increase with number of steps
	
	# how does efficiency increase with no. time steps
#	plot(ri[, c("ntt...dt", "run.time")], type="l", xlab="Number of time steps", ylab="Run time (s)")
	
}

DisplayResponseResults <- function(proj, qsim, 
																	 sel=c(start(qsim), end(qsim)))
{
	par(family="serif")
  evap <- proj$pe
  evap[]<-NA
 # names(resp$qsim)<-paste("dt", resp$run.info[,"dt"])
  DisplayDischargeSelection(proj$groups, qobs=qsim, qsim=proj$qobs, 
                                        proj$rain, evap, sel=sel, 
                                        qmax=max(qsim[])*1000,   # upper bound for discharge display, in mm/hr
                                        ichan=proj$ichan, proj$disp.par, title="")
 

}

DisplayResponseResults2 <- function(resp, 
                                    s=NULL, e=NULL, legend=F)
{
	par(mar=c(1,3,3,3))
  proj<- resp$proj
  evap <- proj$pe
  evap[]<-NA
  qsim <- resp$qsim

  sel <- c(start(qsim), end(qsim))
  if(!is.null(s))
  {
    sel[1]<-as.POSIXct(s)  
  }

	if(!is.null(e))
	{
	  sel[2]<-as.POSIXct(e)	  
	}

	browser()  
#  sel <- as.vector(sel)
  par(family="serif")
  # names(resp$qsim)<-paste("dt", resp$run.info[,"dt"])
  DisplayDischargeSelection(proj$groups, qobs=qsim, qsim=proj$qobs, 
                            proj$rain, evap, sel=sel, 
                            qmax=max(qsim[])*1000,   # upper bound for discharge display, in mm/hr
                            ichan=proj$ichan, proj$disp.par, title="") 
                         #   qsim.lwd=2, qsim.col="green")
  
	# sub <- timeStr #""
	abline(h=pretty(c(0,proj$qmax*1000)),lty=2, col="gray")
	if(legend)
	{
	 xy <- click()
	# xy[1] <- 0

	legend(x=xy, xpd=T, title="Number of groupings", legend=resp$run.info[,"ncuts"], bg="white",
	       fill=get.qobs.cols(NULL), ncol=2)#nrow(resp$run.info))
	
#	 legend(x=xy, xpd=T, title="Time step (mins)", legend=resp$run.info[,"dt"], bg="white",
	  			 #fill=get.qobs.cols(NULL), ncol=2)#nrow(resp$run.info))

	}
}


GetDiffs <- function(range, results)
{
  nres <- length(range)
  ndiff <- nres-1
  diffs <- results[,2:nres]-results[,1:(nres-1)]

  # normalise against  discharges
  diffs.norm <- NULL
  

  # divide by the actual discharge 
  for(icol in 1:ndiff)
  {
    diffs.norm<- cbind(diffs.norm, diffs[,icol]/results[,icol])
  }
  browser()
  diff.from.max<- NULL

  for(icol in 1:nres)
  {
     diff.from.max <- cbind(diff.from.max, abs(results[,icol]-results[,nres]))
  }

  res.tab <- cbind(range[-nres], range[-1], colMeans(abs(diffs)), diff.from.max) 
  # table of results: mean (abs) difference between adjacnet sets and 
  browser()
  print(res.tab)
  
  return(diffs)
}

PlotResponseResults <- function(res, main="", 
                                leg=1:ncol(res), 
                                lty=rep(1, ncol(res)), 
                                lwd=rep(1, ncol(res)))
{
  n <- ncol(res)
  plot.zoo(1000*res, plot.type="single", col=rainbow(n), xlab="", ylab="Specific discharge (mm/hr)",main=main, cex.main=0.75,
           lty = lty, lwd=lwd)
  grid()
  legend(x="topright", fill=rainbow(n), title="Time step (mins)", legend=leg, bg="white", ncol=ceiling(n/5), cex=0.75)
}

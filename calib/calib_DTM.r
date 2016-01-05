# akima: Interpolation of irregularly spaced data
require(akima) 

load.source("dynatopmod.r")

get.calib.dir <- function(proj.ex)
{
	s <- format(proj.ex$run.par$start, "s=%Y-%m-%d")
	e <- format(proj.ex$run.par$end, "e=%Y-%m-%d")
	return(paste0(s, e, collapse=","))	
}


get.calib.fn <- function(proj.ex, par.min.max)
{
  ngroup <- nrow(proj.ex$groups)
  par.nm <- paste0(names(par.min.max), "=", par.min.max, collapse=",")
  tm <- format(Sys.time(), "%Y-%m-%d %H.%M.%S")
  s <- format(proj.ex$run.par$start, "s=%Y-%m-%d")
  e <- format(proj.ex$run.par$end, "e=%Y-%m-%d")
  
  return(paste0(s, ",", e, "/", par.nm, ",", tm, ".dat"))
  
}

build.calib.set <- function(groups, params)
{
	par.sets <- expand.grid(params)
	nms <- names(params) #'c("vchan", "vof", "td", "srz0", "m", "ln_t0", "srz_max")
	# create a table bifg enough to hold the permuations of parameters
	gr <- as.numeric(as.matrix(groups[2,nms]))
	 
	# params must be named
	n.set <- nrow(par.sets)
	
	calib.set <- data.frame(matrix(rep(gr, n.set), nrow=n.set, byrow=T))
	colnames(calib.set)<- nms
	calib.set[,colnames(par.sets)] <- par.sets
	
	return(calib.set)
}

          

plot.resp <- function(res, eff="NSE", var="m", ...)
{
  par(mgp=c(2,1,0))
  par(family="serif")
  par("mar"=c(4,3.5,3.5,3.5))

  rsr.col<- "brown"
  
  plot(res[,c(var, eff)], col="black", pch=20, cex=0.5, yaxt="ny")  #, ylim=c(0.85,0.89))
    
  axis(2, at=pretty(range(round(res[,eff],2), na.rm=T)))
  grid()
  # determine maximum bounding line (pareto front)
  u.bound <- sapply(sort(unique(res[,var])),
    function(x)
    {
      # values at this x
      vals <- res[which(res[,var]==x),]    
      max.eff <- max(vals[,eff], na.rm=T)  
      return(c(x, max.eff))
       
    }
    )
  u.bound <- t(u.bound)
  lines(u.bound, lwd=1, lty=1, col="slategray")
  max.eff <- u.bound[which.max(u.bound[,2]),2] 
  max.eff.val <- u.bound[which.max(u.bound[,2]),1] 
  abline(v=max.eff.val, col="red")
  abline(h=max.eff, col="red")
  mtext(side=4, at=max.eff, signif(max.eff,2), cex=0.75, line=1)
  
  par(las=2)
  mtext(side=3, at=max.eff.val, signif(max.eff.val,3), cex=0.75, line=1)
  
 # par(new=T)
#  plot(res[,c(var, "RSR")], col=rsr.col, pch=".", cex=3, axes=F, ann=F)
#  axis(side=4)
  
 # mtext("RSR", side=4, line=2, col=rsr.col)

}


# par.min.max.2 <- list("m"=c(0, 0.001), "ln_t0"=c(5.9,6.5))
# #resp <-calib.DTM(proj, par.min, par.max, n.run=20, fn=fn)

calib.names <- function(){return(list("ME", "MAE", "MSE", "RMSE", 
							 "NRMSE %", "PBIAS %", "RSR", "rSD", "NSE", "mNSE", 
							 "rNSE", "d", "md", "rd", "cp", "r", "R2", "bR2", "KGE", "VE", "logNSE", 
							 "wb", "qtot", "ovf", "ae", "ntt", "dt", "vchan", 
							 "vof", "td", "srz0", "m", "ln_t0", "srz_max"))} # "date", "tm", 


# apply the given result set to the project groups
apply.set <- function(disc, res, which.n, apply.to=1:nrow(disc$groups))
{
  nms <- intersect(names(disc$groups), colnames(res))
  disc$groups[apply.to,nms]<-res[which.n, nms]  
  return(disc)
  
}

par.names <- function(){return(c("vchan", "td", "srz0", "m", "ln_t0", "srz_max", "dz_drain", "ks_drain"))}
nms <- function(){c(names(run.gof(1:2,2:3)),
         "wb","ttp", "qtot","ovf","ae","ntt","dt", par.names())}


# Run all calibration sets for given project. load any existing results in the 
run.calib <- function(proj, disc=proj$disc[[1]], calib.dir=disc$dir, calib.set, 
										 apply.to = 1:nrow(disc$groups))
{	
	if(!file.exists(calib.dir)){
		dir.create(calib.dir, recursive=T)
	}
	fn <- tempfile("calib_", calib.dir, ".dat")

	calib.dat<- load.calib.dat(calib.dir)	
	
	all <- rbind(calib.set, calib.dat[,colnames(calib.set)])
	# build complete set and remove duplicates
	which <- which(!duplicated(all, fromLast=T)[1:nrow(calib.set)])
	# find sets unique to those already run
	# unrun <- setdiff(calib.set, res[,colnames(calib.set)])
	cat("Found ", length(which), " new calibration sets", "\n")
	
	
}



# apply parameters from given results, run again, and return 
run.sets <- function(proj, disc=proj$disc[[1]], calib.set, 
                     apply.to = 1:nrow(disc$groups), 
                     ichan=disc$ichan,
										 fn=NULL)
{
  res<- NULL
	try(messages.off())
  conn<-stdout()

  dn <- getwd()
  if(!is.null(fn))
  {
  	dn <- dirname(fn)
  	if(!file.exists(dn)){
  		dir.create(dn, recursive=T)
  	}

  	res <- NULL
  	if(file.exists(fn))
  	{
  		try(res <- read.table(fn, header=T, sep="\t"))
      nms <- names(res)
  	}
    if(is.null(res))
  	{
      nms <- colnames(disc$groups)[-(1:9)] # colnames(calib.set)
  		# build an initial set
  		res<-data.frame(t(rep(NA, length(nms))))
  		res<-res[-1,]  		
    }
    # keep only non-NA columsn
    res <- res[,apply(res, MARGIN=2, FUN=function(col){!any(is.na(col))})]
#     if(ncol(res)<length(nms()))
#     {
#       pos <- which(nms()=="vchan")
#       res <- cbind(res[,1:(pos-1)], calib.set$vchan[1], res[,pos:ncol(res)])
#     }

    # all names. need to pick these up from calib results. this ensure data is in correct format
    colnames(res) <- nms	
  	# write back to disk and ensure headings included
#  	write.table(res, sep="\t", file=fn, row.names=F, quote=F)
#  	conn <- file(fn, open="a")
#  	on.exit(close(conn))
  }  
  # if previous results haven't these names then ignore

  all.nms <- intersect(names(calib.set), nms)

#  inm<- which(nms %in% colnames(calib.set))

  all <- rbind(calib.set[,all.nms], res[,all.nms])
  # build complete set and remove duplicates
  which <- which(!duplicated(all, fromLast=T)[1:nrow(calib.set)])
  # find sets unique to those already run
 # unrun <- setdiff(calib.set, res[,colnames(calib.set)])

  cat("Found ", length(which), " new calibration set(s)", "\n")
 
  if(interactive())  
  {
 #   pb <- winProgressBar( "Running calibration...", min=1, max=length(which))
    
  #  on.exit(close(pb), add=T)
  }
  i.calib <- 1
  for(i in which)
  {
    lab <- paste0(colnames(calib.set), "=", calib.set[i,], collapse=",")
    if(interactive())  
    {    
#      setWinProgressBar(pb, value=i.calib, label=lab, title = paste("calibration ", i.calib, " of ", length(which)))
    }

    disc <- apply.set(disc, calib.set, i, apply.to=apply.to)
    
    cm <- disc$groups[2,nms]  #-(1:9)]#   signif(as.numeric(disc$groups[apply.to[1],par.names()]),4)
  #  print(cm)   
    proj$disp.par$title.main <- lab

    # run it!
    run <- run.proj(proj, disc=disc)

    # calculate efficiencies (note order!!!)
    effs <- run.gof(run$qsim, run$qobs, digits=3)
    ttp <- time.at.peak(run$qsim)
    run.info <- c(effs, "wb"=run$wb, "ttp"=ttp, "qmax"=max(run$qsim), "qtot"=sum(run$qsim)*proj$dt, 
                  "ovf"=run$ovf, "ea"=sum(run$evap[,"ae"])*proj$dt,
                  "ntt"=proj$ntt, "dt"=proj$dt, cm)
  #  run.info <- signif(unlist(run.info),4)
  #  conn <- file(fn, open="a")
 #   on.exit(close(conn))
 #  	cat(paste(run.info, "\t"), "\n", file=conn)
#    close(conn)
  	res <- rbind(res, run.info)
 #   colnames(res)<-nms()
    cat("Writing results to file ", fn, "\n")
    write.table(res, sep="\t", file=fn, row.names=F, quote=F)

    i.calib <- i.calib +1 
  # 	flush(conn)	   
    # save 
#     if(nrow(res) > 50 & i.calib %% save.int == 0)
#     {
#     	jpeg(save.fn, width=1024, height=1024)
#     	results.summary(res)
#     	dev.off()
#     }
  }
#  close(pb)
#  close(conn)
  return(data.frame(res))
} 
get.calib.fn <- function(calib.dir, calib.set)
{
  suff <- paste0(names(calib.set), "=", calib.set, collapse="_")

# output file for results, always a single one per run
  fn <- file.path(calib.dir, paste0("calib_", suff, ".dat"))

  return(fn)
}

# vary the parameters for just one HRU in turn
sensitivity.HRU <- function(calib.set, proj, disc=proj$disc[[1]], benchmark=NULL, dn=NULL)
{
# 	if(is.null(benchmark))
# 	{
# 		# compute benchmark efficiency
# 		run <- run.proj(proj, disc=disc)
# 		benchmark <- gof.run(run)$NSE
# 	}
	resp <- list()
	ihrus <- setdiff(1:nrow(disc$groups), disc$ichan)
	resp <- lapply(ihrus,
	
		function(ihru)
		{
			fn <- NULL
			if(!is.null(dn))
			{
				fn <- file.path(dn, paste("calibHRU_", ihru, ".dat"))
			}
			run.sets(proj, disc, calib.set, apply.to=ihru, fn=fn)	
		}
	)
	names(resp) <- disc$groups[ihrus,"id"]
	return(resp)
}

  
# calibration run consists of individual files 
collect.calib.results <- function(calib.dir, fn=NULL, rebuild=F)
{
  fns <- dir(calib.dir, "*.dat", full.names=T)
  calib.dat <- NULL     
  
  if(!rebuild && !is.null(fn) && file.exists(fn))
  {
    message(paste0("Loading calibration results from ", fn))
    calib.dat <- read.table(fn, sep="\t", header=T)
  #  return(calib.dat)
    # exclude results already generated
    nms <- apply(calib.dat, MARGIN=1, 
                 function(x)get.calib.fn(calib.dir, x[-1]))
  
    fns <- setdiff(fns, nms)
  }
  

  cat("Loading ", length(fns), " results\n")
  if(length(fns)>0)
  {
    if(interactive())
    {
    pb <- txtProgressBar(max=length(fns), style=3)
    on.exit(close(pb))
    }
 
    i <- 1
    cat("\nCollecting ", length(fns), " results\n")
    for(fni in fns)
    {
       dat <- read.table(fni, sep="\t", header=T)
      
#        names(dat) <- c(names(dat[1:3]), "ttp", names(dat[4:16]))
     # try({
#            #dat <- dat[,c("NSE", names(calib.set))]
#           if(ncol(dat)>16)
#           {
             dat$ttp <- as.POSIXct(dat$ttp, origin="1970-01-01")
#             names(dat)[4]<- "ttp"
#             names(dat)[5]<- "wb"
#             write.table(dat, file=fni, sep="\t", quote=F)
#             
             calib.dat <- rbind(calib.dat, dat)
#           }
          # }, silent=T)
      i <- i + 1
      if(interactive())
      {
        setTxtProgressBar(pb, i)
      }
    }
    
    
    
    igood <- which(calib.dat$NSE > 0)
    
    calib.dat <- calib.dat[igood,]
    if(!is.null(fn))
    {
      message(paste0("Saving collated results to ", fn))
      # file.path(data.dir, "calib.dat")
      write.table(calib.dat, file=fn, sep="\t", quote=F, row.names=F)
    }
    return(calib.dat)
  }
  return(NULL)
}





calib.nms <- function()
{
	return(c("NSE",	"R2", "wb",	"qtot",	"ovf",	"ae",	"ntt",	"dt", "vchan",	"vof",	"td",	"srz0",	"m",	"ln_t0", "srz_max"))
}
require(rgl)
calib.show <- function(fn,...)
{
  calib.dat <- read.dat(fn)
  
  good<- which(apply(calib.dat, MARGIN=1, function(x)all(!is.na(x))))
  calib.dat<- calib.dat[good,]
#	calib.dat$logNSE <- log(calib)
 # calib.dat <- cbind(calib.dat, "NSE2"=(calib.dat$NSE + calib.dat$logNSE)/2)
#	colnames(calib.dat) <- c("eff", "r2", "pbias", "rsr", "wb","qtot", "ovf", "ntt", "dt", "vchan", "vof", "td", "srz0", "m", "ln_t0", "srz_max")
	try(results.summary(calib.dat,...))
  return(calib.dat)
}

read.dat <- function(fn, nms = calib.nms())
{
	cat("reading file ", fn, "\n")
	calib.dat <- NULL
	try(calib.dat<- read.table(fn, header=F, sep="\t", skip=1))
		
	if(is.null(calib.dat))
	{
		try(calib.dat<- read.csv(fn, sep="\t", fill=T))
		
	}
	if(is.data.frame(calib.dat))
	{
		# remove na colums
		good.col <- which(apply(calib.dat, MARGIN=2, function(x)all(!is.na(x))))
		calib.dat <- calib.dat[,good.col]
		good.row <- which(apply(calib.dat, MARGIN=1, function(x)all(!is.na(x))))
		calib.dat <- calib.dat[good.row,]
		
		if(ncol(calib.dat)==length(nms))
		{
				colnames(calib.dat)<-nms
		}
			return(calib.dat)
	}
}

load.source("filled.contours.r")

# load everything in the given directory, discarding duff files and dupplicates
load.calib.dir <- function(dn, ...)
{
	dat <- NULL
	fns <- dir(dn, "*.dat", full.names=T)
	dats <- lapply(fns, read.dat)

	
	dats <- delete.NULLs(dats)
	dat <- do.call(rbind, dats)
	dat <- unique(dat)
	if(length(dat)>0)
	{
		message(paste(nrow(dat), " lines"))
		print(dat[which.max(dat$NSE),])
	}
			
	return(dat)

}

# column names for result summaries
names <- list(x="m", 
							y="ln_t0", 
							eff="NSE", 
							colour="NSE")




# <-c("eff", "wb", "qtot", "ovf", "td", "x", "m", "ln_t0", "srz_max")
results.summary <- function(calib.dat,nx=40,ny=nx, thresh=NA, scale=F, 
														nms = names,    # column name to use for colouring facets
                            nlev=15,
														...)
{
  
	if(!is.na(thresh))
	{
  	calib.dat <- calib.dat[which(calib.dat[,names$eff]>=thresh),]
  	calib.dat <- calib.dat[which(calib.dat[,names$colour]>=thresh),]
	}
  open.devices(2)
	#layout(matrix(c(1,2,3,4),ncol=2,byrow=T))
	#plot(100*calib.dat[,"ovf"]/calib.dat[,"qtot"], calib.dat[,"eff"])

	# carpet plot of x=m, y=ln_t0 with response var z=eff
	x <- calib.dat[,nms$x]
	y <- calib.dat[,nms$y]
	z <- calib.dat[,nms$eff]
	z2 <- calib.dat[,nms$colour]
	x0 <- seq(min(x), max(x), length = nx)
	y0 <- seq(min(y), max(y), length = ny)
  
	# regular gridded mesh  using cubic spline interpolation
	calib.mesh <- interp(x, y, z, xo=x0,yo=y0, duplicate="strip")
  calib.mesh2 <- interp(x, y, z2, xo=x0,yo=y0, duplicate="strip")
  
	x0 <- calib.mesh$x
	y0 <- calib.mesh$y
	z0 <- calib.mesh$z
  z02  <- calib.mesh2$z
  if(scale)
  {
    x0 <- as.vector(scale(x0))
    y0 <- as.vector(scale(y0))
  }

	#not.na <- which(!is.na(z02))
#	z02[which(is.na(z02))]<-0
#	activate.device(1)
	#image(x=x0, y=y0, z=z02, col=tim.colors(20))
	#  z0[which(z0<thresh)] <- NA
#	filled.contour(x0, y0, z02, col=tim.colors(30))

	activate.device(1)
	filled.contour(x0, y0, z0)
activate.device(2)
plot.new()	
  par(cex.axis=0.9, cex.axis=0.75, cex.lab=1, font.lab=2)
#  layout(matrix(c(1,2,3,3),ncol=2))
  layout(matrix(c(1,2),ncol=2))
#browser()
	calib.interp <- data.frame(expand.grid(x0, y0), as.vector(z0))
  colnames(calib.interp) <- nms[c("x", "y", "eff")]
 # calib.interp2 <- calib.interp[which(calib.interp[,nms$"y"]>14 & calib.interp[,nms$"y"]<=15.2),]
  #& calib.dat$m>0.0063),]  #,]

  plot.resp(calib.interp, eff=names$eff, names$y)
  #plot.resp(calib.dat, eff=names$eff, names$y)

 #calib.interp<- calib.interp[which(calib.interp$"NSE">0.85),]
	calib.interp <- data.frame(expand.grid(x0, y0), as.vector(z02))
	colnames(calib.interp) <- nms[c("x", "y", "colour")]

	plot.new()

	plot.resp(calib.interp, eff=names$colour, names$x)
#plot.resp(calib.dat, eff=names$eff, names$x)
  return(NULL)
  par(new=T)
  par(fig=c(0.5,1,0,1))
  cols<-colorRampPalette(c("red", "green"))

#filled.contour2(x, y, z, xlab=names$x, 
#                ylab=names$y,col=tim.colors(nlev+6), nlevels=nlev)
 # z02 <- round(z02, 4)
#  z0 <- round(z0, 4)

#  cols <- tim.colors(nlev)
   levels <- pretty(range(z0, na.rm=T), nlev)
   filled.contour2(x0, y0, z0, xlab=names$x, 
                   ylab=names$y, col=tim.colors(length(levels)+1), levels=levels)#,
                   
                  # nlevels=nlev)
  cols <- tim.colors(length(levels)+1)
#  pt.cols <- col[cut(z02, breaks=nlev, lab=F)]
#  plot(calib.interp[,c(names$x, names$y)], col=pt.cols)
	par(fig=c(0.5,1,0,1))
	par(new=T)
	levels <- pretty(range(z0, na.rm=T), nlev)	
  contour(x0,y0,z0, levels=levels, add=T)
  #        nlevels=nlev,add=T)
  
}

vis.nms <- function(...)
{
  nms <- data.frame(x="m",
                   y="ln_t0",
                   eff="NSE",col="NSE")
  nms <- merge.lists(nms, list(...))
 return(nms) 
}
# 
results.summary.2 <- function(calib.dat,nx=40,ny=nx, thresh=NA, scale=F, 
                            nms = data.frame(x="m",
                                       y="ln_t0",
                                       eff="NSE",col="NSE"),                             
                            nlevs=10,
                            zexp=0.5,
                            col=tim.colors(64),
                       #     col.set <- rev(topo.colors(3,alpha=0.75))
                            ...)
{
  nms <- merge.lists(nms, list(...))
  if(!is.na(thresh))
  {
    calib.dat <- calib.dat[which(calib.dat[,nms$eff]>=thresh),]
  }
  #	dev.set(3)
  #layout(matrix(c(1,2,3,4),ncol=2,byrow=T))
  #plot(100*calib.dat[,"ovf"]/calib.dat[,"qtot"], calib.dat[,"eff"])
  
  # carpet plot of x=m, y=ln_t0 with response var z=eff
  x <- calib.dat[,nms$x]
  y <- calib.dat[,nms$y]
  z <- calib.dat[,nms$eff]
  # values for colour
  z2 <- calib.dat[,nms$col]
  x0 <- seq(min(x), max(x), length = nx)
  y0 <- seq(min(y), max(y), length = ny)
  
  # regular gridded mesh  using cubic spline interpolation
  calib.mesh <- interp(x, y, z, xo=x0,yo=y0, duplicate="strip")
  calib.mesh2 <- interp(x, y, z2, xo=x0,yo=y0, duplicate="strip")
  
  x0 <- calib.mesh$x
  y0 <- calib.mesh$y
  z0 <- calib.mesh$z
  z02  <- calib.mesh2$z
  if(scale)
  {
    x0 <- as.vector(scale(x0))
    y0 <- as.vector(scale(y0))
  }
 	cols <- matrix(fields::drape.color(z02, col=col)$color.index, nrow=(length(x0)-1), ncol=(length(y0)-1))
 	cols <- cbind(cols[,1], cols)
 	cols <- rbind(cols[1,], cols)
 
 
#  cols <- tim.colors(nlevs)[cut(z, breaks=nlevs, labels=F)]
  
 # cols[]<- col.set[1]  #"gray"

  par(family="serif")

#  cols[which(accept.fit(calib.dat, gof.good))]<-col.set[2]
#  cols[which(accept.fit(calib.dat, gof.very.good))]<-col.set[3]

 # cols[which.max(res$NSE)]<-"red"
 # par(new=T)
  par(mar=c(1,1,1,1))
  pmat <- persp(x0,y0,z0,phi=60,box=T, axes=T, border=NA, theta=45, col=NA, exp=0.8, 
                xlab="", ylab="", zlab="", 
                ticktype = "detailed")
  xyt<-trans3d(x,y,z, pmat)
 
 #	cols <- terrain.colors(nlevs)[cut(z2, breaks=nlevs, labels=F)]
  points(xyt$x,xyt$y,col=cols,pch=16, cex=1)
  
#  legend(fill=rev(col.set), title="Fit", legend=c("very good", "good", "acceptable"), x="bottomleft")

  # scatterplot3d(x,y,z, color=cols, xlab=nms$x, ylab=nms$y, zlab=nms$eff, box=F, pch=16)
	#points3d(x,y,z, color=cols, size=5, alpha=0.75)
#  rgl.surface(x0,y0,z0,color=cols, box=F)
par3d(family="serif")
 	persp3d(x0,y0,z0,xlab=nms$x, ylab=nms$y, zlab=nms$eff,
           color=cols, 
           box=F)
  axes3d(labels=F)
# mtext3d(nms$x, edge="x-", line=1)
# mtext3d(nms$y, edge="y-", line=1)
# mtext3d(nms$eff, edge="z+", line=1)
  aspect3d(1,1,zexp)
  return(calib.dat)
}

# create a raster from the calibration results
results.raster <- function(calib.dat, par.min.max=NULL, nx=40, ny=nx, 
                           x.nm = "m", y.nm = "ln_t0",
                           z.nm="NSE")
{
  x <- calib.dat[,x.nm]
  y <- calib.dat[,y.nm]
  
  if(is.null(par.min.max))
  {
    xmn <- min(x); xmx <- max(x)
    ymn <- min(y); ymx <- max(y)
  }
  else
  {
    xmn <- par.min.max$m[1]
    xmx <- par.min.max$m[2]
    ymn <- par.min.max$ln_t0[1]
    ymx <- par.min.max$ln_t0[2]    
  }
  
  x0 <- seq(min(x), max(x), length = nx)
  y0 <- seq(min(y), max(y), length = ny)
  
  # values for layer layer
  z <- calib.dat[,z.nm]
  
  # regular gridded mesh  using cubic spline interpolation
  calib.mesh <- interp(x, y, z, xo=x0,yo=y0, duplicate="strip")
  x0<- calib.mesh$x # should be same
  y0<- calib.mesh$y
  z0<- calib.mesh$z  
  
  calib.rast <- 
    raster(xmn=xmn, xmx=xmx, 
           ymn=ymn, ymx=ymx, 
           ncol=nx, nrow=ny)
  
  res <- calib.rast

    # values for layer layer
    z <- calib.dat[,z.nm]    
    # regular gridded mesh  using cubic spline interpolation
    calib.mesh <- interp(x, y, z, xo=x0,yo=y0, duplicate="strip")
    x0<- calib.mesh$x # should be same
    y0<- calib.mesh$y
    z0<- calib.mesh$z  
    res <- addLayer(res, raster::setValues(calib.rast, z0))    
  

  return(res)
}

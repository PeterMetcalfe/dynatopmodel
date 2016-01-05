# if any reservoirs specified determine their output using current storage and 
# lookup table calculated from topographic analysis 
update.reservoir.output.flow <- function(groups, flows, stores, reservoirs)
{
	if(length(reservoirs)>0)
	{
		#groups.res <- lapply(names(reservoirs)
		q.res <- sapply(names(reservoirs),
				 function(r.id)
				 {
					 	r <- reservoirs[[r.id]]
					 	
					 	gr.id <- groups[which(groups$id==r.id),"id"]
					 	# look up current storage (deficit) from id
					 	r.store <- stores[which(stores$id==r.id),"sd"]
					 	# use approximation determined in build.reservoirs to estimate 
					 	q <- r(r.store)
										
					 	if(is.na(q))
					 	{
					 		# flood or negative flows
					 		browser()					 		
					 	}		
				 	return(q)
				 }
		)

		flows[which(flows$id %in% names(reservoirs)),"qbf"] <- q.res		
	}
	return(flows)
}

# bernoulli equation for pipeflow at base of reservoir outlet
celerity.reservoir <- function(qin, qbf, params, ichan)
{
  # discharge area
  A <- params["A"]
  c <- qin-qin
  try(c <- -A^2/(qin+qout))
  return(c)
}

# kinematic routing for reservoir
route.kinematic.reservoir <- function(groups,
flowst1,      # fluxes at previous time step (all known)
flows,        # fluxes at current time step (prediction time)
stores,       # current storage
dtt, 
ichan, 
time, 
w=0.5,                      #  weighting factor in numerical scheme
# a function that provdes wave celerities for each groups given fluxes and other parameters
nIter=30) 
{
	lapply(names(reservoirs),
				 function(r.id)
				 {
				 	r <- reservoirs[[r.id]]
				 	
				 	gr.id <- groups[which(groups$id==r.id),"id"]
				 	flowst1 <- flowst1[which(flowst1$id==r.id),]
				 	flows <- flows[which(flows$id==r.id),]
				 	stores <- flowst1[which(stores$id==r.id),]
				 }
	)

	
	
}

res.celerity.func <- function(groups, qin, qout, ...)
{
	# the reservoir info is now maintained in groups
	
	# estimate for local gradient
	q0 <- r(s*0.99)
	q1 <- r(s*1.01)
	ds <- s*0.02
	
	dsdq <- ds/(q1-q0)
	return(dsdq)
	
}

# construct a cubic spline function to approximate reservoir discharges given storage
# use lookup table calculated from topographic analysis
build.reservoirs <- function(reservoirs, groups)
	{
		if(length(reservoirs)>0)
		{ 
			nms <- names(reservoirs)
			#groups.res <- lapply(
			reservoirs <- lapply(names(reservoirs),
											function(r.id)
											{
												r <- reservoirs[[r.id]]
												# find the current storage
												gr.id <- groups[which(groups$id==r.id),"id"]																							
												a<-groups[gr.id,]$area
												# to allow interpolation ensure first storage - discharge entry is for nil flows 
												if(r[1,]>0){r <- rbind(rep(0,ncol(r)), r)}
												# lookup discharge from table of total storage vs discharge
												max.store <- max(r$S)/a
												# convert total storage to specific defict			
												sd <- max.store - r$S/a
												groups[gr.id,]$sd.max <- max.store/a
												
												# function to calculate the specific discharge from the storage deficit
												fun <- approxfun(sd, r$Q/a)
											}
			)
			names(reservoirs) <- nms					
		}

		return(reservoirs)
	}
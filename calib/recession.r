

# analysis of recession parameters from hydrograph
nlen <- 8*4

qobs <-  as.vector(proj.ex$qobs)
sel <- rep(T, length(qobs)-1)
for(idiff in 1:nlen)
{
  
  diff <-diff(qobs, idiff)
  sel <- sel & diff<0
  
  
}

# series with non recession regions removed

qrn<- proj.ex$qobs
qrn[which(!sel)]<- NA

plot(qrn)


# allow user to select a region
s1 <-locator(n=2, "l")

s1 <- as.POSIXlt(s1$x, origin=t.origin, tz="GMT")
rec1<-qrn

rec1<- subset_zoo(rec1, start=s1[1], end=s1[2])

plot.zoo(rec1)

it0 <- which.max(rec1)
q0 <- as.numeric(rec1[it0])
t0 <- index(rec1)[it0]
q <- rec1[it0:length(rec1),]

# time in hours from start of recession period
t <- as.numeric(index(q)-t0)/3600

q <- as.numeric(q)
# solve for m at each point
#m <- t/(1/q-1/q0))

# slope of 1/q should be approx 1/m see Beven 2012 Chapter 6, p214
q1 <- 1/q
#m1 <- (q1[length(q1)]-q1[1])/max(t)

# estimate at each point

m1 <- (q1-1/q0)/t

m <- 1/m1

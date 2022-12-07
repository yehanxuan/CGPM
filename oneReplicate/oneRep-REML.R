
numElem = 1
obsNumLower = 5
obsNumUpper = 10
if (DataType == "easy"){
    tmin = 0
    tmax = 1
    eigenfList = get_eigenfunList_Paul("easy")
    obsCol = get_allObs_Paul(samplesize, M=5,r=3, "easy", scoreType, alpha = 0.6, nmin=2, nmax=10, a=1, b=1, sig=1/4)
    M.set = seq(4, 12, by = 1)
    r.set = 3
    sig.EM = 1
    CVmethod = "like"
} else if (DataType == "prac"){
    tmin = 0
    tmax = 1
    eigenfList = get_eigenfunList_Paul("prac")
    obsCol = get_allObs_Paul(samplesize, M=10, r=5, "prac", scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4)
    M.set = seq(8, 15, by = 1)
    r.set = 5
    sig.EM = 1
    CVmethod = "like"
} else if (DataType == "hybrid"){
    tmin = 0
    tmax = 1
    eigenfList = get_eigenfunList_Paul("hybrid")
    obsCol = get_allObs_Paul(samplesize, M=10, r=10, "hybrid", scoreType, alpha = 0.6, nmin=2, nmax=10, a=1, b=1, sig=1/4)
    M.set = seq(8, 15, by = 1)
    r.set = 5
    sig.EM = 1
}  else if (DataType == "pracNU") {
    tmin = 0
    tmax = 1
    shape1 = 2
    shape2 = 2
    eigenfList = get_eigenfunList_Paul_NU("pracNU", shape1, shape2)
    obsCol = get_allObs_Paul(samplesize, M=10, r=5, "pracNU", scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4, noiseType, shape1, shape2)
    M.set = seq(8, 15, by = 1)
    r.set = 5
    sig.EM = 1
    CVmethod = "like"
} else if (DataType == "easySin") {
    Mdata = 5
    rdata = 3
    pcaTrans = matrix(rnorm(Mdata*rdata), Mdata, rdata)
    pcaTrans = qr.Q(qr(pcaTrans))
    eigenfList = get_eigenfunList_PaulNoSpline("easySin", pcaTrans)
    obsCol = get_allObs_PaulNoSpline(samplesize, M = Mdata, r = rdata, DataType, pcaTrans,
                                     scoreType, alpha = 0.6, nmin = 2, nmax = 10, a=1, b=1, sig = 1/4, noiseType)
    M.set = seq(4, 12, by = 1)
    r.set = 3
    sig.EM = 1
    CVmethod = "like"
} else if (DataType == "pracSin") {
    Mdata = 10
    rdata = 5
    pcaTrans = matrix(rnorm(Mdata*rdata), Mdata, rdata)
    pcaTrans = qr.Q(qr(pcaTrans))
    eigenfList = get_eigenfunList_PaulNoSpline("pracSin", pcaTrans)
    obsCol = get_allObs_PaulNoSpline(samplesize, M=Mdata, r=rdata, DataType, pcaTrans,
                                     scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4, noiseType)
    M.set = seq(13, 23 ,by = 1)
    r.set = 5
    sig.EM = 1
    CVmethod = "like"
} 

## Data processing
newObsCol = obsCol[, -2]
newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
basis.method = "bs"
sl.v = rep(0.5, 10)
max.step = 50
grid.l = seq(0,1,0.01)
grids= seq(0,1,0.002)

nmax<-max(table(newObsCol[,1]))  
L1<-min(as.numeric(newObsCol[,3]))           
L2<-max(as.numeric(newObsCol[,3]))
data.list = fpca.format(newObsCol)
n<-length(data.list)    ##number of subjects
if(n==0){ 
    print("error: no subject has more than one measurements!")
    return(0)
}
data.list.new = data.list
for (i in 1:n){
    cur<-data.list[[i]][[1]]
    temp<-(cur[,2]-L1)/(L2-L1)
    temp[temp<0.00001]<-0.00001
    temp[temp>0.99999]<-0.99999 
    cur[,2]<-temp
    data.list.new[[i]][[1]]<-cur
}




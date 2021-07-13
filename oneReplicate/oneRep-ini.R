fExpNum = 4
pcaCompNum = fExpNum * 2
pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
pcaTrans = qr.Q(qr(pcaTrans))
splitP = c(-1e-9, 2+1e-9)

numSample = 2e2
numElem = 1
obsNumLower = 5
obsNumUpper = 10
pcaKappa = exp(-(0:(pcaCompNum-1))/2)*10
pcaKappaSqrt = sqrt(pcaKappa)
noiseSigma = 1
eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
source("./data/generate_data.R")
obsCol = get_allObs(eigenfList, numSample, pcaKappaSqrt, noiseSigma, 
                    obsNumLower, obsNumUpper)


## Data processing
newObsCol = obsCol[, -2]
newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]

# Change the data to matrix form
#newObsCol = as.matrix(sapply(newObsCol, as.numeric))
#data(easy)
#newObsCol = easy$data

#data(prac)
#newObsCol = prac$data

r.set = c(8)

ini.method = "EM.self"
basis.method = "bs"
sl.v = rep(0.5, 10)
max.step = 50
grid.l = seq(0,1,0.01)
grids= seq(0,1,0.002)

sig.EM = 1
nmax<-max(table(newObsCol[,1]))  
L1<-min(as.numeric(newObsCol[,3]))              ##range of the time 
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




#PACEobsCol = obsCol[order(obsCol$obsID, obsCol$obsT) , ]

#PACEobsCol = PACEobsCol[ ,-2]
#PACEobsCol[, c(2, 3)] = PACEobsCol[, c(3, 2)]
#colnames(PACEobsCol)[c(2,3)] = colnames(PACEobsCol)[c(3,2)]


#DatanewObsCol = list() 
#for (i in 1:length(unique(PACEobsCol[,1])) ){
#    idx = PACEobsCol[ ,1] == i
#    DatanewObsCol[["Lt"]][[i]] = PACEobsCol[idx, 3]
#    DatanewObsCol[["Ly"]][[i]] = PACEobsCol[idx, 2]
#} 



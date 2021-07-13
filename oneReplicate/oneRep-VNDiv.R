# fExpNum = 5
# pcaCompNum = fExpNum * 2
# pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
# pcaTrans = qr.Q(qr(pcaTrans))
# splitP = c(-1e-9, 2+1e-9)
# 
# numSample = 500
# numElem = 1
# obsNumLower = 5
# obsNumUpper = 10
# pcaKappa = exp(-(0:(pcaCompNum-1))/2)*10
# pcaKappaSqrt = sqrt(pcaKappa)
# noiseSigma = 1
# eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
# obsCol = get_allObs(eigenfList, numSample, pcaKappaSqrt, noiseSigma, 
#                     obsNumLower, obsNumUpper)

numElem = 1
obsNumLower = 5
obsNumUpper = 10
if (DataType == "cai"){
    #nKnots = 8
    nKnots = 13
    #M.set = seq(6, 12, by = 2)
    M.set = nKnots
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    K = splineObj$getDoF()
    pcaCompNum = 6
    eigenfList = get_eigenfunList_cai(pcaCompNum)
    obsCol = Generate_Cai(eigenfList, samplesize, pcaCompNum, noiseSigma = 1/4, alpha = 1)
    r.set = pcaCompNum
    sig.EM = 1/4
    #Init = EMInit(obsCol, splineObj,  M.set, r.set, (sig.EM)^2 )
    #XInit = Init[[1]]
    #sig2hat = Init[[2]]
    # InitType = "LS"
} else if (DataType == "Peng"){
    #nKnots = 8
    #nKnots = 13
    # scoreType = "Gaussian"
    nKnots = 18
    M.set = nKnots
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    pcaCompNum = 4
    pcaKappa = c(1, 0.66, 0.52, 0.07)
    pcaKappaSqrt = sqrt(pcaKappa)
    # The true eigenfunction of data can be represented by cubic spline with 20 knots
    eigenfList = get_eigenfunList_exact(order = 4, nknots = 10, pcaCompNum)  
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1/4, obsNumLower, obsNumUpper, scoreType)
    sig.EM = 1/4
    sig2hat = (1/4)^2
    r.set = pcaCompNum
    #Init = EMInit(obsCol, splineObj, M.set, r.set, (sig.EM)^2 )
    #XInit = Init[[1]]
    #sig2hat = Init[[2]]
    #InitType = "LS"
} else if (DataType == "PengHybrid"){
    #scoreType = "Gaussian"
    nKnots = 8
    #nKnots = 15
    M.set = seq(8, 15, by = 1)
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    pcaCompNum = 8 
    pcaKappa = c(1, 0.66, 0.52, 0.07, 9.47e-3, 1.28e-3, 1.74e-4, 2.35e-5 )
    pcaKappaSqrt = sqrt(pcaKappa)
    eigenfList = get_eigenfunList_exact(order = 4, 8, pcaCompNum)
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1/4, obsNumLower, obsNumUpper, scoreType)
    r.set = c(2:6)
    sig2hat = (1/4)^2
    #InitType = "EM"
} else if (DataType == "FourierDense"){
    fExpNum = 5
    pcaCompNum = fExpNum * 2
    a = seq(0,1, length.out = 7)  ### Useless, only create a sequence of length 7
       pcaTrans = Generate_rescale_U(seedJ, a)
    #   pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
    #    pcaTrans = qr.Q(qr(pcaTrans))
    splitP = c(-1e-9, 2+1e-9)
    # numSample = 2e2
    pcaKappa = exp(-(0:(pcaCompNum-1))/2)*10
    pcaKappaSqrt = sqrt(pcaKappa)
    eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    # scoreType = "Gaussian"
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 0.7,
                        obsNumLower, obsNumUpper, scoreType)
    nKnots = 20
    M.EM = nKnots + 2
    r.set = 8
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = (0.7)^2
    # InitType = "LS"
    #InitType = "LOC"
    #InitType = "EM"
} else if (DataType == "FourierOrth"){
    a = seq(0,1, length.out = 7)  ### Useless, only create a sequence of length 7
    pcaTrans = Generate_rescale_U(seedJ, a)
    fExpNum = 5
    pcaCompNum = fExpNum * 2
#    pcaTrans = Generate_pcaTrans(seedJ = 7, pcaCompNum)
    splitP = c(-1e-9, 2+1e-9)
    pcaKappa = c(1, 0.66, 0.52, 0.07, 9.47e-3, 1.28e-3, 1.74e-4, 2.35e-5, 3.18e-6, 4.30e-7 )
    pcaKappaSqrt = sqrt(pcaKappa)
    eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    ## Set score Type, uniform or t-distribution
    #scoreType = "uniform"
    #scoreType = "t"
    # scoreType = "Gaussian"
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1/4, obsNumLower, obsNumUpper, scoreType)
    ## Set knots and ranks
    nKnots = 18
    r.set = 6
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = (1/4)^2
    #InitType = "LS"
    #InitType = "LOC"
    # InitType = "EM"
} else if (DataType == "NonUnif"){
    tmin = 0
    tmax = 1
    knots = Generate_knots(seedJ, 7)
    #knots = sort(runif(9))
    pcaTrans = Generate_rescale_U(seedJ, knots)
    #pcaTrans = Generate_pcaTrans(seedJ = 77, 11)
    #knots = c(0.02, 0.14, 0.17, 0.33, 0.39, 0.5, 0.67, 0.9)
    eigenfList = get_eigenfunList_NonUnif(8, pcaTrans, knots)
    pcaKappaSqrt = sqrt( c(3,2,1, 0.7, 0.1, 0.25, 0.07, 1e-3, 5e-4, 5e-5) )
    # scoreType = "Gaussian"
    obsCol = get_allObs_NonUnif(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1/2, obsNumLower = 7, obsNumUpper = 14,
                                scoreType)
    nKnots = 25
    r.set = 6
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = (1/2)^2
    # InitType = "EM"
} else if (DataType == "Fourier"){
    fExpNum = 5
    pcaCompNum = fExpNum * 2
    tmin = 0
    tmax = 1
    pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
    pcaTrans = qr.Q(qr(pcaTrans))
    splitP = c(-1e-9, 2+1e-9)
    pcaKappa = exp(-(0:(pcaCompNum-1))/2)*10
    pcaKappaSqrt = sqrt(pcaKappa)
    eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1, obsNumLower, obsNumUpper, scoreType)
    nKnots = 18
    M.EM = nKnots + 2
    r.set = 8
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = 1
} else if (DataType == "easy"){
    tmin = 0
    tmax = 1
    eigenfList = get_eigenfunList_Paul("easy")
    obsCol = get_allObs_Paul(samplesize, M=5,r=3, "easy", "Gaussian", alpha = 0.6, nmin=2, nmax=10, a=1, b=1, sig=1/4)
    nKnots = 10
    M.EM = nKnots + 2
    r.set = 3
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = 1/16
} else if (DataType == "prac"){
    tmin = 0
    tmax = 1
    eigenfList = get_eigenfunList_Paul("prac")
    obsCol = get_allObs_Paul(samplesize, M=10, r=5, "prac", "Gaussian", alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4)
    nKnots = 15
    M.EM = nKnots + 2
    r.set = 5
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = 1/16
} else if (DataType == "hybrid"){
    tmin = 0
    tmax = 1
    eigenfList = get_eigenfunList_Paul("hybrid")
    obsCol = get_allObs_Paul(samplesize, M=10, r=10, "hybrid", "Gaussian", alpha = 0.6, nmin=2, nmax=10, a=1, b=1, sig=1/4)
    nKnots = 15
    M.EM = nKnots + 2
    r.set = 5
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = 1/16
}





### Define the spline objective 

splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
K = splineObj$getDoF()


# newObsCol = obsCol[, -2]
# newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
# colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
# 
# ini.method = "EM"
# basis.method = "bs"
# sl.v = rep(0.5, 10)
# max.step = 50
# grid.l = seq(0,1,0.01)
# grids= seq(0,1,0.002)
# #sig.EM = 1/4
# 
# # Init_method = function(newObsCol, r.set, ini.method){
# #     nmax<-max(table(newObsCol[,1]))  
# #     L1<-min(as.numeric(newObsCol[,3]))              ##range of the time 
# #     L2<-max(as.numeric(newObsCol[,3]))
# #     data.list = fpca.format(newObsCol)
# #     n<-length(data.list)    ##number of subjects
# #     if(n==0){ 
# #         print("error: no subject has more than one measurements!")
# #         return(0)
# #     }
# #     data.list.new = data.list
# #     for (i in 1:n){
# #         cur<-data.list[[i]][[1]]
# #         temp<-(cur[,2]-L1)/(L2-L1)
# #         temp[temp<0.00001]<-0.00001
# #         temp[temp>0.99999]<-0.99999 
# #         cur[,2]<-temp
# #         data.list.new[[i]][[1]]<-cur
# #     }
# #     for (k in 1:length(r.set)){
# #         r.c = r.set[k]
# #         print(paste("r=",r.c))
# #     }
# #     
# #     if (ini.method == "loc"){
# #         IniVal = Initial(r.set, ini.method, data.list.new, n, nmax, grid.l, grids, iter.num = 50, basis.EM = "ns", sig.EM)
# #     } else if (ini.method == "EM.self"){
# #         basis.EM="poly"
# #         M.self = 20
# #         IniVal = Initial(r.set, ini.method="EM",data.list.new,n,nmax,grid.l,grids,basis.EM=basis.EM,M.EM=M.self,sig.EM=sig.EM)
# #     }
# #     return(IniVal)
# # }
# 
# nmax<-max(table(newObsCol[,1]))  
# L1<-min(as.numeric(newObsCol[,3]))              ##range of the time 
# L2<-max(as.numeric(newObsCol[,3]))
# data.list = fpca.format(newObsCol)
# n<-length(data.list)    ##number of subjects
# if(n==0){ 
#     print("error: no subject has more than one measurements!")
#     return(0)
# }
# data.list.new = data.list
# for (i in 1:n){
#     cur<-data.list[[i]][[1]]
#     temp<-(cur[,2]-L1)/(L2-L1)
#     temp[temp<0.00001]<-0.00001
#     temp[temp>0.99999]<-0.99999 
#     cur[,2]<-temp
#     data.list.new[[i]][[1]]<-cur
# }
# 
# #select = EM_selection(newObsCol, M.set, r.set, basis.EM = "poly")
# IniVal = Initial(r.set, ini.method ="EM", data.list.new, n, nmax, grid, grids, basis.EM = "poly",
#                  M.EM = nKnots, sig.EM = sig.EM)
# 
# sig2hat = IniVal[[1]]
# covmatrix.ini<-IniVal[[2]]
# eigenf.ini<-IniVal[[3]]
# eigenv.ini<-IniVal[[4]]
# like.ini<-IniVal[[5]]
# 
# eigenfest = t(eigenf.ini)
# basisMat = splineObj$evalSpline(grids)
# UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
# WInit = eigenv.ini
# 
# #decomp = qr(UInit)
# #Q = qr.Q(decomp)
# #R = qr.R(decomp)
# XInit = list(UInit, diag(WInit))


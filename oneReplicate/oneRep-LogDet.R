numElem = 1
obsNumLower = 5
obsNumUpper = 10
if (DataType == "cai"){
    # nKnots = 8
    nKnots = 13 
    M.set = nKnots
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    K = splineObj$getDoF()
    #seq(6, 10, by = 1)
    pcaCompNum = 6
    eigenfList = get_eigenfunList_cai(pcaCompNum)
    obsCol = Generate_Cai(eigenfList, samplesize, pcaCompNum, noiseSigma = 1/4, alpha = 1)
    r.set = pcaCompNum
    #sig.EM = 1/4
    sig2hat = (1/4)^2
   # Init = EMInit(obsCol, splineObj, M.set, r.set, (sig.EM)^2 )
    #XInit = Init[[1]]
    #sig2hat = Init[[2]]
    #InitType = "EM"
} else if (DataType == "Peng"){
    #nKnots = 8
    # nKnots = 13
    #scoreType = "Gaussian"
    nKnots = 18
    M.set = nKnots
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)

    pcaCompNum = 4
    pcaKappa = c(1, 0.66, 0.52, 0.07)
    pcaKappaSqrt = sqrt(pcaKappa)
    # The true eigenfunction of data can be represented by cubic spline with 20 knots
    eigenfList = get_eigenfunList_exact(order = 4, nknots = 8, pcaCompNum)  
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1/4, obsNumLower, obsNumUpper, scoreType)
   # sig.EM = 1/4
    r.set = pcaCompNum
 #   Init = EMInit(obsCol, splineObj, M.set, r.set, (sig.EM)^2 )
  #  XInit = Init[[1]]
  #  sig2hat = Init[[2]]
    sig2hat = (1/4)^2
    #InitType = "EM"
} else if (DataType == "PengHybrid"){
    nKnots = 8
    #nKnots = 15
    scoreType = "Gaussian"
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    
    M.set = seq(8, 15, by = 1)
    pcaCompNum = 8 
    pcaKappa = c(1, 0.66, 0.52, 0.07, 9.47e-3, 1.28e-3, 1.74e-4, 2.35e-5 )
    pcaKappaSqrt = sqrt(pcaKappa)
    eigenfList = get_eigenfunList_exact(order = 4, 8, pcaCompNum)
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1/4, obsNumLower, obsNumUpper, scoretype)
    r.set = c(2:6)
    #sig.EM = 1/4
    sig2hat = (1/4)^2
    #InitType = "EM"
} else if (DataType == "FourierDense"){
    fExpNum = 5
    pcaCompNum = fExpNum * 2
    a = seq(0,1, length.out = 7)  ### Useless, only create a sequence of length 7
    pcaTrans = Generate_rescale_U(seedJ, a)
    #    pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
    #    pcaTrans = qr.Q(qr(pcaTrans))
    splitP = c(-1e-9, 2+1e-9)
    # numSample = 2e2
    pcaKappa = exp(-(0:(pcaCompNum-1))/2)*10
    pcaKappaSqrt = sqrt(pcaKappa)
    eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    #scoreType = "Gaussian"
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 0.7,
                        obsNumLower, obsNumUpper, scoreType)
    nKnots = 20
    M.EM = nKnots + 2
    r.set = 8
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = (0.7)^2
    #InitType = "LS"
    #InitType = "LOC"
    #InitType = "EM"
} else if (DataType == "FourierOrth"){
    fExpNum = 5
    pcaCompNum = fExpNum * 2
    pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
    pcaTrans = qr.Q(qr(pcaTrans))
    
    splitP = c(-1e-9, 2+1e-9)
    pcaKappa = c(1, 0.66, 0.52, 0.07, 9.47e-3, 1.28e-3, 1.74e-4, 2.35e-5, 3.18e-6, 4.30e-7 )
    pcaKappaSqrt = sqrt(pcaKappa)
    eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    ## Set score Type, uniform or t-distribution
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1/4, obsNumLower, obsNumUpper, scoreType, noiseType)
    ## Set knots and ranks
    M.set = c(9:18)
    nKnots = 18
    M.EM = nKnots + 2
    r.set = 4
    select_Method = "penalty"
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = (1/4)^2
    #InitType = "LS"
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
  M.set = seq(10, 20, by = 4)
  nKnots = 18
  M.EM = nKnots + 2
  r.set = 8
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1
} else if (DataType == "easy"){
  tmin = 0
  tmax = 1
  eigenfList = get_eigenfunList_Paul("easy")
  obsCol = get_allObs_Paul(samplesize, M=5,r=3, "easy", scoreType, alpha = 0.6, nmin=2, nmax=10, a=1, b=1, sig=1/4, noiseType)
  nKnots = 10
  M.EM = nKnots + 2
  select_Method = "knots"
  M.set = seq(4, 12, by = 1)  # select by knots without penalty
  r.set = 3
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1/16
} else if (DataType == "prac"){
  tmin = 0
  tmax = 1
  eigenfList = get_eigenfunList_Paul("prac")
  obsCol = get_allObs_Paul(samplesize, M=10, r=5, "prac", scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4, noiseType)
  #nKnots = 15
  nKnots = 8
  M.set = seq(8, 15, by = 1)
  M.EM = nKnots + 2
  select_Method = "knots"
  r.set = 5
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1/16
} else if (DataType == "hybrid"){
  tmin = 0
  tmax = 1
  eigenfList = get_eigenfunList_Paul("hybrid")
  obsCol = get_allObs_Paul(samplesize, M=10, r=10, "hybrid", scoreType, alpha = 0.6, nmin=2, nmax=10, a=1, b=1, sig=1/4)
  nKnots = 15
  M.EM = nKnots + 2
  M.set = seq(8, 15, by = 1)
  r.set = 5
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1/16
} else if (DataType == "easyNU") {
  tmin = 0
  tmax = 1
  shape1 = 2
  shape2 = 2
  eigenfList = get_eigenfunList_Paul_NU("easyNU", shape1, shape2)
  obsCol = get_allObs_Paul(samplesize, M=5,r=3, "easyNU", scoreType, alpha = 0.6, nmin=2, nmax=10, a=1, b=1, sig=1/4, noiseType, shape1, shape2)
  nKnots = 8
  M.EM = nKnots + 2
  
  select_Method = "knots"
  M.set = seq(4, 12, by = 1)  # select by knots without penalty
  r.set = 3
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1/16
} else if (DataType == "pracNU") {
  tmin = 0
  tmax = 1
  shape1 = 2
  shape2 = 2
  eigenfList = get_eigenfunList_Paul_NU("pracNU", shape1, shape2)
  obsCol = get_allObs_Paul(samplesize, M=10, r=5, "pracNU", scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4, noiseType, shape1, shape2)
  #nKnots = 15
  nKnots = 13
  M.set = seq(8, 15, by = 1)
  M.EM = nKnots + 2
  select_Method = "knots"
  r.set = 5
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1/16
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
    nKnots = 8
    M.EM = nKnots + 2
    select_Method = "knots"
    sig2hat = 1/16
} else if (DataType == "pracSin") {
    Mdata = 10
    rdata = 5
    pcaTrans = matrix(rnorm(Mdata*rdata), Mdata, rdata)
    pcaTrans = qr.Q(qr(pcaTrans))
    eigenfList = get_eigenfunList_PaulNoSpline("pracSin", pcaTrans)
    obsCol = get_allObs_PaulNoSpline(samplesize, M=Mdata, r=rdata, DataType, pcaTrans,
                                     scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4, noiseType)
  #  M.set = seq(8, 15, by = 1)
    M.set = seq(10, 20 ,by = 1)
    r.set = 5
    nKnots = 13
    M.EM = nKnots + 2
    select_Method = "knots"
    sig2hat = 1/16
} else if ((DataType == "NonEqualSquare") || (DataType == "NonEqual3/2") || (DataType == "NonEqualRoot")) {
  M = 10
  r = 5
  pcaTrans = matrix(rnorm(M*r), M, r)
  pcaTrans = qr.Q(qr(pcaTrans))
  eigenfList = get_eigenfunList_NonEqual(M, pcaTrans, method = DataType)
  obsCol = get_allObs_PaulNoEqual(samplesize, M, r, eigenfList, pcaTrans, scoreType,
                                  alpha = 0.6, nmin = 2, nmax = 10, a=1, b=1, sig = 1/4, 
                                  noiseType)
  M.set = seq(8, 15, by = 1)
  r.set = 5
  nKnots = 13
  M.EM = nKnots + 2
  select_Method = "knots"
  sig2hat = 1/16
} else if (DataType == "NonEqualRootFix") {
  M = 10
  r = 5
  pcaTrans = diag(M)[, c(3:5, 8:9)]
  eigenfList = get_eigenfunList_NonEqual(M, pcaTrans, method = "NonEqualRoot")
  eigenfList = get_eigenfunList_NonEqual(M, pcaTrans, method = "NonEqualRoot")
  obsCol = get_allObs_PaulNoEqual(samplesize, M, r, eigenfList, pcaTrans, scoreType,
                                  alpha = 0.6, nmin = 2, nmax = 10, a=1, b=1, sig = 1/4, 
                                  noiseType)
  M.set = seq(8, 15, by = 1)
  r.set = 5
  nKnots = 13
  M.EM = nKnots + 2
  select_Method = "knots"
  sig2hat = 1/16
}





### Define the spline objective 











#sig.EM = 1/4

# Init_method = function(newObsCol, r.set, ini.method){
#     nmax<-max(table(newObsCol[,1]))  
#     L1<-min(as.numeric(newObsCol[,3]))              ##range of the time 
#     L2<-max(as.numeric(newObsCol[,3]))
#     data.list = fpca.format(newObsCol)
#     n<-length(data.list)    ##number of subjects
#     if(n==0){ 
#         print("error: no subject has more than one measurements!")
#         return(0)
#     }
#     data.list.new = data.list
#     for (i in 1:n){
#         cur<-data.list[[i]][[1]]
#         temp<-(cur[,2]-L1)/(L2-L1)
#         temp[temp<0.00001]<-0.00001
#         temp[temp>0.99999]<-0.99999 
#         cur[,2]<-temp
#         data.list.new[[i]][[1]]<-cur
#     }
#     for (k in 1:length(r.set)){
#         r.c = r.set[k]
#         print(paste("r=",r.c))
#     }
#     
#     if (ini.method == "loc"){
#         IniVal = Initial(r.set, ini.method, data.list.new, n, nmax, grid.l, grids, iter.num = 50, basis.EM = "ns", sig.EM)
#     } else if (ini.method == "EM.self"){
#         #basis.EM="BSpline"
#         basis.EM = "poly"
#         M.self = 20
#         IniVal = Initial(r.set, ini.method="EM",data.list.new,n,nmax,grid.l,grids,basis.EM=basis.EM,M.EM=M.self,sig.EM=sig.EM)
#     } 
#     return(IniVal)
# }




#sigmaSq = sig2hat

#fList_est = Generate_fList(grids,eigenfest)
#ComputeLoss(fList_est, eigenfList)

#UInit = Q
#WInit = diag((R * WInit )%*%t(R))
#XInit = list(UInit, diag(WInit))

#nSample = length(unique(obsCol$obsID))
#cf = -1
#cvMembership = getCVPartition(nSample, nFold)
#cvParam = list(cvMembership = cvMembership, cf = cf)
#select1 = MLE_selection(obsCol, splineObj, rankSeq, lambdaSeq, controlList1, controlList2, nFold = 10, SInit = NULL, sigmaSq = 1 )
#model1 = MFPCA_EstimateMLE(obsCol, splineObj, rankSeq, 1e-7, controlList1, controlList2, cvParam, 
#                          SInit = XInit, sigmaSq = sig2hat)
#select2 = frobDiverg_selection(obsCol, splineObj , rankSeq, lambdaSeq, controlList1, controlList2, nFold = 10, SInit = NULL, sigmaSq = 1)
#model2 = MFPCA_EstimatefrobDiverg(obsCol, splineObj, rankSeq, 1e-5, controlList1, controlList2, cvParam, 
#                          SInit = XInit, sigmaSq =sig2hat)

#select3 = VNDiv_selection(obsCol, splineObj , rankSeq, lambdaSeq, controlList1, controlList2, nFold = 10, SInit =NULL, sigmaSq = 1)
#model3 = MFPCA_EstimateVNDiv(obsCol, splineObj, rankSeq, 5e-6, controlList1, controlList2, cvParam, 
#                                  SInit = XInit, sigmaSq = sig2hat)


#modelini = Convert(IniVal)
#modelini = c(modelini, 
#             list(tmin = tmin, tmax = tmax,
#                  SFinal = XInit, sigmaSq = sig2hat,
#                  numPCA = rankSeq, numElem = 1,
#                  elemLevels = levels(obsCol$elemID)))
#ComputeLoss(fList_est, eigenfList)
#evaluateLoss(model1,eigenfList)
#evaluateLoss(model2,eigenfList)
#evaluateLoss(model3,eigenfList)
#plotFpcaCompare(model1, eigenfList)
#plotFpcaCompare(model2, eigenfList)
#plotFpcaCompare(model3, eigenfList)
#plotFpcaCompare(modelini, eigenfList)





#plotFpcaCompare = function(pcaModel, trueEigenFList, selK = NULL){
#    tmin = pcaModel$tmin
#    tmax = pcaModel$tmax
#    numPCA = pcaModel$numPCA
#    numElem = pcaModel$numElem
    
#    nSeq = 500
#    tSeq = seq(tmin, tmax, length.out = nSeq)
#    plotData = data.frame();#matrix(0, nSeq * numPCA * numElem, 4)  
#    selR = 1:nSeq
#    for(k in 1:numPCA){
#        for(e in 1:numElem){
#            em = pcaModel$elemLevels[e]
#            fSeq0 = trueEigenFList[[k]][[e]](tSeq)
#            fSeqHat = pcaModel$eigenFunctions[[k]][[e]](tSeq)
#            if(sum(fSeq0*fSeqHat) < 0) fSeqHat = -fSeqHat
#            tmp = data.frame(obsT = tSeq, obsY = fSeq0, 
#                             pcaID =  k, elemID =  em,
#                             curveID = "true", stringsAsFactors =  F)
#            tmp2 = data.frame(obsT = tSeq, obsY = fSeqHat, 
#                              pcaID =  k, elemID =  em,
#                              curveID = "estimate", stringsAsFactors =  F)
#            tmp = rbind(tmp, tmp2)
#            plotData = rbind(plotData, tmp)
#            selR = selR + nSeq
#        }
#    }
#    colnames(plotData) = c("obsT", "obsY", "pcaID", "elemID", "curveID")
#    plotData$elemID = factor(plotData$elemID, levels = pcaModel$elemLevels)
    
    
#    if(!is.null(selK)){
#        plotData = subset(plotData, plotData$pcaID == selK)
        
#        p = ggplot(plotData, aes(obsT, obsY, 
#                                 group = curveID, color = curveID)) +
#            geom_line()
#        p = p + facet_wrap(~elemID)
        
#    }else{
#        p = ggplot(plotData, aes(obsT, obsY, 
#                                 group = curveID, color = curveID)) +
#            geom_line()
#        p = p + facet_wrap(pcaID~elemID)
        
#    }
#    return(p)
#}


# Convert = function(IniVal){
#     sig2hat = IniVal[[1]]
#     covmatrix.ini<-IniVal[[2]]
#     eigenf.ini<-IniVal[[3]]
#     eigenv.ini<-IniVal[[4]]
#     like.ini<-IniVal[[5]]
#     eVector = t(eigenf.ini)
#     eValues = eigenv.ini
#     #tSeq = seq(tmin, tmax, length.out = 200)
#     tSeq = seq(0,1,0.002)
#     eigenFunctions = list()
#     eigenFunctionsDeriv1 = list()
#     eigenFunctionsDeriv2 = list()
#     
#     for (k in 1:length(eValues)){
#         compFunctions = list()
#         compFunctionsDeriv1 = list()
#         compFunctionsDeriv2 = list()
#         eFunSeq = seq(0,1,0.002)
#         eFun = approxfun(tSeq, eVector[k, ])
#         
#         compFunctions = c(compFunctions, list(eFun))
#         
#         eigenFunctions = c(eigenFunctions, list(compFunctions))
#         
#     }
#     
#     model = list(eigenValues = eValues, 
#                  eigenFunctions = eigenFunctions,
#                  eigenFunctionsDeriv1 = eigenFunctionsDeriv1,
#                  eigenFunctionsDeriv2 = eigenFunctionsDeriv2)
#     return(model)
# }

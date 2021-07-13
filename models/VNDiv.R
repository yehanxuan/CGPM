createMFPCAObjVNDiv = function(dataF, splineObj, obsSigma = 1){
    if(is.character(dataF$obsID))
        dataF$obsID = factor(dataF$obsID, ordered = TRUE)
    if(is.factor(dataF$obsID))
        dataF$obsID = as.numeric(dataF$obsID)
    sortI = sort(dataF$obsID, index.return = TRUE)$ix
    obsID = dataF$obsID[sortI]
    
    elemID = dataF$elemID[sortI]
    obsT = dataF$obsT[sortI]
    #### rescale obsT to grids in EM 
    
    grids= seq(0,1,0.002)
    timeindex = floor(obsT * length(grids)) + 1
    obsT = grids[timeindex]
    
    obsY = dataF$obsY[sortI]
    # obsSigma = obsSigma[sortI]
    
    ## The number of points for each obs ID
    obsIDCount = as.numeric(table(obsID)) 
    
    bMatLarge = generate_bMatLarge(obsT, elemID, splineObj)
    
    # smoothness penalty
    Omega = splineObj$get_Omega()
    ## 重新加了不知效果如何
    #Omega = Omega / max(Omega)
    # the number of multivariate components
    elemNum = length(unique(elemID))
    Gamma = Omega
    if(elemNum > 1){
        for(e in 2:elemNum)
            Gamma = bdiag(Gamma, Omega)
    }
    Gamma = as.matrix(Gamma)
    
    optObj = new(mfpcaVNDiv, obsY, bMatLarge, obsSigma, obsIDCount)
    optObj$set_penaltyMatrix(Gamma)
    return(optObj)
}


MFPCA_EstimateVNDiv = function(obsData, splineObj, 
                               optRank, mu,
                               controlList1 = NULL,
                               controlList2 = NULL, 
                               cvParam,
                               SInit = NULL, sigmaSq = 1, eigenfList = NULL){
    tmin = splineObj$getTMin()
    tmax = splineObj$getTMax()
    
    # Check the existence of the column elemID.
    # If not there, this is the univariate FPCA
    if(is.null(obsData$elemID)){
        obsData$elemID = factor(rep(1,length(obsData$obsY)))
    }else{
        obsData$elemID = factor(obsData$elemID)
    }
    
    trainObj = createMFPCAObjVNDiv(obsData, splineObj, sigmaSq)
    
    #controlList = list(alpha = 0.5, tol = 1e-4, iterMax = 100)
    
    # Roughness penalty
    trainObj$set_tuningParameter(mu)
    
    if(is.null(SInit)){
        SInit = MFPCA_Initial(trainObj, optRank, controlList1)
    }
    if(!is.null(cvParam)){
        trainObj$setCVFold(cvParam$cvMembership)
        trainObj$activateCV(cvParam$cf)
    }
    
    estimate = MFPCA_SecondVNDiv(trainObj, optRank, controlList2, SInit, splineObj, eigenfList)
    sigmaSq =  trainObj$get_sigmaSq() 
    SFinal = estimate$SEstimate
    step = estimate$iter
    like = estimate$obj
    gap = estimate$gap 
    Time = estimate$Time
    
    loss = estimate$lossVec
    YVec = estimate$YVec
    
    # convert matrices SFinal to R functions
    numElem = length(levels(obsData$elemID))
    model = MFPCA_convert(SFinal, splineObj, optRank, numElem)
    model = c(model, 
              list(tmin = tmin, tmax = tmax,
                   SFinal = SFinal, sigmaSq = sigmaSq,
                   numPCA = optRank, numElem = numElem,
                   elemLevels = levels(obsData$elemID), YVec = YVec, step = step, like = like, gap = gap, Time = Time, loss = loss))
    
    return(model)
}


MFPCA_SecondVNDiv = function(optObj, optRank, controlList, SInit, splineObj, eigenfList){
    SEstimate = SInit
    iter = 0
    flag = TRUE
    #params = c(0.0000618, 0.9,1)
    params = c(0.1, 0.618, 0.01, 0.1,0)
    gap = 1
    tol = 1e-3
    obj = optObj$objF(SEstimate)
    objVector = c(obj)
    #Time = c(0, 0, 0)
    Time = 0
    YVec = list()
    YVec[[1]] = SInit
    if (!is.null(splineObj)){
        loss = list()
        if (!is.null(eigenfList)){
            loss[[1]] =  ComputeLoss_S(SInit, eigenfList, splineObj )
        }
    }
    while( (gap > tol) & (!is.na(gap>tol) ) ){
        
        
       # iter = iter + 1
    #     start = proc.time()
        Fit = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
        SEstimate = Fit$SFinal
        objVector = c(objVector, Fit$objVec)
        Time = c(Time, Fit$timeVec)
        iter = iter + length(Fit$timeVec)
        YVec = c(YVec, Fit$SVec)
        
        
   #     SEstimate = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
       optObj$updateSigmaSq(SEstimate, params)
        
   #     if (!is.null(splineObj)){
            # loss = c(loss, mean(ComputeLoss_S(SEstimate, eigenfList, splineObj) ) )
    #        loss[[iter+1]] = ComputeLoss_S(SEstimate, eigenfList, splineObj)  
    #    }
    #    Interval = proc.time() - start
        #Time[[iter]] = Interval
    #    Time = rbind(Time, Interval[1:3])
        obj_new = optObj$objF(SEstimate)
     #   objVector = c(objVector, obj_new)
        gap = obj - obj_new
        obj = obj_new
    }
   # Time = as.matrix(Time)
    # TimeMatrix =  apply(Time, 2, cumsum)
    TimeMatrix = cumsum(Time)
    iter = iter + 1
#    SEstimate = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
    #params = c(0.001, 0.9,1)
#    while (iter < 3){
#    optObj$updateSigmaSq(SEstimate, params)
#    SEstimate = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
#       iter = iter + 1
    #optObj$updateSigmaSq(SEstimate, params)
    #SEstimate = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
    #optObj$updateSigmaSq(SEstimate, params)
    #SEstimate = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
    #optObj$updateSigmaSq(SEstimate, params)
    return(list(SEstimate = SEstimate, iter = iter, "obj" = objVector, gap = gap, Time = TimeMatrix, lossVec = loss, YVec = YVec))
}


VNDiv_selection = function(obsCol, splineObj, rankSeq, lambdaSeq, controlList1 = NULL, controlList2, nFold = 10, 
                           sigmaSq, eigenfList = NULL,InitType = NULL, SInit = NULL){
    nSample = length(unique(obsCol$obsID))
    cvMembership = getCVPartition(nSample, nFold)
    L1 = length(rankSeq)
    L2 = length(lambdaSeq)
    K = splineObj$getDoF()
    # ErrorMat = matrix(0, L1, L2)
    ErrorMat = matrix(1e+10, L1, L2)
    for (i in 1:L1){
        if (InitType == "EM"){
            Init = EMInit(obsCol, splineObj, K-2, rankSeq[i], sigmaSq, eigenfList)
            SInit = Init[[1]]
            sigmaSq = Init[[2]]
        } else if (InitType == "LS"){
            Init = Init_LS(obsCol, splineObj, rankSeq[i], pert = 0.01, sigmaSq)
            SInit = Init[[1]]
            sigmaSq = Init[[2]]
        } else if (InitType == "LOC"){
            Init = LOCInit(obsCol, splineObj, rankSeq[i], sigmaSq)
            SInit = Init[[1]]
            sigmaSq = Init[[2]]
        }
        for (j in 1:L2){
            testerror = rep(1e+4, nFold)
            try({
                for (cf in 1:nFold){
                    
                    cvParam = list(cvMembership = cvMembership, cf = cf)
                    model = MFPCA_EstimateVNDiv(obsCol, splineObj, rankSeq[i], lambdaSeq[j], 
                                                controlList1, controlList2, cvParam, SInit = SInit, sigmaSq, eigenfList)
                    ErrorObj = createMFPCAObjVNDiv(obsCol, splineObj, obsSigma = model$sigmaSq)
                    ErrorObj$setCVFold(cvParam$cvMembership)
                    ErrorObj$activateCV(cvParam$cf)
                    ErrorObj$set_tuningParameter(0)
                    testerror[cf] = ErrorObj$outOfBagError(model$SFinal)
                }
            })
            ErrorMat[i, j] = mean(testerror)
        }
    }
    index = which(ErrorMat== min(ErrorMat), arr.ind = TRUE)
    index1 = index[1]
    index2 = index[2]
    opt_rank = rankSeq[index1]
    opt_lambda = lambdaSeq[index2]
    opt_error = ErrorMat[index1, index2]
    return(list("ErrorMat" = ErrorMat, "opt_mu" = opt_lambda, "opt_rank" = opt_rank, "opt_error" = opt_error))
}

MFPCA_VNDiv = function(obsCol, order, nknots, rankSeq, lambdaSeq, controlList1, controlList2,
                       nFold, sigmaSq, eigenfList = NULL, InitType = NULL, SInit = NULL){
    tmin = 0
    tmax = 1
    splineObj = new(orthoSpline, tmin, tmax, order, nknots)
    K = splineObj$getDoF()
    ### Add meanModel here
    meanModel = fitMeanCurve(obsCol, splineObj, lambda = 0)
    obsCol = subtractMeanCurve(meanModel, obsCol)
    #####
    select = VNDiv_selection(obsCol, splineObj, rankSeq, lambdaSeq, controlList1,
                             controlList2, nFold, sigmaSq, eigenfList, InitType = InitType)
    nSample = length(unique(obsCol$obsID))
    cf = -1
    cvMembership = getCVPartition(nSample, nFold)
    cvParam = list(cvMembership = cvMembership, cf = cf)
    if (InitType == "EM"){
        Init = EMInit(obsCol, splineObj, nknots, select$opt_rank, sigmaSq, eigenfList)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
    } else if (InitType == "LS"){
        Init = Init_LS(obsCol, splineObj, select$opt_rank, pert = 0.01, sigmaSq)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
    } else if (InitType == "LOC"){
        Init = LOCInit(obsCol, splineObj, select$opt_rank, sigmaSq)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
    }
    model = MFPCA_EstimateVNDiv(obsCol, splineObj, select$opt_rank, select$opt_mu, 
                                controlList1, controlList2, cvParam, SInit, sigmaSq, eigenfList)
    if (is.null(eigenfList)){
        loss = NULL
    } else {
        loss = evaluateLoss(model, eigenfList)
    }
    opt_lambda = select$opt_mu
    opt_rank = select$opt_rank
    return(list("opt_lambda" = opt_lambda, "opt_rank"= opt_rank, "loss" = loss, "model" = model))
}


oneReplicate_VNDiv = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("./oneReplicate/oneRep-VNDiv.R")

#    fit = MFPCA_VNDiv(obsCol, mOrder, nKnots, r.set, lambdaSeq, controlList1,
#                      controlList2, nFold, sig2hat, InitType = "EM")
    fit = MFPCA_VNDiv(obsCol, mOrder, nKnots, r.set, lambdaSeq, controlList1,
                      controlList2, nFold, sig2hat, eigenfList, InitType)
    opt_lambda = fit$opt_lambda
    opt_rank = fit$opt_rank
    loss = fit$loss
    return(list("loss" = loss, "lambda" = opt_lambda, "rank" = opt_rank))
}


oneReplicateWrap_VNDiv = function(seedJ){
    try({
        eval = oneReplicate_VNDiv(seedJ)
    })
    return(eval)
}


#VNDiv_selection(obsCol, splineObj, rankSeq, lambdaSeq = 1e-10, controlList1 = NULL, controlList2, nFold = 10, SInit = XInit, sigmaSq = 1)

createMFPCAObjPower = function(dataF, splineObj, obsSigma, b){
    if(is.character(dataF$obsID))
        dataF$obsID = factor(dataF$obsID, ordered = TRUE)
    if(is.factor(dataF$obsID))
        dataF$obsID = as.numeric(dataF$obsID)
    sortI = sort(dataF$obsID, index.return = TRUE)$ix
    obsID = dataF$obsID[sortI]
    
    elemID = dataF$elemID[sortI]
    obsT = dataF$obsT[sortI]
    obsY = dataF$obsY[sortI]
    # obsSigma = obsSigma[sortI]
    
    ## The number of points for each obs ID
    obsIDCount = as.numeric(table(obsID)) 
    
    bMatLarge = generate_bMatLarge(obsT, elemID, splineObj)
    
    # smoothness penalty
    Omega = splineObj$get_Omega()
    
    ### 重新加了Omega在这里
    # Omega = Omega / max(Omega)
    # the number of multivariate components
    elemNum = length(unique(elemID))
    Gamma = Omega
    if(elemNum > 1){
        for(e in 2:elemNum)
            Gamma = bdiag(Gamma, Omega)
    }
    Gamma = as.matrix(Gamma)
    
    optObj = new(mfpcaPower, obsY, bMatLarge, obsSigma, obsIDCount, b)
    optObj$set_penaltyMatrix(Gamma)
    return(optObj)
}

MFPCA_EstimatePower = function(obsData, splineObj, b, optRank, mu,
                               controlList1 = NULL,
                               controlList2 = NULL,
                               cvParam, SInit = NULL, sigmaSq = 1, eigenfList){
    tmin = splineObj$getTMin()
    tmax = splineObj$getTMax()
    
    # Check the existence of the column elemID.
    # If not there, this is the univariate FPCA
    if(is.null(obsData$elemID)){
        obsData$elemID = factor(rep(1,length(obsData$obsY)))
    }else{
        obsData$elemID = factor(obsData$elemID)
    }
    
    trainObj = createMFPCAObjPower(obsData, splineObj, sigmaSq, b)
    trainObj$set_tuningParameter(mu)
    if(is.null(SInit)){
        print("Initial values required!")
        break
    }
    
    if(!is.null(cvParam)){
        trainObj$setCVFold(cvParam$cvMembership)
        trainObj$activateCV(cvParam$cf)
    }
    
    estimate = MFPCA_SecondPower(trainObj, optRank, controlList2, SInit, splineObj, eigenfList)
    sigmaSq = trainObj$get_sigmaSq()
    SFinal = estimate$SEstimate 
    step = estimate$iter
    like = estimate$obj
    gap = estimate$gap 
    Time = estimate$Time
    
    YVec = estimate$YVec
    loss = estimate$lossVec
    
    numElem = length(levels(obsData$elemID))
    model = MFPCA_convert(SFinal, splineObj, optRank, numElem)
    model = c(model,
              list(tmin = tmin, tmax = tmax,
                   SFinal = SFinal, sigmaSq = sigmaSq,
                   numPCA = optRank, numElem = numElem,
                   elemLevels = levels(obsData$elemID), YVec = YVec, step = step, like = like, gap = gap, Time = Time, loss = loss))
    return(model)
    
}

MFPCA_SecondPower = function(optObj, optRank, controlList, SInit, splineObj, eigenfList){
    SEstimate = SInit
    iter = 0
    flag = TRUE
    params = c(0.1, 0.618, 0.01, 0.1, 0)
    gap = 1
    tol = 1e-3
    obj = optObj$objF(SEstimate)
    objVector = c(obj)
    
    Time = 0
    YVec = list()
    YVec[[1]] = SInit
    if (!is.null(splineObj)){
        loss = list()
        loss[[1]] =  ComputeLoss_S(SInit, eigenfList, splineObj )
    }
    
    while (gap > tol){
        Fit = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
        SEstimate = Fit$SFinal
        objVector = c(objVector, Fit$objVec)
        Time = c(Time, Fit$timeVec)
        iter = iter + length(Fit$timeVec)
        YVec = c(YVec, Fit$SVec)
        
        #optObj$updateSigmaSq(SEstimate, params)
        
        obj_new = optObj$objF(SEstimate)
        gap = obj - obj_new
        obj = obj_new
    }
    
    TimeMatrix = Time 
    iter = iter + 1
    return(list(SEstimate = SEstimate, iter = iter, "obj" = objVector, gap = gap, Time = TimeMatrix, lossVec = loss, YVec = YVec))
}


Power_selection = function(obsCol, splineObj, b, rankSeq, lambdaSeq, controlList1, controlList2, nFold = 10, sigmaSq, eigenfList, InitType = NULL, SInit = NULL){
    nSample = length(unique(obsCol$obsID))
    # nSample = nlevels(obsCol$obsID)
    cvMembership = getCVPartition(nSample, nFold)
    L1 = length(rankSeq)
    L2 = length(lambdaSeq)
    K = splineObj$getDoF()
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
                    model = MFPCA_EstimatePower(obsCol, splineObj, b,rankSeq[i], lambdaSeq[j],
                                                controlList1, controlList2, cvParam, SInit, sigmaSq, eigenfList)
                    ErrorObj = createMFPCAObjPower(obsCol, splineObj, obsSigma = model$sigmaSq, b)
                    ErrorObj$setCVFold(cvParam$cvMembership)
                    ErrorObj$activateCV(cvParam$cf)
                    ErrorObj$set_tuningParameter(0)
                    testerror[cf] = ErrorObj$outOfBagError(model$SFinal)
                }
            })
            ErrorMat[i, j] = mean(testerror)
        }
    }
    index = which(ErrorMat==min(ErrorMat), arr.ind = TRUE)
    index1 = index[1]
    index2 = index[2]
    opt_rank = rankSeq[index1]
    opt_lambda = lambdaSeq[index2]
    opt_error = ErrorMat[index1, index2]
    return(list("ErrorMat" = ErrorMat, "opt_mu" = opt_lambda, "opt_rank" = opt_rank, "opt_error" = opt_error))
}


MFPCA_Power = function(obsCol, order, nknots, b ,rankSeq, lambdaSeq, controlList1 = NULL, controlList2, 
                       nFold, sigmaSq, eigenfList, InitType = NULL, SInit = NULL){
    tmin = 0
    tmax = 1
    splineObj = new(orthoSpline, tmin, tmax, order, nknots)
    K = splineObj$getDoF()
    ### Add meanModel here
    meanModel = fitMeanCurve(obsCol, splineObj, lambda = 0)
    obsCol = subtractMeanCurve(meanModel, obsCol)
    
    select = Power_selection(obsCol, splineObj, b, rankSeq, lambdaSeq, controlList1,
                             controlList2, nFold, sigmaSq, eigenfList, InitType)
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
    
    model = MFPCA_EstimatePower(obsCol, splineObj, b, select$opt_rank, select$opt_mu,
                                controlList1, controlList2, cvParam, SInit, sigmaSq, eigenfList)
    loss = evaluateLoss(model, eigenfList)
    opt_lambda = select$opt_mu
    opt_rank = select$opt_rank
    return(list("opt_lambda" = opt_lambda, "opt_rank"= opt_rank, "loss" = loss, "model" = model))
}


oneReplicate_Power = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("./oneReplicate/oneRep-Power.R")
    fit = MFPCA_Power(obsCol, mOrder, nKnots, b, r.set, lambdaSeq, controlList1, controlList2, nFold, sig2hat,
                      eigenfList, InitType)
    opt_lambda = fit$opt_lambda
    opt_rank = fit$opt_rank
    loss = fit$loss
    return(list("loss" = loss, "lambda" = opt_lambda, "rank" = opt_rank))
}

oneReplicateWrap_Power = function(seedJ){
    try({
        eval = oneReplicate_Power(seedJ)
    })
}



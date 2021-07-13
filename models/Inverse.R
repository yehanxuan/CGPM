#### Seed function: (b/2)x^2 + 1/x
createMFPCAObjBregman = function(dataF, splineObj, obsSigma, b){
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
    
    optObj = new(mfpcaBreg, obsY, bMatLarge, obsSigma, obsIDCount, b)
    optObj$set_penaltyMatrix(Gamma)
    return(optObj)
}

MFPCA_EstimateBregman = function(obsData, splineObj, b, rank, mu, 
                                 controlList = NULL,
                                 controlList2 = NULL,
                                 cvParam,
                                 SInit = NULL, sigmaSq){
    tmin = splineObj$getTMin()
    tmax = splineObj$getTMax()
    
    # Check the existence of the column elemID.
    # If not there, this is the univariate FPCA
    if(is.null(obsData$elemID)){
        obsData$elemID = factor(rep(1,length(obsData$obsY)))
    }else{
        obsData$elemID = factor(obsData$elemID)
    }
    
    trainObj = createMFPCAObjBregman(obsData, splineObj, sigmaSq, b)
    trainObj$set_tuningParameter(mu)
    if(is.null(SInit)){
        print("Initial values required!")
        break
    }
    if(!is.null(cvParam)){
        trainObj$setCVFold(cvParam$cvMembership)
        trainObj$activateCV(cvParam$cf)
    }
    
    SFinal = MFPCA_SecondBregman(trainObj, rank, controlList2, SInit)
    sigmaSq = trainObj$get_sigmaSq()
    # convert matrices SFinal to R functions
    numElem = length(levels(obsData$elemID))
    model = MFPCA_convert(SFinal, splineObj, rank, numElem)
    model = c(model,
              list(tmin = tmin, tmax = tmax,
                   SFinal = SFinal, sigmaSq = sigmaSq, numPCA = rank, numElem = numElem,
                   elemLevels = levels(obsData$elemID)))
    return(model)
}

MFPCA_SecondBregman = function(optObj, optRank, controlList, SInit){
    SEstimate = SInit
    iter = 1
    flag = TRUE
    params = c(0.1, 0.618, 0.01, 0.1, 0)
    gap = 1
    tol = 1e-3
    obj = optObj$objF(SEstimate)
    
    while(gap > tol){
        SEstimate = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
        optObj$updateSigmaSq(SEstimate, params)
        obj_new = optObj$objF(SEstimate)
        gap = obj_new - obj
        obj = obj_new
    }
    
    return(SEstimate)
}

Bregman_selection = function(obsCol, splineObj, b, rankSeq, lambdaSeq, controlList = NULL,
                             controlList2, nFold = 10, sigmaSq, InitType = NULL, SInit = NULL){
    nSample = length(unique(obsCol$obsID))
    cvMembership = getCVPartition(nSample, nFold)
    L1 = length(rankSeq)
    L2 = length(lambdaSeq)
    K = splineObj$getDoF()
    ErrorMat = matrix(1e+10, L1, L2)
    for (i in 1:L1){
        if (InitType == "EM"){
            Init = EMInit(obsCol, splineObj, K, rankSeq[i], sigmaSq)
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
                model = MFPCA_EstimateBregman(obsCol, splineObj, b, rankSeq[i], lambdaSeq[j], 
                                              controlList1, controlList2, cvParam, SInit, sigmaSq)
                ErrorObj = createMFPCAObjBregman(obsCol, splineObj, obsSigma = model$sigmaSq, b)
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

MFPCA_Bregman = function(obsCol, order, nknots, b, rankSeq, lambdaSeq, controlList1 = NULL, controlList2,
                         nFold, sigmaSq, InitType = NULL, SInit = NULL){
    tmin = 0
    tmax = 1
    splineObj = new(orthoSpline, tmin, tmax, order, nknots)
    K = splineObj$getDoF()
    ### Add meanModel here
    meanModel = fitMeanCurve(obsCol, splineObj, lambda = 0)
    obsCol = subtractMeanCurve(meanModel, obsCol)
    
    select = Bregman_selection(obsCol, splineObj, b, rankSeq, lambdaSeq, controlList1,
                               controlList2, nFold, sigmaSq, InitType)
    nSample = length(unique(obsCol$obsID))
    cf = -1
    cvMembership = getCVPartition(nSample, nFold)
    cvParam = list(cvMembership = cvMembership, cf = cf)
    if (InitType == "EM"){
        Init = EMInit(obsCol, splineObj, nknots, select$opt_rank, sigmaSq)
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
    model = MFPCA_EstimateBregman(obsCol, splineObj, b, select$opt_rank, select$opt_mu, 
                                  controlList1, controlList2, cvParam, SInit, sigmaSq)
    loss = evaluateLoss(model, eigenfList)
    opt_lambda = select$opt_mu
    opt_rank = select$opt_rank
    return(list("opt_lambda" = opt_lambda, "opt_rank"= opt_rank, "loss" = loss, "model" = model))
}


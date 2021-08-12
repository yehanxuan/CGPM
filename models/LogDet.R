
## observation ID, n = 1, .. , N
## observation element ID, each obs has Q vector elements 
## observation obsT(time), obsY (value), obsSigma (uncertainty)
## everything sorted by obsID
## The column obsElem is optional. If not specified, the result is univariate FPCA.

#' Create MFPCA with least square objective function
#' @param dataF the data frame of the observation data. This should have four columns.
#' @param splineObj spline object 
#' @return an object of the class MFPCLoss, 
#' the objective and gradient function can be evaluated.
#' 
createMFPCAObjLS = function(dataF, splineObj){
    #obsID = as.numeric(as.factor(obsID))
    if(is.character(dataF$obsID))
        dataF$obsID = factor(dataF$obsID, ordered = TRUE)
    if(is.factor(dataF$obsID))
        dataF$obsID = as.numeric(dataF$obsID)
    
    sortI = sort(dataF$obsID, index.return = TRUE)$ix
    obsID = dataF$obsID[sortI]
    obsElem = dataF$obsElem[sortI]
    obsT = dataF$obsT[sortI]
    obsY = dataF$obsY[sortI]

    bMat = splineObj$evalSpline(obsT)
    ## The number of points for each obs ID
    obsYList = list()
    bMatList = list()
    
    for(i in unique(obsID)){
        sel = (obsID == i)
        obsYList = c(obsYList, list(obsY[sel]))
        bMatList = c(bMatList, list(bMat[,sel]))
    }
    Omega = splineObj$get_Omega()
    Omega = Omega / max(Omega)
    optObj = new(uFPCAClassv2, obsYList, bMatList)
    optObj$set_penaltyMatrix(Omega)
    return(optObj)
}



#' Title
#'
#' @param dataF 
#' @param splineObj spline object 
#' @param obsSigma 
#'
#' @return
#' @export
#'
#' @examples
createMFPCAObjMLE = function(dataF, splineObj, obsSigma = 1){
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
    
    # the number of multivariate components
    elemNum = length(unique(elemID))
    Gamma = Omega
    if(elemNum > 1){
        for(e in 2:elemNum)
            Gamma = bdiag(Gamma, Omega)
    }
    Gamma = as.matrix(Gamma)
    
    optObj = new(mfpcaMLE, obsY, bMatLarge, obsSigma, obsIDCount)
    optObj$set_penaltyMatrix(Gamma)
    return(optObj)
}

#' @param obsCol the observation dataset
#' @param lambda the tuning parameter
#' @param SInit the initial value
MFPCA_EstimateLS = function(obsCol, splineObj, 
                          optRank,
                          lambda,
                          controlList1,controlList2, cvParam,
                          SInit = NULL){
    # Only one factor, this is the 
    if(is.null(obsCol$elemID))
        obsCol$elemID = factor(rep(1,length(obsCol$obsY)))

    optObj = createMFPCAObjLS(obsCol, splineObj)
    if(is.null(SInit))
        SInit = MFPCA_Initial_LS(optObj, optRank, controlList1 )
    optObj$set_tuningParameter(lambda)
    
    if(!is.null(cvParam)){
        optObj$setCVFold(cvParam$cvMembership)
        optObj$activateCV(cvParam$cf)
    }
    
    estimate = MFPCA_SecondCore(optObj, optRank, controlList2, SInit)
    SFinal = estimate$SEstimate
    step = estimate$iter
    like = estimate$obj
    gap = estimate$gap
    Time = estimate$Time

    numElem = 1
    model = MFPCA_convert(SFinal, splineObj, optRank, numElem)
    model = c(model, 
              list(tmin = tmin, tmax = tmax,
                   SFinal = SFinal, 
                   numPCA = optRank, numElem = numElem,
                   elemLevels = levels(obsCol[,"elemID"])))
    
    return(model)
}


#' @param obsData the data frame for the observation data.
#' @param splineObj the spline object for training.
#' @param optRank Select the rank for model fitting.
#' @param mu2 the tuning parameter for smoothness penalty.
#' @param controlList1 the control parameter for the first round of optimization.
#' @param controlList2 the control parameter for the second round of optimization.
#' @param SInit Initial Value
#' @param sigmaSq the initial value for sigma squared
#' 
#' The useful version is in mFPCA package
MFPCA_EstimateMLE = function(obsData, splineObj, 
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
    
    trainObj = createMFPCAObjMLE(obsData, splineObj, sigmaSq)
    
    #controlList = list(alpha = 0.5, tol = 1e-4, iterMax = 100)
    
    # Roughness penalty
    trainObj$set_tuningParameter(mu)
    
    if(is.null(SInit))
        SInit = MFPCA_Initial(trainObj, optRank, controlList1)
    
    if(!is.null(cvParam)){
        trainObj$setCVFold(cvParam$cvMembership)
        trainObj$activateCV(cvParam$cf)
    }
    
    estimate = MFPCA_SecondMLE(trainObj, optRank, controlList2, SInit, splineObj, eigenfList)
    sigmaSq =  trainObj$get_sigmaSq()
    SFinal = estimate$SEstimate
    step = estimate$iter
    like = estimate$obj
    gap = estimate$gap
    Time = estimate$Time
     
    YVec = estimate$YVec
    loss = estimate$lossVec
    converge = estimate$converge 
    # convert matrices SFinal to R functions
    numElem = length(levels(obsData$elemID))
    model = MFPCA_convert(SFinal, splineObj, optRank, numElem)
    model = c(model, 
              list(tmin = tmin, tmax = tmax,
                   SFinal = SFinal, sigmaSq = sigmaSq,
                   numPCA = optRank, numElem = numElem,
                   elemLevels = levels(obsData$elemID), YVec = YVec, step = step, like = like, gap = gap, Time = Time, loss = loss, converge = converge))
    
    return(model)
}


# Second step estimation on the manifold
MFPCA_SecondMLE = function(optObj, optRank, controlList, SInit, splineObj = NULL, eigenfList = NULL ){
  SEstimate = SInit
  iter = 0
  flag = TRUE
  gap = 1
  tol = 1e-3
  params = c(0.1, 0.618, 0.01, 0.1,0)
  obj = optObj$objF(SEstimate)
  objVector = c(obj)
  #  Time = c(0, 0, 0)
  Time = 0
  YVec = list()
  YVec[[1]] = SInit
  if (!is.null(splineObj)){
    loss = list()
    if (!is.null(eigenfList)){
      loss[[1]] = ComputeLoss_S(SInit, eigenfList, splineObj)
    }
  }
  while(gap > tol){
    Fit = MFPCA_SecondCore(optObj, optRank, controlList, SEstimate)
    SEstimate = Fit$SFinal
    objVector = c(objVector, Fit$objVec)
    Time = c(Time, Fit$timeVec)
    iter = iter + length(Fit$timeVec)
    YVec = c(YVec, Fit$SVec)
    optObj$updateSigmaSq(SEstimate, params)
    
    obj_new = optObj$objF(SEstimate)
    
    gap = obj - obj_new 
    obj = obj_new
    
  }
  TimeMatrix = cumsum(Time)
  iter = iter + 1
  if (gap > tol) {
    converge = 0
  } else {
    converge = 1
  }
  return(list(SEstimate = SEstimate, iter = iter, obj = objVector, gap = gap, Time = TimeMatrix, lossVec = loss, YVec = YVec, converge = converge))
}

# Initial value estimation in the Euclidean space
MFPCA_Initial_LS = function(optObj, optRank, controlList){
    #list(alpha = 1e3, sigma = 0.5, tol = 1e-4, iterMax = 30)
    # total DoF, with component number included.
    splineDF = optObj$get_totalDF()
    
    problem = new(manifoldOpt)
    problem$set_euclidean(splineDF, splineDF)
    # The euclidean loss version
    problem$setObjective(optObj$objF_Euc)
    problem$setGradient(optObj$gradF_Euc)
    problem$update_control(controlList)
    problem$solve()
    SInit = problem$get_optimizer()
    SInit = (SInit + t(SInit)) / 2.0 
    SInit_eigen = eigen(SInit)
    UInit = SInit_eigen$vectors[,1:optRank]
    WInit = SInit_eigen$values[1:optRank]
    mm = min(WInit[WInit>0])
    WInit[WInit <= 0] = mm
    XInit = list(UInit, 
                 diag(WInit))
    return(XInit)
}




#MFPCA_SecondCore = function(optObj, optRank,  controlList, XInit){
    #list(alpha = 1e-2, sigma = 1e-2, tol = 1e-10,
    #        iterMax = 400, iterSubMax = 30)
    
#    splineDF = optObj$get_totalDF()
#    problem = new(manifoldOpt)
#    problem$select_algorithm("sd")
#    problem$set_eigenReg(splineDF, optRank)
    
    # The manifold loss version
#    problem$setObjective(optObj$objF)
#    problem$setGradient(optObj$gradF)
#    problem$initial_point(XInit)
#    problem$update_control(controlList)
#    problem$solve()
#    SFinal = problem$get_optimizer()
    
#    return(SFinal)
#} 



#### Model selection function
#### output: optimal lambda, optimal rank, Error Matrix
MFPCA_Selection = function(obsCol, splineObj, rankSeq, lambdaSeq, controlList1, controlList2, nFold = 10){
    nSample = length(unique(obsCol$obsID))
    cvMembership = getCVPartition(nSample, nFold)
    S1 = length(rankSeq)
    S2 = length(lambdaSeq)
    K = splineObj$getDoF()
    ErrorMat = matrix(0, S1, S2)
    for (i in 1:S1){
        if (InitType == "EM"){
            Init = EMInit(obsCol, splineObj, K, rankSeq[i], sigmaSq, eigenfList)
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
        for (j in 1:S2){
            testerror = rep(0, nFold)
            for (cf in 1:nFold){
                cvParam = list(cvMembership = cvMembership, cf = cf)
                model = MFPCA_EstimateLS(obsCol, splineObj, rankSeq[i], lambdaSeq[j], controlList1, controlList2, cvParam, NULL)
                ## compute the out of bag error
                ErrorObj = createMFPCAObjLS(obsCol, splineObj)
                ErrorObj$setCVFold(cvParam$cvMembership)
                ErrorObj$activateCV(cvParam$cf)
                testerror[cf] = ErrorObj$outOfBagError(model$SFinal)
            }
            ErrorMat[i, j] = mean(testerror)
        }
       
    }
    index = which(ErrorMat== min(ErrorMat), arr.ind = TRUE)
    index1 = index[1]
    index2 = index[2]
    opt_rank = rankSeq[index1]
    opt_lambda = lambdaSeq[index2]
    opt_error = ErrorMat[index1, index2]
    return(list("ErrorMat" = ErrorMat, "opt_lambda" = opt_lambda, "opt_rank" = opt_rank, "opt_error" = opt_error))
}

### MLE selection
MLE_selection = function(obsCol, splineObj, rankSeq, lambdaSeq, controlList1, controlList2, nFold = 10, sigmaSq, eigenfList = NULL, InitType = NULL, SInit = NULL, cvMembership = NULL){
   
    nSample = length(unique(obsCol$obsID))
   # nSample = nlevels(obsCol$obsID)
    if (is.null(cvMembership)){
      cvMembership = getCVPartition(nSample, nFold)  
    }
    
    L1 = length(rankSeq)
    L2 = length(lambdaSeq)
    K = splineObj$getDoF()
    
    #ErrorMat = matrix(0, L1, L2)
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
            testerror = rep(1e+7, nFold)
            
                for (cf in 1:nFold){
                  try({
                    cvParam = list(cvMembership = cvMembership, cf = cf)
                    model = MFPCA_EstimateMLE(obsCol, splineObj, rankSeq[i], lambdaSeq[j], 
                                              controlList1, controlList2, cvParam, SInit, sigmaSq, eigenfList)
                    ## compute the out of bag error
                    ErrorObj = createMFPCAObjMLE(obsCol, splineObj, obsSigma = model$sigmaSq)
                    ErrorObj$setCVFold(cvParam$cvMembership)
                    ErrorObj$activateCV(cvParam$cf)
                    ErrorObj$set_tuningParameter(0)
                    testerror[cf] = ErrorObj$outOfBagError(model$SFinal)
                  })
                }  
            
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


LogDet_select_knots = function(obsCol, order, r.set, M.set, controlList1, controlList2, nFold = 10, sigmaSq, eigenfList = NULL, InitType = NULL, SInit = NULL, cvMembership = NULL){
  tmin = 0
  tmax = 1
  nSample = length(unique(obsCol$obsID))
  # nSample = nlevels(obsCol$obsID)
  if (is.null(cvMembership)){
    cvMembership = getCVPartition(nSample, nFold)
  }
  
  L1 = length(r.set)
  L2 = length(M.set)
  #ErrorMat = matrix(1e+10, L1, L2)
  like.list = NULL
  for (i in 1:L1){
    r = r.set[i]
    M.subset = M.set[M.set >= r]
    M.l = length(M.subset)
    if(M.l==0){               
      print("all M < r.c")
      return(0)
    }
    like.result = numeric(M.l)
    for (j in 1:M.l){
      splineObj = new(orthoSpline, tmin, tmax, order, M.subset[j]-2)
      K = splineObj$getDoF()
      if (InitType == "EM"){
        Init = EMInit(obsCol, splineObj, M.subset[j]-2, r, sigmaSq, eigenfList)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
      } else if (InitType == "LS"){
        Init = Init_LS(obsCol, splineObj, r, pert = 0.01, sigmaSq)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
      } else if (InitType == "LOC"){
        Init = LOCInit(obsCol, splineObj,r, sigmaSq)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
      }
      testerror = rep(1e+7, nFold)
      
        for (cf in 1:nFold){
          try({
          cvParam = list(cvMembership = cvMembership, cf = cf)
          model = MFPCA_EstimateMLE(obsCol, splineObj, r.set[i], 0, 
                                    controlList1, controlList2, cvParam, SInit, sigmaSq, eigenfList)
          ErrorObj = createMFPCAObjMLE(obsCol, splineObj, obsSigma = model$sigmaSq)
          ErrorObj$setCVFold(cvParam$cvMembership)
          ErrorObj$activateCV(cvParam$cf)
          ErrorObj$set_tuningParameter(0)
          testerror[cf] = ErrorObj$outOfBagError(model$SFinal)
          })
        }
      
      like.result[j] = mean(testerror)
    }
    like.list[[i]] = like.result
  }
  
  temp = LogDet_CV(like.list, M.set, r.set)
  index.r = temp[1]
  index.M = temp[2]
  r_opt = r.set[index.r]
  M_opt = M.set[index.M]
  return(list("M_opt" = M_opt, "r_opt" = r_opt, "like_result" = like.list))
}


LogDet_CV = function(result, M.set, r.set){
  cv.mod = matrix(1e+10, length(r.set), length(M.set))
  colnames(cv.mod) = M.set
  rownames(cv.mod) = r.set
  for (k in 1:length(r.set)){
    r.c = r.set[k]
    M.set.c = M.set[M.set >=r.c]
    index.c = sum(M.set < r.c)
    if (length(M.set.c) > 0){
      for (j in 1:length(M.set.c)){
        cv.mod[k, j+index.c] = result[[k]][j]
      }
    }
  }
  
  index.r = 1
  index.M = 1
  for (j in 1:length(M.set)){
    for (k in 1:length(r.set)){
      if (cv.mod[k, j] < cv.mod[index.r, index.M]){
        index.r = k
        index.M = j
      }
    }
  }
  
  temp = c(index.r, index.M)
  names(temp) = c("r", "M")
  return(temp)
}


#### The main function, select the parameter first then output the loss 
MFPCA_leastSquares = function(obsCol, order, nknots, rankSeq, lambdaSeq, controlList1, controlList2, nFold, SInit = NULL){
    tmin = 0
    tmax = 1
    splineObj = new(orthoSpline, tmin, tmax, order, nknots)
    K = splineObj$getDoF()
    
    
    select = MFPCA_Selection(obsCol, splineObj, rankSeq, lambdaSeq, controlList1, controlList2, nFold, SInit = SInit)
    nSample = length(unique(obsCol$obsID))
    cf = -1
    cvMembership = getCVPartition(nSample, nFold)
    cvParam = list(cvMembership = cvMembership, cf = cf)
    model = MFPCA_EstimateLS(obsCol, splineObj, select$opt_rank, select$opt_lambda, controlList1, controlList2, cvParam, NULL)
    loss = evaluateLoss(model,eigenfList)
    opt_lambda = select$opt_lambda
    opt_rank = select$opt_rank
    return(list("opt_lambda" = opt_lambda, "opt_rank"= opt_rank, "loss" = loss))
} 


#### Main function for Rcpp MLE
MFPCA_RcppMLE = function(obsCol, order, nknots, rankSeq, lambdaSeq, controlList1, controlList2, nFold, sigmaSq, eigenfList = NULL, InitType = NULL,  SInit = NULL, cvMembership = NULL){
    tmin = 0
    tmax = 1
    splineObj = new(orthoSpline, tmin, tmax, order, nknots)
    K = splineObj$getDoF()
    ### Add mean here 
     meanModel = fitMeanCurve(obsCol, splineObj, lambda = 0)
     obsCol = subtractMeanCurve(meanModel, obsCol)
 #   dataF = obsCol
 #   if(is.character(dataF$obsID))
 #       dataF$obsID = factor(dataF$obsID, ordered = TRUE)
 #   if(is.factor(dataF$obsID))
 #       dataF$obsID = as.numeric(dataF$obsID)
 #    sortI = sort(dataF$obsID, index.return = TRUE)$ix
 #    obsID = obsCol$obsID[sortI]
     
 #    elemID = obsCol$elemID[sortI]
 #    obsT = obsCol$obsT[sortI]
  #   grids= seq(0,1,0.002)
  #   timeindex = floor(obsT * length(grids)) + 1
  #   obsT = grids[timeindex]
  #   obsY = dataF$obsY[sortI]
     
   #  splineBasis = new(orthoSpline, 0, 1, 4, nknots)
#     B.orth <- t(splineBasis$evalSpline(grids))
  #   B <- B.orth[timeindex, ]
     
   #  theta.zero <- solve(t(B) %*% B, t(B)) %*% obsY
    # obsY <- obsY - B %*% theta.zero
    # obsCol$obsY = obsY
     
    
    select = MLE_selection(obsCol, splineObj, rankSeq, lambdaSeq, controlList1, controlList2, nFold, sigmaSq, eigenfList, InitType = InitType, cvMembership)
    nSample = length(unique(obsCol$obsID))
    cf = -1
    if (is.null(cvMembership)){
            cvMembership = getCVPartition(nSample, nFold)
    }
    
    cvParam = list(cvMembership = cvMembership, cf = cf)
    if (InitType == "EM"){
        Init = EMInit(obsCol, splineObj, nknots, select$opt_rank, sigmaSq, eigenfList)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
    } else if (InitType == "LS"){
        Init = Init_LS(obsCol, splineObj, select$opt_rank, pert = 0.01,  sigmaSq)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
    } else if (InitType == "LOC"){
        Init = LOCInit(obsCol, splineObj, select$opt_rank, sigmaSq)
        SInit = Init[[1]]
        sigmaSq = Init[[2]]
    }
    model = MFPCA_EstimateMLE(obsCol, splineObj, select$opt_rank, select$opt_mu, controlList1, controlList2, 
                              cvParam, SInit, sigmaSq, eigenfList)
    if (!is.null(eigenfList)){
      loss = evaluateLoss(model,eigenfList)
    } else {
      loss = NULL
    }
    opt_lambda = select$opt_mu
    opt_rank = select$opt_rank
    converge = model$converge
    return(list("opt_lambda" = opt_lambda, "opt_rank"= opt_rank, "loss" = loss, "model" = model, "converge" = converge))
}


MFPCA_LogDet = function(obsCol, order, M.set, r.set, controlList1, controlList2, nFold, sigmaSq, eigenfList = NULL,
                        InitType = NULL, SInit = NULL, cvMembership = NULL){
  tmin = 0
  tmax = 1
  
  select = LogDet_select_knots(obsCol, order, r.set, M.set, controlList1, controlList2, nFold = 10,
                               sigmaSq, eigenfList, InitType = InitType, cvMembership)
  nSample = length(unique(obsCol$obsID))
  cf = -1
  if (is.null(cvMembership)){
    cvMembership = getCVPartition(nSample, nFold)
  }
  cvParam = list(cvMembership = cvMembership, cf = cf)
  splineObj = new(orthoSpline, tmin, tmax, order, select$M_opt - 2)
  meanModel = fitMeanCurve(obsCol, splineObj, lambda = 0)
  obsCol = subtractMeanCurve(meanModel, obsCol)
  if (InitType == "EM"){
    Init = EMInit(obsCol, splineObj, select$M_opt - 2, select$r_opt, sigmaSq, eigenfList)
    SInit = Init[[1]]
    sigmaSq = Init[[2]]
  } else if (InitType == "LS"){
    Init = Init_LS(obsCol, splineObj, select$r_opt, pert = 0.01,  sigmaSq)
    SInit = Init[[1]]
    sigmaSq = Init[[2]]
  } else if (InitType == "LOC"){
    Init = LOCInit(obsCol, splineObj, select$r_opt, sigmaSq)
    SInit = Init[[1]]
    sigmaSq = Init[[2]]
  }
  model = MFPCA_EstimateMLE(obsCol, splineObj, select$r_opt, 0, controlList1, controlList2, 
                            cvParam, SInit, sigmaSq, eigenfList)
  if (!is.null(eigenfList)){
    loss = evaluateLoss(model,eigenfList)
  } else {
    loss = NULL
  }
  M_opt = select$M_opt
  r_opt = select$r_opt
  converge = model$converge 
  return(list("opt_lambda" = M_opt, "opt_rank" = r_opt, "loss" = loss, "model" = model, "converge" = converge ))
}



### One replication of mFPCA(our) methods
oneReplicate_mFPCA = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("./code/oneRep-mFPCA.R")
    
    fit = MFPCA_leastSquares(obsCol, mOrder, nKnots, rankSeq, lambdaSeq, controlList1, controlList2, nFold, SInit=XInit)
    opt_lambda = fit$opt_lambda
    opt_rank = fit$opt_rank
    loss = fit$loss
    return(list("loss" = loss, "lambda" = opt_lambda, "rank" = opt_rank))
}

oneReplicateWrap_mFPCA = function(seedJ){
    try({
        eval = oneReplicate_mFPCA(seedJ)
        })
    
    return(eval)
}

### restricted MLE in fpca package
# oneReplicate_MLE = function(seedJ){
#     set.seed(seedJ + repID * 300)
#     source("./code/oneRep-MLE.R")
#     result = fpca.mle(newObsCol, M.set, r.set, ini.method, basis.method, sl.v, max.step, grid.l, grids)
#     grids.new = result$grid
#     M_opt = result$selected_model[[1]]
#     r_opt = result$selected_model[[2]]
#     eigenfest = result$eigenfunctions
#     # True function interpolation
#     eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
#     fList_est = Generate_fList(grids.new, eigenfest)
#     loss = ComputeLoss(fList_est, eigenfList)
#     return(list( "loss" = loss, "M_opt" = M_opt, "rank" = r_opt))
# }
# 
# oneReplicateWrap_MLE = function(seedJ){
#     try({
#        eval = oneReplicate_MLE(seedJ) 
#     })
#     return(eval)
# }

oneReplicate_ini = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("./code/oneRep-ini.R")
#    res = FPCA(DatanewObsCol$Ly, DatanewObsCol$Lt, list(methodSelectK = 8))
    if (ini.method == "loc"){
        IniVal = Initial(r.set, ini.method, data.list.new, n, nmax, grid.l, grids, iter.num = 100, basis.EM = "ns", sig.EM)
    } else if (ini.method == "EM.self"){
        basis.EM="poly"
        M.self = 20
        IniVal = Initial(r.set, ini.method="EM",data.list.new,n,nmax,grid.l,grids,basis.EM=basis.EM,M.EM=M.self,sig.EM=sig.EM)
    } else if (ini.method == "EM"){
        select = EM_selection(newObsCol, M.set, r.set, basis.EM = "poly")
        IniVal = Initial(r.set, ini.method ="EM", data.list.new, n, nmax, grid, grids, basis.EM = "poly",
                         M.EM = select$M_opt, sig.EM = sig.EM)
    }
    
    
    sig2hat = IniVal[[1]]
    covmatrix.ini<-IniVal[[2]]
    eigenf.ini<-IniVal[[3]]
    eigenv.ini<-IniVal[[4]]
    like.ini<-IniVal[[5]]
#    grids.new = res$workGrid
    eigenfest = t(eigenf.ini)
    eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    fList_est = Generate_fList(grids, eigenfest)
    loss = ComputeLoss(fList_est, eigenfList)
    r_opt = r.set
    return(list("loss" = loss, "result2" = 0, "rank" = r_opt))
}

oneReplicateWrap_ini = function(seedJ){
    try({
        eval = oneReplicate_ini(seedJ)
    })
    return(eval)
}

oneReplicate_LogDet = function(seedJ){
#    seedJ = seedJ
    set.seed(seedJ + repID * 300)
    source("./oneReplicate/oneRep-LogDet.R")
    cvMembership = getCVPartition_seed(samplesize, nFold = 10, seedJ)
    if (select_Method == "knots"){
      fit = try(MFPCA_LogDet(obsCol, mOrder, M.set, r.set, controlList1, controlList2, nFold, sig2hat, eigenfList, InitType, cvMembership) )
    } else if (select_Method == "penalty"){
      fit = try(MFPCA_RcppMLE(obsCol, mOrder, nKnots, r.set, lambdaSeq, controlList1, controlList2, nFold, sig2hat, eigenfList, InitType, cvMembership) )
    }
    
    if (inherits(fit, "try-error")){
      converge = 0
      # opt_lambda = fit$opt_lambda
      #opt_rank = fit$opt_rank
      #loss = fit$loss
      opt_lambda = NULL
      opt_rank = NULL
      loss = NULL
    } else {
      converge = fit$converge 
      # converge = 1
      opt_lambda = fit$opt_lambda
      opt_rank = fit$opt_rank
      loss = fit$loss
    }
#    fit = MFPCA_RcppMLE(obsCol, mOrder, nKnots, r.set, lambdaSeq, controlList1, controlList2, nFold, sig2hat, InitType = "EM")
    
    return(list("loss" = loss, "lambda" = opt_lambda, "rank" = opt_rank, "converge" = converge))
}

oneReplicateWrap_LogDet = function(seedJ){
     try({
       eval = oneReplicate_LogDet(seedJ)
    })
    return(eval)
}

Init_LS = function(dataF, splineObj, r, pert = 0.01, sigmaSq){
    #   if(is.character(dataF$obsID))
    #    dataF$obsID = factor(dataF$obsID, ordered = TRUE)
    #if(is.factor(dataF$obsID))
    #    dataF$obsID = as.numeric(dataF$obsID)
    #sortI = sort(dataF$obsID, index.return = TRUE)$ix
    #obsID = dataF$obsID[sortI]
    newObsCol = dataF[, -2]
    newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
    colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
    basis.method = "bs"
    sl.v = rep(0.5, 10)
    max.step = 50
    grid.l = seq(0,1,0.01)
    grids= seq(0,1,0.002)
    
    sig.EM = sqrt( sigmaSq )
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
    
    data.obj = Format.EM(data.list.new, n, nmax, grids)
    
    obj = data.obj
    y <- obj$y
    ### data represented as a single vector (by stacking the observations for different curves)
    
    timeindex <- obj$timeindex
    ### a vector with the index of location (on the grid) of each point in y
    
    curve <- obj$curve
    ### a vector with the curve number corresponding to each value of y
    
    N <- obj$n
    ### number of observations (= n, in our notation)
    
    
    
    #elemID = dataF$elemID[sortI]
    #obsT = dataF$obsT[sortI]
    
    #grids= seq(0,1,0.002)
    #timeindex = floor(obsT * length(grids)) + 1
    #obsT = grids[timeindex]
    
    #obsY = dataF$obsY[sortI]
    # obsSigma = obsSigma[sortI]
    #N = length(unique(obsID))
    K = splineObj$getDoF()
    
    
   # grids= seq(0,1,0.002)
    # R.inv<-R.inverse(K,grids)
    #lmin<-0
    #lmax<-1
    #delta<-(lmax-lmin)/(K-3)
    #knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
    #bs<-apply(matrix(grids),1, BC.orth,knots=knots,R.inv=R.inv)
    #bs<-t(bs)*sqrt(grids[2]-grids[1])
    #B<-bs[obsT,]
    
    # B = t(splineObj$evalSpline(obsT))
    B.orth = t(splineObj$evalSpline(grids))
    B = B.orth[timeindex, ]
    ## The number of points for each obs ID
    #obsIDCount = as.numeric(table(obsID))
    # bMatLarge = generate_bMatLarge(obsT, elemID, splineObj)
    gamma = matrix(0, N, K)
    # for (i in 1:N){
    #     sel = (obsID == i)
    #     X = t(bMat[ , sel])
    #     tempy = obsY[sel]
    #     gamma[i, ] = solve(t(X) %*% X + pert * diag(K), t(X)) %*% tempy
    # }
    
    #obsY = y
    #obsID = curve
    theta.zero = solve(t(B) %*% B, t(B)) %*% y
    y <- y - B%*%theta.zero
    
    alpha <- list(alpha = matrix(1, N, r), alphaprod = array(1, c(N, r, r)))
    theta <- as.matrix(init(y, timeindex = NULL, curve, B, r, pert = 0.01)$theta)
    alpha <- getemalpha(y, curve, theta, B, alpha) 
    sigmaSq <- getsigma(y, curve, B, theta, alpha, sigmaSq)
    alpha <- getemalpha(y, curve, theta, B, alpha, sigma = sigmaSq) 
    
    UInit = theta
    WInit = solve(alpha$Dinv)
#    prc = prcomp(gamma)
#    UInit = as.matrix(prc$rotation[ ,1:r])
    # smoothness penalty
 #   Omega = splineObj$get_Omega()
    XInit = list(UInit, WInit)
    sighat_Init = sigmaSq
    return(list(XInit, sighat_Init))
}







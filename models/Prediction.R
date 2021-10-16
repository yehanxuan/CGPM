MeanModel_GCV = function(obsData, splineObj, lambdaSeq, ySigma = NULL){
    tVec = obsData[, "obsT"]
    yVec = obsData[, "obsY"]
    multiElement = "elemID" %in% colnames(obsData)
    Omega = splineObj$get_Omega()
    Bmat = splineObj$evalSpline(tVec)
    Bmat = t(Bmat)
    BmatW = Bmat
    GCV = rep(0, length(lambdaSeq))
    for (i in 1:length(lambdaSeq)){
        if(!multiElement){
            meanModel = fitMeanCurveSingle(tVec, yVec, splineObj, 
                                           lambdaSeq[i], ySigma)
            meanModel = list(modelList = meanModel,
                             elemLevels = NULL)
        }else{
            elemID = obsData[,"elemID"]
            meanModel = fitMeanCurveMultiple(tVec, yVec, splineObj, 
                                             lambdaSeq[i], elemID, ySigma)
        }
        #beta = meanModel$modelList[[1]]$beta
        
        tmp = t(BmatW) %*% BmatW + lambdaSeq[i] * Omega
        S = BmatW %*% solve(tmp, t(BmatW)) 
        res = meanModel$modelList[[1]]$residualVec
        GCV[i] = 1/length(yVec) * sum(res^2)/( 1 - sum(diag(S))/ length(yVec) )^2
    }
    index = which.min(GCV)
    opt_lambda = lambdaSeq[index]
    return(opt_lambda)
}

Split_data = function(test_data, splineObj){
    # Split the data into two part, one estimate the score, the other 
    # make prediction
    if(is.character(test_data$obsID)){
        test_data$obsID = factor(test_data$obsID, ordered = TRUE) }
    if(is.factor(test_data$obsID)){
        test_data$obsID = as.numeric(test_data$obsID)}
    sortI = sort(test_data$obsID, index.return = TRUE)$ix
    obsID = test_data$obsID[sortI]
    
    elemID = test_data$elemID[sortI]
    obsT = test_data$obsT[sortI]
    obsY = test_data$obsY[sortI]
    obsIDCount = as.numeric(table(obsID))
    
    bMat = splineObj$evalSpline(obsT)
    
    pred_obsYList = list()
    pred_bMatList = list()
    pred_TList= list()
    est_obsYList = list()
    est_bMatList = list()
    est_TList = list()
    
    for (i in unique(obsID)){
        sel = (obsID == i)
        size = length(obsY[sel])
        pred_size = floor(0.5 * size)
        estimate_size = size - pred_size
        
        
        index = which(sel == TRUE)
        pred_ind = sort(sample(index, pred_size))
        estimate_ind = index[!index %in% pred_ind]
        pred_obsYList = c(pred_obsYList, list(obsY[pred_ind]))
        est_obsYList = c(est_obsYList, list(obsY[estimate_ind]))
        pred_bMatList = c(pred_bMatList, list(bMat[ ,pred_ind, drop = F]))
        est_bMatList = c(est_bMatList, list(bMat[,estimate_ind, drop = F]))
        pred_TList = c(pred_TList, list(obsT[pred_ind]))
        est_TList = c(est_TList, list(obsT[estimate_ind]))
    }
    DataList = list("est_obsYList" = est_obsYList, "pred_obsYList" = pred_obsYList, "est_bMatList" = est_bMatList,
                    "pred_bMatList" = pred_bMatList, "est_TList" = est_TList, "pred_TList" = pred_TList)
    return(DataList)
}


Estimate_score = function(model, splineObj, DataList, meanModel = NULL){
    est_obsYList = DataList$est_obsYList
    pred_obsYList = DataList$pred_obsYList
    est_bMatList = DataList$est_bMatList
    pred_bMatList = DataList$pred_bMatList
    est_TList = DataList$est_TList
    pred_TList = DataList$pred_TList
    
    pred_YList = list()
    varList = list()
    for (i in 1:length(est_obsYList)){
        Best = est_bMatList[[i]]
        
        yVec_est = est_obsYList[[i]]
        tVec_est = est_TList[[i]]
        U = model$SFinal[[1]]
        W = model$SFinal[[2]]
        h = t(U)%*%Best
        SigmaMat = t(h)%*%W%*%h + model$sigmaSq * diag(dim(Best)[2])
        tmp1 = W%*%h
        
        if (is.null(model$meanEst)){
            
            if (!is.null(meanModel)){
                score_exp = tmp1%*%solve(SigmaMat, yVec_est - meanModel$modelList[[1]]$meanFunction(tVec_est))
            } else {
                score_exp = tmp1%*%solve(SigmaMat, yVec_est)
            }
            score_cov = W - tmp1 %*% solve(SigmaMat, t(tmp1))
            
            Bpred = pred_bMatList[[i]]
            yVec_true = pred_obsYList[[i]]
            tVec_pred = pred_TList[[i]]
            tmp2 = t(Bpred)%*%U
            # modified, something wrong with the prediction 
            ypred = tmp2%*%score_exp
           # ypred = ( tmp2%*% diag( sqrt(diag(W)))  )%*%score_exp
            if (!is.null(meanModel)){
                ypred = ypred + meanModel$modelList[[1]]$meanFunction(tVec_pred)
            }
            
        } 
        var_pred = diag( tmp2%*%score_cov%*%t(tmp2) )  + model$sigmaSq
        
        pred_YList = c(pred_YList, list(ypred))
        varList = c(varList, list(var_pred))
    }
    
    return( list("pred_YList" = pred_YList, "varList" = varList, "pred_TList" = pred_TList) )
}

###################
##Data Processing##
###################

# Process the real data
myFunctionForApply <- function(x, ...) {
    # Do your processing
    # Let's say it ends up in variable 'ret':
    if (length(x) == 0)
        return(NA)
    return(x)
}



# Produce the astronomical data 
Astro_data = function(file_list){
    numElem = 1
    obsID = 1
    obsMat = data.frame()
    for (i in 1:length(file_list)){
        Table = read_tsv( paste0("./CSP_Photometry_DR2/", file_list[i]), skip = 4, na = "99.900")
        colnames(Table)[1] = "MJD"
        
        obsT = Table$MJD
        obsY_u = Table$u
        obsY_B = Table$B
        obsY_V = Table$V
        obsY_g = Table$g
        obsY_r = Table$r
        
        #obsY_i = Table$i
        
        ## Clean the NA data
        tmp1 = cbind(obsT[!is.na(obsY_u)], obsY_u[!is.na(obsY_u)])
        tmp2 = cbind(obsT[!is.na(obsY_B)], obsY_B[!is.na(obsY_B)])
        tmp3 = cbind(obsT[!is.na(obsY_V)], obsY_V[!is.na(obsY_V)])
        tmp4 = cbind(obsT[!is.na(obsY_g)], obsY_g[!is.na(obsY_g)])
        tmp5 = cbind(obsT[!is.na(obsY_r)], obsY_r[!is.na(obsY_r)])
        
        
        Ind_min1 = myFunctionForApply( which.min(tmp1[, 2]) )
        Ind_min2 = myFunctionForApply( which.min(tmp2[, 2]) )
        Ind_min3 = myFunctionForApply( which.min(tmp3[, 2]) )
        Ind_min4 = myFunctionForApply( which.min(tmp4[, 2]) )
        Ind_min5 = myFunctionForApply( which.min(tmp5[, 2]) )
        
        Ind_min_vec = c(Ind_min1, Ind_min2, Ind_min3, Ind_min4, Ind_min5)
        
        
        if (Ind_min2 > 1){
            maxpoint = tmp2[Ind_min2, ]
            for (i in 1:5){
                assign( paste0("tmp", i), sweep(get(paste0("tmp", i)), 2, maxpoint) )
                #get(paste0("tmp", i))[ ,2:3] = get(paste0("tmp", i))[ ,2:3] - maxpoint
            }
            
            tmp1 = tmp1[ (tmp1[, 1] > -10)&(tmp1[, 1] < 50), ]
            tmp2 = tmp2[ (tmp2[, 1] > -10)&(tmp2[, 1] < 50), ]
            tmp3 = tmp3[ (tmp3[, 1] > -10)&(tmp3[, 1] < 50), ]
            tmp4 = tmp4[ (tmp4[, 1] > -10)&(tmp4[, 1] < 50), ]
            tmp5 = tmp5[ (tmp5[, 1] > -10)&(tmp5[, 1] < 50), ]
            for (i in 1:5){
                # Need to use & 
                if ( (Ind_min_vec[i] > 1) & (!is.na(Ind_min_vec[i]) ) ){
                    Tmp = cbind(obsID, numElem, get(paste0("tmp", i)) )
                    obsMat = rbind(obsMat, Tmp)
                    obsID = obsID + 1
                }
            }
        }
    }
    
    colnames(obsMat) =c("obsID", "elemID", "obsT", "obsY")
    
    nmax = max(table(obsMat[ ,1]))
    L1 = min(as.numeric(obsMat[ ,3]))
    L2 = max(as.numeric(obsMat[ ,3]))
    tmp = (obsMat[ ,3] - L1)/(L2 - L1)
    tmp[tmp<0.00001]<-0.00001
    tmp[tmp>0.99999]<-0.99999 
    
    ##### 不做scaling怎么样
    tmp2 = as.numeric(obsMat[,4])
    
    
    obsCol = obsMat
    obsCol[ ,3] = tmp
    obsCol[ ,4] = -tmp2
    colnames(obsCol) =c("obsID", "elemID", "obsT", "obsY")
    obsCol$elemID = as.factor(obsCol$elemID)
    return(obsCol)
}


MSFE = function(pred_YList, pred_obsYList){
    Error = 0
    for (i in 1:length(pred_obsYList)){
        ypred = pred_obsYList[[i]]
        yVec = pred_YList[[i]]
        err = mean( (ypred - yVec)^2 ) 
        Error = Error + err
    }
    Error = Error/length(pred_obsYList)
    return(Error)
}


MFPCA_pred_new = function(train_origin, mOrder, nKnots, r.set, lambdaSeq, controlList1, controlList2,
                          nFold, sig2hat, method, InitType, cvMembership ){
    tmin = 0
    tmax = 1
    
    if (method == "LogDetPenalty") {
        splineObj = new(orthoSpline, 0, 1, 4, nKnots)
        meanSeq = exp(seq(-10,-2, length.out = 7))
        mean_mu = MeanModel_GCV(train_origin, splineObj, meanSeq)
        meanModel = fitMeanCurve(train_origin, splineObj, lambda = mean_mu)
        train_data = subtractMeanCurve(meanModel, train_origin)
        
        fit = MFPCA_RcppMLE(train_data, mOrder, nKnots, r.set, lambdaSeq, controlList1, controlList2, nFold, sig2hat, eigenfList = NULL,
                            InitType)
        opt_lambda  = fit$opt_lambda
        opt_rank = fit$opt_rank
        model = fit$model
        Time = tail(model$Time, n = 1) 
        
        DataList = Split_data(test_data, splineObj)
        pred_result = Estimate_score(model, splineObj, DataList, meanModel)
        MSE = MSFE(pred_result$pred_YList, DataList$pred_obsYList)
    } else if (method == "EM"){
        fit_EM = MFPCA_EM(train_origin, M.set, r.set, sig2hat)
        opt_lambda = fit_EM$opt_knots 
        opt_rank = fit_EM$opt_rank
        model = fit_EM$model
        Time = tail(fit_EM$Time, n = 1)
        
        splineObj = new(orthoSpline, 0, 1, 4, opt_lambda - 2)
        meanModel = fitMeanCurve(train_origin, splineObj, 0)
        DataList = Split_data(test_data, splineObj) 
        pred_result = Estimate_score(model, splineObj, DataList, meanModel)  
        MSE = MSFE(pred_result$pred_YList, DataList$pred_obsYList)
    } else if (method == "REML"){
        fit_REML = MFPCA_REML(train_origin, M.set, r.set, InitType, sig2hat, 
                              CVmethod = "like")
        
        opt_lambda = fit_REML$opt_knots 
        opt_rank = fit_REML$opt_rank
        model = fit_REML$model
        Time = tail(fit_REML$Time, n = 1)
        
        splineObj = new(orthoSpline, 0, 1, 4, opt_lambda - 2)
        meanModel = fitMeanCurve(train_origin, splineObj, 0)
        DataList = Split_data(test_data, splineObj)  
        pred_result = Estimate_score(model, splineObj, DataList, meanModel)
        MSE = MSFE(pred_result$pred_YList, DataList$pred_obsYList)
    } else if (method == "LogDet") {
        fit_Knots = MFPCA_LogDet(train_origin, mOrder, M.set, r.set, controlList1, controlList2, 
                                 nFold, sig2hat, InitType = InitType)
        
        opt_lambda = fit_Knots$opt_lambda
        opt_rank = fit_Knots$opt_rank
        model = fit_Knots$model
        Time = tail(model$Time, n = 1)
        meanModel = fit_Knots$meanModel 
        
        splineObj = new(orthoSpline, 0, 1, 4, opt_lambda - 2)
        DataList = Split_data(test_data, splineObj)  # Maybe seed should add here 
        pred_result = Estimate_score(model, splineObj, DataList, meanModel)  # Meanmodel missed here 
        MSE = MSFE(pred_result$pred_YList, DataList$pred_obsYList)
    }
    
    
    return(list("MSFE" = MSE, "Time" = Time, "lambda" = opt_lambda, "rank" = opt_rank )) 
}



oneReplicate_pred = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("simuSettings/simuSetting-pred.R")
    source("oneReplicate/oneRep-pred.R")
    samplesize = length(unique(train_origin$obsID))
    cvMembership = getCVPartition_seed(samplesize, nFold = 10, seedJ)
    fit = MFPCA_pred_new(train_origin, mOrder, nKnots, r.set, lambdaSeq, controlList1, controlList2, nFold, sig2hat, method, InitType, cvMembership)
    Error = fit$MSFE
    lambda = fit$lambda
    rank = fit$rank
    Time = fit$Time
    return(list("MSFE" = Error, "Time" = Time, "lambda" = lambda, "rank" = rank))
}


oneReplicateWrap_pred = function(seedJ){
    
    try({
        eval = oneReplicate_pred(seedJ)
    })
    return(eval)
}





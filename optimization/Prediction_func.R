Split_data = function(test_data, splineObj){
    
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
        #obsYList = c(obsYList, list(obsY[sel]))
        #bMatList = c(bMatList, list(bMat[, sel]))
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
        ypred = tmp2%*%score_exp
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


### Predict on the grid interval 
Estimate_score_whole = function(model, splineObj, DataList, meanModel = NULL){
    est_obsYList = DataList$est_obsYList
    pred_obsYList = DataList$pred_obsYList
    est_bMatList = DataList$est_bMatList
    pred_bMatList = DataList$pred_bMatList
    est_TList = DataList$est_TList
    pred_TList = DataList$pred_TList
    #pred_TList = list()
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
            
         # If you want to plot, use these sentence     
           # predT = seq(0.01, 0.99, 0.01)
            #Bpred = splineObj$evalSpline(predT)
            #tVec_pred = predT
            
            Bpred = pred_bMatList[[i]]
            yVec_true = pred_obsYList[[i]]
            predT = pred_TList[[i]]
            tVec_pred = predT
            
            tmp2 = t(Bpred)%*%U
            ypred = tmp2%*%score_exp
            if (!is.null(meanModel)){
                ypred = ypred + meanModel$modelList[[1]]$meanFunction(tVec_pred)
            }
            
        } 
        var_pred = diag( tmp2%*%score_cov%*%t(tmp2) )  + model$sigmaSq
       # pred_TList = c(pred_TList, list(predT))
        pred_YList = c(pred_YList, list(ypred))
        varList = c(varList, list(var_pred))
    }
    
    return( list("pred_YList" = pred_YList, "varList" = varList, "pred_TList" = pred_TList) )
}   



#DataList = Split_data(train_data, splineObj)


#pred_YList = pred_result$pred_YList
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




#### Prediction functions 
plotFpcaCompare = function(pcaModel, trueEigenFList, selK = NULL){
    tmin = pcaModel$tmin
    tmax = pcaModel$tmax
    numPCA = pcaModel$numPCA
    numElem = pcaModel$numElem
    
    nSeq = 500
    tSeq = seq(tmin, tmax, length.out = nSeq)
    plotData = data.frame();#matrix(0, nSeq * numPCA * numElem, 4)  
    selR = 1:nSeq
    for(k in 1:numPCA){
        for(e in 1:numElem){
            em = pcaModel$elemLevels[e]
            fSeq0 = trueEigenFList[[k]][[e]](tSeq)
            fSeqHat = pcaModel$eigenFunctions[[k]][[e]](tSeq)
            if(sum(fSeq0*fSeqHat) < 0) fSeqHat = -fSeqHat
            tmp = data.frame(obsT = tSeq, obsY = fSeq0, 
                             pcaID =  k, elemID =  em,
                             curveID = "true", stringsAsFactors =  F)
            tmp2 = data.frame(obsT = tSeq, obsY = fSeqHat, 
                              pcaID =  k, elemID =  em,
                              curveID = "estimate", stringsAsFactors =  F)
            tmp = rbind(tmp, tmp2)
            plotData = rbind(plotData, tmp)
            selR = selR + nSeq
        }
    }
    colnames(plotData) = c("obsT", "obsY", "pcaID", "elemID", "curveID")
    plotData$elemID = factor(plotData$elemID, levels = pcaModel$elemLevels)
    
    
    if(!is.null(selK)){
        plotData = subset(plotData, plotData$pcaID == selK)
        
        p = ggplot(plotData, aes(obsT, obsY, 
                                 group = curveID, color = curveID)) +
            geom_line()
        p = p + facet_wrap(~elemID)
        
    }else{
        p = ggplot(plotData, aes(obsT, obsY, 
                                 group = curveID, color = curveID)) +
            geom_line()
        p = p + facet_wrap(pcaID~elemID)
        
    }
    return(p)
}


Convert = function(IniVal){
    sig2hat = IniVal[[1]]
    covmatrix.ini<-IniVal[[2]]
    eigenf.ini<-IniVal[[3]]
    eigenv.ini<-IniVal[[4]]
    like.ini<-IniVal[[5]]
    eVector = t(eigenf.ini)
    eValues = eigenv.ini
    #tSeq = seq(tmin, tmax, length.out = 200)
    tSeq = seq(0,1,0.002)
    eigenFunctions = list()
    eigenFunctionsDeriv1 = list()
    eigenFunctionsDeriv2 = list()
    
    for (k in 1:length(eValues)){
        compFunctions = list()
        compFunctionsDeriv1 = list()
        compFunctionsDeriv2 = list()
        eFunSeq = seq(0,1,0.002)
        eFun = approxfun(tSeq, eVector[k, ])
        
        compFunctions = c(compFunctions, list(eFun))
        
        eigenFunctions = c(eigenFunctions, list(compFunctions))
        
    }
    
    model = list(eigenValues = eValues, 
                 eigenFunctions = eigenFunctions,
                 eigenFunctionsDeriv1 = eigenFunctionsDeriv1,
                 eigenFunctionsDeriv2 = eigenFunctionsDeriv2)
    return(model)
}

plotRealEstimate = function(pcaModel, selK = NULL){
    tmin = pcaModel$tmin
    tmax = pcaModel$tmax
    numPCA = pcaModel$numPCA
    numElem = pcaModel$numElem
    
    nSeq = 500
    tSeq = seq(tmin, tmax, length.out = nSeq)
    plotData = data.frame();#matrix(0, nSeq * numPCA * numElem, 4)  
    selR = 1:nSeq
    for(k in 1:numPCA){
        for(e in 1:numElem){
            em = pcaModel$elemLevels[e]
            #fSeq0 = trueEigenFList[[k]][[e]](tSeq)
            fSeqHat = pcaModel$eigenFunctions[[k]][[e]](tSeq)
            #if(sum(fSeq0*fSeqHat) < 0) fSeqHat = -fSeqHat
            #tmp = data.frame(obsT = tSeq, obsY = fSeq0, 
            #                 pcaID =  k, elemID =  em,
            #                 curveID = "true", stringsAsFactors =  F)
            tmp = data.frame(obsT = tSeq, obsY = fSeqHat, 
                             pcaID =  k, elemID =  em,
                             curveID = "estimate", stringsAsFactors =  F)
            #tmp = rbind(tmp, tmp2)
            plotData = rbind(plotData, tmp)
            selR = selR + nSeq
        }
    }
    colnames(plotData) = c("obsT", "obsY", "pcaID", "elemID", "curveID")
    plotData$elemID = factor(plotData$elemID, levels = pcaModel$elemLevels)
    
    
    if(!is.null(selK)){
        plotData = subset(plotData, plotData$pcaID == selK)
        
        p = ggplot(plotData, aes(obsT, obsY, 
                                 group = curveID, color = curveID)) +
            geom_line()
        p = p + facet_wrap(~elemID)
        
    }else{
        p = ggplot(plotData, aes(obsT, obsY, 
                                 group = curveID, color = curveID)) +
            geom_line()
        p = p + facet_wrap(pcaID~elemID)
        
    }
    return(p)
}

######
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





###### 

EM_estimate = function(train_data, rank, splineObj, ini.method){

newObsCol = train_data[, -2]
newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
    

#ini.method = "EM.self"
basis.method = "bs"
sl.v = rep(0.5, 10)
max.step = 50
grid.l = seq(0,1,0.01)
grids= seq(0,1,0.002)
sig.EM = 1

IniVal = Init_method(newObsCol, rank, ini.method)
sig2hat = IniVal[[1]]
covmatrix.ini<-IniVal[[2]]
eigenf.ini<-IniVal[[3]]
eigenv.ini<-IniVal[[4]]
like.ini<-IniVal[[5]]

eigenfest = t(eigenf.ini)
#obsMat$elemID  = as.factor(obsMat$elemID)

K = splineObj$getDoF()
basisMat = splineObj$evalSpline(grids)

UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
WInit = eigenv.ini
XInit = list(UInit, diag(WInit))
modelini = Convert(IniVal)
modelini = c(modelini, 
             list(tmin = tmin, tmax = tmax,
                  SFinal = XInit, sigmaSq = sig2hat,
                  numPCA = rank, numElem = 1,
                  elemLevels = levels(train_data$elemID) ))

return(modelini)
}



Newton_estimate = function(train_data, M.set, r.set, splineObj){
    tmin = 0
    tmax = 1
    newObsCol = train_data[, -2]
    newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
    colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
    
    
    ini.method = "EM"
    basis.method = "bs"
    sl.v = rep(0.5, 10)
    max.step = 50
    grid.l = seq(0,1,0.01)
    grids= seq(0,1,0.002)
    sig.EM = 1
    
    iter.EM.num = 50   ## bug in the code of fpca package code
    result = fpca.mle(newObsCol, M.set, r.set, ini.method , basis.method, sl.v, max.step, grid.l, grids)
    grids.new = result$grid
    M_opt = result$selected_model[[1]]
    r_opt = result$selected_model[[2]]
    eigenfest = result$eigenfunctions
    K = splineObj$getDoF()
    basisMat = splineObj$evalSpline(grids.new)
    
    UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
    WInit = result$eigenvalues
    XInit = list(UInit, diag(WInit))
    modelMLE = Convert_Newton(result)
    modelMLE = c(modelMLE, list(tmin = tmin, tmax = tmax,
                                SFinal = XInit, sigmaSq = result$error_var,
                                numPCA = r_opt, numElem = 1,
                                elemLevels = levels(train_data$elemID), M_opt = M_opt) )
}



Convert_Newton = function(result){
    eigenfest = result$eigenfunctions
    eVector = eigenfest
    eValues = result$eigenvalues
    grids.new = result$grid
    #tSeq = seq(tmin, tmax, length.out = 200)
    tSeq = grids.new
    eigenFunctions = list()
    eigenFunctionsDeriv1 = list()
    eigenFunctionsDeriv2 = list()
    
    for (k in 1:length(eValues)){
        compFunctions = list()
        compFunctionsDeriv1 = list()
        compFunctionsDeriv2 = list()
        eFunSeq = seq(0,1,0.002)
        eFun = approxfun(tSeq, eVector[k, ])
        
        compFunctions = c(compFunctions, list(eFun))
        
        eigenFunctions = c(eigenFunctions, list(compFunctions))
        
    }
    
    model = list(eigenValues = eValues, 
                 eigenFunctions = eigenFunctions,
                 eigenFunctionsDeriv1 = eigenFunctionsDeriv1,
                 eigenFunctionsDeriv2 = eigenFunctionsDeriv2)
    return(model)
}

######## Prediction wrap

training = function(train_data, splineObj, rank,mu, method = "VNDiv", is_mean = FALSE, mean_mu = 0){
#    if ( is_mean == TRUE ){
#        meanModel = fitMeanCurve(train_data, splineObj, lambda = mean_mu)
#        train_data = subtractMeanCurve(meanModel, train_data)
#    }
    
    nSample = length(unique(train_data$obsID))
    cf = -1
    cvMembership = getCVPartition(nSample, nFold = 10)
    cvParam = list(cvMembership = cvMembership, cf = cf)
    controlList1 = list(alpha = 0.618, tol = 1e-6, iterMax = 200, sigma = 0.05,verbose = 0)
    #controlList2 = list(alpha = 1e-1, tol = 1e-7, sigma = 1e-3, beta = 0.618,
    #                    iterMax = 500, verbose = 0)
    
    controlList2 = list(alpha = 0.1, tol = 1e-7, sigma = 5e-2, beta = 0.618,
                        iterMax = 2000, verbose = 1)
    
    modelini = EM_estimate(train_data, rank, splineObj)
    XInit = modelini$SFinal
    sig2hat = modelini$sigmaSq
    if (method == "VNDiv"){
        #select = VNDiv_selection(train_data, splineObj, 4, lambdaSeq, controlList1, controlList2, nFold = 10, SInit = NULL, sigmaSq = 0.01)
        model = MFPCA_EstimateVNDiv(train_data, splineObj, rank, mu, controlList1, controlList2, cvParam,
                                    SInit = XInit, sigmaSq = sig2hat)
    } else if (method == "frobDiverg"){
        
        model = MFPCA_EstimatefrobDiverg(train_data, splineObj, rank, mu, controlList1, controlList2, cvParam,
                                         SInit = XInit, sigmaSq = sig2hat)
    } else if (method == "LogDet"){
        #select = MLE_selection(train_data, splineObj, 4, lambdaSeq, controlList1, controlList2, nFold = 10, SInit = NULL, sigmaSq = 0.01)
        model = MFPCA_EstimateMLE(train_data, splineObj, rank, mu, controlList1, controlList2, cvParam, SInit = XInit, sigmaSq = sighat)
    } else if (method == "EM"){
        model = modelini
    }
    
    return(model)
}

Prediction_Wrap = function(DataList, model, splineObj, meanModel){
    #DataList = Split_data(test_data, splineObj)
    pred_result = Estimate_score(model, splineObj, DataList, meanModel)
    MSFE = MSFE(pred_result$pred_YList, DataList$pred_obsYList)
    return(list("pred_result" = pred_result, "MSFE" = MSFE))
}

MakePlotData = function(X, Time, se, Curve){
    plotData = cbind(0, se, Time, X)
    plotData = as.data.frame(plotData)
    colnames(plotData) = c("Curve", "se", "Time", "X")
    plotData$Curve = Curve
    return(plotData)
}


######## Real data processing



myFunctionForApply <- function(x, ...) {
    # Do your processing
    # Let's say it ends up in variable 'ret':
    if (length(x) == 0)
        return(NA)
    return(x)
}


###### EM algorithm
Init_method = function(newObsCol, M = NULL, r.set, ini.method){
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
    for (k in 1:length(r.set)){
        r.c = r.set[k]
        print(paste("r=",r.c))
    }
    
    if (ini.method == "loc"){
        IniVal = Initial(r.set, ini.method, data.list.new, n, nmax, grid.l, grids, iter.num = 50, basis.EM = "ns", sig.EM)
    } else if (ini.method == "EM.self"){
        basis.EM="poly"
        M.self = 20
        IniVal = Initial(r.set, ini.method="EM",data.list.new,n,nmax,grid.l,grids,basis.EM=basis.EM,M.EM=M.self,sig.EM=sig.EM)
    } else if (ini.method == "EM"){
        basis.EM = "poly"
        IniVal = Initial(r.set, ini.method = "EM",data.list.new,n,nmax,grid.l,grids,basis.EM=basis.EM,M.EM=M, sig.EM=sig.EM)
    }
    return(IniVal)
}



#########
MFPCA_pred = function(train_data, order,nknots, rank, lambdaSeq, controlList1, controlList2, nFold, SInit, sigmaSq = 0.1){
    
    tmin = 0
    tmax = 1
    splineObj = new(orthoSpline, tmin, tmax, order, nknots)
    nSample = length(unique(train_data$obsID))
    cf = -1
    cvMembership = getCVPartition(nSample, nFold = 10)
    cvParam = list(cvMembership = cvMembership, cf = cf)
    
    if (method == "LogDet"){
        select = MLE_selection(train_data, splineObj, rank, lambdaSeq, controlList1, controlList2,
                          nFold, SInit =SInit, sigmaSq = sigmaSq)
        model = MFPCA_EstimateMLE(train_data, splineObj, rank, select$opt_mu, controlList1, controlList2, cvParam, SInit = SInit,  sigmaSq)
        opt_mu = select$opt_mu
    }else if (method == "VNDiv"){
        select = VNDiv_selection(train_data, splineObj, rank, lambdaSeq, controlList1, controlList2,
                                   nFold, SInit = SInit, sigmaSq)
        model = MFPCA_EstimateVNDiv(train_data, splineObj, rank, select$opt_mu, controlList1, controlList2,
                                        cvParam, SInit, sigmaSq)
        opt_mu = select$opt_mu
    } else if (method == "frobDiverg"){
            select = frobDiverg_selection(train_data, splineObj, rank, lambdaSeq, controlList1, controlList2, nFold, 
                                          SInit = XInit, sigmaSq = sig2hat)
            model = MFPCA_EstimatefrobDiverg(train_data, splineObj, rank, select$opt_mu, controlList1, controlList2,
                                             cvParam, XInit, sig2hat)
            opt_mu = select$opt_mu
    } else if (method == "EM"){
            model = modelini
            opt_mu = 0
    } else if (method == "Newton"){
        model = Newton_estimate(train_data, M.set, r.set, splineObj)
        opt_mu = 2
    } else if (method == "loc"){
        model = modelini
        opt_mu = 1
    }
    pred_result = Estimate_score(model, splineObj, DataList, meanModel)
    MSFE = MSFE(pred_result$pred_YList, DataList$pred_obsYList)
    return(list("MSFE" = MSFE, "lambda" = opt_mu )) 
}

MFPCA_pred_new = function(train_data, mOrder, nKnots, rank, lambdaSeq, controlList1, controlList2,
                          nFold, sig2hat, method, InitType, cvMembership ){
    tmin = 0
    tmax = 1
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    if (method == "LogDet") {
        fit = MFPCA_RcppMLE(train_data, mOrder, nKnots, rank, lambdaSeq, controlList1, controlList2, nFold, sig2hat, eigenfList = NULL,
                            InitType)
        opt_lambda  = fit$opt_lambda
        opt_rank = fit$opt_rank
        model = fit$model
    } else if (method == "VNDiv"){
        fit = MFPCA_VNDiv(train_data, mOrder, nKnots, rank, lambdaSeq, controlList1, controlList2, nFold, sig2hat, eigenfList = NULL,
                          InitType)
        opt_lambda = fit$opt_lambda
        opt_rank = fit$opt_rank
        model = fit$model
    } else if (method == "frobDiverg") {
        fit = MFPCA_frobDiverg(train_data, mOrder, nKnots, rank, lambdaSeq, controlList1, controlList2, nFold, sig2hat,
                               eigenfList = NULL, InitType)
        opt_lambda = fit$opt_lambda
        opt_rank = fit$opt_rank
        model = fit$model
    } else if (method == "EM"){
        fit = MFPCA_EM(train_data, M.set, rank, sig2hat, splineObj)
        opt_lambda = fit$opt_knots 
        opt_rank = fit$opt_rank
        model = fit$model
    } else if (method == "REML"){
        fit = MFPCA_REML(train_data, M.set, rank, ini.method = "EM", sig2hat, splineObj, eigenfList = NULL,CVmethod = "like", cvMembership = cvMembership)
        opt_lambda = fit$opt_knots
        opt_rank = fit$opt_rank
        model = fit$model
    } else if (method == "LOC"){
        fit = MFPCA_LOC(train_data, rank, splineObj)
        opt_lambda = 0 
        opt_rank = rank 
        model = fit 
    } else if (method == "LogDetKnots") {
        fit = MFPCA_LogDet(train_data, mOrder, M.set, rank, controlList1, controlList2, 
                           nFold, sig2hat, eigenfList = NULL, InitType)
        opt_lambda = fit$opt_lambda
        opt_rank = fit$opt_rank
        model = fit$model 
    }
    
    if (method == "LogDetKnots"){
        splineObj_new = new(orthoSpline, 0, 1, mOrder, opt_lambda-2)
        pred_result = Estimate_score(model, splineObj_new, DataList, meanModel)
    }
    pred_result = Estimate_score(model, splineObj, DataList, meanModel)
    MSFE = MSFE(pred_result$pred_YList, DataList$pred_obsYList)
    return(list("MSFE" = MSFE, "lambda" = opt_lambda, "rank" = opt_rank )) 
}






oneReplicate_pred = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("simuSettings/simuSetting-pred.R")
    source("oneReplicate/oneRep-pred.R")
    samplesize = length(unique(train_data$obsID))
    cvMembership = getCVPartition_seed(samplesize, nFold = 10, seedJ)
    # fit = MFPCA_pred(train_data, mOrder, nKnots, rank, lambdaSeq, controlList1, controlList2, nFold = 10, SInit = XInit, sigmaSq = sig2hat)
    fit = MFPCA_pred_new(train_data, mOrder, nKnots, rank, lambdaSeq, controlList1, controlList2, nFold, sig2hat, method, InitType, cvMembership)
    MSFE = fit$MSFE
    lambda = fit$lambda
    rank = fit$rank
    return(list("MSFE" = MSFE, "lambda" = lambda, "rank" = rank))
}


oneReplicateWrap_pred = function(seedJ){
    
    try({
        eval = oneReplicate_pred(seedJ)
    })
    return(eval)
}


oneReplicateWrap_pred_new = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("./code/oneRep-pred.R")
    fit = MFPCA_pred
}



####### Add a new function here to select the model of EM
ini_estimate  = function(train_data, rank, splineObj, ini.method, M.set = NULL){
    newObsCol = train_data[, -2]
    newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
    colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
    
    
    #ini.method = "EM.self"
    basis.method = "bs"
    sl.v = rep(0.5, 10)
    max.step = 50
    grid.l = seq(0,1,0.01)
    grids= seq(0,1,0.002)
    sig.EM = 1
    
    
    if ( (ini.method == "loc") || (ini.method == "EM.self")  ){
    IniVal = Init_method(newObsCol, M.set, rank, ini.method)
    sig2hat = IniVal[[1]]
    covmatrix.ini<-IniVal[[2]]
    eigenf.ini<-IniVal[[3]]
    eigenv.ini<-IniVal[[4]]
    like.ini<-IniVal[[5]]
    meanEst = IniVal[[7]]
    meanModel = NULL
    }  else if (ini.method == "EM"){
        if (is.null(M.set)) {
            print("warning: M.set required")
            break
        } 
        select = EM_selection(newObsCol, M.set, rank,  basis.EM = "poly")
        IniVal = Init_method(newObsCol, select$M_opt, rank, ini.method )
        sig2hat = IniVal[[1]]
        covmatrix.ini<-IniVal[[2]]
        eigenf.ini<-IniVal[[3]]
        eigenv.ini<-IniVal[[4]]
        like.ini<-IniVal[[5]]
        meanEst = IniVal[[7]]
        
        
        
        elemLevels = 1
        modelList = list()
        k = 0
        for (e in elemLevels){
            k = k + 1
            fitSeq = approxfun(newObsCol$obsT, meanEst)
            modelS = list(tmin = 0, tmax = 1, meanFunction = fitSeq)
            modelList = c(modelList, list(modelS))
        }
        meanModel = list(elemLevels = elemLevels,
                         modelList = modelList)
        
        
        
    }
    eigenfest = t(eigenf.ini)
    #obsMat$elemID  = as.factor(obsMat$elemID)
    
    K = splineObj$getDoF()
    basisMat = splineObj$evalSpline(grids)
    
    UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
    WInit = eigenv.ini
    XInit = list(UInit, diag(WInit))
    modelini = Convert(IniVal)
    modelini = c(modelini, 
                 list(tmin = tmin, tmax = tmax,
                      SFinal = XInit, sigmaSq = sig2hat,
                      numPCA = rank, numElem = 1,
                      elemLevels = levels(train_data$elemID), meanModel = meanModel))
    
    return(modelini)
}


EM_selection = function(newObsCol, M.set,r, basis.EM = "poly"){
    
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
    
    M.set = M.set[M.set >= r]
    M.l = length(M.set)
    if(M.l==0){              ##if all M in the set M.set are less than r, return zero  
        print("all M<r")
        return(0)
    }
    
    eigen.result<-array(0,dim=c(1,r,M.l))   ##estimated eigenvalues for three methods(Newton,ini, ini/EM) under different M
    eigenf.result<-NULL                     ##estimated eigenfunctions 
    sig.result<-matrix(0,1,M.l)             ##estimated error sd. for three methods under differen M
    like.result<-matrix(0,1,M.l)
    
    names.m <- c("EM")
    
    for (i in 1:M.l){
        M.EM = M.set[i]
        temp.EM = EM(data.list.new, n, nmax, grids, M.EM, iter.num = 50, r, basis.EM = basis.EM, sig.EM)
        #temp = Initial(r, ini.method = "EM", data.list.new, n, nmax, grid.l, grids, M.EM, basis.EM = basis.EM)
        EMeigenvec.est<-temp.EM[[1]]*sqrt(length(grids))
        EMeigenval.est<-temp.EM[[2]]
        EMsigma.est<-temp.EM[[3]]
        covmatrix.EM<-EMeigenvec.est%*%diag(EMeigenval.est)%*%t(EMeigenvec.est)
        like.EM<-loglike.cov.all(covmatrix.EM,EMsigma.est,data.list.new ,n)
        
        sig2hat<-EMsigma.est^2
        covmatrix.ini<-covmatrix.EM
        eigenf.ini<-EMeigenvec.est
        eigenv.ini<-EMeigenval.est
        like.ini<-like.EM
        ### -2*loglike
        like.result[, i] = like.EM
    }

    index = which.min(like.result)
    M_opt = M.set[index]
    return(list("M_opt" = M_opt, "like_result" = like.result  ) )
}


loc_meanModel = function(train_data){
    newObsCol = train_data[, -2]
    newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
    colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
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
    
    temp<-LocLin.mean(data.list.new,n,nmax, grids)  ##subtract estimated mean
    #data.list.new<-temp[[1]]
    fitmu<-temp[[2]]
    
    elemLevels = 1
    modelList = list()
    k = 0
    for (e in elemLevels){
        k = k + 1
        fitSeq = approxfun(grids, fitmu)
        modelS = list(tmin = 0, tmax = 1, meanFunction = fitSeq)
        modelList = c(modelList, list(modelS))
    }
    meanModel = list(elemLevels = elemLevels,
                     modelList = modelList)
    return(meanModel)
}


LOC_PACE = function(train_origin, rank, splineObj){
    
    sortI = sort(train_origin$obsID, index.return = TRUE)$ix
    obsID = train_origin$obsID[sortI]
    elemID = train_origin$elemID[sortI]
    obsT = train_origin$obsT[sortI]
    obsY = train_origin$obsY[sortI]
    obsIDCount = as.numeric(table(obsID))
    
    Ly = list()
    Lt = list()
    for (i in unique(obsID)){
        sel = (obsID == i)
        Ly = c(Ly, list(obsY[sel]))
        Lt = c(Lt, list(obsT[sel]))
    }
    
    
    
    #Estimate = FPCA(Ly, Lt, list(methodBwMu = 'CV', methodBwCov = 'CV', useBinnedCov = FALSE, nRegGrid = 100, maxK = rank))
    #Estimate = FPCA(Ly, Lt, list(userBwMu = 0.5, maxK = rank))
    Estimate = FPCA(Ly, Lt, list(userBwCov = 0.27,  FVEthreshold = 1, maxK = rank) ) 
    sig2hat = Estimate$sigma2
    sig2hat = Estimate$sigma2
    
    BwMu = Estimate$bwMu
    BwCov = Estimate$bwCov
    eigenfest = t(Estimate$phi)
    #obsMat$elemID  = as.factor(obsMat$elemID)
    
    K = splineObj$getDoF()
    basisMat = splineObj$evalSpline(Estimate$workGrid)
    
    UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
    WInit = Estimate$lambda
    XInit = list(UInit, diag(WInit))
    
    
    model = Convert_loc(Estimate)
    model = c(model, list(tmin = tmin, tmax = tmax,
                          SFinal = XInit, sigmaSq = sig2hat,
                          numPCA = rank, numElem = 1,
                          elemLevels = levels(train_origin$elemID), BwMu = BwMu, BwCov = BwCov))
    return(model)
}


Convert_loc = function(Estimate){
    sig2hat = Estimate$sigma2
    eigenf.loc = Estimate$phi
    eigenv.loc = Estimate$lambda
    eVector = t(eigenf.loc)
    eValues = eigenv.loc
    tSeq = Estimate$workGrid
    eigenFunctions = list()
    eigenFunctionsDeriv1 = list()
    eigenFunctionsDeriv2 = list()
    mean = Estimate$mu
    
    elemLevels = 1
    modelList = list()
    k = 0
    for (e in elemLevels){
        k = k + 1
        fitSeq = approxfun(tSeq, mean)
        modelS = list(tmin = 0, tmax = 1, meanFunction = fitSeq)
        modelList = c(modelList, list(modelS))
    }
    meanModel = list(elemLevels = elemLevels,
                     modelList = modelList)
    
    for (k in 1:length(eValues)){
        compFunctions = list()
        compFunctionsDeriv1 = list()
        compFunctionsDeriv2 = list()
        eFunSeq = Estimate$workGrid
        eFun = approxfun(tSeq, eVector[k, ])
        
        compFunctions = c(compFunctions, list(eFun))
        
        eigenFunctions = c(eigenFunctions, list(compFunctions))
        
    }
    
    model = list(eigenValues = eValues, 
                 eigenFunctions = eigenFunctions,
                 eigenFunctionsDeriv1 = eigenFunctionsDeriv1,
                 eigenFunctionsDeriv2 = eigenFunctionsDeriv2,
                 meanModel = meanModel)
    return(model)
}


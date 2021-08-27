
# getEigenFunction = function( pcaID, elemID, splitP, fExpNum, pcaTrans){
#     tSeq = seq(0, 2, length.out = 20000)
#     tSeq = tSeq[-c(1, length(tSeq))]
#     denseFourior = get_fourierBasisMat(tSeq, fExpNum)%*% pcaTrans
#     elemCut = as.numeric(cut(tSeq, breaks = splitP))
#     selElem = (elemID == elemCut)
#     eigenf = denseFourior[selElem, pcaID]
#     diffP = splitP[elemID+1] - splitP[elemID]
#     tSeq = (tSeq[selElem] - splitP[elemID]) / diffP
#     result = cbind(tSeq, eigenf * sqrt(diffP))
#     colnames(result) = c("tSeq", "eigenf")
#     efun = approxfun(result[,1], result[,2], rule = 2)
#     return(efun)
# }
# 
# getEigenFHat = function(trueF,splineObj, DEigenVectorHat,
#                         pcaID, elemID){
#     denseSpline = splineObj$evalSpline(trueF[,"tSeq"])
#     splineDF = nrow(splineObj$get_Omega())
#     selR = 1:splineDF + splineDF*(elemID-1)
#     phiHat = denseSpline %*% DEigenVectorHat[selR, pcaID]
#     ss = sign(sum(phiHat * trueF[,"eigenf"]))
#     result = cbind(trueF[,"tSeq"], ss * phiHat)
#     colnames(result) = c("tSeq", "eigenf")
#     return(result)
# }


evaluateLoss = function(model, eigenfList){
    pcaCompNum = model$numPCA
    numElem = model$numElem
    tmin = model$tmin
    tmax = model$tmax
    tSeq = seq(tmin, tmax, length.out = 5e3)
    deltaT = tSeq[2] - tSeq[1]
    lossMat = matrix(0, pcaCompNum, numElem)
    for(pcaID in 1:pcaCompNum){
        for(elemID in 1:numElem){
            hatF = model$eigenFunctions[[pcaID]][[elemID]](tSeq)
            trueF = eigenfList[[pcaID]][[elemID]](tSeq)
            ss = sign(sum(hatF * trueF))
            diffF = hatF - trueF * ss
            lossMat[pcaID, elemID] = sum(diffF^2) * deltaT
        }
    }
    return(sqrt(lossMat))
}




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

compute = function(resultT){
    res = Reduce("+", resultT)/length(resultT)
    resultS2 = lapply(resultT, function(tmp) tmp^2)
    res2 = Reduce("+", resultS2)/length(resultT)
    resSD = sqrt(res2 - res^2)/sqrt(length(resultT))
    return(list(res = res, resSD = resSD))
}


Generate_fList = function(grids, eigenfunctions){
    fList = list()
    numPCA = nrow(eigenfunctions) 
    for (pcaID in 1:numPCA){
        tmp = list()
        efun = approxfun(grids, eigenfunctions[pcaID, ], rule = 2)
        tmp = c(tmp, list(efun))
        fList = c(fList, list(tmp))
    }
    return(fList)
}

ComputeLoss = function(fList_est , fList_true){
    #pcaCompNum = length(fList_est)
  pcaCompNum = length(fList_true)
    lossMat = matrix(0, pcaCompNum, 1)
    tSeq = seq(0, 1, length.out = 5e3)
    deltaT = tSeq[2] - tSeq[1]
    numElem = 1
    for (pcaID in 1:pcaCompNum){
        hatF = fList_est[[pcaID]][[numElem]](tSeq)
        trueF = fList_true[[pcaID]][[numElem]](tSeq)
        ss = sign(sum(hatF * trueF))
        diffF = hatF - trueF * ss
        lossMat[pcaID, numElem] = sum(diffF^2) * deltaT
    }
    return(sqrt(lossMat))
}


StoF = function(SHat, splineObj){
    splineDF = splineObj$getDoF()
    tmin = splineObj$getTMin()
    tmax = splineObj$getTMax()
    
  #  eiDec = eigen(SHat[[2]])
  #  eVector = SHat[[1]] %*% eiDec$vectors
  #  eValues = eiDec$values
    
    # tSeq = seq(tmin, tmax, length.out = 200)
    tSeq = seq(0,1,0.002)
    bMatSeq = t(splineObj$evalSpline(tSeq))
    
    # eigenmat = bMatSeq %*% eVector
    eigenmat = bMatSeq %*% SHat[[1]]
    
    eigenFunctions = list()
    numPCA = dim(SHat[[2]])[1]
    numElem = 1
    
    Gamma = eigenmat %*% (SHat[[2]] %*% t(eigenmat))
    Gamma.svd  = eigen(Gamma, symmetric=TRUE)
    eigenval.est = Gamma.svd$values[1: numPCA]/length(tSeq)
    # Add drop = F 8/27/2021
    eigenvec.est = Gamma.svd$vectors[ , 1:numPCA, drop = F]
    
        for (k in 1:numPCA){
        selR = 1:splineDF
        compFunctions = list()
        
            for (e in 1:numElem){
            # eFunSeq = bMatSeq %*% eVector[selR, k] ### 这个地方应该没有正则化
            eFunSeq = eigenvec.est[ , k, drop = F] * sqrt(length(tSeq))
            ###  Add something here to make it
            #Gamma = eFunSeq %*% (eValues * t(eFunSeq))
            eFun = approxfun(tSeq, eFunSeq)
            compFunctions = c(compFunctions, list(eFun))
            selR = selR + splineDF
        }
        
        eigenFunctions = c(eigenFunctions, list(compFunctions))
    }
    return(eigenFunctions)
}


ComputeLoss_S = function(S, fList_true, splineObj){
    fList_est = StoF(S, splineObj)
    pcaCompNum = length(fList_est)
    lossMat = matrix(0, pcaCompNum, 1)
    tSeq = seq(0, 1, length.out = 5e3)
    deltaT = tSeq[2] - tSeq[1]
    numElem = 1
    for (pcaID in 1:pcaCompNum){
        hatF = fList_est[[pcaID]][[numElem]](tSeq)
        trueF = fList_true[[pcaID]][[numElem]](tSeq)
        ss = sign(sum(hatF * trueF))
        diffF = hatF - trueF * ss
        lossMat[pcaID, numElem] = sum(diffF^2) * deltaT
    }
    return(sqrt(lossMat))
}

ComputeLoss_SVec = function(YVec, fList_true, splineObj){
  LossList = list()
  for (i in 1:length(YVec)){
    S = YVec[[i]]
    LossList[[i]] = ComputeLoss_S(S, fList_true, splineObj)
  }
  return(LossList)
}

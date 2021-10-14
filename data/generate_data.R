# Get the simulation data for the first settings


library(tibble)
library(mvtnorm)
#' Repeat get_oneObs to get multiple observation groups
#' 
#' @param numSample (int) total number of observation groups
get_allObs = function(eigenFList,
                      numSample, pcaKappaSqrt, noiseSigma, 
                      obsNumLower, obsNumUpper, scoreType, noiseType = "Gaussian"){
    # pre allocation of memory to speed up code
    numElem = length(eigenFList[[1]])  # True, How many variate 
    obsMat = matrix(0, numSample*obsNumUpper*numElem, 4)
    startI = 1
    for(obsID in 1:numSample){
        tmp = get_oneObs(eigenFList, pcaKappaSqrt, noiseSigma, 
                         obsNumLower, obsNumUpper, scoreType, noiseType)
        tmp = cbind(obsID, tmp)
        endI = startI + nrow(tmp) - 1
        #cat(endI, nrow(obsMat),"\n")
        obsMat[startI:endI,] =  tmp
        startI = endI + 1
    }
    obsMat = data.frame(obsMat)
    colnames(obsMat) = c("obsID", "elemID", "obsT", "obsY")
    obsMat = obsMat[1:(endI-1),]
    obsMat$elemID = factor(obsMat$elem)
    ## Add a modification here, introducing a factor here
    obsMat$obsID = factor(obsMat$obsID)
    return(obsMat)
}


## Generate one group of observation
#' @param eigenFList (list) the list of list of eigenfunctions
#' @param pcaKappaSqrt FCPA scores
#' @param noiseSigma (double) the standard deviation of obs noise
#' @param obsNumLower (int) the lower bound of obs number
#' @param obsNumUpper (int) the upper bound of obs number
get_oneObs = function(eigenFList,
                      pcaKappaSqrt, noiseSigma, 
                      obsNumLower, obsNumUpper, scoreType, noiseType = "Gaussian"){

    numElem = length(eigenFList[[1]])
    numPCA = length(eigenFList)
    if (scoreType == "Gaussian"){
        score = rnorm(numPCA, sd = pcaKappaSqrt)
    } else if (scoreType == "uniform"){
        score = runif(numPCA, min = -sqrt(3), max = sqrt(3))
        score = score*pcaKappaSqrt
    } else if (scoreType == "t"){
        dfv = 3
        rtvar = dfv/(dfv-2)
        score = rt(numPCA, df = dfv)
        score = score/sqrt(rtvar)
        score = score * pcaKappaSqrt
    }
    
    #score = rnorm(numPCA, sd = pcaKappaSqrt)
#    dfv = 6
#    rtvar = dfv / (dfv - 2)
#    score = rt(numPCA, df = dfv)
    #score = rtvar / sqrt(rtvar)
#    score = score * pcaKappaSqrt
    result = numeric(0)
    for(e in 1:numElem){
        elemObsNum = sample(obsNumLower:obsNumUpper, 1)
        
        obsT = runif(elemObsNum, 0, 1)
        if (noiseType == "Gaussian"){
            obsY = rnorm(elemObsNum, sd = noiseSigma)
        } else if (noiseType == "t") {
            obsY = rt(elemObsNum, df = 3)*noiseSigma/1.725
        } else if (noiseType == "uniform") {
            obsY = runif(elemObsNum, min = -sqrt(3), max = sqrt(3))*noiseSigma
        }
        
        for(k in 1:numPCA){
            obsY = obsY + score[k] * eigenFList[[k]][[e]](obsT)
        }
        tmp = cbind(e, obsT, obsY)
        result = rbind(result, tmp)
    }
    
    return(result)
}


get_eigenfunList = function(pcaTrans, fExpNum, splitP){
    tSeq = seq(0, 2, length.out = 30000)
    tSeq = tSeq[-c(1, length(tSeq))]
    denseFourior = get_fourierBasisMat(tSeq, fExpNum)%*% pcaTrans

    numPCA = ncol(pcaTrans)
    numElem = length(splitP) - 1
    
    elemCut = as.numeric(cut(tSeq, breaks = splitP)) # All the factor is 1
    result = list()
    for(pcaID in 1:numPCA){
        tmp = list()
        for(elemID in 1:numElem){
            selElem = (elemID == elemCut)
            diffP = splitP[elemID + 1] - splitP[elemID]
            eigenT = (tSeq[selElem] - splitP[elemID]) / diffP
            eigenF = denseFourior[selElem, pcaID] * sqrt(diffP)
            efun = approxfun(eigenT, eigenF, rule = 2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}



#' The Fourier Basis defined on [0,2]
#' 
#' @param  obsT (vector) the observation time points inside [0,2]
#' @param fExpNum (int) the number of fourier bases
get_fourierBasisMat = function(obsT, fExpNum){
    ## orthonormal basis constructed from fourior basis
    basisMat = matrix(0, nrow = length(obsT), ncol = fExpNum*2)
    for(i in 1:fExpNum){
        basisMat[,2*i-1] = cos(i*obsT*pi) 
        basisMat[,2*i] = sin(i*obsT*pi) 
    }     
    ## t(denseFourior) %*% denseFourior * 2 * pi /1000 = Identity matrix
    return(basisMat)
}


# 
# ## Each row is one FPCA component
# ## Each column is one vector component
# par(mfrow = c(pcaCompNum, numComp), mar = rep(1,4))
# denseFourior = get_fourierBasisMat(tSeq, fExpNum)
# tmpM = numeric()
# for(i in 1:pcaCompNum){
#     tmp = numeric()
#     for(j in 1:numComp){
#         selR = 1:6 + 6*(j-1)
#         tmp = c(tmp, denseFourior %*% gVector[selR,i])
#         plot(tSeq, denseFourior %*% gVector[selR,i], type = "l")
#     }
#     tmpM = cbind(tmpM, tmp)
# }
# t(tmpM) %*% tmpM * 2*pi /length(tSeq)

# This function transfer some data set without elemID to dataframe with factorized elemID 
# this is only the univariate version
TransferData = function(obsCol){
    obsCol = as.data.frame(obsCol)
    obsCol = add_column(obsCol, "elemID" = rep(1, nrow(obsCol)), .after = "obsID")
    obsCol$elemID = factor(obsCol$elemID)
    return(obsCol)
}


####Generate data frame
### Generate our settings similar to Tony Cai
Generate_Cai = function(eigenfList, numSample, numPCA, noiseSigma = 1/4, alpha = 2){
   numElem = 1
   obsMat = matrix(0, numSample*obsNumUpper*numElem, 4)
   startI = 1
   for (obsID in 1:numSample){
       tmp = get_oneObs_Cai(eigenfList, noiseSigma, numPCA, obsNumLower, obsNumUpper, alpha)
       tmp = cbind(obsID, tmp)
       endI = startI + nrow(tmp) - 1
       obsMat[startI:endI,] = tmp
       startI = endI + 1
   }
   obsMat = data.frame(obsMat)
   colnames(obsMat) = c("obsID", "elemID", "obsT", "obsY")
   obsMat = obsMat[1:(endI-1),]
   obsMat$elemID = factor(obsMat$elem)
   ## Add a modification here, introducing a factor here
   obsMat$obsID = factor(obsMat$obsID)
   return(obsMat)
}

get_oneObs_Cai = function(eigenfList, noiseSigma = 1/4, numPCA, obsNumLower, obsNumUpper, alpha){
    score = rnorm(numPCA, 0, 1)
    result = numeric(0)
    elemObsNum = sample(obsNumLower:obsNumUpper, 1)
    obsT = runif(elemObsNum, 0, 1)
    obsY = rnorm(elemObsNum, sd = noiseSigma)
    weight = c(1:numPCA)
    zeta = (-1)^(weight + 1)*( weight^(-alpha) ) *2
    # Phi = matrix(1, elemObsNum, numPCA)
    #for (k in 1:numPCA){
        # Phi[, k] = sqrt(2)*cos(k*pi*obsT)
    #    phi[, k] 
    #}
    e = 1
    for(k in 1:numPCA){
        obsY = obsY + zeta[k]*score[k]*eigenfList[[k]][[e]](obsT)  # The noisesigma should be small to keep the signal noise ratio
    }
    
     # we are the univariate data
    tmp = cbind(e, obsT, obsY)
    result = rbind(result, tmp)
    return(result)
}

get_eigenfunList_cai = function(numPCA){
    tSeq = seq(0, 1, length.out = 30000)
    tSeq = tSeq[-c(1, length(tSeq))]
    numElem = 1
    result = list()
    for (pcaID in 1:numPCA){
        tmp = list()
        for (elemID in 1:numElem){
            eigenT = tSeq
            eigenF = sqrt(2)*cos(pcaID*pi*tSeq)
            efun = approxfun(eigenT, eigenF, rule = 2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}


#### This part of code can give the eigenfunction that can be exactly represented by B-spline basis
library(FDABasics)
get_eigenfunList_exact = function(order = 4, nknots, numPCA){
    tSeq = seq(0, 1, length.out = 30000)
    tSeq = tSeq[-c(1, length(tSeq))]
    splineObj = new(orthoSpline, tmin = 0, tmax = 1, order, nknots)
    K = splineObj$getDoF()
    Basis = t(splineObj$evalSpline(tSeq) )
    ## We need to define a QR decomposition here
    M = matrix(rnorm(K*numPCA), K, numPCA)
    Q = qr.Q(qr(M))
    
    
    temp = t(Basis)%*%Basis*(tSeq[2] - tSeq[1])
    R = t(chol(temp))
    
    if (ncol(Basis) != nrow(Q)){
        print("error in dimension between Basis and Q")
    } else {
        eigenFunc = t( solve(R, t(Basis)) )%*%Q
    }
    
    
    if (numPCA > K){
        print("K should be larger than r")
        break
    }
    numElem = 1
    result = list()
    for (pcaID in 1:numPCA){
        tmp = list()
        for (elemID in 1:numElem){
            eigenT  = tSeq
            eigenF = eigenFunc[, pcaID]
            efun = approxfun(eigenT, eigenF, rule = 2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}

Generate_pcaTrans = function(seedJ, pcaCompNum){
    old <- .Random.seed
    set.seed(seedJ)
    pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
    pcaTrans = qr.Q(qr(pcaTrans))
    .Random.seed <<- old
    return(pcaTrans)
}


Generate_knots = function(seedJ, IntKnots){
    old <- .Random.seed
    set.seed(seedJ)
    knots = sort(runif(IntKnots))
    .Random.seed <<- old
    return(knots)
}

Generate_rescale_U = function(seedJ, knots){
    tmin = 0
    tmax = 1
    mOrder = 4
    splineObj_generate = new(orthoSpline, tmin, tmax, mOrder, length(knots) + 1)
    Gamma = splineObj_generate$get_Omega()
    old <- .Random.seed
    set.seed(seedJ)
    A = rmvnorm(dim(Gamma)[1], mean = rep(0,dim(Gamma)[1]), sigma = solve(Gamma + 0.01*diag(dim(Gamma)[1])))
    A = t(A)
    U = qr.Q(qr(A))
    pcaTrans =  U[, sort(diag(t(U) %*% Gamma %*% U), index.return = TRUE)$ix]
    .Random.seed <<- old
    return(pcaTrans)
}


get_eigenfunList_NonUnif = function(df.M, pcaTrans, knots){
    lmin = 0
    lmax = 1
    tSeq = seq(0, 1, length.out = 30000)
    tSeq = tSeq[-c(1, length(tSeq))]
    delta<-(lmax-lmin)/(df.M-3)
    #knots<-c(seq(lmin,lmax/2,by=delta/2), seq(lmax/2+delta, lmax, by = delta))
    #Boundary.knots = c(lmax+delta,lmax+2*delta)
    #df = length(knots) + length(Boundary.knots)
    ## length(knots) = 18
    #grids = seq(0, 1, length.out = 500)
    Basis = bs(tSeq, knots =  knots, degree = 3, Boundary.knots = c(0, 1))
    temp = t(Basis)%*%Basis*(tSeq[2] - tSeq[1])
    R = t(chol(temp))
    Basis = t( solve(R, t(Basis)) )%*% pcaTrans
    
    eigenFunc = Basis
    #df = dim(Basis)[2]
    numPCA = dim(Basis)[2]
    numElem = 1
    result = list()
    for (pcaID in 1:numPCA){
        tmp = list()
        for (elemID in 1:numElem){
            eigenT = tSeq
            eigenF = eigenFunc[ ,pcaID]
            efun = approxfun(eigenT, eigenF, rule = 2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}

get_oneObs_NonUnif = function(eigenFList, pcaKappaSqrt, noiseSigma, obsNumLower, 
                              obsNumUpper, scoreType){
    numElem = length(eigenFList[[1]])
    numPCA = length(eigenFList)
    if (scoreType == "Gaussian"){
        score = rnorm(numPCA, sd = pcaKappaSqrt)
    } else if (scoreType == "uniform"){
        score = runif(numPCA, min = -sqrt(3), max = sqrt(3))
        score = score*pcaKappaSqrt
    } else if (scoreType == "t"){
        dfv = 3
        rtvar = dfv/(dfv-2)
        score = rt(numPCA, df = dfv)
        score = score/sqrt(rtvar)
        score = score * pcaKappaSqrt
    }
    
    result = numeric(0)
    for(e in 1:numElem){
        elemObsNum = sample(obsNumLower:obsNumUpper, 1)
        
        obsT = runif(elemObsNum, 0, 1)
        obsY = rnorm(elemObsNum, sd = noiseSigma)
        for(k in 1:numPCA){
            obsY = obsY + score[k] * eigenFList[[k]][[e]](obsT)
        }
        tmp = cbind(e, obsT, obsY)
        result = rbind(result, tmp)
    }
    
    return(result)
}

get_allObs_NonUnif = function(eigenFList, numSample, pcaKappaSqrt, 
                              noiseSigma, obsNumLower, obsNumUpper, scoreType){
    numElem = 1
    obsMat = matrix(0, numSample*obsNumUpper*numElem, 4)
    startI = 1
    for(obsID in 1:numSample){
        tmp = get_oneObs_NonUnif(eigenFList, pcaKappaSqrt, noiseSigma, 
                                 obsNumLower, obsNumUpper, scoreType)
        tmp = cbind(obsID, tmp)
        endI = startI + nrow(tmp) - 1
        #cat(endI, nrow(obsMat),"\n")
        obsMat[startI:endI,] =  tmp
        startI = endI + 1
    }
    obsMat = data.frame(obsMat)
    colnames(obsMat) = c("obsID", "elemID", "obsT", "obsY")
    obsMat = obsMat[1:(endI-1),]
    obsMat$elemID = factor(obsMat$elem)
    ## Add a modification here, introducing a factor here
    obsMat$obsID = factor(obsMat$obsID)
    return(obsMat)
}


get_eigenfunList_NonUnif = function(df.M, pcaTrans, knots){
    lmin = 0
    lmax = 1
    tSeq = seq(0, 1, length.out = 30000)
    tSeq = tSeq[-c(1, length(tSeq))]
    #delta<-(lmax-lmin)/(df.M-3)
    #knots<-c(seq(lmin,lmax/2,by=delta/2), seq(lmax/2+delta, lmax, by = delta))
    #Boundary.knots = c(lmax+delta,lmax+2*delta)
    #df = length(knots) + length(Boundary.knots)
    ## length(knots) = 18
    #grids = seq(0, 1, length.out = 500)
    #Boudary.knots = c(0, 1)
    Basis = bs(tSeq, knots =  knots, degree = 3)
    temp = t(Basis)%*%Basis*(tSeq[2] - tSeq[1])
    R = t(chol(temp))
    Basis = t( solve(R, t(Basis)) )%*% pcaTrans
    
    eigenFunc = Basis
    #df = dim(Basis)[2]
    numPCA = dim(Basis)[2]
    numElem = 1
    result = list()
    for (pcaID in 1:numPCA){
        tmp = list()
        for (elemID in 1:numElem){
            eigenT = tSeq
            eigenF = eigenFunc[ ,pcaID]
            efun = approxfun(eigenT, eigenF, rule = 2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}

get_oneObs_NonUnif = function(eigenFList, pcaKappaSqrt, noiseSigma, obsNumLower, 
                              obsNumUpper, scoreType){
    numElem = length(eigenFList[[1]])
    numPCA = length(eigenFList)
    if (scoreType == "Gaussian"){
        score = rnorm(numPCA, sd = pcaKappaSqrt)
    } else if (scoreType == "uniform"){
        score = runif(numPCA, min = -sqrt(3), max = sqrt(3))
        score = score*pcaKappaSqrt
    } else if (scoreType == "t"){
        dfv = 3
        rtvar = dfv/(dfv-2)
        score = rt(numPCA, df = dfv)
        score = score/sqrt(rtvar)
        score = score * pcaKappaSqrt
    }
    
    result = numeric(0)
    for(e in 1:numElem){
        elemObsNum = sample(obsNumLower:obsNumUpper, 1)
        
        obsT = runif(elemObsNum, 0, 1)
        obsY = rnorm(elemObsNum, sd = noiseSigma)
        for(k in 1:numPCA){
            obsY = obsY + score[k] * eigenFList[[k]][[e]](obsT)
        }
        tmp = cbind(e, obsT, obsY)
        result = rbind(result, tmp)
    }
    
    return(result)
}

get_allObs_NonUnif = function(eigenFList, numSample, pcaKappaSqrt, 
                              noiseSigma, obsNumLower, obsNumUpper, scoreType){
    numElem = 1
    obsMat = matrix(0, numSample*obsNumUpper*numElem, 4)
    startI = 1
    for(obsID in 1:numSample){
        tmp = get_oneObs_NonUnif(eigenFList, pcaKappaSqrt, noiseSigma, 
                                 obsNumLower, obsNumUpper, scoreType)
        tmp = cbind(obsID, tmp)
        endI = startI + nrow(tmp) - 1
        #cat(endI, nrow(obsMat),"\n")
        obsMat[startI:endI,] =  tmp
        startI = endI + 1
    }
    obsMat = data.frame(obsMat)
    colnames(obsMat) = c("obsID", "elemID", "obsT", "obsY")
    obsMat = obsMat[1:(endI-1),]
    obsMat$elemID = factor(obsMat$elem)
    ## Add a modification here, introducing a factor here
    obsMat$obsID = factor(obsMat$obsID)
    return(obsMat)
}


get_eigenfunList_Paul = function(Type, seedU = 7){
    grids= seq(0,1,0.002)
    numElem = 1
    result = list()
    if (Type == "easy"){
        data(easy)
        eigenF = easy$eigenfunctions
        numPCA = nrow(eigenF)
    } else if (Type == "prac"){
        data(prac)
        eigenF = prac$eigenfunctions
        numPCA = nrow(eigenF)
    } else if (Type == "hybrid"){
        data(prac)
        eigenF_prac = prac$eigenfunctions
        R.inv = R.inverse(M, grids)
        PhiMat = apply(matrix(grids), MARGIN = 1, Basis.Poly, M = M, R.inv = R.inv)
        Ut = eigenF_prac%*%t(PhiMat)/length(grids)
        U = t(Ut)
        U = qr.Q(qr(U))
        UComp = Ucomplete(seedU, U)
        Ubasis = cbind(U, UComp)
        eigenF = t(Ubasis)%*%PhiMat
        numPCA = nrow(eigenF)
    }
    
    for (pcaID in 1:numPCA){
        tmp = list()
        for (elemID in 1:numElem){
            eigenT = grids
            efun = approxfun(eigenT, eigenF[pcaID, ], rule=2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}


get_eigenfunList_Paul_NU = function(Type, shape1 = 2, shape2 = 2) {
    grids= seq(0,1,0.002)
    numElem = 1
    result = list()
    if (Type == "easyNU") {
        data(easy)
        eigenF = easy$eigenfunctions
        numPCA = nrow(eigenF)
    } else if (Type == "pracNU"){
        data(prac)
        eigenF = prac$eigenfunctions
        numPCA = nrow(eigenF)
    }
    
    NewGrids = qbeta(grids, shape1, shape2)
    tmpMat = c()
    for (pcaID in 1:numPCA) {
        efun_NU = approxfun(NewGrids, eigenF[pcaID, ], rule = 2)
        tmpMat = rbind(tmpMat, efun_NU(grids))
    }
    
    Norm = tmpMat%*%t(tmpMat)*(grids[2] - grids[1])
    R <- t(chol(Norm))
    BS = solve(R, tmpMat)
    
    
    for (pcaID in 1:numPCA) {
        tmp = list()
        for (elemID in 1:numElem) {
            efun = approxfun(grids, BS[pcaID, ], rule = 2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}

get_allObs_Paul = function(numSample, M, r, Type, ScoreType, alpha, nmin, nmax, a, b, sig, noiseType = "Gaussian", shape1 = 2, shape2 = 2){
    numElem = 1
    #obsMat = matrix(0, numSample*obsNumUpper*numElem, 4)
    obsMat = matrix(0, numSample*nmax*numElem, 4)
    startI = 1
    for (obsID in 1:numSample){
        tmp = get_oneObs_Paul(M, r, Type, ScoreType, alpha, nmin, nmax, a, b, sig, noiseType = noiseType, shape1, shape2)
        tmp = cbind(obsID, tmp)
        endI = startI + nrow(tmp) - 1
        obsMat[startI:endI, ] = tmp 
        startI = endI + 1
    }
    obsMat = data.frame(obsMat)
    colnames(obsMat) = c("obsID", "elemID", "obsT", "obsY")
    obsMat = obsMat[1:(endI-1),]
    obsMat$elemID = factor(obsMat$elem)
    ## Add a modification here, introducing a factor here
    obsMat$obsID = factor(obsMat$obsID)
    return(obsMat)
}





get_oneObs_Paul = function(M, r, Type, ScoreType, alpha = 0.6, nmin, nmax, a=1,
                           b=1, sig, seedU = 7, noiseType = "Gaussian", shape1 = 2, shape2 = 2){
    result = numeric(0)
    L.v<-MeasL(nmin,nmax,a,b)
    t.v<-L.v
    if (noiseType == "Gaussian") {
        e.v<-rnorm(length(L.v))*sig
    } else if (noiseType == "t") {
        e.v <- rt(length(L.v), df = 3)*sig/1.725
    } else if (noiseType == "uniform") {
        e.v <- runif(length(L.v), min = -sqrt(3), max = sqrt(3))*sig
    }
    
    
    e = 1
    result.all = RandomCurve_Paul(t.v, M, r, Type, ScoreType, alpha, seedU, shape1, shape2)
    
    obsY = result.all[[1]]
    #Ksi.all = result.all[[2]]
    
    obsY = obsY[1:length(L.v)] + e.v
    obsT = t.v
    tmp = cbind(e, obsT, obsY)
    result = rbind(result, tmp)
    return(result)
}


RandomCurve_Paul = function(t.v, M, r, Type, ScoreType, alpha, seedU = 7, shape1 = 2, shape2 = 2){
    grids= seq(0,1,0.002)
    if (Type == "easy"){
        data(easy)
        eigenF = easy$eigenfunctions
        FList = list()
        for (i in 1:nrow(eigenF)){
            efun = approxfun(grids, eigenF[i, ], rule = 2)
            FList = c(FList, list(efun))
        }
    } else if (Type == "prac"){
        data(prac)
        eigenF = prac$eigenfunctions
        FList = list()
        for (i in 1:nrow(eigenF)){
            efun = approxfun(grids, eigenF[i, ], rule = 2)
            FList = c(FList, list(efun))
        }
    } else if (Type == "hybrid"){
        data(prac)
        eigenF_prac = prac$eigenfunctions
        R.inv = R.inverse(M, grids)
        PhiMat = apply(matrix(grids), MARGIN = 1, Basis.Poly, M = M, R.inv = R.inv)
        Ut = eigenF_prac%*%t(PhiMat)/length(grids)
        U = t(Ut)
        U = qr.Q(qr(U))
        UComp = Ucomplete(seedU, U)
        Ubasis = cbind(U, UComp)
        
        eigenF = t(Ubasis)%*%PhiMat
        FList = list()
        for (i in 1:nrow(eigenF)){
            efun = approxfun(grids, eigenF[i, ], rule = 2)
            FList = c(FList, list(efun))
        }
    } else if (Type == "easyNU") {
        data(easy)
        eigenF = easy$eigenfunctions
        FList = list()
        NewGrids = qbeta(grids, shape1, shape2 )
        tmpMat = c()
        for (i in 1:nrow(eigenF)) {
            efun_NU = approxfun(NewGrids, eigenF[i, ], rule = 2)
            tmpMat = rbind(tmpMat, efun_NU(grids))
        }
        
        Norm = tmpMat%*%t(tmpMat)*(grids[2] - grids[1])
        R <- t(chol(Norm))
        BS = solve(R, tmpMat)
        
        for (i in 1:nrow(eigenF)){
            efun = approxfun(grids, BS[i, ], rule = 2)
            FList = c(FList, list(efun))
        }
    } else if (Type == "pracNU") {
        data(prac)
        eigenF = prac$eigenfunctions
        FList = list()
        NewGrids = qbeta(grids, shape1, shape2 )
        tmpMat = c()
        for (i in 1:nrow(eigenF)) {
            efun_NU = approxfun(NewGrids, eigenF[i, ], rule = 2)
            tmpMat = rbind(tmpMat, efun_NU(grids))
        }
        
        Norm = tmpMat%*%t(tmpMat)*(grids[2] - grids[1])
        R <- t(chol(Norm))
        BS = solve(R, tmpMat)
        
        for (i in 1:nrow(eigenF)){
            efun = approxfun(grids, BS[i, ], rule = 2)
            FList = c(FList, list(efun))
        }
    }
    
    result = numeric(length(t.v))
    eigen.f = NULL
    for (j in 1:length(FList)){
        eigen.f = rbind(eigen.f, FList[[j]](t.v))
    }
    
    if ( (Type == "easy")||(Type == "prac") || (Type == "easyNU") ||(Type == "pracNU")){
        eigen.v = matrix(Eigenv(r, method = "algebra", alpha, beta = 1, r1 = 3, r2 = 7, gamma = 2), r, length(t.v))
    } else if (Type == "hybrid"){
        eigen.v = matrix(Eigenv(r, method = "hybrid", alpha, beta = 1, r1 = 3, r2 = 7, gamma = 2), r1+r2, length(t.v))
    }
    
    if (ScoreType == "Gaussian"){
        ksi.v = rnorm(r)
    } else if (ScoreType == "t"){
        df = 3
        ksi.v = rt(r, df)/sqrt(df/(df-2))
    }
    ksi.m = matrix(ksi.v, r, length(t.v))
    
    temp = sqrt(eigen.v)*eigen.f*ksi.m
    result<-apply(temp,2,sum)
    
    return(list(result, FList))
}

Ucomplete = function(seedU, U){
    old <- .Random.seed
    set.seed(seedU)
    M = nrow(U)
    PU = U%*%t(U)
    X2 = matrix(rnorm(M*(M-5)), M, (M-5) )
    tmp = (diag(M) - PU) %*% X2
    UComp = qr.Q(qr(tmp))
    .Random.seed <<- old
    return(UComp)
}

get_eigenfunList_PaulNoSpline = function(Type, pcaTrans, R.inv = NULL, grid = NULL){
    grids= seq(0,1,0.002)
    numElem = 1
    result = list()
    M <-nrow(pcaTrans)
    r <- ncol(pcaTrans)
    if ( (Type == "easySin") || (Type == "pracSin")){
    # grid = NULL? R.inv = NULL
        eigenf.method = "sin"
        eigen.f<-apply(matrix(grids),MARGIN=1,Eigenf,basis.method=eigenf.method,B= pcaTrans,R.inv=R.inv,grid=grid)
        numPCA = nrow(eigen.f)
    } else if ( (Type == "easySpike")|| (Type == "pracSpike")){
        eigenf.method = "spike"  ## R.inv required!
        a<-c(50,70,160)
        b<-c(600,300,40)
        c<-c(1,1,10)  # As in Basis.Spike function 
        R.inv = Spike.orthM(a, b, c, grids)
        eigen.f<-apply(matrix(grids),MARGIN=1,Eigenf,basis.method=eigenf.method,B= pcaTrans,R.inv=R.inv,grid=grid)
    } 
    
    for (pcaID in 1:numPCA){
        tmp = list()
        for (elemID in 1:numElem){
            eigenT = grids
            efun = approxfun(eigenT, eigen.f[pcaID, ], rule=2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}


get_allObs_PaulNoSpline = function(numSample, M, r, Type, pcaTrans, ScoreType, alpha, nmin, nmax, a, b, sig, noiseType = "Gaussian"){
    numElem = 1
    #obsMat = matrix(0, numSample*obsNumUpper*numElem, 4)
    obsMat = matrix(0, numSample*nmax*numElem, 4)
    startI = 1
    for (obsID in 1:numSample){
        tmp = get_oneObs_PaulNoSpline(M, r, Type, pcaTrans, ScoreType, alpha, nmin, nmax, a, b, sig, noiseType = noiseType)
        tmp = cbind(obsID, tmp)
        endI = startI + nrow(tmp) - 1
        obsMat[startI:endI, ] = tmp 
        startI = endI + 1
    }
    obsMat = data.frame(obsMat)
    colnames(obsMat) = c("obsID", "elemID", "obsT", "obsY")
    obsMat = obsMat[1:(endI-1),]
    obsMat$elemID = factor(obsMat$elem)
    ## Add a modification here, introducing a factor here
    obsMat$obsID = factor(obsMat$obsID)
    return(obsMat)
}

# generate one sparsely observed realization with measurement error; 
# Not exactly represented by spline function 
get_oneObs_PaulNoSpline = function(M, r, Type, Q, ScoreType, alpha = 0.6, nmin, nmax, a = 1,
                                   b = 1, sig, noiseType = "Gaussian") {
    result = numeric(0)
    L.v<-MeasL(nmin,nmax,a,b)
    t.v<-L.v
    if (noiseType == "Gaussian") {
        e.v<-rnorm(length(L.v))*sig
    } else if (noiseType == "t") {
        e.v <- rt(length(L.v), df = 3)*sig/1.725
    } else if (noiseType == "uniform") {
        e.v <- runif(length(L.v), min = -sqrt(3), max = sqrt(3))*sig
    }
    
    if ( (dim(Q)[1] != M) || (dim(Q)[2] != r) ) {
        cat("dimension of Transform matrix is different!")
        break()
    }
    
    e = 1
    result.all = RandomCurve_PaulNoSpline(t.v, r, Type, Q, ScoreType, eigenv.method = "algebra", alpha)
    
    obsY = result.all
    obsY = obsY[1:length(L.v)] + e.v
    obsT = t.v
    tmp = cbind(e, obsT, obsY)
    result = rbind(result, tmp)
    return(result)
}

RandomCurve_PaulNoSpline = function(t.v, r, Type, B, ScoreType, eigenv.method="algebra", alpha, 
                                    beta=1, r1=3, r2=7, gamma=2,R.inv=NULL,grid=NULL) {
    ##t.v--vector of evaluation points 
    ##r--number of nonzero eigenvalues
    ##eigenf.method: the way to generate eigenfunctions; one of "sin", "poly",
    ##eigenv.method: the way to generate eigenvalues; one of "algebra", "geometric" and "hyperbolic"
    ##return: random curve evaluated at t.v: length(t.v) by 1; and pca scores: r by 1
    result<-numeric(length(t.v))
    if (Type == "easySin"){
        eigenf.method = "sin"
        eigen.f<-apply(matrix(t.v),MARGIN=1,Eigenf,basis.method=eigenf.method,B=B,R.inv=R.inv,grid=grid)
    } else if (Type == "easySpike") {
        eigenf.method = "spike"
        eigen.f<-apply(matrix(t.v),MARGIN=1,Eigenf,basis.method=eigenf.method,B=B,R.inv=R.inv,grid=grid)
    } else if (Type == "pracSin") {
        eigenf.method = "sin"
        eigen.f<-apply(matrix(t.v),MARGIN=1,Eigenf,basis.method=eigenf.method,B=B,R.inv=R.inv,grid=grid)
    } else if (Type == "pracSpike") {
        eigenf.method = "spike"
        eigen.f<-apply(matrix(t.v),MARGIN=1,Eigenf,basis.method=eigenf.method,B=B,R.inv=R.inv,grid=grid)
    }
   
    eigen.v<-matrix(Eigenv(r,eigenv.method,alpha,beta,r1,r2,gamma),r,length(t.v))
    
    if (ScoreType == "Gaussian"){
        ksi.v = rnorm(r)
    } else if (ScoreType == "t"){
        df = 3
        ksi.v = rt(r, df)/sqrt(df/(df-2))
    } else if (ScoreType == "uniform") {
        ksi.v = runif(r, min = -sqrt(3), max = sqrt(3))
    }
    
    ksi.m = matrix(ksi.v, r, length(t.v))
    temp = sqrt(eigen.v)*eigen.f*ksi.m
    result<-apply(temp,2,sum)
    return(result)
}


get_eigenfunList_NonEqual = function(df.M, pcaTrans, method){
    lmin = 0
    lmax = 1
    grids= seq(0,1,0.002)
    tSeq = grids[-c(1, length(grids))]
    delta<-(lmax-lmin)/(df.M-3)
    knots_old = seq(0, 1, delta)[-c(1,length(seq(0, 1, delta)))]
    if (method == "NonEqualSquare") {
        knots = knots_old^2
    } else if (method == "NonEqual3/2"){
        knots = knots_old^(3/2)
    } else if (method == "NonEqualRoot") {
        knots = knots_old^(1/2)
    }
    Basis = bs(tSeq, knots = knots, degree = 3, intercept = TRUE, Boundary.knots = c(0, 1))
    temp = t(Basis)%*%Basis*(grids[2] - grids[1])
    R = t(chol(temp))
    
    eigenfunctions = t( solve(R, t(Basis)) ) %*% pcaTrans
    
    numPCA = dim(eigenfunctions)[2]
    numElem = 1
    result = list()
    for (pcaID in 1:numPCA) {
        tmp = list()
        for (elemID in 1:numElem) {
            eigenT = tSeq
            eigenF = eigenfunctions[, pcaID]
            efun = approxfun(eigenT, eigenF, rule = 2)
            tmp = c(tmp, list(efun))
        }
        result = c(result, list(tmp))
    }
    return(result)
}


RandomCurve_PaulNonEqual = function(t.v, r, eigenFList, ScoreType, eigenv.method="algebra", alpha, 
                                    beta=1, r1=3, r2=7, gamma=2,R.inv=NULL,grid=NULL) {
    ##t.v--vector of evaluation points 
    ##r--number of nonzero eigenvalues
    ##eigenf.method: the way to generate eigenfunctions; one of "sin", "poly",
    ##eigenv.method: the way to generate eigenvalues; one of "algebra", "geometric" and "hyperbolic"
    ##return: random curve evaluated at t.v: length(t.v) by 1; and pca scores: r by 1
    result<-numeric(length(t.v))
    e = 1
    eigen.f = c()
    for (k in 1:r) {
        eigen.f = rbind(eigen.f, eigenFList[[k]][[e]](t.v))
    }
    
    
    eigen.v<-matrix(Eigenv(r,eigenv.method,alpha,beta,r1,r2,gamma),r,length(t.v))
    
    if (ScoreType == "Gaussian"){
        ksi.v = rnorm(r)
    } else if (ScoreType == "t"){
        df = 3
        ksi.v = rt(r, df)/sqrt(df/(df-2))
    } else if (ScoreType == "uniform") {
        ksi.v = runif(r, min = -sqrt(3), max = sqrt(3))
    }
    
    ksi.m = matrix(ksi.v, r, length(t.v))
    temp = sqrt(eigen.v)*eigen.f*ksi.m
    result<-apply(temp,2,sum)
    return(result)
}

get_oneObs_PaulNoEqual = function(eigenFList, ScoreType, alpha = 0.6, nmin, nmax, a = 1,
                                  b = 1, sig, noiseType = "Gaussian") {
    result = numeric(0)
    L.v<-MeasL(nmin,nmax,a,b)
    t.v<-L.v
    if (noiseType == "Gaussian") {
        e.v<-rnorm(length(L.v))*sig
    } else if (noiseType == "t") {
        e.v <- rt(length(L.v), df = 3)*sig/1.725
    } else if (noiseType == "uniform") {
        e.v <- runif(length(L.v), min = -sqrt(3), max = sqrt(3))*sig
    }
    
    e = 1
    result.all = RandomCurve_PaulNonEqual(t.v, r, eigenFList, ScoreType, 
                                          eigenv.method = "algebra", alpha)
    
    obsY = result.all
    obsY = obsY[1:length(L.v)] + e.v
    obsT = t.v
    tmp = cbind(e, obsT, obsY)
    result = rbind(result, tmp)
    return(result)
}

get_allObs_PaulNoEqual = function(numSample, M, r, eigenFList, pcaTrans, ScoreType, alpha, nmin, nmax,
                                  a, b, sig, noiseType = "Gaussian"){
    numElem = 1
    obsMat = matrix(0, numSample*nmax*numElem, 4)
    startI = 1
    for (obsID in 1:numSample) {
        tmp = get_oneObs_PaulNoEqual(eigenFList, ScoreType, alpha, nmin, nmax,
                                     a, b, sig, noiseType = noiseType)
        tmp = cbind(obsID, tmp)
        endI = startI + nrow(tmp) - 1
        obsMat[startI:endI, ] = tmp 
        startI = endI + 1
    }
    obsMat = data.frame(obsMat)
    colnames(obsMat) = c("obsID", "elemID", "obsT", "obsY")
    obsMat = obsMat[1:(endI-1),]
    obsMat$elemID = factor(obsMat$elem)
    ## Add a modification here, introducing a factor here
    obsMat$obsID = factor(obsMat$obsID)
    return(obsMat)
}

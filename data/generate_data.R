# Get the simulation data for the first settings


library(tibble)
#' Repeat get_oneObs to get multiple observation groups
#' 
#' @param numSample (int) total number of observation groups
get_allObs = function(eigenFList,
                      numSample, pcaKappaSqrt, noiseSigma, 
                      obsNumLower, obsNumUpper, scoreType){
    # pre allocation of memory to speed up code
    numElem = length(eigenFList[[1]])  # True, How many variate 
    obsMat = matrix(0, numSample*obsNumUpper*numElem, 4)
    startI = 1
    for(obsID in 1:numSample){
        tmp = get_oneObs(eigenFList, pcaKappaSqrt, noiseSigma, 
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


## Generate one group of observation
#' @param eigenFList (list) the list of list of eigenfunctions
#' @param pcaKappaSqrt FCPA scores
#' @param noiseSigma (double) the standard deviation of obs noise
#' @param obsNumLower (int) the lower bound of obs number
#' @param obsNumUpper (int) the upper bound of obs number
get_oneObs = function(eigenFList,
                      pcaKappaSqrt, noiseSigma, 
                      obsNumLower, obsNumUpper, scoreType){

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
        obsY = rnorm(elemObsNum, sd = noiseSigma)
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


get_allObs_Paul = function(numSample, M, r, Type, ScoreType, alpha, nmin, nmax, a, b, sig){
    numElem = 1
    obsMat = matrix(0, numSample*obsNumUpper*numElem, 4)
    startI = 1
    for (obsID in 1:numSample){
        tmp = get_oneObs_Paul(M, r, Type, ScoreType, alpha, nmin, nmax, a, b, sig)
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
                           b=1, sig, seedU = 7){
    result = numeric(0)
    L.v<-MeasL(nmin,nmax,a,b)
    t.v<-L.v
    e.v<-rnorm(length(L.v))*sig
    
    e = 1
    result.all = RandomCurve_Paul(t.v, M, r, Type, ScoreType, alpha, seedU)
    
    obsY = result.all[[1]]
    #Ksi.all = result.all[[2]]
    
    obsY = obsY[1:length(L.v)] + e.v
    obsT = t.v
    tmp = cbind(e, obsT, obsY)
    result = rbind(result, tmp)
    return(result)
}


RandomCurve_Paul = function(t.v, M, r, Type, ScoreType, alpha, seedU = 7){
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
    }
    
    result = numeric(length(t.v))
    eigen.f = NULL
    for (j in 1:length(FList)){
        eigen.f = rbind(eigen.f, FList[[j]](t.v))
    }
    
    if ( (Type == "easy")||(Type == "prac") ){
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

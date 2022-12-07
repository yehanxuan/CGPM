# Get the simulation data for the first settings
library(tibble)
library(mvtnorm)
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
        r1 = 3
        r2 = 7
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
    
    if (Type == "hybridSin") {
        eigenf.method = "sin"
        eigen.f <- apply(matrix(grids),MARGIN=1,Eigenf,basis.method=eigenf.method,B= pcaTrans,R.inv=R.inv,grid=grid)
        numPCA = nrow(eigen.f)
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
    if (Type == "hybridSin") {
        result.all = RandomCurve_PaulNoSpline(t.v, r, Type, Q, ScoreType, eigenv.method = "hybrid", alpha)
    } else {
        result.all = RandomCurve_PaulNoSpline(t.v, r, Type, Q, ScoreType, eigenv.method = "algebra", alpha)
    }
    
    
    obsY = result.all
    obsY = obsY[1:length(L.v)] + e.v
    obsT = t.v
    tmp = cbind(e, obsT, obsY)
    result = rbind(result, tmp)
    return(result)
}

RandomCurve_PaulNoSpline = function(t.v, r, Type, B, ScoreType, eigenv.method, alpha,
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
    } else if (Type == "hybridSin") {
        eigenf.method = "sin"
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










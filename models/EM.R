####functions for EM
####5-22-07
################################EM
### use the function "ns" in library(splines)

##formatting data : create data object suitable for implementation of EM 
library(splines)
Format.EM <- function(data.list,n,nmax,grid){
    data.f<-Format.data(data.list,n,nmax)
    Obs<-data.f[[1]]
    T<-data.f[[2]]
    N<-data.f[[3]]
    data<-TranMtoV(Obs,T,N)
    y<-data[,2]
    t<-data[,3]
    
    
    timeindex = floor(t*length(grid))+1
    result = list(y=y,curve=data[,1],n=n,timeindex=timeindex, t=t, N = N, T = T, Obs = Obs)  # add a t here (modified)
    return(result)
}


###EM
EM<-function(data.list,n,nmax,grids,M.EM,iter.num,r.EM,basis.EM,sig.EM, eigenfList = NULL, InitType = NULL){
    ## sig.EM is standard error
    
    # R.inv<-R.inverse(M.EM,grids)
    ## Change here 
    splineBasis = new(orthoSpline, 0, 1, 4, M.EM - 2)
    R.inv = Spline.inverse(splineBasis, grids)
    ## (a) formatting data for use in the EM routine
    data.obj <- Format.EM(data.list,n,nmax,grids) 
    ## (b) EM estimation of eigenfunctions (thetand PC scores (alpha)
    EMest <- fpcaEM(data.obj,k=r.EM,df=M.EM, grid = grids, maxit = iter.num, tol = 0.001, pert = 0.01, sigma = sig.EM^2,basis.EM, eigenfList, R.inv, InitType)
    
    ####Add one line to compute the mean function of EM
    B = EMest$B
    theta.zero = EMest$theta.zero
    meanEst = B %*% theta.zero
    
    ####
    if(basis.EM=="ns"){
        
        B.basis <- cbind(1, ns(grids, df = M.EM)) ### spline basis (with constant) used in EM
        B.orthobasis <- qr.Q(qr(B.basis)) ### orthonormalized spline basis used in EM
    }
    
    if (basis.EM == "BSpline"){
        splineBasis = new(orthoSpline, 0, 1, 3, M.EM)
        B.basis <- t(splineBasis$evalSpline(grids))
        B.orthobasis = qr.Q(qr(B.basis))
        
    }
    
    if(basis.EM=="poly"){
        
    #    lmin<-0
    #    lmax<-1 
    #    delta<-(lmax-lmin)/(M.EM-3)
        # M.EM degree of freedom 
    #    knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
    #    bs<-apply(matrix(grids),1, BC.orth,knots=knots,R.inv=R.inv)
    #    B.orthobasis<-t(bs)*sqrt(grids[2]-grids[1])
        
        ## How about use our basis ?
        splineBasis = new(orthoSpline, 0, 1, 4, M.EM - 2)
        B.basis = t(splineBasis$evalSpline(grids))
        B.orthobasis = B.basis
    } 
    
    EMeigenmat <- B.orthobasis %*% EMest$theta ### (unscaled) eigenfunctions in original time scale
    ## add 
    if (r.EM == 1){
        EMeigen.est <- EMeigenmat * as.numeric( 1/sqrt(t(EMeigenmat)%*%EMeigenmat ) )
    } else {
        EMeigen.est <- EMeigenmat%*% diag(1/sqrt(diag(t(EMeigenmat)%*%EMeigenmat) )) ###  normalization
    }
    ## (c) Crude ''estimated'' of eigenvalues from normalization
    EMDinv <- EMest$alpha[[3]]  ## estimate of D^{-1} from EM algorithm
    ##EMgamma <- EMeigenmat %*% solve(EMDinv) %*% t(EMeigenmat)
    EMgamma <- EMeigenmat %*% solve(EMDinv, t(EMeigenmat))
    ### estimate of Gamma
    
    ####EMgamma.svd <- svd(EMgamma)  ## svd of Gamma (* length(grid))
    ####EMeigenvec.est <- EMgamma.svd$u[,1:r.EM]
    ### first r.EM eigenvectors of Gamma
    #####EMeigenval.est <- EMgamma.svd$d[1:r.EM]/length(grids)
    
    
    EMgamma.svd <- eigen(EMgamma,symmetric=TRUE)     ###use eigen with symmetric=TRUE
    ### first r.EM eigenvectors of Gamma
    EMeigenval.est <- EMgamma.svd$values[1:r.EM]/length(grids)
    ### first r.EM eigenvectors of Gamma
    EMeigenvec.est <- EMgamma.svd$vectors[,1:r.EM]
    ### estimated eigenvalues
    EMsigma.est = sqrt(EMest$sigma)   ##estimated sigma
    
    step = EMest$step
    obj = EMest$obj
    Time = EMest$Time
    loss  =EMest$lossVec
    converge = EMest$converge 
    ##result
    result<-list(EMeigenvec.est,EMeigenval.est,EMsigma.est, meanEst, step, obj, Time, loss, converge)
    return(result)
    
}


calclike <- function(y, sigma, Dinv, theta, B, curve){  
    
    ### calculates the log likelihood of data	
    
    like <- 0	
    N <- length(table(curve))	
    for(i in 1:N){       
        X <- B[curve == i, ] %*% theta		
        C.my <- X %*% solve(Dinv, t(X)) + sigma * diag(dim(X)[1])
        C.my<-(C.my+t(C.my))/2 
        like <- like - t(y[curve == i]) %*% solve(C.my, y[curve == i])/2 - sum(logb(eigen(C.my,symmetric=TRUE, only.values=TRUE)$value))/2	
    }
    #     print(paste("Like =", like) 
    return(like)
}


fpcaEM <- function(obj, k = 2, df = 5, grid = seq(0.01, 1, length = 100), maxit = 50, tol = 0.001, pert = 0.01, sigma = 1, basis.method="ns", eigenfList = NULL, R.inv=NULL, InitType = NULL){
    
    ## computes the MLEs of alpha (PC score), theta (eigenfunctions represented in spline basis)
    
    y <- obj$y                 
    ### data represented as a single vector (by stacking the observations for different curves)
    
    timeindex <- obj$timeindex             
    ### a vector with the index of location (on the grid) of each point in y 
    
    curve <- obj$curve                   
    ### a vector with the curve number corresponding to each value of y 
    
    N <- obj$n     
    t = obj$t
    N_LOC = obj$N
    T_LOC = obj$T 
    Obs_LOC = obj$Obs
    data<-TranMtoV(Obs_LOC,T_LOC,N_LOC)
    ### number of observations (= n, in our notation)
    
    if(basis.method=="ns"){
        #  print("ns")
        B <- cbind(1, ns(grid, df = df))      
        ### spline basis evaluated at the grid locations, 
        #### df = d.f. for the splines (= M in our notation); grid = a vector of possible time points
        
        
        B.orth <- qr.Q(qr(B))  #### orthonormalizing the columns of B
        B <- B.orth[timeindex, ]  ### evaluating the basis at the observed times
    } 
    
    
    
    if (basis.method == "BSpline"){
        splineBasis = new(orthoSpline, 0, 1, 4, df)
        B <- t(splineBasis$evalSpline(grid))
        B.orth = qr.Q(qr(B)) 
        B <- B.orth[timeindex, ]
    }
    # Only use poly basis 
    basis.method = "poly"
    
    if(basis.method=="poly"){
        #  print("poly")
        ###R.inv<-R.inverse(df)
     #   lmin<-0
     #   lmax<-1 
     #   delta<-(lmax-lmin)/(df-3)
     #   knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
     #   bs<-apply(matrix(grid),1, BC.orth,knots=knots,R.inv=R.inv)
     #   bs<-t(bs)*sqrt(grid[2]-grid[1])
     #   B<-bs[timeindex,]
        splineBasis = new(orthoSpline, 0, 1, 4, df-2)
        B.orth <- t(splineBasis$evalSpline(grid))
        B <- B.orth[timeindex, ]
    } 
    
    R.old <- 0    
    #R.new <- 1  
    R.new = 1e10
    #ind <- 1  
    ind <- 0
    theta.zero <- solve(t(B) %*% B, t(B)) %*% y 
    y <- y - B %*% theta.zero 
    
    alpha <- list(alpha = matrix(1, N, k), alphaprod = array(1, c(N, k, k)))   
    ### k = dimension of alpha (= r, in our notation)
    
    theta <- as.matrix(init(y, timeindex, curve, B, k, pert)$theta)  
    alpha <- getemalpha(y, curve, theta, B, alpha)  
    sigma <- getsigma(y, curve, B, theta, alpha, sigma)
    # Add another sentence here
    alpha <- getemalpha(y, curve, theta, B, alpha, sigma = sigma) 
    ## Add likelihood here, we want to monitor the likelihood
    
    if ( isTRUE(InitType == "LOC") ){
        grid.l=seq(0,1,0.01)
        hmu.cv<-h.select(t,y,method="cv")
        fitmu<-sm.regression(t,y,h=hmu.cv,poly.index=1,eval.points=t)$estimate 
        conver<-HatGOff(fitmu,data,N_LOC)  ##convert to off diagonal pair forma
        Cova<-conver[[2]]
        CovT<-conver[[3]]
        hcov.cv<-h.select(t(CovT),Cova,method="cv")
        hsig.cv = h.select(t,(y-fitmu)^2,method = "cv")
        lin.result<-LocLinEst.new(Obs_LOC,T_LOC,N_LOC,grid.l,hmu.cv,hcov=hcov.cv,hsig.cv)
        muhat<-lin.result[[1]]
        covmatrix<-lin.result[[2]]
        sig2hat<-lin.result[[4]]/(1-0)
        ###project on a finer grid
        timeindex_LOC<-floor(grids*length(grid.l))+1
        timeindex_LOC[length(grids)]<-length(grid.l)
        covmatrix<-covmatrix[timeindex_LOC,timeindex_LOC]
        
        ##get eigenfunctions and eigenvalues 
        eigen.c<-EigenC_LOC(covmatrix,grids)
        eigenf.c<-eigen.c[[1]][,1:k]
        eigenv.c<-eigen.c[[2]][1:k]
        
        eigenfest_LOC = t(eigenf.c)
        
        tmp1 = matrix(rep(sqrt(eigenv.c), N), nrow = N, byrow=TRUE)
        prod = array(1, c(N, k, k))
        for (i in 1:N){
          prod[i, ,] = sqrt(eigenv.c)%*%t(sqrt(eigenv.c))
        }
        alpha <- list(alpha = tmp1, alphaprod = prod) 
    
        theta = t(eigenfest_LOC%*%B.orth/dim(B.orth)[1])
        alpha <- getemalpha(y, curve, theta, B, alpha) 
        sigma = sig2hat ## square here 
        
          
    }
    
    loss = list()
    if (!is.null(eigenfList)){
      EMeigenvec = Generate_EM_eigen(basis.method = "poly", alpha, theta, k, df, grid, R.inv)
      loss[[1]] = computeLoss_EM(EMeigenvec, grid, eigenfList)
    }
    
    #loss = computeLoss_EM(EMeigenvec, grid, eigenfList)
    
    objF = c( -2 * calclike(y, sigma, alpha$Dinv, theta, B, curve)/N )
   # Time = c(0,0,0)
    Time = 0
    #while(abs(R.old - R.new)/R.new > tol & (ind < maxit)){  
    while(abs(R.old - R.new)> tol & (ind < maxit)){    
        start = as.vector(proc.time()[1])
        ind <- ind + 1
       # sigma <- getsigma(y, curve, B, theta, alpha, sigma)   
        theta <- getemtheta(y, curve, alpha, B, theta) 
        alpha <- getemalpha(y, curve, theta, B, alpha, sigma = sigma) 
        sigma <- getsigma(y, curve, B, theta, alpha, sigma) 
        
        #### If we want to make the coparison fair, we should make change 
        
        
        ##### Original code 
        R.old <- R.new 
        #R.new <- sum((y - (B %*% theta * alpha$alpha[curve, ]) %*% rep(1, k))^2)
        #         print(R.new) 
        #R.new
        
        calclike(y, sigma, alpha$Dinv, theta, B, curve) 
        EMeigenvec = Generate_EM_eigen(basis.method = "poly",alpha, theta, k, df, grid, R.inv)
        #loss = c( loss, computeLoss_EM(EMeigenvec, grid, eigenfList))
        if (!is.null(eigenfList)){
          loss[[ind+1]] = computeLoss_EM(EMeigenvec, grid, eigenfList)
        }
        
        
        Interval = as.vector(proc.time()[1]) - start
        
        R.new =  -2 * calclike(y, sigma, alpha$Dinv, theta, B, curve)/N 
        objF = c(objF, -2 * calclike(y, sigma, alpha$Dinv, theta, B, curve)/N )
       # Time = rbind(Time, Interval[1:3])
        #Time = rbind(Time, Interval)
        Time = c(Time, Interval)
        tol.c = abs(R.old - R.new)
    }  
    # Time = as.matrix(Time)
    # TimeMatrix = apply(Time, 2, cumsum)
    TimeMatrix = cumsum(Time)
    step = ind + 1
    temp <- svd(theta)  
    ### add iterations number and likelihood 
    if (tol.c < tol) {
      converge = 1
    } else {
      converge = 0
    }
    result = list(alpha = alpha, theta = theta, B = B, theta.zero = theta.zero, sigma = sigma, step = step, obj = objF, Time = TimeMatrix, lossVec = loss, converge = converge)
    return(result)
}


getemalpha <-function(y, curve, theta, B, alphaobj, sigma = 1){ 
    
    ### computes EM update for alpha (E-step of EM)
    
    alpha <- alphaobj$alpha  
    alphaprod <- alphaobj$alphaprod
    n <- length(table(curve))
    N <- dim(alpha)[1]  
    if(dim(alpha)[2] > 1) {
        alphasq <- apply(alphaprod, 1, diag)
        Dinv <- N * diag(as.vector(1/(alphasq %*% rep(1, N))))  
    }  
    else 
        Dinv <- as.matrix(1/mean(alphaprod))  
    for(i in 1:n) {  
        X <- B[curve == i, ] %*% theta 
        Calpha <- solve(sigma * Dinv + t(X) %*% X)  
        
        alpha[i, ] <- Calpha %*% t(X) %*% y[curve == i]              
        ###   hat{alpha_i} (eq. (27))
        
        alphaprod[i,  ,  ] <- sigma * Calpha + alpha[i,  ] %*% t(alpha[i,  ])      
        ### hat{alpha_i alapha_i^T}  (eq. (28))
    }  
    result = list(alpha = alpha, alphaprod = alphaprod, Dinv = Dinv)
    return(result)
}


getemtheta <- function(y, curve, alphaobj, B, theta, tol = 0.0001){ 
    
    ### computes EM update for theta (M step)
    
    q <- dim(B)[2]
    R.old <- 1
    R.new <- 0
    alpha <- alphaobj$alpha
    alphaprod <- alphaobj$alphaprod
    
    k <- dim(alpha)[2]              
    ### = r, in our notation
    
    while(abs(R.old - R.new)/R.new > tol) {	
        for(j in 1:k) {
            ind <- rep(1, k)   
            ind[j] <- 0       
            tempy <- alpha[curve, j] * y - ((B %*% theta) * alphaprod[curve,  ,  j]) %*% ind   
            tempX <- B * sqrt(alphaprod[curve, j, j])      
            theta[, j] <- solve(t(tempX) %*% tempX, t(B)) %*% tempy		   
        }    
        R.old <- R.new  
        R.new <- sum((y - ((B %*% theta) * alpha[curve,  ]) %*% rep(1, k))^2)
    }	    
    return(theta)
}


getsigma <- function(y, curve, B, theta, alpha, sigma){	
    
    ### EM update for sigma^2 (M step)
    
    tempsigma <- 0  
    Dinv <- alpha$Dinv   
    N <- dim(alpha$alpha)[1]
    fit <- NULL
    for(i in 1:N){   
        X <- B[curve == i,  ] %*% theta  
        Calpha <- solve(Dinv + t(X) %*% X/sigma)  
        fit <- c(fit, B[curve == i,  ] %*% theta %*% alpha$alpha[i,])
        tempsigma <- tempsigma + sum(diag(X %*% Calpha %*% t(X)))
    }	
    sigma <- (sum((y - fit)^2) + tempsigma)/length(y)
    #     print(paste("sigma =", sigma))	
    return(sigma)
}


init <- function(y, timeindex, curve, B, k, pert = 0){     
    
    ### initial values of theta and gamma
    
    tab <- table(curve)  
    N <- length(tab) 
    s <- c(0, cumsum(tab)) 
    q <- dim(B)[2]   
    gamma <- matrix(0, N, q) 
    for(i in 1:N) {
        X <- B[(s[i] + 1):s[i + 1],  ]  
        tempy <- y[(s[i] + 1):s[i + 1]]     
        gamma[i,  ] <- solve(t(X) %*% X + pert * diag(q), t(X)) %*% tempy
    }	
    theta <- prcomp(gamma)$rotation[, 1:k]
    result = list(theta = theta, gamma = gamma)
    return(result)
}


EM_selection = function(newObsCol, M.set, r.set, basis.EM="poly", nFold = 10, eigenfList = NULL, InitType = NULL, cvMembership = NULL, grids= seq(0,1,0.002)){
    
    nmax<-max(table(newObsCol[,1]))  
    L1<-min(as.numeric(newObsCol[,3]))              ##range of the time 
    L2<-max(as.numeric(newObsCol[,3]))
    data.list = fpca.format(newObsCol)
    n<-length(data.list)    ##number of subjects
    if (is.null(cvMembership)){
      cvMembership = getCVPartition(n, nFold)
    }
    
    # max( unlist( lapply(data.list, function(x) dim(x[[1]])[1] ) ) )
    
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
    
    result  = NULL
    
    for (k in 1:length(r.set)){
        r.c = r.set[k]
        M.set = M.set[M.set >= r.c]
        M.l = length(M.set)
        if(M.l==0){              ##if all M in the set M.set are less than r, return zero  
            print("all M < r.c")
            return(0)
        }
        
        eigen.result<-array(0,dim=c(1,r.c ,M.l))   ##estimated eigenvalues for three methods(Newton,ini, ini/EM) under different M
        eigenf.result<-NULL                     ##estimated eigenfunctions 
        sig.result<-matrix(0,1,M.l)             ##estimated error sd. for three methods under differen M
        like.result<-matrix(0,1,M.l)
        
        names.m <- c("EM")
        
        ErrorVec = rep(1e10, M.l)
        for (i in 1:M.l){
            M.EM = M.set[i]
            testerror = rep(1e7, nFold)
            for (cf in 1:nFold){
              try({
                trainIndex = c(cvMembership != cf)
                testIndex = c(cvMembership == cf)
                train.data.list = data.list.new[trainIndex]
                test.data.list = data.list.new[testIndex]
                ntest = length(test.data.list)
                ### Do we need to modify nmax here?? 
                nmax_train = max( unlist( lapply(train.data.list, function(x) dim(x[[1]])[1] ) ) )
                temp.EM = EM(train.data.list, length(train.data.list), nmax_train, grids, M.EM, iter.num = 50, r.c , basis.EM = basis.EM, sig.EM, eigenfList, InitType)
                EMeigenvec.est<-temp.EM[[1]]*sqrt(length(grids))
                EMeigenval.est<-temp.EM[[2]]
                EMsigma.est<-temp.EM[[3]]
                if (r.c == 1){
                  covmatrix.EM = EMeigenval.est*EMeigenvec.est%*%t(EMeigenvec.est)
                } else {
                  covmatrix.EM<-EMeigenvec.est%*%diag(EMeigenval.est)%*%t(EMeigenvec.est)
                }
                test.like.EM<-loglike.cov.all(covmatrix.EM,EMsigma.est, test.data.list, ntest) ### mean(-2*loglikelihood) in optimization.R
                testerror[cf] = test.like.EM
              })
            }
            ErrorVec[i] = mean(testerror)
            #temp = Initial(r, ini.method = "EM", data.list.new, n, nmax, grid.l, grids, M.EM, basis.EM = basis.EM)
            # sig2hat<-EMsigma.est^2
            #  covmatrix.ini<-covmatrix.EM
            #  eigenf.ini<-EMeigenvec.est
            #  eigenv.ini<-EMeigenval.est
            #  like.ini<-like.EM
            #  ### -2*loglike
            #  like.result[, i] = like.EM
        }
        result[[k]] = ErrorVec
    }
    
    temp = EM.CV(result, M.set, r.set)
    index.r = temp[1]
    index.M = temp[2]
    r_opt = r.set[index.r]
    M_opt = M.set[index.M]
    #index = which.min(ErrorVec)
    #M_opt = M.set[index]
    return(list("M_opt" = M_opt, "r_opt" = r_opt, "like_result" = ErrorVec  ) )
}


EM.CV = function(result, M.set, r.set){
    cv.mod = matrix(1e+10, length(r.set), length(M.set))
    colnames(cv.mod) = M.set
    rownames(cv.mod) = r.set
    
    for (k in 1:length(r.set)){
        r.c = r.set[k]
        M.set.c = M.set[M.set>=r.c]
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
            if (cv.mod[k, j] < cv.mod[index.r, index.M] ){
                index.r = k
                index.M = j
            }
        }
    }
    
    temp = c(index.r, index.M)
    names(temp) = c("r", "M")
    return(temp)
}


MFPCA_EM = function(obsCol, M.set, r.set, sig.EM, splineObj = NULL, eigenfList = NULL,InitType = NULL){
    tmin = 0
    tmax = 1
    newObsCol = obsCol[, -2]
    newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
    colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
    
    basis.method = "bs"
    basis.EM = "poly"
    sl.v = rep(0.5, 10)
    max.step = 50
    grid.l = seq(0,1,0.01)
    grids= seq(0,1,0.002)
    # sig.EM = 1/4
    nmax<-max(table(newObsCol[,1]))
    L1<-min(as.numeric(newObsCol[,3]))              ##range of the time
    L2<-max(as.numeric(newObsCol[,3]))
    data.list = fpca.format(newObsCol)
    n<-length(data.list)    ##number of subjects
    if(n==0){
        print("error: no subject has more than one measurements!") #/Users/hanxuan/Documents/mFPCA/mFPCA/R/MultivariateFpcaTraining.R
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
   #select = EM_selection(newObsCol, M.set, r.set, basis.EM = "poly", nFold = 10)
   select = EM_selection(newObsCol, M.set, r.set, basis.EM = "poly", nFold = 10, eigenfList, InitType)
   IniVal = Initial(select$r_opt, ini.method = "EM", data.list.new, n, nmax, grid.l, grids,
                    M.EM = select$M_opt, iter.num = 50, basis.EM = "poly", sig.EM = sig.EM, eigenfList, InitType)
   sig2hat = IniVal[[1]]
   covmatrix.ini<-IniVal[[2]]
   eigenf.ini<-IniVal[[3]]
   eigenv.ini<-IniVal[[4]]
   like.ini<-IniVal[[5]]
   #    grids.new = res$workGrid
   step  = IniVal[[8]]
   obj = IniVal[[9]]
   
   eigenfest = t(eigenf.ini)
   
   if (!is.null(splineObj)){
     K = splineObj$getDoF()
     basisMat = splineObj$evalSpline(grids)
     UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
     WInit = eigenv.ini
     XInit = list(UInit, diag(WInit))
   } else {
     XInit = NULL
   }
   
   # eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
   fList_est = Generate_fList(grids, eigenfest)
   if (is.null(eigenfList)){
       loss = NULL
   } else {
       loss = ComputeLoss(fList_est, eigenfList)
   }
   
   r_opt = select$r_opt
   M_opt = select$M_opt
   
   model = list(eigenFunctions = fList_est)
   model = c(model, list( tmin = tmin, tmax = tmax, SFinal = XInit,
   sigmaSq = sig2hat, numPCA = r_opt, numElem = 1, elemLevels = "1"
   ))
   
   if (!is.null(eigenfList)){
     loss = evaluateLoss(model, eigenfList)
   }
   
   
   return(list("opt_knots" = M_opt, "opt_rank" = r_opt, "loss" = loss, "model" = model, "step" = step, "obj" = obj))
}


oneReplicate_EM = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("./oneReplicate/oneRep-EM.R")
    cvMembership = getCVPartition_seed(samplesize, nFold = 10, seedJ)
    select = EM_selection(newObsCol, M.set, r.set, basis.EM = "poly", nFold = 10, eigenfList, InitType, cvMembership)
    IniVal = Initial(select$r_opt, ini.method ="EM", data.list.new, n, nmax, grid.l, grids,M.EM = select$M_opt, iter.num = 50,basis.EM = "poly",
                    sig.EM = sig.EM, eigenfList, InitType)
    
    
    sig2hat = IniVal[[1]]
    covmatrix.ini<-IniVal[[2]]
    eigenf.ini<-IniVal[[3]]
    eigenv.ini<-IniVal[[4]]
    like.ini<-IniVal[[5]]
    #    grids.new = res$workGrid
    step = IniVal[[8]]
    obj = IniVal[[9]]
    converge = IniVal[[12]]
    eigenfest = t(eigenf.ini)
    # eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    fList_est = Generate_fList(grids, eigenfest)
    loss = ComputeLoss(fList_est, eigenfList)
    r_opt = select$r_opt
    M_opt = select$M_opt
    return(list("loss" = loss, "result2" = M_opt, "rank" = r_opt, "converge" = converge))
}


oneReplicateWrap_EM = function(seedJ){
    try({
        eval = oneReplicate_EM(seedJ)
    })
    return(eval)
}


EMInit = function(obsCol, splineObj, nKnots, r, sigmaSq, eigenfList = NULL){
    nmax = max(table(obsCol[ ,1]))
    L1 = min(as.numeric(obsCol[ ,3]))
    L2 = max(as.numeric(obsCol[ ,3]))
    tmp = (obsCol[ ,3] - L1)/(L2 - L1)
    tmp[tmp<0.00001]<-0.00001
    tmp[tmp>0.99999]<-0.99999 
    obsCol[, 3] = tmp
    ## Data processing for EM algorithm
    newObsCol = obsCol[, -2]
    newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
    colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
    
    ini.method = "EM"
    basis.method = "bs"
    basis.EM = "poly"
    sl.v = rep(0.5, 10)
    max.step = 50
    grid.l = seq(0,1,0.01)
    grids= seq(0,1,0.002)
    iter.num = 50
    
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
    #select = EM_selection(newObsCol, M.set,r.set, basis.EM = "poly")
    M.EM = nKnots+2
    splineBasis = new(orthoSpline, 0, 1, 4, M.EM - 2)
    R.inv = Spline.inverse(splineBasis, grids)
    data.obj <- Format.EM(data.list.new, n,nmax,grids) 
    EMest <- fpcaEM(data.obj,k=r,df=M.EM, grid = grids, maxit = iter.num, tol = 0.001, pert = 0.01, sigma = sigmaSq,basis.EM, eigenfList, R.inv)
  #  IniVal = Initial(r, ini.method ="EM", data.list.new, n, nmax, grid.l, grids, M.EM, iter.num = 50, basis.EM = "poly",
  #                    sig.EM = sqrt(sigmaSq), eigenfList)
  
    UInit = EMest$theta
    Dinv = EMest$alpha[[3]]
    WInit = solve(Dinv)
    sig2hat = EMest$sigma
  #  sig2hat = IniVal[[1]]
  #  covmatrix.ini<-IniVal[[2]]
  #  eigenf.ini<-IniVal[[3]]
  #  eigenv.ini<-IniVal[[4]]
  #  like.ini<-IniVal[[5]]
    
  #  eigenfest = t(eigenf.ini)
  #  basisMat = splineObj$evalSpline(grids)
  #  UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
  #  WInit = eigenv.ini
    if (r == 1){
        XInit = list(UInit, as.matrix(WInit))
    } else {
    #XInit = list(UInit, diag(WInit))
      XInit = list(UInit, WInit)
    }
    return(list(XInit, sig2hat))
}




Generate_EM_eigen = function(basis.method,alpha, theta, k, df, grids, R.inv){
    if(basis.method=="poly"){
        
     #   lmin<-0
     #   lmax<-1
     #   delta<-(lmax-lmin)/(df-3)
        # M.EM degree of freedom
      #  knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
      #  bs<-apply(matrix(grids),1, BC.orth,knots=knots,R.inv=R.inv)
      #  B.orthobasis<-t(bs)*sqrt(grids[2]-grids[1])
        splineBasis = new(orthoSpline, 0, 1, 4, df - 2)
        B.basis = t(splineBasis$evalSpline(grids))
        B.orthobasis = B.basis
        
    }
    
    EMeigenmat <- B.orthobasis %*% theta
    
    
    
    # EMeigen.est = B.orthobasis %*% theta
    EMDinv = alpha[[3]]
    EMgamma <- EMeigenmat %*% solve(EMDinv, t(EMeigenmat))
    EMgamma.svd <- eigen(EMgamma,symmetric=TRUE)
    EMeigenval.est <- EMgamma.svd$values[1:k]/length(grids)
    ### first r.EM eigenvectors of Gamma
    EMeigenvec.est <- EMgamma.svd$vectors[,1:k]
    ### estimated eigenvalues
    #EMsigma.est = sqrt(sigma)   ##estimated sigma
    
    return(EMeigenvec.est)
}

computeLoss_EM = function(eigenvec.est, grids, eigenfList){
    eigenf = eigenvec.est * sqrt(length(grids))
    eigenfest = t(eigenf)
    fList_est = Generate_fList(grids, eigenfest)
#    loss = mean( ComputeLoss(fList_est, eigenfList) )
    ## Output loss of each component 
    loss = ComputeLoss(fList_est, eigenfList)
    return(loss)
}




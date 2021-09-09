##################################
######02-06-08: R package "fpca"
######R functions for MLE of FPCA to longitudinal data
##name: fpca.mle
##prupose: use Newton's method to fit MLE for FPCA, as well as do model selection by 
##         the approximate CV score

fpca.mle<-function(data.m, M.set,r.set,ini.method="EM", basis.method="bs",sl.v=rep(0.5,10),max.step=50,
                   grid.l=seq(0,1,0.01),grids=seq(0,1,0.002), eigenfList = NULL, CVmethod = "score"){
    ##para:data.m: the data matrix with three columns: column 1: subject ID, column 2: observation, column 3: measurement time
    ##M.set--# of basis functions; r.set: dimension of the process
    ##ini.method: initial method for Newton, one of "EM","loc"; 
    ##basis.method: basis functions to use for Newton, one of "bs" (cubic Bsplines), "ns" (nature splines)
    ##sl.v: shrinkage steps for Newton, e.g., sl.v<-c(0.5,0.5,0.5) means the first three steps are 0.5
    ##max.step: max number of iterations of Newton
    ##grid.l: grids for loc; grids: denser grids for EM/Newton;
    ##return: a list of (i) the selected model; (ii) the corresponding eigenfunctions; (iii)eigenvalues
    ##        (iv) error variance; (v) fitted mean (by local linear); (vi) the grid where the evaluation takes place
    
    ##(i)set some parameters in fpca.fit
    tol<-1e-3                                        #tolerance level to determine convergence in Newton 
    cond.tol<-1e+10                                  #tolerance to determine singularity in Newton 
    if(length(sl.v)<max.step){
        sl.v<-c(sl.v,rep(1,max.step-length(sl.v)))
    }else{
        sl.v<-sl.v[1:max.step]
    }
    
    if(basis.method=="bs"){
        basis.method<-"poly"
    }
    
    M.self<-20
    basis.EM<-"poly"                                  ##basis for EM
    
    if(ini.method=="EM"){
        no.EM<-FALSE
    }else{
        no.EM<-TRUE
    }
    
    ##initial value of sig in EM 
    
    
    nmax<-max(table(data.m[,1]))   ##maximum number of measurements per subject  
    L1<-min(as.numeric(data.m[,3]))              ##range of the time 
    L2<-max(as.numeric(data.m[,3]))
    
    ##(ii) format to data.list
    data.list<-fpca.format(data.m)
    n<-length(data.list)    ##number of subjects
    if(n==0){ 
        print("error: no subject has more than one measurements!")
        return(0)
    }
    
    ##(iii)rescale time onto [0,1]
    data.list.new<-data.list
    for (i in 1:n){
        cur<-data.list[[i]][[1]]
        temp<-(cur[,2]-L1)/(L2-L1)
        temp[temp<0.00001]<-0.00001
        temp[temp>0.99999]<-0.99999 
        cur[,2]<-temp
        data.list.new[[i]][[1]]<-cur
    }
    
    ## (iv)estimate mean and substract
    #library(sm) 
    
    ## 去掉LocLin.mean 试一试
  #  temp<-LocLin.mean(data.list.new,n,nmax, grids)  ##subtract estimated mean
  #  data.list.new<-temp[[1]]
  #  fitmu<-temp[[2]]
    
    ## (v)apply fpca.fit on all combinations of M and r in M.set and r.set and rescale back
    result<-NULL
    
    iter.EM.num = 50 ### some bugs there 
    
    
    #### Add grids.new before fpca.fit
    grids.new<-grids*(L2-L1)+L1
    for (k in 1:length(r.set)){
        r.c<-r.set[k]
        print(paste("r=",r.c))
        result.c<-fpca.fit(M.set,r.c,data.list.new,n,nmax,grid.l,grids,iter.EM.num,sig.EM,ini.method,basis.method,
                           sl.v,max.step,tol,cond.tol,M.self,no.EM,basis.EM, eigenfList, L1, L2)
        
        if(is.vector(result.c[[1]])){
            print("warning: see warning code")
        } 
        
        ##rescale back
        grids.new<-grids*(L2-L1)+L1
        
        temp<-result.c$eigenvec
        M.set.u<-M.set[M.set>=r.c]   ##only when M>=r, there is a result
        if(length(M.set.u)==0){
            temp<-NULL
        }else{
            for(j in 1:length(M.set.u)){   
                temp[[j]]<-temp[[j]]/sqrt(L2-L1)
            }
        }
        
        result.c$eigenvec<-temp
        result.c$eigenval<-result.c$eigenval*(L2-L1)
        result.c$grids<-grids.new  
        
        
        result[[k]]<-result.c
    }  
    ##(vi) model selection
    if (CVmethod == "like"){
      mod.sele = REML_CV(result, M.set, r.set)
    } else if (CVmethod == "score"){
      mod.sele<-fpca.cv(result,M.set,r.set,tol=tol)
    }
    
    temp<-mod.sele[[3]]
    index.r<-temp[1] 
    index.M<-temp[2]
    cv.result<-mod.sele[[1]]
    con.result<-mod.sele[[2]]
    
    ##(vii) return the selected model
    result.sele<-result[[index.r]]
    eigenf<-result.sele$eigenvec
    eigenf.sele<-eigenf[[index.M]][,,1]
    
    eigenv<-result.sele$eigenval
    eigenv.sele<-eigenv[1,,index.M]
    
    sig<-result.sele$sig
    sig.sele<-sig[1,index.M]
    
    ## Add step here
    step = result.sele$step[index.M]
    like = result.sele$like[[index.M]]
    Time = result.sele$Time[[index.M]]
    loss = result.sele$loss[[index.M]]
    
    temp.model<-c(M.set[index.M],r.set[index.r])
    names(temp.model)<-c("M","r")
    if (r.set == 1){
        eigenf.sele = as.matrix(t(eigenf.sele))
    }
    rownames(eigenf.sele)<-paste("eigenfunction",1:r.set[index.r])
    names(eigenv.sele)<-paste("eigenvalue",1:r.set[index.r])
    #"fitted_mean"=fitmu
    temp<-list("selected_model"=temp.model,"eigenfunctions"=eigenf.sele,"eigenvalues"=eigenv.sele,"error_var"=sig.sele^2,"grid"=grids.new,"cv_scores"=cv.result,"converge"=con.result, "step" = step, "like" = like,
               "Time" = Time, "loss" = loss)
    
    return(temp)  
}


fpca.fit<-function(M.set,r,data.list,n,nmax,grid.l=seq(0,1,0.01),grids=seq(0,1,0.002),
                   iter.num=50,sig.EM=1,ini.method, basis.method="poly",sl.v,max.step=50,tol=1e-3,cond.tol=1e+10,
                   M.self=20,no.EM=FALSE,basis.EM="poly", eigenfList = NULL, L1, L2){
    
    ##M.set--# of basis functions; r: dimension of the process
    ##data.list: data of the fpca.format format 
    ##n: sample size; nmax: maximum number of measurements per subject;
    ##grid.l: grids for loc; grids: denser grids for EM and Newton;
    ##iter.num: max number of iteration of EM; sig.EM: initial values of sig (erros sd.) for EM
    ##ini.method: initial method for Newton, one of "EM","loc","EM.self" (meaning using a fixed M.self); 
    ##basis.method: basis functions to use for Newton, one of "poly" (cubic Bspline), "ns" (nature spline)
    ##sl.v: shrinkage steps for Newton, e.g., sl.v<-c(0.5,0.5,0.5) means the first three steps are 0.5
    ##max.step: max number of iterations of Newton
    ##tol: tolerance level to determine convergence in Newton in terms of the l_2 norm of the gradient 
    ##cond.tol: tolerance to determine singularity in Newton in terms of condition number of the Hessian matrix
    ##M.self: dimension to use for EM.self; 
    ##no.EM: do EM (F) or not (T);
    ##basis.EM: basis to use for EM, one of "poly" or "ns"
    grid.news = grids*(L2-L1)+L1
    
    M.set<-M.set[M.set>=r]    ## only fit Newton for those M that is >= r.  
    M.l<-length(M.set)
    if(M.l==0){              ##if all M in the set M.set are less than r, return zero  
        print("all M<r")
        return(0)
    }
    
    ##results for return
    eigen.result<-array(0,dim=c(3,r,M.l))   ##estimated eigenvalues for three methods(Newton,ini, ini/EM) under different M
    eigenf.result<-NULL                     ##estimated eigenfunctions 
    sig.result<-matrix(0,3,M.l)             ##estimated error sd. for three methods under differen M
    like.result<-matrix(0,3,M.l)            ##likelihood for three methods under different M  
    cv.result <-numeric(M.l)                ##approximate CV score for Newton under different M
    converge.result<-numeric(M.l)           ##whether Newton converges under different M
    step.result<-numeric(M.l)               ## number of iterations for Newton to converge 
    cond.result<-numeric(M.l)               ## the condition number of the final Hessian in Newton
    
    like.list = list()    ## the likelihood list for each step 
    Time.list = list()
    loss.list = list()
    if(!no.EM){                             ## if do EM 
        names.m<-c("newton","ini","EM")  
    }else{
        names.m<-c("newton","ini",ini.method)   ## if do not do EM
    }
    
    dimnames(eigen.result)<-list(names.m,NULL,NULL)
    rownames(sig.result)<-names.m
    rownames(like.result)<-names.m
    
    ###initial estimates by loc  
    if(ini.method=="loc"){
        IniVal<-try(Initial(r,ini.method,data.list,n,nmax,grid.l,grids))
        if(inherits(IniVal, "try-error")){
            print("warning: error in initial step by loc")
            return(-2)
        } 
        
        sig2hat<-IniVal[[1]]
        covmatrix.ini<-IniVal[[2]] 
        eigenf.ini<-IniVal[[3]]
        eigenv.ini<-IniVal[[4]]
        like.ini<-IniVal[[5]]
    }
    
    ###initial values from EM with M=M.self
    if(ini.method=="EM.self"){ 
        IniVal<-try(Initial(r,ini.method="EM",data.list,n,nmax,grid.l,grids,basis.EM=basis.EM,M.EM=M.self,sig.EM=sig.EM))
        if(inherits(IniVal, "try-error")){
            print("warning: error in initial step by EM.self")
            return(-3)
        }
        
        sig2hat<-IniVal[[1]]
        covmatrix.ini<-IniVal[[2]] 
        eigenf.ini<-IniVal[[3]]
        eigenv.ini<-IniVal[[4]]
        like.ini<-IniVal[[5]]
    }
    
    for(i in 1:M.l){
        #print(paste("M=", M.set[i]))
        ###################### EM
        ##(i)
        if(!no.EM){
            M.EM<-M.set[i]      ##number of basis use
            temp.EM<-try(EM(data.list,n,nmax,grids,M.EM,iter.num,r,basis.EM=basis.EM,sig.EM, eigenfList))
            if(inherits(temp.EM, "try-error")){
                print("warning: error in initial step by EM")
                return(-1)
            }
            
            EMeigenvec.est<-temp.EM[[1]]*sqrt(length(grids))
            EMeigenval.est<-temp.EM[[2]]
            EMsigma.est<-temp.EM[[3]]
            
            if (r == 1){
                covmatrix.EM <- EMeigenval.est * EMeigenvec.est%*%t(EMeigenvec.est)
            } else {
                covmatrix.EM <- EMeigenvec.est%*%diag(EMeigenval.est)%*%t(EMeigenvec.est)
            }
            
            ##(iii)likelihood of EM result
            like.EM<-loglike.cov.all(covmatrix.EM,EMsigma.est,data.list,n)
        } else if (ini.method == "LS"){  ###  I add it by myself
            M.EM  = M.set[i]
            temp.LS <- try(REML_LS(data.list, n, nmax, grids, M.EM, iter.num, r, basis.EM=basis.EM, sig.EM))
            if(inherits(temp.LS, "try-error")){
                print("warning:error in initial step by LS")
                return(-1)
            }
            
            EMeigenvec.est <- temp.LS[[1]]*sqrt(length(grids))
            EMeigenval.est <- temp.LS[[2]]
            EMsigma.est <- temp.LS[[3]]
            
            if (r == 1){
                covmatrix.EM <- EMeigenval.est * EMeigenvec.est%*%t(EMeigenvec.est)
            } else {
                covmatrix.EM <- EMeigenvec.est%*%diag(EMeigenval.est)%*%t(EMeigenvec.est)
            }
            like.EM = loglike.cov.all(covmatrix.EM, EMsigma.est, data.list, n)
        }else{ 
            ##use other initial value in the place of EM
            EMeigenvec.est<-eigenf.ini
            EMeigenval.est<-eigenv.ini
            EMsigma.est<-sqrt(sig2hat)
            like.EM<-like.ini 
        }
        
        
        
        ##########################newton's method
        #####(0): Initial value for Newton's method
        
        if(ini.method=="EM"){
            sig2hat<-EMsigma.est^2
            covmatrix.ini<-covmatrix.EM
            eigenf.ini<-EMeigenvec.est
            eigenv.ini<-EMeigenval.est
            like.ini<-like.EM
        } else if (ini.method == "LS"){
            sig2hat<-EMsigma.est^2
            covmatrix.ini<-covmatrix.EM
            eigenf.ini<-EMeigenvec.est
            eigenv.ini<-EMeigenval.est
            like.ini<-like.EM
        }
        
        ### compute the first step loss
        fList_est = Generate_fList(grids, t(eigenf.ini))
       # loss.ini = mean(ComputeLoss(fList_est, eigenfList))
        loss.ini = list()
        if (!is.null(eigenfList)){
          loss.ini[[1]] = ComputeLoss(fList_est, eigenfList)
        }
        
        ####(i) projection on basis functions to get inital values of B and Lambda
        ##
        M<-M.set[i]   ##number of basis to use
        temp<-try(Proj.Ini(covmatrix.ini,M,r,basis.method,grids))
        if(inherits(temp, "try-error")){
            print("warning: error in Newton step due to projection")
            return(-4)
        }
        B.ini<-temp[[1]]
        Lambda.ini<-temp[[2]]        
        if(any(Lambda.ini<0)) {  ##should >0
            print("warning: error in Newton step: initial covariance not p.d.!")
            return(-4)
        }
        
        if(sig2hat>0){
            sig.ini<-sqrt(sig2hat)
        }else{
            sig.ini<-sqrt(Lambda.ini[r])
        }
        
        #### (ii)get auxilary things
        if(basis.method=="poly"){ ##cubic B-spline with equal spaced knots
            R.inv<-R.inverse(M,grids)
       #     splineBasis = new(orthoSpline, 0, 1, 4, M - 2)
       #     R.inv = Spline.inverse(splineBasis, grids)
            phi.aux<-apply(matrix(1:n),MARGIN=1,Phi.aux,M=M,basis.method=basis.method,data.list=data.list,R.inv=R.inv)
      #      phi.aux = apply(matrix(1:n), MARGIN=1, Spline.aux, basis.method = basis.method, data.list = data.list, splineBasis = splineBasis, R.inv = R.inv)
        }
        
        if(basis.method=="ns"){ ##nature spline
            bs.u<-NS.orth(grids,df=M)/sqrt(grids[2]-grids[1])
            phi.aux<-apply(matrix(1:n),MARGIN=1,Phi.aux,M=M,basis.method=basis.method,data.list=data.list,R.inv=bs.u, grid=grids)
        }
        
        
        ### Add one row here
        if (r == 1){
            B.ini = as.matrix(B.ini)
        }
        
        #### (iv) Newton iterations
        newton.result<-Newton.New(B.ini,phi.aux,sig.ini,Lambda.ini,data.list,n,sl.v,max.step,tol,cond.tol, eigenfList, R.inv, L1, L2)
        step<-newton.result[[7]]
        like<-newton.result[[1]][1:step]
        #like = newton.result[[1]][1:(step+1)]
        B.up<-newton.result[[2]][,,step]
        Lambda.up<-newton.result[[3]][step,]
        sig.up<-newton.result[[4]][step]
        gradB<-newton.result[[5]][,,step]
        gradL<-newton.result[[6]][step,] 
        cond<-newton.result[[8]]
        error.c<-newton.result[[9]]
        Time = newton.result[[10]]
        loss = newton.result[[11]]
       # loss = c(loss.ini, loss)
      #  loss = list(loss.ini, loss)
        loss = c(loss.ini, loss)
        ##check convergence
        # converge<-max(abs(gradB),abs(gradL))          ##converged? if zero
        ## Maybe we need to modify here, 09/09/2021
        converge = newton.result[[12]] # tol.c 
        if (r == 1){
            B.up = as.matrix(B.up)
        }
        #### (v) cross-validation score
        if(error.c){
            print("warning: cv not computed due to errors in the Newton iteration steps")
            cv.score<-(-99)
        }else{
            cv.score<-try(CV(B.up,phi.aux,sig.up,Lambda.up,data.list,n))
            # How do they compute the cv.score
            error.c<-inherits(cv.score, "try-error")
            if(error.c){
                print("warning: cv not compuated due to errors in the CV step")
                cv.score<-(-99)
            }
        }
        
        
        
        ## (vi)estimated eigenfunctions/values 
        ##eigenfunctions
        eigenv.up<-Lambda.up[order(-Lambda.up)]
        
        ## I add an extra drop = F here
        B.f<-B.up[,order(-Lambda.up), drop = F]
        if(basis.method=="poly"){
            eigenf.up<-t(apply(matrix(grids),1,Eigenf,basis.method=basis.method,B=B.f,R.inv=R.inv))
          #  splineBasis = new(orthoSpline, 0, 1, 4, M - 2)
           # eigenf.up = t(apply( matrix(grids), 1, EigenSpline, basis.method = basis.method, B = B.f, splineBasis = splineBasis, R.inv = R.inv))
        }
        if(basis.method=="ns"){
            eigenf.up<-t(apply(matrix(grids),1,Eigenf,basis.method=basis.method,B=B.f,R.inv=bs.u,grid=grids))
        }
        
        ##
        dim.c<-c(r,length(grids),3)
        eigen.est<-array(0,dim.c)
        eigen.est[,,1]<-t(eigenf.up)        ##Newton
        eigen.est[,,2]<-t(eigenf.ini)       ##Initial
        eigen.est[,,3]<-t(EMeigenvec.est)   ##EM/initial   
        dimnames(eigen.est)<-list(NULL,NULL,names.m)
        eigenf.result[[i]]<-eigen.est
        
        ##eigenvalues, sig, like, cv, converge, step, cond
        eigen.temp<-rbind(eigenv.up,eigenv.ini,EMeigenval.est)  
        eigen.result[,,i]<-eigen.temp
        sig.result[,i]<-c(sig.up,sig.ini,EMsigma.est)
        like.result[,i]<-c(like[step],like.ini,like.EM)
        cv.result[i]<-cv.score
        converge.result[i]<-converge
        #step.result[i]<-step
        step.result[i]<-step + 1 ##Add one step here 
        cond.result[i]<-cond
        like.list[[i]] = c(like.ini, like) 
        Time.list[[i]] = Time
        loss.list[[i]] = loss
    }
    
    ####(vii) return result
    result<-list("eigenval"=eigen.result,"sig"=sig.result,"-2*loglike"=like.result,
                 "cv"=cv.result,"converge"=converge.result,"step"=step.result,"condition#"=cond.result,
                 "eigenvec"=eigenf.result, "like" = like.list, "Time" = Time.list, "loss" = loss.list)
    return(result)
}


##name: fpca.cv
##purpose: model selection based on the results of a set of models
fpca.cv<-function(result.new,M.set,r.set,tol=1e-3){
    cv.mod<-matrix(-99,length(r.set),length(M.set))
    colnames(cv.mod)<-M.set
    rownames(cv.mod)<-r.set
    
    cond.mod<-cv.mod
    cv.sele<-cv.mod+1e+10
    
    for (k in 1:length(r.set)){
        r.c<-r.set[k]
        M.set.c<-M.set[M.set>=r.c]
        index.c<-sum(M.set<r.c)
        
        if(length(M.set.c)>0){
            for(j in 1:length(M.set.c)){
                cv.mod[k,j+index.c]<-result.new[[k]]$cv[j]
                cond.mod[k,j+index.c]<-result.new[[k]]$converge[j]
                if(cv.mod[k,j+index.c]!=(-99)&&cond.mod[k,j+index.c]<tol) 
                    cv.sele[k,j+index.c]<-cv.mod[k,j+index.c] 
            }
        } 
        
    }
    
    index.r<-1
    index.M<-1
    
    for (j in 1:length(M.set)){
        for(k in 1:length(r.set)){
            if(cv.sele[k,j]<cv.sele[index.r,index.M]){
                index.r<-k
                index.M<-j
            }
        }
    } 
    
    ##
    temp<-c(index.r,index.M)
    names(temp)<-c("r","M")
    result<-list("cv"=cv.mod,"converge"=cond.mod,"selected model"=temp)   
    return(result)
}


REML_CV = function(result, M.set, r.set, tol=1e-3){
  cv.mod<-matrix(-99,length(r.set),length(M.set))
  colnames(cv.mod)<-M.set
  rownames(cv.mod)<-r.set
  cond.mod <- cv.mod
  
  
  cv.like = cv.mod + 1e+10
  
  for (k in 1:length(r.set)){
    r.c = r.set[k]
    M.set.c = M.set[M.set >= r.c]
    index.c = sum(M.set<r.c)
    
    if (length(M.set.c)>0){
      for (j in 1:length(M.set.c)){
        cv.mod[k, j+index.c] = result[[k]]$cv[j]
        cond.mod[k, j+index.c] = result[[k]]$converge[j]
        ## Modify, we do not care about the CV ,09/09/2021
        #cv.mod[k, j+index.c]!=(-99) &&
        if (  cond.mod[k, j+index.c] < tol){
          cv.like[k, j+index.c] <- result[[k]]$`-2*loglike`["newton",][j]
        }
      }
    }
  }
  
  index.r = 1
  index.M = 1
  for (j in 1:length(M.set)){
    for (k in 1:length(r.set)){
      if (cv.like[k,j] < cv.like[index.r, index.M]){
        index.r = k
        index.M = j
      }
    }
  }
  
  temp <- c(index.r, index.M)
  names(temp)<-c("r","M")
  selection = list("cv" = cv.mod, "converge" = cond.mod , "select model" = temp)
  return(selection)
}


###name: fpca.format
##purpose: format the data into data.list as the input of fpca.fit and also exclude subjects with only one measurement 
fpca.format<-function(data.m){
    ##para: data.m: the data matrix with three columns: column 1: subject ID, column 2: observation, column 3: measurement time
    ##return: a list of n components where n is the number of subjects
    ##The ith component is a matrix of two columns: the first column is the measurements of the ith subject,
    ##and the second column is the corresponding times of measurements of the ith subject
    
    ID<-unique(data.m[,1]) 
    n<-length(ID)
    data.list<-NULL
    
    count<-1
    for(i in 1:n){
        N.c<-sum(data.m[,1]==ID[i]) 
        cur<-data.m[data.m[,1]==ID[i],]
        
        if(N.c>1){
            Obs.c<-as.numeric(cur[,2])
            T.c<-as.numeric(cur[,3])
            temp<-matrix(cbind(Obs.c,T.c),nrow=N.c, ncol=2)
            data.list[[count]]<-list(temp,NULL) 
            count<-count+1
        }
        
    }
    
    return(data.list)
    
}


MFPCA_REML = function(obsCol, M.set, r.set, ini.method, sig.EM = 1, splineObj = NULL, eigenfList = NULL, CVmethod = "like"){
    tmin = 0
    tmax = 1
    newObsCol = obsCol[, -2]
    newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
    colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]
    
    basis.method = "bs"
    sl.v = rep(0.5, 10)
    max.step = 50
    grid.l = seq(0,1,0.01)
    grids= seq(0,1,0.002)
    #sig.EM = 1
    
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
    
    result = fpca.mle(newObsCol, M.set, r.set, ini.method, basis.method, sl.v, max.step, grid.l,
             grids, eigenfList, CVmethod)
    grids.new = result$grid
    M_opt = result$selected_model[[1]]
    r_opt = result$selected_model[[2]]
    eigenfest = result$eigenfunctions
    
    if (!is.null(splineObj)){
      K = splineObj$getDoF()
      basisMat = splineObj$evalSpline(grids.new)
      UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
      WInit = result$eigenvalues
      XInit = list(UInit, diag(WInit))
    } else {
      XInit = NULL
    }
    
    
    fList_est = Generate_fList(grids.new, eigenfest)
    sigmaSq = result$error_var
    # True function interpolation    # eigenfList = get_eigenfunList(pcaTrans, fExpNum
    model = list(eigenFunctions = fList_est)
    model = c(model, list(tmin = tmin, tmax = tmax, SFinal = XInit, sigmaSq = sigmaSq,
    numPCA = r_opt, numElem = 1, elemLevels = "1" ))
    if (is.null(eigenfList)){
        loss = NULL
    } else {
        loss = evaluateLoss(model, eigenfList)
    }
    step = result$step
    return(list("opt_knots" = M_opt, "opt_rank" = r_opt, "loss" = loss, "model" = model, "step" = step))
}


REML_selection = function(newObsCol, M.set, r.set, ini.method, nFold=10, sig.EM=1, eigenfList = NULL,cvMembership = NULL){
  basis.method = "bs"
  sl.v = rep(0.5, 10)
  max.step = 50
  grid.l = seq(0,1,0.01)
  grids= seq(0,1,0.002)
  #sig.EM = 1
  
  L1 = length(r.set)
  L2 = length(M.set)
  ErrorMat = matrix(1e+10, L1, L2)
  for (i in 1:L1){
    for (j in 1:L2){
      testerror = rep(1e+7, nFold)
      for (cf in 1:nFold){
        try({
          trainIndex = which(cf == cvMembership)
          test_newObsCol =  newObsCol[newObsCol$obsID %in% trainIndex, ]
          train_newObsCol = newObsCol[!(newObsCol$obsID %in% trainIndex), ] 
          #cvParam = list(cvMembership = cvMembership, cf = cf)
          train_result = fpca.mle(train_newObsCol, M.set[j], r.set[i], ini.method, basis.method, sl.v, max.step, grid.l,
                                  grids, eigenfList)
          eigenf.train = train_result$eigenfunctions
          eigenv.train = train_result$eigenvalues
          covTrain = t(eigenf.train)%*%diag(eigenv.train)%*%eigenf.train
          train_result$error_var
          test_data.list = fpca.format(test_newObsCol)
          ntest = length(test_data.list)
          testerror[cf] = loglike.cov.all(covTrain, sqrt(train_result$error_var), test_data.list, ntest)
        })
      }
      ErrorMat[i, j] = mean(testerror)
    }
  }
  index = which(ErrorMat== min(ErrorMat, na.rm = TRUE), arr.ind = TRUE)
  index1 = index[1]
  index2 = index[2]
  r_opt = r.set[index1]
  M_opt = M.set[index2]
  opt_error = ErrorMat[index1, index2]
  return(list("ErrorMat" = ErrorMat, "M_opt" = M_opt, "r_opt" = r_opt, "opt_error" = opt_error))
}


oneReplicate_REML = function(seedJ){
    set.seed(seedJ + repID * 300)
    source("./oneReplicate/oneRep-REML.R")
   # cvMembership = getCVPartition_seed(n, nFold = 10, seedJ)
  #  select = REML_selection(newObsCol, M.set,r.set,ini.method,nFold,sig.EM,eigenfList, cvMembership)
  #  M_opt = select$M_opt
  #  r_opt = select$r_opt
    
    result = fpca.mle(newObsCol, M.set, r.set, ini.method, basis.method, sl.v, max.step, grid.l, grids, eigenfList, CVmethod = "like")
  #  result = fpca.mle(newObsCol, M_opt, r_opt,ini.method, basis.method, sl.v, max.step,
  #                   grid.l, grids, eigenfList, CVmethod = "score")
    grids.new = result$grid
    
    M_opt = result$selected_model[[1]]
    r_opt = result$selected_model[[2]]
    
    conv = result$converge[as.character(r_opt), as.character(M_opt)]
    if (conv < 1e-3){
      converge = 1
    } else {
      converge = 0 
    }
    eigenfest = result$eigenfunctions
    # True function interpolation
    # eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    fList_est = Generate_fList(grids.new, eigenfest)
    loss = ComputeLoss(fList_est, eigenfList)
    return(list( "loss" = loss, "M_opt" = M_opt, "rank" = r_opt, "converge" = converge))
}


oneReplicateWrap_REML = function(seedJ){
    try({
        eval = oneReplicate_REML(seedJ) 
    })
    return(eval)
}


fpcaLS = function(obj, k, df, grid = seq(0.01, 1, length = 100), maxit = 50, tol = 0.001, pert = 0.01,
                  sigma, basis.method, R.inv=NULL){
    ## computes the MLEs of alpha (PC score), theta (eigenfunctions represented in spline basis)
    
    y <- obj$y                 
    ### data represented as a single vector (by stacking the observations for different curves)
    
    timeindex <- obj$timeindex             
    ### a vector with the index of location (on the grid) of each point in y 
    
    curve <- obj$curve                   
    ### a vector with the curve number corresponding to each value of y 
    
    N <- obj$n                           
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
        splineBasis = new(orthoSpline, 0, 1, 3, df)
        B <- t(splineBasis$evalSpline(grid))
        B.orth = qr.Q(qr(B)) 
        B <- B.orth[timeindex, ]
    }
    if(basis.method=="poly"){
        #  print("poly")
        ###R.inv<-R.inverse(df)

        #lmin<-0
        #lmax<-1
        #delta<-(lmax-lmin)/(df-3)
        #knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
        #bs<-apply(matrix(grid),1, BC.orth,knots=knots,R.inv=R.inv)
        #bs<-t(bs)*sqrt(grid[2]-grid[1])
        #B<-bs[timeindex,]
        splineBasis = new(orthoSpline, 0, 1, 4, df - 2)
        B.orth <- t(splineBasis$evalSpline(grid))
        B <- B.orth[timeindex, ]
    } 
    
    theta.zero <- solve(t(B) %*% B, t(B)) %*% y 
    y <- y - B %*% theta.zero 
    
    alpha <- list(alpha = matrix(1, N, k), alphaprod = array(1, c(N, k, k)))   
    ### k = dimension of alpha (= r, in our notation)
    theta <- as.matrix(init(y, timeindex, curve, B, k, pert)$theta)  
    alpha <- getemalpha(y, curve, theta, B, alpha)  
    sigma <- getsigma(y, curve, B, theta, alpha, sigma)
    alpha <- getemalpha(y, curve, theta, B, alpha, sigma = sigma)
    
    result = list(alpha = alpha, theta = theta, B = B, theta.zero = theta.zero, sigma = sigma)
    return(result)
}


REML_LS = function(data.list, n, nmax, grids, M.EM, iter.num, r.EM, basis.EM, sig.EM){
    ## sig.EM is standard error
    
    ## Change to our spline Basis
    #
    # R.inv<-R.inverse(M.EM,grids)
    splineBasis  = new(orthoSpline, 0, 1, 4, M.EM - 2)
    R.inv = Spline.inverse(splineBasis, grids)
    ## (a) formatting data for use in the EM routine
    data.obj <- Format.EM(data.list,n,nmax,grids) 
    ## (b) EM estimation of eigenfunctions (theta) and PC scores (alpha)
    LSest = fpcaLS(data.obj, k = r.EM, df = M.EM, grid = grids, maxit = iter.num, tol = 0.001, pert = 0.01,
                   sigma = sig.EM^2, basis.EM, R.inv)
    B = LSest$B
    theta.zero = LSest$theta.zero
    meanEst = B %*% theta.zero
    
    if(basis.EM=="ns"){
        
        B.basis <- cbind(1, ns(grids, df = M.EM)) ### spline basis (with constant) used in EM
        B.orthobasis <- qr.Q(qr(B.basis)) ### orthonormalized spline basis used in EM
    }
    
    if(basis.EM=="poly"){
        
        #lmin<-0
        #lmax<-1
        #delta<-(lmax-lmin)/(M.EM-3)
        # M.EM degree of freedom 
        #knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
        #bs<-apply(matrix(grids),1, BC.orth,knots=knots,R.inv=R.inv)
        #B.orthobasis<-t(bs)*sqrt(grids[2]-grids[1])
        splineBasis = new(orthoSpline, 0, 1, 4, M.EM - 2)
        B.basis <- t(splineBasis$evalSpline(grids))
        B.orthobasis = B.basis
        
    } 
    EMeigenmat <- B.orthobasis %*% LSest$theta ### (unscaled) eigenfunctions in original time scale
    if (r.EM == 1){
        EMeigen.est <- EMeigenmat * as.numeric( 1/sqrt(t(EMeigenmat)%*%EMeigenmat ) )
    } else {
        EMeigen.est <- EMeigenmat%*% diag(1/sqrt(diag(t(EMeigenmat)%*%EMeigenmat) )) ###  normalization
    }
    ## (c) Crude ''estimated'' of eigenvalues from normalization
    EMDinv <- LSest$alpha[[3]]  ## estimate of D^{-1} from EM algorithm
    ##EMgamma <- EMeigenmat %*% solve(EMDinv) %*% t(EMeigenmat)
    EMgamma <- EMeigenmat %*% solve(EMDinv, t(EMeigenmat))
    EMgamma.svd <- eigen(EMgamma,symmetric=TRUE)     ###use eigen with symmetric=TRUE
    ### first r.EM eigenvectors of Gamma
    EMeigenval.est <- EMgamma.svd$values[1:r.EM]/length(grids)
    ### first r.EM eigenvectors of Gamma
    EMeigenvec.est <- EMgamma.svd$vectors[,1:r.EM]
    ### estimated eigenvalues
    EMsigma.est = sqrt(LSest$sigma)   ##estimated sigma
    ## result
    result = list(EMeigenvec.est, EMeigenval.est, EMsigma.est, meanEst)
    return(result)
}

require("sm")
require("splines")


Generate_REML_eigen = function(basis.method, grids, B.up, Lambda.up, R.inv, L1, L2){
    
    # eigenv.up<-Lambda.up[order(-Lambda.up)]
    #if (r == 1){
    #    B.up = as.matrix(B.up)
    #}
    
    B.f<-B.up[,order(-Lambda.up), drop = F]
    if(basis.method=="poly"){
       # M = nrow(B.f)
        eigenf.up<-t(apply(matrix(grids),1,Eigenf,basis.method=basis.method,B=B.f,R.inv=R.inv))
       # splineBasis = new(orthoSpline, 0, 1, 4, M - 2)
      #  eigenf.up = t(apply( matrix(grids), 1, EigenSpline, basis.method = basis.method, B = B.f, splineBasis = splineBasis, R.inv = R.inv))
    }
    eigen.est = t(eigenf.up)  #eigenvec[[index.M]][,,1] eigenfunctions
    eigenfunctions = t(eigenf.up)
    eigenfunctions = eigenfunctions/sqrt(L2 - L1)
    return(eigenfunctions)
}

computeLoss_REML = function(eigenfunctions, eigenfList, grids.new){
    fList_est_REML = Generate_fList(grids.new, eigenfunctions)
    #loss_REML = mean( ComputeLoss(fList_est_REML, eigenfList) )
    loss_REML =  ComputeLoss(fList_est_REML, eigenfList) 
    return(loss_REML)
}

EigenSpline = function(t, B, basis.method, splineBasis, R.inv){
    M<-nrow(B)
    r<-ncol(B)
    result<-numeric(r)
    
    if (basis.method == "poly"){
       # splineBasis = new(orthoSpline, 0, 1, 4, M - 2)
        phi.c = splineBasis$evalSpline(t) # t here ?
        phi.c = R.inv %*% phi.c
    }
    result<-t(B)%*%matrix(phi.c)
    return(result)
}

Spline.aux = function(basis.method, data.list, i, splineBasis, R.inv){
    data.c = data.list[[i]]
    t.v = data.c[[1]][, 2]
    y.v = data.c[[1]][, 1]
    phi.c <- BSpline.Evl(basis.method, t.v, splineBasis, R.inv)
    t.phi.c<-t(phi.c)
    psi.c<-phi.c%*%t.phi.c
    d.c<-phi.c%*%y.v%*%t(y.v)%*%t.phi.c
    return(list(phi.c,psi.c,d.c))
}

Basis.spline.orth = function(grids, splineBasis){
    bs = apply(matrix(grids), 1, Basis.Spline, splineBasis = splineBasis)
    temp = bs%*%t(bs) * (grids[2] - grids[1])
    R = t(chol(temp))
    result = solve(R, bs)
    return(list(result, R))
}

Spline.inverse = function(splineBasis, grids){
    R <- Basis.spline.orth(grids, splineBasis)[[2]]
    result <- solve(R)
    return(result)
}

### Not normalize version
Basis.Spline = function(t, splineBasis){
    bs.t <- splineBasis$evalSpline(t)
    return(bs.t)
}

Basis.Spline.t = function(t, splineBasis, R.inv){
    bs.t <- splineBasis$evalSpline(t)
    result <- R.inv%*%bs.t
    return(result)
}

BSpline.Evl = function( basis.method, t.v, splineBasis, R.inv){
    result <- numeric(length(t.v))
    if (basis.method == "poly"){
        result = apply(matrix(t.v), MARGIN = 1, Basis.Spline.t, splineBasis = splineBasis, R.inv)
    }
    return(result)
}

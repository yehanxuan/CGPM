#M.set = c(10, 11)
#r.set = c(8)

ini.method = "EM"
#ini.method = "LS"
#ini.method = "loc"
basis.method = "bs"
sl.v = rep(0.5, 10)
max.step = 50
#grid.l = seq(0,1,0.005)
grid.l = seq(0,1,0.01) # The default value
grids= seq(0,1,0.002)

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
    pcaCompNum = min(length(fList_est), length(fList_true))
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





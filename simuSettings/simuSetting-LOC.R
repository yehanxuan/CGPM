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


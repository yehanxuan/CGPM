



#------ rdatasimu1 ------
source("./utils/evaluate.R")
allMethod = c("LOC", "EM", "REML", "LogDet", "VNDiv", "frobDiverg")
InitType = "EM"
samplesize = 200 
nn = 200
scoreType = "Gaussian"
DataType = "Fourier"
collectResults_Fourier = function(samplesize, DataType, allMethod, scoreType, InitType){
    colSel = 1
    finalMean = matrix(0, 8, 6)
    finalSD = matrix(0, 8, 6)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        try({
            if ((m == 1) | (m == 2)){
                fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
                                  scoreType, "-", "EM", ".RData")
            } else {
                
                fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
                                  scoreType, "-", InitType, ".RData")
            }
        })
        
    
        try({
            load(fileList)
            resultT = compute(Final$MSE)
            res = resultT$res
            resSD = resultT$resSD
            finalMean[, colSel] = res
            finalSD[, colSel] = resSD
        })
        colSel = colSel + 1
    }
    return(list(finalMean, finalSD))
}

EMdata = collectResults_Fourier(samplesize, DataType, allMethod, scoreType, InitType = "EM")
finalMean = EMdata[[1]]
finalSD = EMdata[[2]]

allMethod = c("REML", "LogDet", "VNDiv", "frobDiverg")

collectSubmatrix = function(samplesize, DataType, allMethod, scoreType, InitType){
    colSel = 1
    finalMean2 = matrix(0, 8, 4)
    finalSD2 = matrix(0, 8, 4)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        try({
            fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
                              scoreType, "-", "LS", ".RData")
        })
        
        try({
            load(fileList)
            resultT = compute(Final$MSE)
            res = resultT$res
            resSD = resultT$resSD
            finalMean2[, colSel] = res
            finalSD2[, colSel] = resSD
        })
        
        colSel = colSel + 1
    }
    return(list(finalMean2, finalSD2))
}

LSdata = collectSubmatrix(samplesize, DataType, allMethod, scoreType, InitType = "LS")
finalMean2 = LSdata[[1]]
finalSD2 = LSdata[[2]]

finalMean_comb = cbind(finalMean, finalMean2)
finalSD_comb = cbind(finalSD, finalSD2)

outputNewFourier = function(finalMean, finalSD, nn){
    rk = 3
     cat("\\hline\n")
    nrow(finalMean)
    cat("\\multicolumn{1}{c}{$n =", nn, "\\quad ",
        "$", "}",  "\\\\\n\\hline\n", sep = "")
    for (rowI in 1:nrow(finalMean)){
       cat(paste0("$ \\| \\hat{\\psi}_{", rowI, "} - \\psi_{0" ,rowI, "} \\|_{L_2}$" ))
        for (colJ in 1:ncol(finalMean)){
            cat(" & ")
            cat(format(round( finalMean[rowI, colJ],rk), nsmall = rk))
            cat("(",format(round(finalSD[rowI, colJ],rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
    cat("\\hline\n")
}

outputNewFourier(finalMean_comb, finalSD_comb, nn)

#------ rdatasimu2 ------
samplesize = 500
nn = 500
allMethod = c("LOC", "EM", "REML", "LogDet", "VNDiv", "frobDiverg")
EMdata = collectResults_Fourier(samplesize, DataType, allMethod, scoreType, InitType = "EM")
finalMean = EMdata[[1]]
finalSD = EMdata[[2]]

allMethod = c("REML", "LogDet", "VNDiv", "frobDiverg")

LSdata = collectSubmatrix(samplesize, DataType, allMethod, scoreType, InitType = "LS")
finalMean2 = LSdata[[1]]
finalSD2 = LSdata[[2]]

finalMean_comb = cbind(finalMean, finalMean2)
finalSD_comb = cbind(finalSD, finalSD2)

outputNewFourier(finalMean_comb, finalSD_comb, nn)

#------ rdatasimut ------
allMethod = c("REML", "EM", "LOC", "frobDiverg", "LogDet", "VNDiv")
#
samplesize = 200 
nn = 200

EMdata = collectResults_Fourier(samplesize, DataType, allMethod, scoreType = "t", InitType = "EM")
finalMean = EMdata[[1]]
finalSD = EMdata[[2]]
outputNewFourier(finalMean, finalSD, nn)

samplesize = 500 
nn = 500

EMdata = collectResults_Fourier(samplesize, DataType, allMethod, scoreType = "t", InitType = "EM")
finalMean = EMdata[[1]]
finalSD = EMdata[[2]]
outputNewFourier(finalMean, finalSD, nn)

# ------ rdatasimuU ------

allMethod = c("REML", "EM", "LOC", "frobDiverg", "LogDet", "VNDiv")
#
samplesize = 200 
nn = 200

EMdata = collectResults_Fourier(samplesize, DataType, allMethod, scoreType = "uniform", InitType = "EM")
finalMean = EMdata[[1]]
finalSD = EMdata[[2]]
outputNewFourier(finalMean, finalSD, nn)

samplesize = 500 
nn = 500

EMdata = collectResults_Fourier(samplesize, DataType, allMethod, scoreType = "uniform", InitType = "EM")
finalMean = EMdata[[1]]
finalSD = EMdata[[2]]
outputNewFourier(finalMean, finalSD, nn)

#---- penalty ------

#### The role of penalty 
allMethod = c("frobDiverg", "LogDet", "VNDiv") 
samplesize = 500

nKnots = 35
PenaltySeq = c("N", "Y")

collectResults = function(samplesize, DataType, allMethod){
    colSel = 1
    finalMean = matrix(0, 8, 6)
    finalSD = matrix(0, 8, 6)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        for (p in 1:length(PenaltySeq)){
            Penalty = PenaltySeq[p]
            fileList = paste0("./penaltyPower/method-", method, "-", DataType, "-", samplesize,"-", nKnots, "-", scoreType,
                              "-", "LS" , "-", Penalty, ".RData")
            load(fileList)
            resultT = compute(Final$MSE)
            res = resultT$res
            resSD = resultT$resSD
            finalMean[, colSel] = res
            finalSD[, colSel] = resSD
            colSel = colSel + 1
        }
    }
    return(list(finalMean, finalSD))
}

outputNew = function(finalMean, finalSD, nn, InitType){
    rk = 3
    cat("\\hline\n")
    cat("\\multicolumn{7}{l}{$n =", nn, "\\quad ",
        "$", InitType,"}",  "\\\\\n\\hline\n", sep = "")
    nrow(finalMean)
    for (rowI in 1:nrow(finalMean)){
        cat(paste0("$ \\| \\hat{\\psi}_{", rowI, "} - \\psi_{0" ,rowI, "} \\|_{L_2}$" ))
        for (colJ in 1:6){
            cat(" & ")
            cat(format(round( finalMean[rowI, colJ],rk), nsmall = rk))
            cat("(",format(round(finalSD[rowI, colJ],rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
    cat("\\hline\n")
}

penaltyData = collectResults(samplesize, DataType, allMethod)
finalMean = penaltyData[[1]]
finalSD = penaltyData[[2]]

outputNew(finalMean, finalSD, nn= 500, InitType = "LS")

#------ realpred ------

allMethod = c("REML", "EM", "LOC", "frobDiverg", "LogDet", "VNDiv")
compSeq = c(2,3,4)
collectResults_Pred = function(allMethod, compSeq){
    colSel = 1
    finalMean = matrix(0, 3, 6)
    finalSD = matrix(0, 3, 6)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        for (j in 1:length(compSeq)){
            if (m == 1){
                fileList = paste0("predata/method-", method, "-", compSeq[j], "-", "LS",
                                  ".RData")
            } else {
                fileList = paste0("predata/method-", method, "-", compSeq[j], "-", "LS",
                                  ".RData")
            }
            
            load(fileList)
            resultT = compute(Final$MSE)
            res = resultT$res
            resSD = resultT$resSD
            finalMean[j, m] = res
            finalSD[j, m] = resSD
        }
       
    }
    return(list(finalMean, finalSD))
}

outputNew_pred = function(finalMean, finalSD){
    rk = 3
    cat("\\hline\n")
    cat("\\multicolumn{7}{l}{No.FPCs}",  "\\\\\n\\hline\n", sep = "")
    nrow(finalMean)
    for (rowI in 1:nrow(finalMean)){
        cat(paste0("$R=", compSeq[rowI], "$"))
        for (colJ in 1:ncol(finalMean)){
            cat(" & ")
            cat(format(round( finalMean[rowI, colJ]*100,rk), nsmall = rk))
            cat("(",format(round(finalSD[rowI, colJ]*100 ,rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
}

PredData = collectResults_Pred(allMethod, compSeq)
finalMean = PredData[[1]]
finalSD = PredData[[2]]
outputNew_pred(finalMean, finalSD)



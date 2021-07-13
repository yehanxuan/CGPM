
rm(list = ls())
library(snowfall)
library(rOptManifold)
library(mFPCA)
library(sm)
library(fpca)
library(fdapace)
library(mvtnorm)
library(FDABasics)

args <- commandArgs(trailingOnly = TRUE)
method = args[1]
DataType = args[2]
repID = as.numeric(args[3])
samplesize = as.numeric(args[4])
scoreType = args[5]
InitType = args[6]

    
#scoreType = NULL


method = "LogDet"
DataType = "prac"
repID = 2
samplesize = 500
scoreType = "Gaussian"
InitType = "EM"
    
source("./data/generate_data.R")
source("./models/LogDet.R")
#devtools::reload("../mFPCA/")
source("./utils/evaluate.R")
source("./optimization/functions_Optimization.R")
source("./models/LOC.R")
source("./models/EM.R")
source("./models/REML.R")
source("./optimization/functions_GenData.R")
source("./models/VNDiv.R")
source("./models/frobDiverg.R")
source("./models/Power.R")
source("./utils/CV.R")

source(paste0("./simuSettings/simuSetting-", method, ".R"))
if ((DataType == "NonUnif")||(DataType == "FourierOrth")|| (DataType == "FourierDense")){
    seedJ = 2
}
source(paste0("./oneReplicate/oneRep-", method, ".R"))





nCPUS = 8
maxIter = 8

result1 = list()
result2 = list()
result3 = list()
if (method == "mFPCA"){
    for (i in 1:ceiling(maxIter/nCPUS)){
        print(i)
        sfInit(parallel = TRUE, cpus = nCPUS)
        sfExportAll()
        sfLibrary(rOptManifold)
        sfLibrary(mFPCA)
        sfLibrary(mvtnorm)
        sBegin = (i-1)*nCPUS +1
        sEnd = min(i*nCPUS, maxIter)
        seedSeq = seq(sBegin, sEnd, by = 1)
        tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_mFPCA)
        sfStop()
        tmp = Filter(function(x) length(x) >= 3, tmp)
        tmp1 = lapply(tmp, function(x) x[[1]])
        tmp2 = lapply(tmp, function(x) x[[2]])
        tmp3 = lapply(tmp, function(x) x[[3]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
    }
} else if (method == "REML"){
    for (i in 1:ceiling(maxIter/nCPUS)){
        print(i)
        sfInit(parallel = TRUE, cpus = nCPUS)
        
        sfExportAll()
        sfLibrary(rOptManifold)
        sfLibrary(mFPCA)
        sfLibrary(sm)
        sfLibrary(fpca)
        sfLibrary(splines)
        sfLibrary(mvtnorm)
        sfLibrary(FDABasics)
        sBegin = (i-1)*nCPUS +1
        sEnd = min(i*nCPUS, maxIter)
        seedSeq = seq(sBegin, sEnd, by = 1)
        tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_REML)
        sfStop()
        tmp = Filter(function(x) length(x) >= 3, tmp)
        tmp1 = lapply(tmp, function(x) x[[1]])
        tmp2 = lapply(tmp, function(x) x[[2]])
        tmp3 = lapply(tmp, function(x) x[[3]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
    }
} else if (method =="EM") {
    for (i in 1:ceiling(maxIter/nCPUS)){
        print(i)
        sfInit(parallel = TRUE, cpus = nCPUS)
        sfExportAll()
        sfLibrary(rOptManifold)
        sfLibrary(mFPCA)
        sfLibrary(sm)
        sfLibrary(splines)
        sfLibrary(mvtnorm)
        sfLibrary(FDABasics)
#        sfLibrary(fdapace)
        sBegin = (i-1)*nCPUS +1
        sEnd = min(i*nCPUS, maxIter)
        seedSeq = seq(sBegin, sEnd, by = 1)
        tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_EM)
        sfStop()
        tmp = Filter(function(x) length(x) >= 3, tmp)
        tmp1 = lapply(tmp, function(x) x[[1]])
        tmp2 = lapply(tmp, function(x) x[[2]])
        tmp3 = lapply(tmp, function(x) x[[3]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
    } 
} else if (method == "LOC") {
    for (i in 1:ceiling(maxIter/nCPUS)){
        print(i)
        sfInit(parallel = TRUE, cpus = nCPUS)
        sfExportAll()
        sfLibrary(rOptManifold)
        sfLibrary(mFPCA)
        sfLibrary(sm)
        sfLibrary(splines)
        sfLibrary(mvtnorm)
        sfLibrary(FDABasics)
        #        sfLibrary(fdapace)
        sBegin = (i-1)*nCPUS +1
        sEnd = min(i*nCPUS, maxIter)
        seedSeq = seq(sBegin, sEnd, by = 1)
        tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_LOC)
        sfStop()
        tmp = Filter(function(x) length(x) >= 3, tmp)
        tmp1 = lapply(tmp, function(x) x[[1]])
        tmp2 = lapply(tmp, function(x) x[[2]])
        tmp3 = lapply(tmp, function(x) x[[3]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
    } 
} else if (method == "LogDet"){
    for (i in 1:ceiling(maxIter/nCPUS)){
        print(i)
        sfInit(parallel = TRUE, cpus = nCPUS)
        sfExportAll()
        sfLibrary(rOptManifold)
        sfLibrary(mFPCA)
        sfLibrary(sm)
        sfLibrary(fpca)
        sfLibrary(splines)
        sfLibrary(mvtnorm)
        sfLibrary(FDABasics)
        #        sfLibrary(fdapace)
        sBegin = (i-1)*nCPUS +1
        sEnd = min(i*nCPUS, maxIter)
        seedSeq = seq(sBegin, sEnd, by = 1)
        tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_LogDet)
        sfStop()
        tmp = Filter(function(x) length(x) >= 3, tmp)
        tmp1 = lapply(tmp, function(x) x[[1]])
        tmp2 = lapply(tmp, function(x) x[[2]])
        tmp3 = lapply(tmp, function(x) x[[3]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
    }
} else if ((method == "VNDiv")||(method == "frobDiverg") ){
    for (i in 1:ceiling(maxIter/nCPUS)){
        print(i)
        sfInit(parallel = TRUE, cpus = nCPUS)
        sfExportAll()
        sfLibrary(rOptManifold)
        sfLibrary(mFPCA)
        sfLibrary(sm)
        sfLibrary(fpca)
        sfLibrary(splines)
        sfLibrary(mvtnorm)
        sfLibrary(FDABasics)
        sBegin = (i-1)*nCPUS +1
        sEnd = min(i*nCPUS, maxIter)
        seedSeq = seq(sBegin, sEnd, by = 1)
        if (method == "VNDiv"){
            tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_VNDiv)
        } else if (method == "frobDiverg"){
            tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_frobDiverg)
        }
        sfStop()
        tmp = Filter(function(x) length(x) >= 3, tmp)
        tmp1 = lapply(tmp, function(x) x[[1]])
        tmp2 = lapply(tmp, function(x) x[[2]])
        tmp3 = lapply(tmp, function(x) x[[3]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
    }
} else if (method == "Power"){
    for (i in 1:ceiling(maxIter/nCPUS)){
        print(i)
        sfInit(parallel = TRUE, cpus = nCPUS)
        sfExportAll()
        sfLibrary(rOptManifold)
        sfLibrary(mFPCA)
        sfLibrary(sm)
        sfLibrary(fpca)
        sfLibrary(splines)
        sfLibrary(mvtnorm)
        sfLibrary(FDABasics)
        #        sfLibrary(fdapace)
        sBegin = (i-1)*nCPUS +1
        sEnd = min(i*nCPUS, maxIter)
        seedSeq = seq(sBegin, sEnd, by = 1)
        tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_Power)
        sfStop()
        tmp = Filter(function(x) length(x) >= 3, tmp)
        tmp1 = lapply(tmp, function(x) x[[1]])
        tmp2 = lapply(tmp, function(x) x[[2]])
        tmp3 = lapply(tmp, function(x) x[[3]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
    }
}



Final = list("MSE" = result1, "lambda" = result2, "rank" = result3)
#Final = list("rank" = rank, "lambda" = lambda, "Err" = Err)
#Err = compute(result1)
#Err 

#lambda = result2
#table(unlist(lambda))

#rank = result3
#table(unlist(rank))

if (method == "mFPCA"){
    savepath = paste0("./data/method-", method, "nKnots-", nKnots, ".RData")
} else if (method == "REML"){
    if ( !is.null(scoreType) ){
        if (!is.null(InitType)){
            savepath = paste0("./data/method-", method , "-", DataType, "-", samplesize, "-", scoreType, "-", InitType, ".RData")
        } else {
            savepath = paste0("./data/method-", method , "-", DataType, "-", samplesize, "-", scoreType, ".RData")
        }
    } else {
    savepath = paste0("./data/method-", method , "-", DataType, "-", samplesize, ".RData") 
    }
} else if (method == "EM"){
    if ( !is.null(scoreType)){
        if (!is.null(scoreType)){
            savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", scoreType, "-", InitType, ".RData")
        } else {
        savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", scoreType, ".RData")
        } 
    } else {
        savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, ".RData")
    }
} else if (method == "LOC"){
    if (!is.null(scoreType)){
        savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", scoreType, "-", InitType, ".RData")
    } else {
        savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, ".RData")
    }
} else if (method == "LogDet"){
    if (!is.null(InitType)){
        if (!is.null(scoreType)){
            savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", nKnots, "-", scoreType, "-", InitType, ".RData")
        } else {
            savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", nKnots, "-", InitType, ".RData")
        }
    } else {
        savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", nKnots, ".RData")
    }
} else if (method == "VNDiv"){
    if (!is.null(InitType)){
        if (!is.null(scoreType)){
            savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", nKnots, "-", scoreType, "-", InitType, ".RData")
        } else {
            savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", nKnots, "-", InitType, ".RData")
        }
    } else {
        savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize,"-", nKnots ,".RData")
    }
} else if (method == "frobDiverg" ){
    if (!is.null(InitType)){
        if (!is.null(scoreType)){
            savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", nKnots, "-", scoreType, "-", InitType, ".RData")
        } else {
            savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", nKnots, "-", InitType, ".RData")
        }
    } else {
        savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize,"-", nKnots ,".RData")
    }
} else if (method == "Power"){
    savepath = paste0("./data/method-", method, "-", DataType, "-", samplesize, "-", nKnots, "-", scoreType, "-", InitType, ".RData")
}


try({
    save(Final, file = savepath)
})


Err = compute(result1)
Err 

lambda = result2
table(unlist(lambda))

rank = result3
table(unlist(rank))


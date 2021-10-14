
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
noiseType = args[7]

    
#scoreType = NULL


method = "REML"
DataType = "pracSin"
repID = 2
samplesize = 500
scoreType = "Gaussian"
InitType = "EM"
noiseType = "Gaussian"
    
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
result4 = list()
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
        tmp4 = lapply(tmp, function(x) x[[4]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
        result4 = c(result4, tmp4)
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
        tmp4 = lapply(tmp, function(x) x[[4]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
        result4 = c(result4, tmp4)
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
        tmp4 = lapply(tmp, function(x) x[[4]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
        result4 = c(result4, tmp4)
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
        tmp4 = lapply(tmp, function(x) x[[4]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
        result4 = c(result4, tmp4)
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



Final = list("MSE" = result1, "lambda" = result2, "rank" = result3, "converge" = result4)
#Final = list("rank" = rank, "lambda" = lambda, "Err" = Err)
#Err = compute(result1)
#Err 
Final$MSE[which(Final$converge != 0 )]
compute(Final$MSE[which(Final$converge != 0 )])
#lambda = result2
#table(unlist(lambda))

#rank = result3
#table(unlist(rank))
folder = DataType

if (file.exists(folder)) {
    cat("The folder already exists")
} else {
    dir.create(folder)
}


if (method == "REML") {
    savepath = paste0("./", DataType, "/method-", method, "-", DataType, "-", samplesize, "-",  scoreType, "-", InitType, "-", noiseType, ".RData")
} else if (method == "EM") {
    savepath = paste0("./", DataType, "/method-", method, "-", DataType, "-", samplesize, "-",  scoreType, "-", InitType, "-", noiseType, ".RData")
} else if (method == "LogDet") {
    savepath = paste0("./", DataType, "/method-", method, "-", DataType, "-", samplesize, "-",  scoreType, "-", InitType, "-", noiseType,
                      "-", nKnots,"-", select_Method, ".RData")
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


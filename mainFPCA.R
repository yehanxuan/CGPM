rm(list = ls())
library(snowfall)
library(rOptManifold)
library(mFPCA)
library(sm)
#library(fpca) # fpca is not available for R version 4.2.1
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

method = "LogDet"
DataType = "easy"
repID = 2
samplesize = 100
scoreType = "Gaussian"
InitType = "EM"
noiseType = "Gaussian"
    
source("./data/generate_data.R")
source("./models/LogDet.R")
source("./utils/evaluate.R")
source("./models/LOC.R")
source("./models/EM.R")
source("./models/REML.R")
source("./optimization/functions_Optimization.R")
source("./optimization/functions_GenData.R")
source("./utils/CV.R")

source(paste0("./simuSettings/simuSetting-", method, ".R"))
source(paste0("./oneReplicate/oneRep-", method, ".R"))

nCPUS = 8
maxIter = 8
result1 = list()
result2 = list()
result3 = list()
result4 = list()
result5 = list()
result6 = list()
result7 = list()
if (method == "REML"){
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
        tmp5 = lapply(tmp, function(x) x[[5]])
        tmp6 = lapply(tmp, function(x) x[[6]])
        tmp7 = lapply(tmp, function(x) x[[7]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
        result4 = c(result4, tmp4)
        result5 = c(result5, tmp5)
        result6 = c(result6, tmp6)
        result7 = c(result7, tmp7)
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
        tmp5 = lapply(tmp, function(x) x[[5]])
        tmp6 = lapply(tmp, function(x) x[[6]])
        tmp7 = lapply(tmp, function(x) x[[7]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
        result4 = c(result4, tmp4)
        result5 = c(result5, tmp5)
        result6 = c(result6, tmp6)
        result7 = c(result7, tmp7)
    } 
} else if (method == "LogDet"){
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
        tmp5 = lapply(tmp, function(x) x[[5]])
        tmp6 = lapply(tmp, function(x) x[[6]])
        tmp7 = lapply(tmp, function(x) x[[7]])
        result1 = c(result1, tmp1)
        result2 = c(result2, tmp2)
        result3 = c(result3, tmp3)
        result4 = c(result4, tmp4)
        result5 = c(result5, tmp5)
        result6 = c(result6, tmp6)
        result7 = c(result7, tmp7)
    }
}



Final = list("MSE" = result1, "lambda" = result2, "rank" = result3, "converge" = result4, "runTime" = result5,
             "sysTime" = result6, "eigenValues" = result7)


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



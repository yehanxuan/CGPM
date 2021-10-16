rm(list = ls())

library(snowfall)
library(sm)
library(readr)
require(readxl)
require(tidyverse)
library(fpca)
library(fdapace)
library(FDABasics)
library(mFPCA)
library(rOptManifold)


args <- commandArgs(trailingOnly = TRUE)
method = args[1]
repID = as.numeric(args[2]) 
rank = as.numeric(args[3])
InitType = args[4]


#method="LogDetPenalty"
#repID=2
#r.set=3
#InitType = "EM"


source("./models/Prediction.R")
file_list = list.files(path = "./CSP_Photometry_DR2/")
obsCol = Astro_data(file_list)




source("./data/generate_data.R")
source("./models/LogDet.R")
source("./models/EM.R")
source("./models/REML.R")
source("./models/LOC.R")

source("./utils/CV.R")
source("./utils/evaluate.R")

source("./optimization/functions_Optimization.R")
source("./optimization/functions_LocLin.R")
source("./optimization/functions_GenData.R")


source("./simuSettings/simuSetting-pred.R")
source("./oneReplicate/oneRep-pred.R")










savepath = paste0("./predata/method-", method, "-", rank, "-", InitType, ".RData") 




nCPUS = 10
maxIter = 10



result1 = list()
result2 = list()
result3 = list()
result4 = list()
for (i in 1:ceiling(maxIter/nCPUS)) {
    print(i)
    sfInit(parallel = TRUE, cpus = nCPUS)
    sfExportAll()
    sfLibrary(rOptManifold)
    sfLibrary(mFPCA)
    sfLibrary(FDABasics)
    sfLibrary(fdapace)
    sfLibrary(sm)
    sfLibrary(readxl)
    sfLibrary(readr)
    sfLibrary(sm)
    sfLibrary(splines)
    sBegin = (i-1)*nCPUS + 1
    sEnd = min(i*nCPUS, maxIter)
    seedSeq = seq(sBegin, sEnd, by = 1)
    tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_pred)
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

Final = list("MSFE" = result1, "Time" = result2, "lambda" = result3, "rank" = result4)

try({
    save(Final, file = savepath)
})

Err = compute(result1)
Err

Time = result2
#table(unlist(lambda))

lambda = result3
table(unlist(lambda))

rank = result4
table(unlist(rank))









library(readr)
require(readxl)
require(tidyverse)
library(fpca)
library(FDABasics)
library(mFPCA)
library(rOptManifold)
library(fdapace)
#Table = read_tsv("./CSP_Photometry_DR2/SN2004dtopt+nir_photo.dat", skip = 4, na = "99.900")
#obsY_i = Table$i
#obsY_Y = Table$Y
#obsY_J =  Table$J
#obsY_H = Table$H


file_list = list.files(path = "./CSP_Photometry_DR2/")
source("./utils/load_data.R")

obsCol = Astro_data(file_list)
source("./data/generate_data.R")
source("./models/leastSquares.R")
#devtools::reload("../mFPCA/")
source("./utils/evaluate.R")
source("./optimization/functions_Optimization.R")
source("./optimization/functions_LocLin.R")
source("./optimization/functions_EM.R")
source("./optimization/fpca_func.R")
source("./optimization/functions_GenData.R")
#source("./code/functions_initial.R")
source("./optimization/prediction_func.R")
source("./models/VNDiv.R")
source("./models/frobDiverg.R")


newObsCol = obsCol[, -2]
newObsCol[ , c(2,3)] <- newObsCol[ ,c(3,2)]
colnames(newObsCol)[c(2,3)] = colnames(newObsCol)[c(3,2)]



rankSeq = 4
ini.method = "EM.self"
basis.method = "bs"
sl.v = rep(0.5, 10)
max.step = 50
grid.l = seq(0,1,0.01)
grids= seq(0,1,0.002)
sig.EM = 1
IniVal = Init_method(newObsCol, M = NULL, rankSeq, ini.method = "EM.self")

sig2hat = IniVal[[1]]
covmatrix.ini<-IniVal[[2]]
eigenf.ini<-IniVal[[3]]
eigenv.ini<-IniVal[[4]]
like.ini<-IniVal[[5]]

eigenfest = t(eigenf.ini)
#obsMat$elemID  = as.factor(obsMat$elemID)
splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
K = splineObj$getDoF()
basisMat = splineObj$evalSpline(grids)

UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
WInit = eigenv.ini

#decomp = qr(UInit)
#Q = qr.Q(decomp)
#R = qr.R(decomp)
XInit = list(UInit, diag(WInit))

source("./code/simuSetting-RcppMLE.R")

modelini = Convert(IniVal)
modelini = c(modelini, 
             list(tmin = tmin, tmax = tmax,
                  SFinal = XInit, sigmaSq = sig2hat,
                  numPCA = rankSeq, numElem = 1,
                  elemLevels = levels(obsCol$elemID) ))


#plot(grids, eigenfest[3,], type = "l")









rankSeq = 4
#select = (obsCol, splineObj, rankSeq, 0, controlList1,
#                       controlList2, nFold = 10, SInit = NULL, sigmaSq = 0.1)
nSample = length(unique(obsCol$obsID))
cf = -1
cvMembership = getCVPartition(nSample, nFold)
cvParam = list(cvMembership = cvMembership, cf = cf)


meanModel = fitMeanCurve(obsCol, splineObj, lambda = 1e-3)

obsCol_c = subtractMeanCurve(meanModel, obsCol)

#obsCol_d = obsCol_c[obsCol_c$obsID < 8]
#obsCol_c$obsID = as.factor( obsCol_c$obsID )
#p1 = ggplot(obsCol_c, aes(x = obsCol_c$obsT, y = obsCol_c$obsY, color = factor(obsCol_c$obsID))) + geom_point()
#p1

plotMeanModel(meanModel)


rankSeq  = 3
lambdaSeq = exp(seq(-20,-12, length.out = 5))
nFold = 10
sig2hat = 1

select1 = VNDiv_selection(obsCol_c, splineObj, rankSeq, lambdaSeq, controlList1, controlList2, nFold, sigmaSq = sig2hat, eigenfList = NULL, InitType = "EM")
model1 = MFPCA_EstimateVNDiv(obsCol_c,splineObj, optRank, mu, controlList1, controlList2, cvParam, SInit = XInit, sigmaSq, eigenfList)

model1 = MFPCA_EstimateVNDiv(obsCol_c, splineObj, 4, 1e-5, controlList1, controlList2, cvParam,
                             SInit = XInit, sigmaSq = sig2hat)


select2 =frobDiverg_selection( obsCol_c, splineObj, rankSeq, lambdaSeq, controlList1, controlList2, nFold = 10, SInit = NULL, sigmaSq = sig2hat)
model2 = MFPCA_EstimatefrobDiverg(obsCol_c, splineObj, 4, 1e-7, controlList1, controlList2, cvParam,
                                  SInit = XInit, sigmaSq = sig2hat)



controlList2 = list(alpha = 0.1, tol = 1e-6, sigma = 5e-2, beta = 0.618,
                   iterMax = 2000, verbose = 1)
rankSeq = c(4)
lambdaSeq = exp(seq(-17,-10, length.out = 5))
select3 = MLE_selection(obsCol_c, splineObj, rankSeq, lambdaSeq, controlList1, controlList2, nFold = 10, SInit = XInit, sigmaSq = sig2hat)
model3 = MFPCA_EstimateMLE(obsCol_c, splineObj, 4, 4e-5, controlList1, controlList2, cvParam,
                           SInit = XInit, sigmaSq = sig2hat)

plotRealEstimate(modelini)
plotRealEstimate(model1)
plotRealEstimate(model2)
plotRealEstimate(model3)

ComparePlot(modelini, model1, model2, model3)

ComparePlot = function(Model1, Model2, Model3, Model4, Model5, Model6, selK = NULL){
    tmin = Model1$tmin
    tmax = Model1$tmax
    numPCA = Model1$numPCA
    numElem = Model1$numElem
    nSeq = 500
    tSeq = seq(0.0001, 0.9999, length.out = nSeq)
    plotData = data.frame();#matrix(0, nSeq * numPCA * numElem, 4)  
    selR = 1:nSeq
    
    for(k in 1:numPCA){
        for(e in 1:numElem){
            em = Model1$elemLevels[e]
            #fSeq0 = trueEigenFList[[k]][[e]](tSeq)
            fSeqHat1 = Model1$eigenFunctions[[k]][[e]](tSeq)
            fSeqHat2 = Model2$eigenFunctions[[k]][[e]](tSeq)
            fSeqHat3 = Model3$eigenFunctions[[k]][[e]](tSeq)
            fSeqHat4 = Model4$eigenFunctions[[k]][[e]](tSeq)
            fSeqHat5 = Model5$eigenFunctions[[k]][[e]](tSeq)
            fSeqHat6 = Model6$eigenFunctions[[k]][[e]](tSeq)
            
            if(sum(fSeqHat1*fSeqHat2) < 0) fSeqHat2 = -fSeqHat2
            if(sum(fSeqHat1*fSeqHat3) < 0) fSeqHat3 = -fSeqHat3
            if(sum(fSeqHat1*fSeqHat4) < 0) fSeqHat4 = -fSeqHat4
            if(sum(fSeqHat1*fSeqHat5) < 0) fSeqHat5 = -fSeqHat5
            if(sum(fSeqHat1*fSeqHat6) < 0) fSeqHat6 = -fSeqHat6
            
            tmp1 = data.frame(obsT = tSeq, obsY = fSeqHat1, 
                             pcaID =  k, elemID =  em,
                             curveID = "EM", stringsAsFactors =  F)
            tmp2 = data.frame(obsT = tSeq, obsY = fSeqHat2, 
                              pcaID =  k, elemID =  em,
                              curveID = "VN", stringsAsFactors =  F)
            tmp3 = data.frame(obsT = tSeq, obsY = fSeqHat3, 
                              pcaID =  k, elemID =  em,
                              curveID = "Frob", stringsAsFactors =  F)
            tmp4 = data.frame(obsT = tSeq, obsY = fSeqHat4, 
                              pcaID =  k, elemID =  em,
                              curveID = "LogDet", stringsAsFactors =  F)
            tmp5 = data.frame(obsT = tSeq, obsY = fSeqHat5, 
                              pcaID =  k, elemID =  em,
                              curveID = "REML", stringsAsFactors =  F)
            tmp6 = data.frame(obsT = tSeq, obsY = fSeqHat6, 
                              pcaID =  k, elemID =  em,
                              curveID = "LOC", stringsAsFactors =  F)
            tmp = rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
            plotData = rbind(plotData, tmp)
            selR = selR + nSeq
        }
    }
    colnames(plotData) = c("obsT", "obsY", "pcaID", "elemID", "curveID")
    plotData$elemID = factor(plotData$elemID, levels = Model1$elemLevels)
    
    
    if(!is.null(selK)){
        plotData = subset(plotData, plotData$pcaID == selK)
        
        p = ggplot(plotData, aes(obsT, obsY, 
                                 group = curveID, color = curveID)) +
            geom_line()
        p = p + facet_wrap(~elemID)
        
    }else{
       # p = ggplot(plotData, aes(obsT, obsY, 
    #                             group = curveID, color = curveID)) +
     #       geom_line( aes(linetype = curveID) ) + labs( color = "Method")
        
        p = ggplot(plotData) +
            geom_line( aes(obsT, obsY, linetype = curveID, group = curveID, color = curveID), size = 0.7 ) + 
            scale_colour_manual(values = c("red",'orange', "green", "#009E73", "Blue", "purple"))
        #+ labs( color = "Method")
       # to_string <- as_labeller( multi_line = FALSE)
        #labeller = to_string
        p = p + facet_wrap(~pcaID, scales= "free") + theme_bw() +xlab("Scaled Time") + ylab("Scaled magnitude") +
            theme(legend.position="top") + labs(linetype = "Method", color = "Method") 
        
    }
    return(p)
}



###### Plot of 10 lines 
library(ggplot2)
obs10  =  obsCol[obsCol$obsID < 11, ]

p1 = ggplot(obs10, aes(obsT, obsY, group = obs10$obsID, color = obs10$obsID)) + geom_line()
p1





elemLevels = meanModel$elemLevels
numPoints = 200

plotData = data.frame()
eI = 0
for(eName in elemLevels){
    eI = eI + 1
    tmin = meanModel$modelList[[eI]]$tmin
    tmax = meanModel$modelList[[eI]]$tmax
    tSeq = seq(tmin, tmax, length.out = numPoints)
    ySeq = meanModel$modelList[[eI]]$meanFunction(tSeq)
    tmp = data.frame(eName, tSeq, ySeq)
    plotData = rbind(plotData, tmp)
}


plotData$obsID = 6
plotData[, c(1,2,3,4)] = plotData[, c(4,1,2,3)]

colnames(plotData) = c("obsID", "elemID", "obsT", "obsY")
#plotData$elemID = factor(plotData$elemID, levels = meanModel$elemLevels)

#### Draw only one supernova 

obsCol = obsCol[ obsCol$obsID < 6 , ] #### Draw only one supernova 


####
file_list_new = file_list[70]
obsCol= Astro_data(file_list_new)

obsCol$obsID = as.numeric(obsCol$obsID)
plotData$obsID = as.numeric(plotData$obsID)
obsCol = rbind(obsCol, plotData)
#q = ggplot(plotData, aes(obsT, obsY)) +
#    geom_line() + xlab("Scaled time") + ylab("Reversed magnitude") + theme_bw()

obsCol$obsID = as.factor(obsCol$obsID)
p2 = ggplot(obsCol, aes(obsT, obsY, group = obsID, color = obsID)) + geom_point(data = filter(obsCol, obsCol$obsID !=6))  +  theme_bw()+
    xlab("Scaled time") + ylab("Scaled magnitude") + scale_colour_manual(labels = c("u", "B", "V", "g", "r", "mean"), values = c("red",'orange', "green", "Blue", "purple", "black"))+
    theme(legend.position="top") + labs(color = "bands")  + geom_line(data = filter(obsCol, obsCol$obsID ==6 ) )
p2



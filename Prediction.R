#### Prediction task
library(readr)
require(readxl)
require(tidyverse)
library(fpca)
library(FDABasics)
library(mFPCA)
library(rOptManifold)
#Table = read_tsv("./CSP_Photometry_DR2/SN2004dtopt+nir_photo.dat", skip = 4, na = "99.900")
#obsY_i = Table$i
#obsY_Y = Table$Y
#obsY_J =  Table$J
#obsY_H = Table$H


file_list = list.files(path = "./CSP_Photometry_DR2/")
obsCol = Astro_data(file_list)

source("./optimization/Prediction_func.R")
nSample = length(unique(obsCol$obsID))
Membership = getCVPartition(nSample, 5)

tmp = (Membership == 1)
train_Index = which(tmp == FALSE)
test_Index = which(tmp == TRUE)
### 80% training data and 20% testing data
train_origin = obsCol[obsCol$obsID %in% train_Index, ] 
test_data = obsCol[obsCol$obsID %in% test_Index, ]
#train_origin = train_data

tmin = 0
tmax = 1
mOrder = 4
nKnots = 30
### Estimate mean function
splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)

meanSeq = exp(seq(-10,-2, length.out = 7))
mean_mu = MeanModel_GCV(train_origin, splineObj, meanSeq)
#meanModel = fitMeanCurve(train_origin, splineObj, lambda = mean_mu)

meanModel = fitMeanCurve(train_origin, splineObj, lambda = mean_mu)
train_data = subtractMeanCurve(meanModel, train_origin)

plotMeanModel_AS = function(meanModel, addPlot = FALSE){
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
    colnames(plotData) = c("elemID", "obsT", "obsY")
    plotData$elemID = factor(plotData$elemID, levels = meanModel$elemLevels)
    
    if(addPlot){
        q = geom_line(aes(obsT, obsY), 
                      data = plotData, color = "black") 
    }else{
        q = ggplot(plotData, aes(obsT, obsY)) +
            geom_line() + xlab("Scaled time") + ylab("Reversed magnitude") + theme_bw()
    }
    return(q)
}

plotMeanModel_AS(meanModel)



ini.method = "EM"
M.set = c(15:25)


#meanModel = loc_meanModel(train_origin)
#train_data = subtractMeanCurve(meanModel, train_origin)
modelini = ini_estimate(train_data, 4, splineObj, ini.method = "EM", M.set)




XInit = modelini$SFinal
sig2hat = modelini$sigmaSq
plotRealEstimate(modelini)
modelini$sigmaSq

nSample = length(unique(train_data$obsID))
cf = -1
cvMembership = getCVPartition(nSample, nFold = 10)
cvParam = list(cvMembership = cvMembership, cf = cf)
controlList1 = list(alpha = 0.618, tol = 1e-6, iterMax = 200, sigma = 0.05,verbose = 0)
#controlList2 = list(alpha = 1e-1, tol = 1e-7, sigma = 1e-3, beta = 0.618,
#                    iterMax = 500, verbose = 0)

controlList2 = list(alpha = 0.1, tol = 1e-7, sigma = 5e-2, beta = 0.618,
                    iterMax = 2000, verbose = 1)
lambdaSeq = exp(seq(-17,-10, length.out = 7))


#model1 = training_Wrap(train_data, splineObj, 4, 5e-6, method = "VNDiv", is_mean = FALSE)
model1 = MFPCA_EstimateVNDiv(train_data, splineObj, 4, 5e-6, controlList1, controlList2, cvParam,
                             SInit = XInit, sigmaSq = sig2hat)


plotRealEstimate(model1)
model1$sigmaSq

model2 = MFPCA_EstimatefrobDiverg(train_data, splineObj, 4, 7e-8, controlList1, controlList2, cvParam,
                         SInit =XInit, sigmaSq = sig2hat)

plotRealEstimate(model2)
model2$sigmaSq

rankSeq = c(3,4)
select3 = MLE_selection(train_data, splineObj, 4, 1e-5, controlList1, controlList2, nFold = 10, SInit = XInit, sigmaSq = sig2hat)
model3 = MFPCA_EstimateMLE(train_data, splineObj, 4, 7e-6, controlList1, controlList2, cvParam, SInit = XInit, sigmaSq = sig2hat)
plotRealEstimate(model3)
model3$sigmaSq
#SFinal = model$SFinal
#sigmaSq = model$sigmaSq



M.set = c(15: 22)
r.set = c(4)

model4 = Newton_estimate(train_data, M.set, r.set, splineObj)


plotRealEstimate(model4)
model4$sigmaSq

library(fdapace)
model5 = LOC_PACE(train_data, 4 , splineObj)
plotRealEstimate(model5)
model5$sigmaSq



ComparePlot(modelini, model1, model2, model3, model4, model5)
################



DataList = Split_data(test_data, splineObj)

#DataList = list("est_obsYList" = est_obsYList, "pred_obsYList" = pred_obsYList, "est_bMatList" = est_bMatList,
#                "pred_bMatList" = pred_bMatList, "est_TList" = est_TList, "pred_TList" = pred_TList)

#############
pred_resultini = Estimate_score_whole(modelini, splineObj, DataList, meanModel = meanModel)
pred_result1 =  Estimate_score_whole(model1, splineObj, DataList, meanModel = meanModel)
pred_result2 = Estimate_score_whole(model2, splineObj, DataList, meanModel = meanModel)
pred_result3 = Estimate_score_whole(model3, splineObj, DataList, meanModel = meanModel)
pred_result4 = Estimate_score_whole(model4, splineObj, DataList, meanModel = meanModel)
pred_result5 = Estimate_score(model5, splineObj, DataList, meanModel = meanModel)

MSFE(pred_resultini$pred_YList, DataList$pred_obsYList)
MSFE(pred_result1$pred_YList, DataList$pred_obsYList)
MSFE(pred_result2$pred_YList, DataList$pred_obsYList)
MSFE(pred_result3$pred_YList, DataList$pred_obsYList)
MSFE(pred_result4$pred_YList, DataList$pred_obsYList)
MSFE(pred_result5$pred_YList, DataList$pred_obsYList)


X3 = pred_result3$pred_YList[[1]]
Time3 = pred_result3$pred_TList[[1]]
se3 = sqrt(pred_result3$varList[[1]])

plotData = MakePlotData(X3,Time3, se3, "LogDet")
#sqrt(pred_result$varList)
#plotData = cbind(0, se3, Time3, X3)
#plotData = as.data.frame(plotData)
#colnames(plotData) =  c("Curve", "se","Time",'X')
#plotData$Curve = "LogDet"
names(plotData)

real_data = test_data[test_data$obsID == 4, ]
real_data[, 2] = 0
real_data$obsID = "True"
names(real_data) <- names(plotData)

plotData = rbind(plotData, real_data)

########### Use red and orange to distinguish est and test data 
Est_Time = DataList$est_TList[[1]]
Est_X = DataList$est_obsYList[[1]]
Est_se = 0
plotData_Est = MakePlotData(Est_X, Est_Time, Est_se, "obs")

Pred_Time = DataList$pred_TList[[1]]
Pred_X = DataList$pred_obsYList[[1]]
Pred_se = 0
plotData_pred = MakePlotData(Pred_X, Pred_Time, Pred_se, "held-out")

plotData = rbind(plotData, plotData_Est, plotData_pred)

X1 = pred_result1$pred_YList[[1]]
Time1 = pred_result1$pred_TList[[1]]
se1 = sqrt(pred_result1$varList[[1]])
plotData1 = MakePlotData(X1, Time1, se1, "VNDiv")
#plotData1 = cbind(0, se1, Time1, X1)
#plotData1 = as.data.frame(plotData1)
#colnames(plotData1) =  c("Curve", "se","Time",'X')
#plotData1$Curve = "VNDiv"
plotData = rbind(plotData, plotData1)

X2 = pred_result2$pred_YList[[1]]
Time2 = pred_result2$pred_TList[[1]]
se2 = sqrt(pred_result2$varList[[1]])
plotData2 = MakePlotData(X2, Time2, se2, "frob")
plotData = rbind(plotData, plotData2)

X_ini = pred_resultini$pred_YList[[1]]
Time_ini = pred_resultini$pred_TList[[1]]
se_ini = sqrt(pred_resultini$varList[[1]])
plotData_ini = MakePlotData(X_ini, Time_ini, se_ini, "EM")
#plotData_ini = as.data.frame(plotData2)
#colnames(plotData_ini) =  c("Curve", "se","Time",'X')
#plotData_ini$Curve = "EM"
plotData = rbind(plotData, plotData_ini)


plotData$Curve = as.factor(plotData$Curve)

plot = ggplot(plotData, aes(Time, X, group = Curve, color = Curve)) + geom_point(data = filter(plotData, plotData$Curve !="LogDet")) + geom_line(data = filter(plotData, plotData$Curve =="LogDet"))+ 
    geom_ribbon(data = filter(plotData, plotData$Curve =="LogDet"), aes(x=Time, ymin = X - 1.96*se, ymax = X + 1.96*se ),alpha = I(1/5) ,size = 0.7, linetype = 0)
plot = plot + theme_bw() + scale_linetype_manual(values=c("dotdash","dotted","solid", "longdash", "3313"))
plot = plot +  scale_colour_manual(values = c("orange",'blue', "red")) +
         scale_fill_manual(values =  c("blue",'orange', "#009E73", "red", "purple") ) + 
         theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
               legend.position = c(0.1, 0.14), legend.title = element_blank()) + xlab("Scaled time") + ylab("Scaled maginitude") 
#+
#        scale_fill_manual(values =  c("orange",'blue', "red") )
plot
# plot = ggplot(plotData, aes(Time, X, group = Curve, color = Curve)) + geom_line(aes(linetype = Curve)) + 
#     geom_ribbon(aes(x=Time, ymin = X - 1.96*se, ymax = X + 1.96*se,linetype = Curve,fill = Curve ), alpha=I(1/7), size = 1)
# plot = plot + theme_bw() + scale_linetype_manual(values=c("dotdash","dotted","solid", "longdash", "3313"))
# plot = plot +  scale_colour_manual(values = c("blue",'orange', "#009E73", "red", "purple")) +
#     scale_fill_manual(values =  c("blue",'orange', "#009E73", "red", "purple") ) + 
#     theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
#           legend.position = c(0.1, 0.12), legend.title = element_blank()) + xlab("Scaled time") + ylab("Reversed maginitude") 
# plot



#### We could try on package fdapace 

library(fdapace)
#Ly = train_origin$obsY
#Lt = train_origin$obsT



modelini = LOC_PACE(train_origin, rank, splineObj)



Convert_loc = function(Estimate){
    sig2hat = Estimate$sigma2
    eigenf.loc = Estimate$phi
    eigenv.loc = Estimate$lambda
    eVector = t(eigenf.loc)
    eValues = eigenv.loc
    tSeq = Estimate$workGrid
    eigenFunctions = list()
    eigenFunctionsDeriv1 = list()
    eigenFunctionsDeriv2 = list()
    mean = Estimate$mu
    
    elemLevels = 1
    modelList = list()
    k = 0
    for (e in elemLevels){
        k = k + 1
        fitSeq = approxfun(tSeq, mean)
        modelS = list(tmin = 0, tmax = 1, meanFunction = fitSeq)
        modelList = c(modelList, list(modelS))
    }
    meanModel = list(elemLevels = elemLevels,
                     modelList = modelList)
    
    for (k in 1:length(eValues)){
        compFunctions = list()
        compFunctionsDeriv1 = list()
        compFunctionsDeriv2 = list()
        eFunSeq = Estimate$workGrid
        eFun = approxfun(tSeq, eVector[k, ])
        
        compFunctions = c(compFunctions, list(eFun))
        
        eigenFunctions = c(eigenFunctions, list(compFunctions))
        
    }
    
    model = list(eigenValues = eValues, 
                 eigenFunctions = eigenFunctions,
                 eigenFunctionsDeriv1 = eigenFunctionsDeriv1,
                 eigenFunctionsDeriv2 = eigenFunctionsDeriv2,
                 meanModel = meanModel)
    return(model)
}






#modelini = Convert_loc(Estimate)
#modelini = c(modelini, list(tmin = tmin, tmax = tmax,
#                      SFinal = XInit, sigmaSq = sig2hat,
#                      numPCA = rank, numElem = 1,
#                      elemLevels = levels(train_origin$elemID), BwMu = BwMu, BwCov = BwCov))

#Estimate = FPCA(Ly, Lt, list(methodBwMu = 'CV', methodBwCov = 'CV', useBinnedCov=FALSE, useBW1SE = TRUE, nRegGrid = 100, maxK = rank))
library(fdapace)

sortI = sort(train_origin$obsID, index.return = TRUE)$ix
obsID = train_origin$obsID[sortI]
elemID = train_origin$elemID[sortI]
obsT = train_origin$obsT[sortI]
obsY = train_origin$obsY[sortI]
obsIDCount = as.numeric(table(obsID))

meanSeq = exp(seq(-10,-2, length.out = 7))
mean_mu = MeanModel_GCV(train_origin, splineObj, meanSeq)
#meanModel = fitMeanCurve(train_origin, splineObj, lambda = mean_mu)

meanModel = fitMeanCurve(train_origin, splineObj, lambda = mean_mu)
train_data = subtractMeanCurve(meanModel, train_origin)

rank = 2
#Estimate = FPCA(Ly, Lt, list(userMu = userMu, methodBwCov = 'CV', useBinnedCov=FALSE, maxK = rank))
Estimate = FPCA(Ly, Lt, list(userMu = userMu, maxK = rank))
#Estimate = FPCA(Ly, Lt, list(userBwCov = 0.12, maxK = rank))
sig2hat = Estimate$sigma2
#Estimate$mu
BwMu = Estimate$bwMu
BwCov = Estimate$bwCov
eigenfest = t(Estimate$phi)
#obsMat$elemID  = as.factor(obsMat$elemID)

K = splineObj$getDoF()
basisMat = splineObj$evalSpline(Estimate$workGrid)

UInit = t(eigenfest%*%t(basisMat)/dim(basisMat)[2])
WInit = Estimate$lambda
XInit = list(UInit, diag(WInit))

modelini = Convert_loc(Estimate)
modelini = c(modelini, list(tmin = tmin, tmax = tmax,
                      SFinal = XInit, sigmaSq = sig2hat,
                      numPCA = Estimate$selectK, numElem = 1,
                      elemLevels = levels(train_origin$elemID), BwMu = BwMu, BwCov = BwCov))

meanModel = modelini$meanModel

plotRealEstimate(modelini)
modelini$sigmaSq

DataList = Split_data(test_data, splineObj)
pred_resultini = Estimate_score(modelini, splineObj, DataList, meanModel = meanModel)
MSFE(pred_resultini$pred_YList, DataList$pred_obsYList)


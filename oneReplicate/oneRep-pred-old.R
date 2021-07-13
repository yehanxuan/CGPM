
obsCol = Astro_data(file_list)

nSample = length(unique(obsCol$obsID))
Membership = getCVPartition(nSample, 5)

tmp = (Membership == 1)
train_Index = which(tmp == FALSE)
test_Index = which(tmp == TRUE)
### 80% training data and 20% testing data
train_data = obsCol[obsCol$obsID %in% train_Index, ] 
test_data = obsCol[obsCol$obsID %in% test_Index, ]
train_origin = train_data


if ( (method == "LogDet") || (method =="VNDiv") ||(method == "frobDiverg") ){
meanSeq = exp(seq(-10,-2, length.out = 7))
mean_mu = MeanModel_GCV(train_data, splineObj, meanSeq)
meanModel = fitMeanCurve(train_data, splineObj, lambda = mean_mu)
train_data = subtractMeanCurve(meanModel, train_data)

modelini = ini_estimate(train_data, rank, splineObj, ini.method = "EM", M.set)

} else if ( (method == "loc")){
#meanModel = loc_meanModel(train_data)
#meanModel = NULL
#modelini = ini_estimate(train_data, rank, splineObj, ini.method = "loc", M.set)
modelini = LOC_PACE(train_origin, rank, splineObj)
meanModel = modelini$meanModel
} else if ((method == "Newton") ){
#meanModel = loc_meanModel(train_data)    
modelini = ini_estimate(train_data, rank, splineObj, ini.method = "EM", M.set) 
#meanModel = modelini$meanModel
} else if ( method =="EM"){
modelini =  ini_estimate(train_data, rank, splineObj, ini.method = "EM", M.set)
meanModel = modelini$meanModel
}

DataList = Split_data(test_data, splineObj)

#modelini = EM_estimate(train_data, rank, splineObj, ini.method)
#modelini = ini_estimate(train_data, rank, splineObj, ini.method, M.set)
XInit = modelini$SFinal
sig2hat = modelini$sigmaSq



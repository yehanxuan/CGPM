
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

meanSeq = exp(seq(-10,-2, length.out = 7))
mean_mu = MeanModel_GCV(train_data, splineObj, meanSeq)
meanModel = fitMeanCurve(train_data, splineObj, lambda = mean_mu)
train_data = subtractMeanCurve(meanModel, train_data)

sig2hat = 1
sig.EM = 1



DataList = Split_data(test_data, splineObj)

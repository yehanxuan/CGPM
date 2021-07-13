#### Boxplot #####

Sum_Loss = function(MSE, comp){
    loss = lapply(MSE, function(x) mean(x[1:comp]))
    return(loss)
}


library(ggplot2)
source("utils/evaluate.R")
# create a data frame
variaty=rep(LETTERS[1:7], each=40)
treatment=rep(c("high","low"),each=20)
note=seq(1:280)+sample(1:150, 280, replace=T)
data=data.frame(variety, treatment ,  note)
####################
# grouped boxplot
ggplot(data, aes(x=variaty, y=note, fill=treatment)) + 
    geom_boxplot()

methods = rep(c("REML", "LogDet", "VNDiv", "Frob"), each = 3*maxIter)
Init = rep(c("EM", "LS","LOC"), each = maxIter)
AveLoss = unlist(Sum_Loss(Final$MSE, comp = 4))

CreateBoxDataFrame = function(methodsName, InitName, AveLoss){
    L = length(AveLoss)
    methods = rep(methodsName, L)
    Init = rep(InitName, L)
    Frame = data.frame(methods, Init, AveLoss)
    return(Frame)
}


load("Boxdata.Rdata")
#### REML 

#load("dataServer/method-REML-EM-FourierOrth-200-Gaussian.RData")
#load("dataServer/method-REML-EM-FourierOrth-200-uniform.RData")
#load("dataServer/method-REML-EM-FourierOrth-200-t.RData")

#load("dataServer/method-REML-Fourier-200-Gaussian-EM.RData")
#load("dataServer/method-REML-Fourier-200-uniform-EM.RData")

load("dataServer/method-REML-Fourier-200-t-EM.RData")

#load("dataServer/method-REML-EM-Fourier-200-t-EM.RData")
compute(Final$MSE)
AveLoss_REML_EM = unlist(Sum_Loss(Final$MSE, 4))
Frame_REML_EM = CreateBoxDataFrame("REML", "EM", AveLoss_REML_EM)
#load("data/method-REML-loc-FourierOrth-200-Gaussian.RData")
#load("dataServer/method-REML-loc-FourierOrth-200.RData")
#load("dataServer/method-REML-loc-FourierOrth-200-uniform.RData")
#load("dataServer/method-REML-loc-FourierOrth-200-t.RData")

#load("dataServer/method-REML-Fourier-200-Gaussian-LOC.RData")
#load("dataServer/method-REML-Fourier-200-uniform-LOC.RData")
load("dataServer/method-REML-Fourier-200-t-LOC.RData")

compute(Final$MSE)
AveLoss_REML_LOC = unlist(Sum_Loss(Final$MSE, 4))
Frame_REML_LOC = CreateBoxDataFrame("REML", "LOC", AveLoss_REML_LOC)

#load("dataServer/method-REML-LS-FourierOrth-200.RData")
#load("dataServer/method-REML-LS-FourierOrth-200-uniform.RData")
#load("dataServer/method-REML-LS-FourierOrth-200-t.RData")

#load("dataServer/method-REML-Fourier-200-Gaussian-LS.RData")
#load("dataServer/method-REML-Fourier-200-uniform-LS.RData")
load("dataServer/method-REML-Fourier-200-t-LS.RData")

compute(Final$MSE)
AveLoss_REML_LS = unlist(Sum_Loss(Final$MSE, 4))
Frame_REML_LS = CreateBoxDataFrame("REML", "LS", AveLoss_REML_LS)

Frame_REML = rbind(Frame_REML_EM, Frame_REML_LOC, Frame_REML_LS)
#### LogDet 
#load("data/method-LogDet-FourierOrth-200-18-EM-Gaussian.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-EM-Gaussian.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-EM-uniform.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-EM-t.RData")

#load("dataServer/method-LogDet-Fourier-200-18-EM-Gaussian.RData")
load("dataServer/method-LogDet-Fourier-200-18-EM-t.RData")
#load("dataServer/method-LogDet-Fourier-200-18-EM-uniform.RData")

compute(Final$MSE)
AveLoss_LogDet_EM = unlist(Sum_Loss(Final$MSE, 4))
Frame_LogDet_EM = CreateBoxDataFrame("LogDet", "EM", AveLoss_LogDet_EM)

#load("data/method-LogDet-FourierOrth-200-18-LOC-Gaussian.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-LOC-Gaussian.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-LOC-uniform.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-LOC-t.RData")

#load("dataServer/method-LogDet-Fourier-200-18-LOC-Gaussian.RData")
load("dataServer/method-LogDet-Fourier-200-18-LOC-t.RData")
#load("dataServer/method-LogDet-Fourier-200-18-LOC-uniform.RData")

compute(Final$MSE)
AveLoss_LogDet_LOC = unlist(Sum_Loss(Final$MSE, 4))
Frame_LogDet_LOC = CreateBoxDataFrame("LogDet", "LOC", AveLoss_LogDet_LOC)

#load("data/method-LogDet-FourierOrth-200-18-LS-Gaussian.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-LS-Gaussian.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-LS-uniform.RData")
#load("dataServer/method-LogDet-FourierOrth-200-18-EM-t.RData")

#load("dataServer/method-LogDet-Fourier-200-18-LS-Gaussian.RData")
load("dataServer/method-LogDet-Fourier-200-18-LS-t.RData")
#load("dataServer/method-LogDet-Fourier-200-18-LS-uniform.RData")

compute(Final$MSE)
AveLoss_LogDet_LS = unlist(Sum_Loss(Final$MSE, 4))
Frame_LogDet_LS = CreateBoxDataFrame("LogDet", "LS", AveLoss_LogDet_LS)

Frame_LogDet = rbind(Frame_LogDet_EM, Frame_LogDet_LOC, Frame_LogDet_LS)
### VNDiv
#load("data/method-VNDiv-FourierOrth-200-18-EM-Gaussian.RData")
#load("dataServer/method-VNDiv-FourierOrth-200-18-EM-Gaussian.RData")
#load("dataServer/method-VNDiv-FourierOrth-200-18-EM-uniform.RData")
#load("dataServer/method-VNDiv-FourierOrth-200-18-EM-t.RData")

#load("dataServer/method-VNDiv-Fourier-200-18-EM-uniform.RData")
load("dataServer/method-VNDiv-Fourier-200-18-EM-t.RData")
#load("dataServer/method-VNDiv-Fourier-200-18-EM-Gaussian.RData")
compute(Final$MSE)
AveLoss_VN_EM = unlist(Sum_Loss(Final$MSE, 4))
Frame_VN_EM = CreateBoxDataFrame("VN", "EM", AveLoss_VN_EM)

#load("data/method-VNDiv-FourierOrth-200-18-LOC-Gaussian.RData")
#load("dataServer/method-VNDiv-FourierOrth-200-18-LOC-Gaussian.RData")
#load("dataServer/method-VNDiv-FourierOrth-200-18-LOC-uniform.RData")
#load("dataServer/method-VNDiv-FourierOrth-200-18-LOC-t.RData")


#load("dataServer/method-VNDiv-Fourier-200-18-LS-Gaussian.RData")
#load("dataServer/method-VNDiv-Fourier-200-18-LOC-uniform.RData")
load("dataServer/method-VNDiv-Fourier-200-18-LOC-t.RData")
compute(Final$MSE)
AveLoss_VN_LOC = unlist(Sum_Loss(Final$MSE, 4))
Frame_VN_LOC = CreateBoxDataFrame("VN", "LOC", AveLoss_VN_LOC)

#load("dataServer/method-VNDiv-FourierOrth-200-18-LS-Gaussian.RData")
#load("dataServer/method-VNDiv-FourierOrth-200-18-LS-uniform.RData")
#load("dataServer/method-VNDiv-FourierOrth-200-18-LS-t.RData")



#load("dataServer/method-VNDiv-Fourier-200-18-LS-Gaussian.RData")
#load("dataServer/method-VNDiv-Fourier-200-18-LS-uniform.RData")
load("dataServer/method-VNDiv-Fourier-200-18-LS-t.RData")
compute(Final$MSE)
AveLoss_VN_LS = unlist(Sum_Loss(Final$MSE, 4))
Frame_VN_LS = CreateBoxDataFrame("VN", "LS", AveLoss_VN_LS)

Frame_VN = rbind(Frame_VN_EM, Frame_VN_LOC, Frame_VN_LS)
### Frob
#load("data/method-frobDiverg-FourierOrth-200-18-EM-Gaussian.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-EM-Gaussian.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-EM-uniform.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-EM-t.RData")

#load("dataServer/method-frobDiverg-Fourier-200-18-EM-Gaussian.RData")
#load("dataServer/method-frobDiverg-Fourier-200-18-EM-uniform.RData")
load("dataServer/method-frobDiverg-Fourier-200-18-EM-t.RData")
compute(Final$MSE)
AveLoss_Frob_EM = unlist(Sum_Loss(Final$MSE, 4))
Frame_Frob_EM = CreateBoxDataFrame("Frob", "EM", AveLoss_Frob_EM)
#load("data/method-frobDiverg-FourierOrth-200-18-LOC-Gaussian.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-LOC-Gaussian.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-LOC-uniform.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-LOC-t.RData")

#load("dataServer/method-frobDiverg-Fourier-200-18-LS-Gaussian.RData")
#load("dataServer/method-frobDiverg-Fourier-200-18-LOC-uniform.RData")
load("dataServer/method-frobDiverg-Fourier-200-18-LOC-t.RData")
compute(Final$MSE)
AveLoss_Frob_LOC = unlist(Sum_Loss(Final$MSE, 4))
Frame_Frob_LOC = CreateBoxDataFrame("Frob", "LOC", AveLoss_Frob_LOC)

#load("data/method-frobDiverg-FourierOrth-200-18-LS-Gaussian.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-LS-Gaussian.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-LS-uniform.RData")
#load("dataServer/method-frobDiverg-FourierOrth-200-18-LS-t.RData")

#load("dataServer/method-frobDiverg-Fourier-200-18-LS-Gaussian.RData")
#load("dataServer/method-frobDiverg-Fourier-200-18-LS-uniform.RData")
load("dataServer/method-frobDiverg-Fourier-200-18-LS-t.RData")
compute(Final$MSE)
AveLoss_Frob_LS = unlist(Sum_Loss(Final$MSE , 4))
Frame_Frob_LS = CreateBoxDataFrame("Frob", "LS", AveLoss_Frob_LS)

Frame_Frob = rbind(Frame_Frob_EM, Frame_Frob_LOC, Frame_Frob_LS)
#Boxdata = data.frame(methods, Init, AveLoss)


### load EM data 
#load("dataServer/method-EM-FourierOrth-200-Gaussian.RData")
#load("dataServer/method-EM-FourierOrth-200-uniform.RData")
#load("dataServer/method-EM-FourierOrth-200-t.RData")

#load("dataServer/method-EM-Fourier-200-Gaussian.RData")
#load("dataServer/method-EM-Fourier-200-uniform.RData")
load("dataServer/method-EM-Fourier-200-t.RData")
compute(Final$MSE)
AveLoss_EM = unlist(Sum_Loss(Final$MSE, 4))
Frame_EM = CreateBoxDataFrame("EM", "LS", AveLoss_EM)


#### Note

Boxdata_Gaussian = rbind(Frame_EM, Frame_REML, Frame_LogDet, Frame_VN, Frame_Frob)
Boxdata_Gaussian$DataType = "Gaussian"








Boxdata_t = rbind(Frame_EM, Frame_REML, Frame_LogDet, Frame_VN, Frame_Frob) 
Boxdata_t$DataType = "t"


Boxdata_unif = rbind(Frame_EM, Frame_REML, Frame_LogDet, Frame_VN, Frame_Frob)
Boxdata_unif$DataType = "uniform"

#Boxdata = rbind(Frame_EM, Frame_REML, Frame_LogDet, Frame_VN, Frame_Frob) 

Boxdata = rbind(Boxdata_Gaussian, Boxdata_t, Boxdata_unif)

colnames(Boxdata)[2] <- "Initial"

save(Boxdata, file = "Boxdata.Rdata")
#Boxdata = rbind( Frame_REML, Frame_LogDet, Frame_VN, Frame_Frob) 
    
ggplot(Boxdata, aes(x = methods, y = AveLoss, fill=Initial)) + geom_boxplot() + 
    ylab("Average error") + theme(legend.key.size = unit(1, 'cm')) + ylim(0.1, 1.65) + 
    facet_wrap(~DataType) + theme(legend.position="top") 
    

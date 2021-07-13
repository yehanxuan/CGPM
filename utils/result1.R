for(i in 1:4){
    if(i==1) cat("PE")
    if(i==3) cat("\\# of selected predictors")
    for(j in 1:6){
        if(i==1 || i==3){
            cat(" &",1,sep= "")
            
        }else{
            cat(" & (",0.1 ,")",sep= "")
            
        }
    }
    cat("\\\\ \n ")
}


cat("\\hline\n")
cat("\\multicolumn{7}{l}{$n =", nn, "\\quad p = ", pp ,
    "$}  \\\\\n\\hline\n", sep = "")


allMethod = c("REML", "EM", "LOC", "LogDet", "Frob", "VN")
finalMean = matrix(0, 6, length(allMethod))
finalSD = matrix(0, 6, length(allMethod))
colSel = 1

for (m in 1:length(allMethod)){
    fileList = "./dataServer/method-EM-FourierOrth-200-Gaussian.RData"
    fileList = "./dataServer/method-LogDet-FourierOrth-200-18-EM-Gaussian.RData"
    load(fileList)
    
    resultT = compute(Final$MSE)
    res = resultT$res
    resSD = resultT$resSD
    
    
    finalMean[, colSel] = res
    finalSD[, colSel] = resSD
    colSel = colSel + 1
}


allMethod = c("LogDet", "frobDiverg", "VNDiv")
DataType = "Fourier"
samplesize = 500
scoreType = "Gaussian"
InitType = "LS"
PenaltySeq = c("N", "Y")
nKnots = 35


load("penaltyPower/method-frobDiverg-Fourier-200-30-Gaussian-LS-N.RData")
collectResults = function(samplesize, DataType, allMethod){
    colSel = 1
    finalMean = matrix(0, 8, 6)
    finalSD = matrix(0, 8, 6)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        for (p in 1:length(PenaltySeq)){
            Penalty = PenaltySeq[p]
            fileList = paste0("./penaltyPower/method-", method, "-", DataType, "-", samplesize,"-", nKnots, "-", scoreType,
                              "-", InitType, "-", Penalty, ".RData")
            load(fileList)
            resultT = compute(Final$MSE)
            res = resultT$res
            resSD = resultT$resSD
            finalMean[, colSel] = res
            finalSD[, colSel] = resSD
            colSel = colSel + 1
        }
    }
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



#### Output of Power 
str = '\\begin{table}
\\centering
\\begin{tabular}{c|rrrr}
\\hline
\\multicolumn{2}{r}{0.8} & 0.9 & 1.1 & 1.2 \\\\
\\hline \n'

writeLines(str)
####

str = '
\\hline
\\end{tabular}	
\\caption{Simulation results the Power method under different choice of $\\alpha$ on the first setting. The four subtables correspond to various sample sizes and number of covariates. The numbers in the parentheses are the standard errors. The columns correspond to the estimation methods.}	
\\end{table}'

writeLines(str)

bSeq = c("0.8", "0.9", "1.1", "1.2")
DataType = "Fourier"
samplesize = 200
nn = 200
scoreType
InitType = "LS"
collectResults_Power = function(samplesize, DataType, bSeq){
    colSel = 1
    finalMean = matrix(0, 8, length(bSeq))
    finalSD = matrix(0, 8, length(bSeq))
    for (m in 1:length(bSeq)){
        b = bSeq[m]
        fileList = paste0("./PowerData/method-", "Power" , "-", DataType, "-", samplesize,"-", "18", "-", scoreType,
                          "-", InitType, "-", b, ".RData")
        load(fileList)
        resultT = compute(Final$MSE)
        res = resultT$res
        resSD = resultT$resSD
        finalMean[, colSel] = res
        finalSD[, colSel] = resSD
        colSel = colSel + 1
    }
}



outputNewPower = function(finalMean, finalSD, nn, InitType){
    rk = 3
    cat("\\hline\n")
    cat("\\multicolumn{5}{l}{$n =", nn, "\\quad ",
        "$", InitType,"}",  "\\\\\n\\hline\n", sep = "")
    nrow(finalMean)
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



#### The role of penalty 
allMethod = c("frobDiverg", "LogDet", "VNDiv") 
samplesize = 500
nKnots = 35
InitType
PenaltySeq

collectResults = function(samplesize, DataType, allMethod){
    colSel = 1
    finalMean = matrix(0, 8, 6)
    finalSD = matrix(0, 8, 6)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        for (p in 1:length(PenaltySeq)){
            Penalty = PenaltySeq[p]
            fileList = paste0("./penaltyPower/method-", method, "-", DataType, "-", samplesize,"-", nKnots, "-", scoreType,
                              "-", InitType, "-", Penalty, ".RData")
            load(fileList)
            resultT = compute(Final$MSE)
            res = resultT$res
            resSD = resultT$resSD
            finalMean[, colSel] = res
            finalSD[, colSel] = resSD
            colSel = colSel + 1
        }
    }
}



#### Load Fourier data 
# REML
load("dataServer/method-REML-Fourier-200-Gaussian-LS.RData")
load("dataServer/method-REML-Fourier-200-Gaussian-LOC.RData")
# EM 
load("dataServer/method-EM-Fourier-200-t.RData")
load("dataServer/method-EM-Fourier-200-uniform.RData")
#LOC
load("dataServer/method-LOC-Fourier-200-t.RData")
## LogDet 
load("dataServer/method-LogDet-Fourier-200-18-LS-Gaussian.RData")
## frobDiverg
load("dataServer/method-frobDiverg-Fourier-200-18-LS-Gaussian.RData")
##
load("dataServer/method-VNDiv-Fourier-200-18-LS-Gaussian.RData")

allMethod = c("REML", "EM", "LOC", "frobDiverg", "LogDet", "VNDiv")

allMethod = c("LOC", "EM", "REML", "LogDet", "VNDiv", "frobDiverg")
InitType = "EM"
samplesize = 500 
nn = 500
scoreType = "Gaussian"
DataType
collectResults_Fourier = function(samplesize, DataType, allMethod){
    colSel = 1
    finalMean = matrix(0, 8, 6)
    finalSD = matrix(0, 8, 6)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        #if (m == 1){
        #    fileList = paste0("./dataServer/method-", method, "-", DataType, "-", samplesize, "-", scoreType,
        #                      "-", InitType, ".RData")
        #} else if (m %in% c(2,3)){
        #    fileList = paste0("./dataServer/method-", method, "-", DataType, "-", samplesize, "-", scoreType,
        #                     ".RData")
        #} else {
        #    fileList = paste0("./dataServer/method-", method, "-", DataType, "-", samplesize,"-", 18, "-", InitType,"-", scoreType,
        #                      ".RData")
        #}
        try({
            if ((m == 1) | (m == 2)){
                fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
                                  scoreType, "-", "EM", ".RData")
            } else {
                
                fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
                                  scoreType, "-", InitType, ".RData")
                
                #fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
                #                 scoreType, "-", InitType, ".RData")
            }
        })
        
        
        #fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
        #                 scoreType, "-", InitType, ".RData")
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
}

allMethod = c("REML", "LogDet", "VNDiv", "frobDiverg")
collectSubmatrix = function(){
    colSel = 1
    finalMean2 = matrix(0, 8, 4)
    finalSD2 = matrix(0, 8, 4)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        #if (m == 1){
        #    fileList = paste0("./dataServer/method-", method, "-", DataType, "-", samplesize, "-", scoreType,
        #                      "-", InitType, ".RData")
        #} else if (m %in% c(2,3)){
        #    fileList = paste0("./dataServer/method-", method, "-", DataType, "-", samplesize, "-", scoreType,
        #                     ".RData")
        #} else {
        #    fileList = paste0("./dataServer/method-", method, "-", DataType, "-", samplesize,"-", 18, "-", InitType,"-", scoreType,
        #                      ".RData")
        #}
        try({
            fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
                              scoreType, "-", "LS", ".RData")
        })
        
        
        #fileList = paste0("./NewFourier/method-", method, "-", DataType, "-", samplesize, "-",
        #                 scoreType, "-", InitType, ".RData")
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
}


finalMean_comb = cbind(finalMean, finalMean2)
finalSD_comb = cbind(finalSD, finalSD2)

outputNewFourier = function(finalMean, finalSD, nn, InitType){
    rk = 3
    cat("\\hline\n")
    cat("\\multicolumn{7}{l}{$n =", nn, "\\quad ",
        "$", InitType,"}",  "\\\\\n\\hline\n", sep = "")
    nrow(finalMean)
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

outputNewFourier(finalMean_comb, finalSD_comb, 200, "LS")




#######  Outout the table of Prediction performance 
allMethod = c("REML", "EM", "LOC", "frobDiverg", "LogDet", "VNDiv")
load("predata/method-frobDiverg-3-LS.RData")
compSeq = c(2,3,4)
collectResults_Pred = function(){
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
        #colSel = colSel + 1
    }
}

outputNew_pred = function(){
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

######### Simu_Peng and Paul 

source("./utils/evaluate.R")
allMethod = c("LOC", "EM", "REML", "LogDet", "VNDiv", "frobDiverg")
allMethod = c("REML", "LogDet", "VNDiv", "frobDiverg")
InitType = "EM"
samplesize = 100
nn = 100
scoreType = "Gaussian"
DataType = "prac"

#nKnots = 15

collectResults_Paul = function(samplesize, DataType, allMethod, scoreType, InitType){
    colSel = 1
    finalMean = matrix(0, 5, 6)
    finalSD = matrix(0, 5, 6)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        #try({
            if ((m == 1) | (m == 2)){
                fileList = paste0("./Pauldata/method-", method, "-", DataType, "-", samplesize, "-",
                                  scoreType, "-", "EM", ".RData")
            } else if ((m == 5) | (m == 6)){
                fileList = paste0("./Pauldata/method-", method, "-", DataType, "-", samplesize, "-",
                                  scoreType, "-", InitType, "-", 15, ".RData")
            } else {
                
                fileList = paste0("./Pauldata/method-", method, "-", DataType, "-", samplesize, "-",
                                  scoreType, "-", InitType, ".RData")
            }
      #  })
        
        
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

EMdata = collectResults_Paul(samplesize, DataType, allMethod, scoreType, InitType = "EM")
finalMean = EMdata[[1]]
finalSD = EMdata[[2]]

allMethod = c("REML", "LogDet", "VNDiv", "frobDiverg")
simu_nknots = 8
DataType = "easy"
samplesize = 50
collectSubmatrix = function(samplesize, DataType, allMethod, scoreType, InitType, simu_nknots){
    colSel = 1
    finalMean2 = matrix(0, 3, 4)
    finalSD2 = matrix(0, 3, 4)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
       # try({
        if (m == 1){
            fileList = paste0("./Pauldata/method-", method, "-", DataType, "-", samplesize, "-",
                              scoreType, "-", InitType, ".RData")
        } else { 
            fileList = paste0("./Pauldata/method-", method, "-", DataType, "-", samplesize, "-",
                                                       scoreType, "-", InitType, "-", simu_nknots, ".RData")
        }
            
        
        #    if ((m == 3) | (m == 4) ){
        #        fileList = paste0("./Pauldata/method-", method, "-", DataType, "-", samplesize, "-",
        #                          scoreType, "-", InitType, "-", 15, ".RData")
        #    } else {
        #        fileList = paste0("./Pauldata/method-", method, "-", DataType, "-", samplesize, "-",
        #                          scoreType, "-", InitType, ".RData")
        #    }
            
            
            
     #   })
        
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

EMdata = collectSubmatrix(samplesize, DataType, allMethod, scoreType, InitType = "EM", simu_nknots = 8)
finalMean = EMdata[[1]]
finalSD = EMdata[[2]]

LOCdata = collectSubmatrix(samplesize, DataType, allMethod, scoreType, InitType = "LOC", simu_nknots = 8)
finalMean2 = LOCdata[[1]]
finalSD2 = LOCdata[[2]]

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


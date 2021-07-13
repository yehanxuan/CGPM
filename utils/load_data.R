###################
##Data Processing##
###################

# Process the real data
myFunctionForApply <- function(x, ...) {
    # Do your processing
    # Let's say it ends up in variable 'ret':
    if (length(x) == 0)
        return(NA)
    return(x)
}



# Produce the astronomical data 
Astro_data = function(file_list){
    numElem = 1
    obsID = 1
    obsMat = data.frame()
    for (i in 1:length(file_list)){
        Table = read_tsv( paste0("./CSP_Photometry_DR2/", file_list[i]), skip = 4, na = "99.900")
        colnames(Table)[1] = "MJD"
        
        obsT = Table$MJD
        obsY_u = Table$u
        obsY_B = Table$B
        obsY_V = Table$V
        obsY_g = Table$g
        obsY_r = Table$r
        
        #obsY_i = Table$i
        
        ## Clean the NA data
        tmp1 = cbind(obsT[!is.na(obsY_u)], obsY_u[!is.na(obsY_u)])
        tmp2 = cbind(obsT[!is.na(obsY_B)], obsY_B[!is.na(obsY_B)])
        tmp3 = cbind(obsT[!is.na(obsY_V)], obsY_V[!is.na(obsY_V)])
        tmp4 = cbind(obsT[!is.na(obsY_g)], obsY_g[!is.na(obsY_g)])
        tmp5 = cbind(obsT[!is.na(obsY_r)], obsY_r[!is.na(obsY_r)])
        
        
        Ind_min1 = myFunctionForApply( which.min(tmp1[, 2]) )
        Ind_min2 = myFunctionForApply( which.min(tmp2[, 2]) )
        Ind_min3 = myFunctionForApply( which.min(tmp3[, 2]) )
        Ind_min4 = myFunctionForApply( which.min(tmp4[, 2]) )
        Ind_min5 = myFunctionForApply( which.min(tmp5[, 2]) )
        
        Ind_min_vec = c(Ind_min1, Ind_min2, Ind_min3, Ind_min4, Ind_min5)
        
        
        if (Ind_min2 > 1){
            maxpoint = tmp2[Ind_min2, ]
            for (i in 1:5){
                assign( paste0("tmp", i), sweep(get(paste0("tmp", i)), 2, maxpoint) )
                #get(paste0("tmp", i))[ ,2:3] = get(paste0("tmp", i))[ ,2:3] - maxpoint
            }
            
            tmp1 = tmp1[ (tmp1[, 1] > -10)&(tmp1[, 1] < 50), ]
            tmp2 = tmp2[ (tmp2[, 1] > -10)&(tmp2[, 1] < 50), ]
            tmp3 = tmp3[ (tmp3[, 1] > -10)&(tmp3[, 1] < 50), ]
            tmp4 = tmp4[ (tmp4[, 1] > -10)&(tmp4[, 1] < 50), ]
            tmp5 = tmp5[ (tmp5[, 1] > -10)&(tmp5[, 1] < 50), ]
            for (i in 1:5){
                # Need to use & 
                if ( (Ind_min_vec[i] > 1) & (!is.na(Ind_min_vec[i]) ) ){
                    Tmp = cbind(obsID, numElem, get(paste0("tmp", i)) )
                    obsMat = rbind(obsMat, Tmp)
                    obsID = obsID + 1
                }
            }
        }
    }
    
    colnames(obsMat) =c("obsID", "elemID", "obsT", "obsY")
    
    nmax = max(table(obsMat[ ,1]))
    L1 = min(as.numeric(obsMat[ ,3]))
    L2 = max(as.numeric(obsMat[ ,3]))
    tmp = (obsMat[ ,3] - L1)/(L2 - L1)
    tmp[tmp<0.00001]<-0.00001
    tmp[tmp>0.99999]<-0.99999 
    
    ##### 不做scaling怎么样
    tmp2 = as.numeric(obsMat[,4])
 #   S1 = min(as.numeric(obsMat[ ,4]))
  #  S2 = max(as.numeric(obsMat[ ,4]))
  #  tmp2 = (obsMat[ ,4] - S1)/(S2 - S1)
  #  tmp2[tmp2<0.00001]<-0.00001
  #  tmp2[tmp2>0.99999]<-0.99999 
  #  tmp2 = tmp2 * 5 - 2.5
        #tmp2 = tmp2 - 0.5
    
    obsCol = obsMat
    obsCol[ ,3] = tmp
    obsCol[ ,4] = -tmp2
    colnames(obsCol) =c("obsID", "elemID", "obsT", "obsY")
    obsCol$elemID = as.factor(obsCol$elemID)
    return(obsCol)
}

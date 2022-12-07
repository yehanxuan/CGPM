#### CV Function#### Get CV partition
getCVPartition = function(nSample, nFold){
    nPerFold = ceiling(nSample / nFold)
    cvMembership = rep(1:nFold, nPerFold)
    shuffleV = sample(1:nSample, size = nSample, replace = FALSE)
    cvMembership = cvMembership[shuffleV]
    return(cvMembership)
}

getCVPartition_seed = function(nSample, nFold, seed = 77){
    nPerFold = ceiling(nSample / nFold)
    cvMembership = rep(1:nFold, nPerFold)
    old <- .Random.seed
    set.seed(seed)
    shuffleV = sample(1:nSample, size = nSample, replace = FALSE)
    cvMembership = cvMembership[shuffleV]
    .Random.seed <<- old
    return(cvMembership)
}

numElem = 1
obsNumLower = 5
obsNumUpper = 10

if (DataType == "Fourier"){
    #b = 0.9
    #b = 1.1
    fExpNum = 5
    pcaCompNum = fExpNum * 2
    tmin = 0
    tmax = 1
    pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
    pcaTrans = qr.Q(qr(pcaTrans))
    splitP = c(-1e-9, 2+1e-9)
    pcaKappa = exp(-(0:(pcaCompNum-1))/2)*10
    pcaKappaSqrt = sqrt(pcaKappa)
    eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
    obsCol = get_allObs(eigenfList, samplesize, pcaKappaSqrt, noiseSigma = 1, obsNumLower, obsNumUpper, scoreType)
    nKnots = 18
    M.EM = nKnots + 2
    r.set = 8
    splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
    sig2hat = 1
}
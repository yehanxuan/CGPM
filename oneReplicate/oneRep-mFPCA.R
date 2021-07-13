


fExpNum = 4
pcaCompNum = fExpNum * 2
pcaTrans = matrix(rnorm(pcaCompNum^2), pcaCompNum, pcaCompNum)
pcaTrans = qr.Q(qr(pcaTrans))
splitP = c(-1e-9, 2+1e-9)

numSample = 2e3
numElem = 1
obsNumLower = 5
obsNumUpper = 10
pcaKappa = exp(-(0:(pcaCompNum-1))/2)*10
pcaKappaSqrt = sqrt(pcaKappa)
noiseSigma = 0.1
eigenfList = get_eigenfunList(pcaTrans, fExpNum, splitP)
obsCol = get_allObs(eigenfList, numSample, pcaKappaSqrt, noiseSigma, 
                    obsNumLower, obsNumUpper)

### Define the spline objective 

splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
K = splineObj$getDoF()

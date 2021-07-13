tmin = 0
tmax = 1
mOrder = 4
#nKnots = 10
#rankSeq = 8
#optRank = 5
#splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
#K = splineObj$getDoF()
nFold = 10
lambdaSeq = exp(seq(-22,-9, length.out = 5))
#lambdaSeq = 0
#lambdaSeq = 1e-10
#muSeq = seq(0,0.01, length.out = 5)
#M.set = c(15:30)

controlList1 = list(alpha = 0.618, tol = 1e-6, iterMax = 200, sigma = 0.05,verbose = 0)
#controlList2 = list(alpha = 1e-1, tol = 1e-7, sigma = 1e-3, beta = 0.618,
#                    iterMax = 500, verbose = 0)

controlList2 = list(alpha = 0.1, tol = 1e-5, sigma = 5e-2, beta = 0.618,
                    iterMax = 2000, verbose = 1)









tmin = 0
tmax = 1
mOrder = 4
#nKnots = 10
#rankSeq = rev(3:6)
## fix the rank to be 8
#rankSeq = 8
nFold = 10
lambdaSeq = exp(seq(-22,-12, length.out = 7))
#lambdaSeq = 0
#lambdaSeq = 1e-7

# M.set = c(15:25)
controlList1 = list(alpha = 1, tol = 1e-6, iterMax = 200, sigma = 0.05,verbose = 0)
#controlList2 = list(alpha = 1e-1, tol = 1e-7, sigma = 1e-3, beta = 0.618,
#                    iterMax = 500, verbose = 0)

controlList2 = list(alpha = 0.1, tol = 1e-5, sigma = 5e-2, beta = 0.618,
                    iterMax = 1000, verbose = 1)
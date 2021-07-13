tmin = 0
tmax = 1
mOrder = 4
nFold = 10
lambdaSeq = exp(seq(-22,-12, length.out = 5))

controlList1 = list(alpha = 0.618, tol = 1e-6, iterMax = 200, sigma = 0.05,verbose = 0)
controlList2 = list(alpha = 0.1, tol = 1e-5, sigma = 5e-2, beta = 0.618,
                    iterMax = 2000, verbose = 1)

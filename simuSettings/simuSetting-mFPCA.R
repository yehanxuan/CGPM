



tmin = 0
tmax = 1
rankSeq = 8
lambdaSeq = rev(exp(seq(0,2, length.out = 2)))
mOrder = 4
nKnots = 10
nFold = 10

controlList1 = list(alpha = 2, tol = 1e-7, iterMax = 10, sigma = 0.1,verbose = 0)
controlList2 = list(alpha = 0.1, tol = 1e-7, sigma = 1e-4, beta = 0.618,
                    iterMax = 500, verbose = 1)





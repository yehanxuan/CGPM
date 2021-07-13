tmin = 0
tmax = 1
mOrder = 4
nKnots = 30
nFold = 10
### Estimate mean function
splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)

# rank = 4

controlList1 = list(alpha = 0.618, tol = 1e-6, iterMax = 200, sigma = 0.05,verbose = 0)
#controlList2 = list(alpha = 1e-1, tol = 1e-7, sigma = 1e-3, beta = 0.618,
#                    iterMax = 500, verbose = 0)

controlList2 = list(alpha = 0.1, tol = 1e-5, sigma = 5e-2, beta = 0.618,
                    iterMax = 2000, verbose = 0)

#lambdaSeq = exp(seq(-17,-10, length.out = 7))
lambdaSeq = exp(seq(-20,-12, length.out = 7))

#ini.method = "EM"
basis.method = "bs"
sl.v = rep(0.5, 10)
max.step = 50
grid.l = seq(0,1,0.01)
grids= seq(0,1,0.002)
sig.EM = 1

M.set = c(15:30)

tmin = 0
tmax = 1
mOrder = 4
nKnots = 20
nFold = 10


# rank = 4

controlList1 = list(alpha = 0.618, tol = 1e-6, iterMax = 200, sigma = 0.05,verbose = 0)


controlList2 = list(alpha = 0.1, tol = 1e-5, sigma = 5e-2, beta = 0.618,
                    iterMax = 2000, verbose = 0)

lambdaSeq = exp(seq(-20,-12, length.out = 7))

#ini.method = "EM"
basis.method = "bs"
sl.v = rep(0.5, 10)
max.step = 50
grid.l = seq(0,1,0.01)
grids= seq(0,1,0.002)
sig.EM = 1

#M.set = c(10:20)
M.set = 15
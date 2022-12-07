numElem = 1
obsNumLower = 5
obsNumUpper = 10
if (DataType == "easy"){
  tmin = 0
  tmax = 1
  eigenfList = get_eigenfunList_Paul("easy")
  obsCol = get_allObs_Paul(samplesize, M=5,r=3, "easy", scoreType, alpha = 0.6, nmin=2, nmax=10, a=1, b=1, sig=1/4, noiseType)
  nKnots = 10
  M.EM = nKnots + 2
  select_Method = "knots"
  M.set = seq(4, 12, by = 1)  # select by knots without penalty
  r.set = 3
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1/16
} else if (DataType == "prac"){
  tmin = 0
  tmax = 1
  eigenfList = get_eigenfunList_Paul("prac")
  obsCol = get_allObs_Paul(samplesize, M=10, r=5, "prac", scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4, noiseType)
  #nKnots = 15
  nKnots = 8
  M.set = seq(8, 15, by = 1)
  M.EM = nKnots + 2
  select_Method = "knots"
  r.set = 5
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1/16
}  else if (DataType == "pracNU") {
  tmin = 0
  tmax = 1
  shape1 = 2
  shape2 = 2
  eigenfList = get_eigenfunList_Paul_NU("pracNU", shape1, shape2)
  obsCol = get_allObs_Paul(samplesize, M=10, r=5, "pracNU", scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4, noiseType, shape1, shape2)
  nKnots = 13
  M.set = seq(8, 15, by = 1)
  M.EM = nKnots + 2
  select_Method = "knots"
  r.set = 5
  splineObj = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  sig2hat = 1/16
} else if (DataType == "easySin") {
    Mdata = 5
    rdata = 3
    pcaTrans = matrix(rnorm(Mdata*rdata), Mdata, rdata)
    pcaTrans = qr.Q(qr(pcaTrans))
    eigenfList = get_eigenfunList_PaulNoSpline("easySin", pcaTrans)
    obsCol = get_allObs_PaulNoSpline(samplesize, M = Mdata, r = rdata, DataType, pcaTrans,
                                     scoreType, alpha = 0.6, nmin = 2, nmax = 10, a=1, b=1, sig = 1/4, noiseType)
    M.set = seq(4, 12, by = 1)
    r.set = 3
    nKnots = 8
    M.EM = nKnots + 2
    select_Method = "knots"
    sig2hat = 1/16
} else if (DataType == "pracSin") {
    Mdata = 10
    rdata = 5
    pcaTrans = matrix(rnorm(Mdata*rdata), Mdata, rdata)
    pcaTrans = qr.Q(qr(pcaTrans))
    eigenfList = get_eigenfunList_PaulNoSpline("pracSin", pcaTrans)
    obsCol = get_allObs_PaulNoSpline(samplesize, M=Mdata, r=rdata, DataType, pcaTrans,
                                     scoreType, alpha = 0.6, nmin=2,nmax=10,a=1,b=1,sig=1/4, noiseType)
    M.set = seq(10, 20 ,by = 1)
    r.set = 5
    nKnots = 13
    M.EM = nKnots + 2
    select_Method = "knots"
    sig2hat = 1/16
} 





require(h5r)

source("tinyTestHarness.R")

TH <- TestHarness()
h5 <- H5File("test-read-points.h5", 'w')

TH('create h5 attribute', {
  createH5Attribute(h5, "e", 10:1)
  all(getH5Attribute(h5, "e")[] == 10:1)
})

TH('create h5 attribute 2', {
  strings <- c("ACGTACGT", "GGGGGGGGGGGGGGGG", "CCGCGCGCG")
  createH5Attribute(h5, "d", strings)
  all(getH5Attribute(h5, "d")[] == strings)
})

indta <- cbind(c(' xx ', '$$$$$$$$$$$', '||||'),
               c(' jjjj', '"$"', '"""""""'))

TH('create dataset', {
  d <- createH5Dataset(h5, "d1", indta)
  all(d[] == indta)
})

TH('write data', {
  ndta <- cbind("yyy", "xxx")
  writeH5Data(d, ndta, c(1,1), dim(ndta))
  indta[1,] <- ndta
  all(getH5Dataset(h5, "d1")[] == indta)
})

TH('read points', {
  d1 <- createH5Dataset(h5, "d3", runif(100000))
  p <- readPoints(d1, ss <- sample(1:length(d1), size = 1000, replace = T))
  all(p == d1[ss])
})

TH(action = 'print')
TH(action = 'throw')

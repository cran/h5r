require(h5r)

source("tinyTestHarness.R")

##
## Make a new TestHarness.
##
TH <- TestHarness()

fileName <- "testwrite.h5"

TH("file create", { h5 <- H5File(fileName, 'w'); TRUE })
TH("file remove", { fn <- h5@fileName; rm(h5) ; gc() ; file.remove(fn); TRUE })

h5 <- H5File(fileName, 'w')

TH("group create 1", { g1 <- createH5Group( h5, "grp1" ); TRUE })
TH("group create 2", {createH5Group( g1, "grp1-1" ); TRUE })
TH("group create 3", {createH5Group( g1, "grp1-2" ); TRUE })

TH("group exists 1", all(names(listH5Contents(h5)) == c(".", "grp1", "grp1/grp1-1", "grp1/grp1-2")))
TH("group delete 1", deleteH5Obj(h5, "grp1/grp1-2"))
TH("group exists 2", all(names(listH5Contents(h5)) == c(".", "grp1", "grp1/grp1-1")))

TH("dataset exists 1", {
  d <- createH5Dataset(g1, "d1", dims = c(10, 10), dType = "integer")
  all(dim(d) == c(10, 10))
})

TH("dataset write 1", {
  indta <- as.integer(outer(1:10, 1:10))
  writeH5Data(d, data = indta, offsets = as.integer(c(1, 1)), extents = as.integer(c(10, 10)))
  all(d[] == indta)

  mdta <- rbind(as.integer(c(rep(1, 10), rep(2, 10))))
  writeH5Data(d, data = mdta, offsets = as.integer(c(1, 1)), extents = as.integer(c(2, 10)))
  indta[1:20] <- mdta

  ## the by-column ordering vs. by-row.
  all(indta == t(rbind(d[1:2,], d[3:10,])))
})

TH("string dataset create 1", {
  d <- createH5Dataset(g1, "d2", z <- as.character(runif(1000)))
  (all(d[1:10] == z[1:10]) &&
   all(z[1:10]==unlist(read1DSlabs(d, 1:10, rep(1, 10)))))
})

TH("string dataset create 2", {
  d <- createH5Dataset(g1, "d3", z <- cbind(as.character(runif(1000)), as.character(runif(1000))))
  (all(d[1:2,1] == z[1:2, 1]) && all(z[1:2, 1] == readSlab(d, c(1,1), c(2, 1))))
})

TH("attribute creation 1", {
  atr <- createH5Attribute(g1, "jim", 20:1)
  atr@name == "jim"
})
TH("attribute fetch 1", {
  all(getH5Attribute(g1, "jim")[] == 20:1)
})
TH("attribute deletion 1", {
  deleteH5Attribute(g1, "jim")
})

TH("attribute creation 2", {
  atr <- createH5Attribute(g1, "jim", as.character(20:1))
  atr@name == "jim"
})
TH("attribute fetch 2", {
  all(getH5Attribute(g1, "jim")[] == as.character(20:1))
})
TH("attribute deletion 2", {
  deleteH5Attribute(g1, "jim")
})

gt <- createH5Group(h5, "grp-type")

TH("dataset type create 1", {
  d <- createH5Dataset(gt, "d-type-1", matrix(rnorm(200, 100, sd = 10), ncol = 10), dType = "integer",
                       chunk = c(10, 10))
  dd <- getH5Dataset(gt, "d-type-1")
  all(dim(dd) == c(20, 10)) && storage.mode(dd[]) == 'integer'
})

TH("dataset type create 2", {
  d <- createH5Dataset(gt, "d-type-2", rnorm(200, 100, sd = 10), dType = "integer",
                       dim = c(10, 20), chunk = c(10, 10))
  dd <- getH5Dataset(gt, "d-type-2")
  all(dim(dd) == c(10, 20)) && storage.mode(dd[]) == 'integer'
})

TH("dataset type 3", {
  d <- createH5Dataset(gt, "d-type-3", rnorm(200, 100, sd = 10), dType = "character",
                       dim = c(10, 20), chunk = c(10, 10))
  dd <- getH5Dataset(gt, "d-type-3")
  all(dim(dd) == c(10, 20)) && storage.mode(dd[]) == 'character'
})

TH(action = 'print')
TH(action = 'throw')


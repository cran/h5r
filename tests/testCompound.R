require(h5r)

source("tinyTestHarness.R")

##
## Make a new TestHarness.
##
TH <- TestHarness()

file <- system.file("h5_files", "compound.h5", package = 'h5r')
f <- H5File(file)

dsets <- c("NS", "S10", "SVLEN") 

for (dset in dsets) {
  d <- getH5Dataset(f, dset)
  TH(paste(dset, "read", sep = "-"), !is.null(d[]))
  TH(paste(dset, "data.frame", sep = "-"),
     nrow(as.data.frame(d[1:10])) == 10)
  TH(paste(dset, "data.frame 2", sep = "-"),
     all(as.data.frame(d[1:10]) == as.data.frame(d[10:1])[10:1,]))
}

TH(action = 'print')
TH(action = 'throw')


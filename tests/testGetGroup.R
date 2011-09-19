require(h5r)

source('tinyTestHarness.R')
TH <- TestHarness()

file <- system.file("h5_files", "ex_1.h5", package = 'h5r')
f <- H5File(file)

TH('get group', {
  ## Not quite sure how to really do this right, I don't expect the
  ## memory to always go back to the original state.
  g = gc()
  v <- replicate(1000, {
    getH5Group(f, "group_1")
  })
  rm(v)
  gc()
  TRUE
})

TH(action = 'print')
TH(action = 'throw')

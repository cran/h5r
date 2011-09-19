##
## This isn't really a test of the h5r package - so probably best to
## just comment it out. Mostly, this is about me allocating too small
## chunks of memory and then running out.
##

##
## Currently, it seems as if I hold on to too much memory - this
## corresponds to me not cleaning something up in HDF5 because
## Valgrind says I'm fine.
##

## Essentially, this code reproduces the spirit of this post:
## https://stat.ethz.ch/pipermail/r-devel/2011-September/062025.html
## Which is related to this false bug:
## https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14611

## the thrust is that with these small allocations you don't allow the
## Kernel to reclaim memory because the memory becomes "fragmented" -
## it is not 100% certain to me that R is faultless on this score. I'm
## finding the forced call of malloc_trim to be necessary in order to
## actually return the memory to the OS.
require(h5r)

showPS <- function() system(paste('ps -eo pid,vsz,%mem | grep', Sys.getpid()))
gcl <- function() { lapply(1:10, gc, verbose = F)[[10]] }

gc()

## 7/2/2013 - malloc_trim commented out for CRAN submission.


showPS()
m <- .Call("h5R_allocate_gig")
b <- 'bar' # from the post, "blocking the memory"
rm(m)
gcl()
showPS()
# h5r:::.mallocTrim()
showPS()

m <- sapply(1:1000, function(a) {
  .Call("h5R_allocate_meg")
})
b <- 'bar' # from the post, "blocking the memory"
rm(m)
gcl()
showPS()
# h5r:::.mallocTrim()
showPS()

m <- sapply(1:100000, function(a) {
  .Call("h5R_allocate_k")
})
b <- 'bar' # from the post, "blocking the memory"
rm(m)
gcl()
showPS()
# h5r:::.mallocTrim()
showPS()

m <- sapply(1:1000000, function(a) {
  .Call("h5R_allocate_k")
})
b <- 'bar' # from the post, "blocking the memory"
rm(m)
gcl()
showPS()
# h5r:::.mallocTrim()
showPS()

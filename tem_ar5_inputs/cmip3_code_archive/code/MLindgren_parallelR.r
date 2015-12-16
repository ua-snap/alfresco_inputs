RNGkind("L'Ecuyer-CMRG")
set.seed(<something>)
## start M workers
s <- .Random.seed
for (i in 1:M) {
s <- nextRNGStream(s)
# send s to worker i as .Random.seed
}
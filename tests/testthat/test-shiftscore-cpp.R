ssr <- function(x, y) {
  stopifnot(dim(x) == dim(y))
  stopifnot(length(x) > 1)
  ex <- cumsum(x) / sum(x)
  ey <- cumsum(y) / sum(y)
  c(
    sum(ex - ey),
    sum(pmax(ex - ey, 0)),
    sum(pmin(ex - ey, 0)),
    sum(abs(ex - ey))
  ) / (length(x) - 1)
}

test_that("fast-score works", {
  x <- matrix(c(10,4:1,rep(0,5)), ncol = 1)
  y <- matrix(c(rep(0,5),1:4,10), ncol = 1)
  ssf <- c(ShiftScoreFast(x, y, sum(x), sum(y)))
  expect_equal(ssf, ssr(x,y))
  
  # skipped locations
  x <- matrix(c(10,4:1,rep(0,7)), ncol = 1)
  y <- matrix(c(rep(0,7),1:4,10), ncol = 1)
  
  ssf <- c(ShiftScoreFast(x, y, sum(x), sum(y)))
  expect_equal(ssf, ssr(x,y))
})


test_that("shift-score works", {
  x <- matrix(c(10,4:1,rep(0,5)), ncol = 1)
  y <- matrix(c(rep(0,5),1:4,10), ncol = 1)
  
  ssf <- c(ShiftScore(x, y, 0L, 100L))
  expect_equal(ssf, c(ssr(x,y),1,1))
  # skipped locations
  x <- matrix(c(10,4:1,rep(0,7)), ncol = 1)
  y <- matrix(c(rep(0,7),1:4,10), ncol = 1)
  
  ssf <- c(ShiftScore(x, y, 0L, 100L))
  expect_equal(ssf, c(ssr(x,y),1,1))
  
  # small p-value
  ssf <- c(ShiftScore(x,y,1L,1000L))
  expect_equal(ssf, c(ssr(x,y),0,0))
  
  # positives and negatives
  x <- matrix(c(8,4,3,2,0,0,6, rep(0,5)))
  y <- matrix(c(5,4,8,1,0,4, rep(0,6)))
  
  ss <- ShiftScore(x, y, 0L, 1000L)[1:4]
  expect_equal(ss, ssr(x,y))
})


test_that("allTheShiftScores works", {
  a <- c(10,4,3,2,1)
  da <- tibble::tibble(
    scores = c(a, rev(a), a, rev(a), a, a),
    distances = c(0:4, 5:9, 0:4, 7:11, 0:4, c(0:3,5)),
    samp = rep(c(0,1,0,1,0,1), each = 5),
    fhash = rep(letters[1:3], each = 10))
  
  t1 <- with(da, 
             allTheShiftScores(fhash, distances, scores, samp, 1L, 1000L, 3L))
  expect_equal(t1[1,], colSums(t1[2:3,]))
  expect_equal(t1[4,], t1[2,] - t1[3,])
  expect_true(all(t1[5:6,3] > 0.1))
  x <- matrix(c(10,4:1,rep(0,5)), ncol = 1)
  y <- matrix(c(rep(0,5),1:4,10), ncol = 1)
  expect_equal(t1[1:4,1], ssr(y,x))
})

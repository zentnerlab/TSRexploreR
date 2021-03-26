test_that("fast-score works", {
  x <- c(10,4,3,2,1,0,0,0,0,0)
  y <- rev(x)
  fx <- x / sum(x)
  fy <- y / sum(y)
  
  k <- 10L
  pos <- seq_along(x)
  ss <- sum((fy - fx) * pos / k)
  
  expect_equal(ShiftScoreFast(x, y, k, sum(x), sum(y), pos), ss)
  
  # skipped locations
  pos <- c(1:5, 6:10 + 2)
  ss <- sum((fy - fx) * pos / k)
  expect_equal(ShiftScoreFast(x, y, k, sum(x), sum(y), pos), ss)
  expect_equal(ShiftScoreFast(x, y, k, sum(x), sum(y), pos - 1), ss)
})


test_that("shift-score works", {
  x <- Matrix::sparseMatrix(1:5, rep(1,5), x = c(10,4,3,2,1), dims = c(10,1))
  y <- Matrix::sparseMatrix(6:10, rep(1,5), x = rev(c(10,4,3,2,1)), dims = c(10,1))
  fx <- x / sum(x)
  fy <- y / sum(y)
  
  k <- 10L
  pos <- 1:10
  ss <- sum((fy - fx) * pos / k)
  ssp <- sum(pmax((fy - fx) * pos / k, 0))
  ssn <- sum(pmin((fy - fx) * pos / k, 0))
  cpp_ss <- ShiftScore(x, y, 0L, 100L)[,,drop = TRUE]
  expect_equal(length(cpp_ss), 4)
  expect_equal(cpp_ss, c(ss,0,ssp,ssn))
  
  # skipped locations
  x <- Matrix::sparseMatrix(1:5, rep(1,5), x = c(10,4,3,2,1), dims = c(12,1))
  y <- Matrix::sparseMatrix(8:12, rep(1,5), x = rev(c(10,4,3,2,1)), dims = c(12,1))
  fx <- x / sum(x)
  fy <- y / sum(y)
  
  k <- 12L
  pos <- 1:12
  ss <- sum((fy - fx) * pos / k)
  ssp <- sum(pmax((fy - fx) * pos / k, 0))
  ssn <- sum(pmin((fy - fx) * pos / k, 0))
  cpp_ss <- ShiftScore(x, y, 0L, 100L)[,,drop = TRUE]
  expect_equal(length(cpp_ss), 4)
  expect_equal(cpp_ss, c(ss,0,ssp,ssn))
  
  
  
})

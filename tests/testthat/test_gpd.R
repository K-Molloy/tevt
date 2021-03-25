test_that("output of rgpd", {
  set.seed(1)
  expect_equal(rgpd(6, 1:6, 0.5, 0.01), c(1.154527, 2.233249, 3.427128, 5.208487, 5.112751, 7.156477), tolerance=1e-3)
  expect_equal(rgpd(6, 1, 0.5, 0.01), c(2.468417, 1.543512, 1.498398, 1.031899, 1.115453, 1.097225), tolerance=1e-3)
  expect_error(rgpd(1, 1:8, 1:5, 0))
  expect_error(rgpd(10, 1:8, 1, 0.01))
})

test_that("output of dgpd", {
  set.seed(1)
  expect_equal(dgpd(2:4, 1, 0.5, 0.01), c(0.270652877, 0.038077000, 0.005560804), tolerance=1e-3)
  expect_equal(dgpd(2, -2:1, 0.5, 0.01), c(0.0008418422, 0.0055608042, 0.0380770002, 0.2706528769), tolerance=1e-3)
})

test_that("output of pgpd", {
  set.seed(1)
  expect_equal(pgpd(2:4, 1, 0.5, 0.01), c(0.8619670, 0.9802000, 0.9970528), tolerance=1e-3)
})

test_that("output of qgpd", {
  set.seed(1)
  expect_equal(qgpd(seq(0.9, 0.6, -0.1), 2, 0.5, 0.01), c(3.164650, 2.811230, 2.605625, 2.460251), tolerance=1e-3)
})


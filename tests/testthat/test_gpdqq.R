
test_that("gpd qq plot", {
  data=rexp(1000)
  expect_error(gpd.qq(data,1000,1,0), NA)
})


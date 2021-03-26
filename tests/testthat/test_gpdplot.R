
test_that("gpd accessory plots", {
  data(danish)
  out = gpd.fit(danish, 10)
  expect_error(gpd.plot(out, pick=1), NA)
  expect_error(gpd.plot(out, pick=2), NA)
  expect_error(gpd.plot(out, pick=3), NA)
  #expect_error(gpd.plot(out, pick=4), NA)
})


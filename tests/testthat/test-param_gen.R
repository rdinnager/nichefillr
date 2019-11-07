test_that("parameter generation works", {
  example_parms <- example_parms("simple_peaks_2")
  expect_equal(length(example_parms), 3)
})

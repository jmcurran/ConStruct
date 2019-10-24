test_that("maxLikelihood",{
  inFile = system.file("extdata", "infile.txt", package = "ConStruct", mustWork = TRUE)
  data.f = readData(inFile, missing = 0)
  results = maxLikelihood(data.f, resolution = 100)

  expect_true(abs(results$fst.max - 0.209714) < 1e-6)
  expect_true(abs(results$fst - 0.04462213) < 1e-6)
  expect_true(abs(results$G - 0.01240437) < 1e-6)
})


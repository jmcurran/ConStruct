test_that("numLoci", {
  inFile = system.file("extdata", "infile.txt", package = "ConStruct", mustWork = TRUE)
  r = readData(inFile)
  expect_equal(r$numLoci, 12)
})

test_that("numAllelesPerLocus", {
  inFile = system.file("extdata", "infile.txt", package = "ConStruct", mustWork = TRUE)
  r = readData(inFile)
  expect_equal(r$numAllelesPerLocus, rep(8, 12))
})

test_that("Counts", {
  inFile = system.file("extdata", "infile.txt", package = "ConStruct", mustWork = TRUE)
  r = readData(inFile)
  Counts = do.call(rbind, r$Counts)
  dimnames(Counts) = NULL

  testCounts = matrix(c(24,1,33,31,37,120,39,87,89,20,23,58,14,34,26,85,26,21,29,26,17,60,40,32,34,47,108,14,32,65,32,22,
                      85,116,47,64,44,22,17,77,90,26,115,41,19,40,41,42,32,68,63,29,16,17,24,42,18,79,65,44,47,44,58,65,
                       9,67,29,21,70,14,101,46,70,66,35,49,53,52,51,144,71,45,28,83,135,118,60,50,137,32,81,17,31,26,55,31),
                    nrow = 12, ncol = 8)

  expect_equal(Counts, testCounts)
})







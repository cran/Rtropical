library(Rtropical)
library(testthat)
data(apicomplexa)
treevecs <- as.matrix(apicomplexa, parallel = TRUE)

# Test for tropical PCA by polytope
pca_fit <- troppca.poly(treevecs)
point_projections <- t(apply(treevecs, 1, tropproj.poly, tconv = t(pca_fit$pc)))
test_that("projection on polytope agrees", {
  expect_equal(point_projections, pca_fit$projection)
})

# Test for tropical PCA by linear space
pca_fit <- troppca.linsp(treevecs)
point_projections <- tropproj.linsp(treevecs, pca_fit$pc)
test_that("projection on linear space agrees", {
  expect_equal(point_projections, pca_fit$projection)
})

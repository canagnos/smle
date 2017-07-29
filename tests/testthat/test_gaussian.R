context("Streaming multivariate Gaussian tests")

n = 1000
outGen <- genAbrupt(n=n, muW=n/2, p=2, seed=2)
params <- list()
etas = rep(NA, n)
for (i in 1:n){
  params <- streamGaussian(outGen$X[i,], params, alpha=0.001) 
  etas[i] <- params$eta
}

etas_after_change = etas[(n/2+1):(n/2+10)]
etas_before_change = etas[(n/2-10):(n/2)]
mean(etas_after_change)

test_that("Learning rate increases immediately after change", {
  expect_gt(mean(diff(etas_after_change)), 0)
  expect_gt(mean(etas_after_change), mean(etas_before_change))
})


test_that("Learning rate stabilises immediately after change", {
  expect_gt(var(etas_after_change), var(etas_before_change))
})

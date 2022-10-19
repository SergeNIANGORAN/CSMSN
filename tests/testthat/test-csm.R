test_that("CSMSN", {

  #run CSMSN
  n.subject = 50
  n.center = 6
  sigma2c = 1
  sigma2r = 4

  # generated data
  res <- data.generebis(n.subject, n.center,sigma2c,sigma2r)

  # test csm bayesian ATAYAMA
  CONSTRUCT_BAYES(res)

})

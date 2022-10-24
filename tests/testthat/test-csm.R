test_that("CSMSN", {

  #run CSMSN
  n.subject = 50
  n.center = 6
  sigma2c = 1
  sigma2r = 4

  # generated data
  res <- data.generebis(n.subject, n.center,sigma2c,sigma2r)

  # Testing HATAYAMA CSM method using Bayesian Finite Mixture Model
  CONSTRUCT_BAYES(res)

  # Testing DESMET CSM method using Mixed Linear Model
   CONSTRUCT_ZDI(res)

  # Testing Student T-test CSM method
   CONSTRUCT_TestMOY(res)

  # Testing Distance CSM method
   CONSTRUCT_DISTANCE(res)

  # Testing the Master CSM program function
  MASTER_CSM_MOY_GLOBAL_SIMS()

})

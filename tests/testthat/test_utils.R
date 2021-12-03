test_that(".assertPackagesAvailable() works properly", {
 
  testfunc <- function(...) .assertPackagesAvailable(...)
  
  # checks
  expect_error(testfunc(1L))
  expect_error(testfunc("test", "error"))
  expect_error(testfunc("test", c(TRUE, FALSE)))
  
  expect_true(testfunc("base"))
  expect_true(testfunc("githubuser/base"))
  expect_true(testfunc(c("base", "methods")))
  expect_error(testfunc(c("error", "error2")), "BiocManager")
  expect_error(testfunc("error1", suggestInstallation = FALSE), "installed.\n$")
  
  rm(testfunc)
})
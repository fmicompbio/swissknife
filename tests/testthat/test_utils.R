test_that(".assertPackagesAvailable() works properly", {
 
  testfunc <- function(...) .assertPackagesAvailable(...)
  
  # checks
  expect_error(testfunc(1L))
  expect_error(testfunc("test", "error"))
  expect_error(testfunc("test", c(TRUE, FALSE)))
  expect_error(testfunc(c("error", "error2")), "BiocManager")
  expect_error(testfunc("error1", bioc = FALSE), "install.packages")
  
  expect_true(testfunc("base"))
  expect_true(testfunc(c("base", "methods")))

  rm(testfunc)
})
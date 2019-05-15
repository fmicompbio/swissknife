context("col2hex")

test_that("col2hex() works properly", {
 
  # checks
  expect_error(col2hex(col="nonExistentColor"))
  expect_error(col2hex(col="red", alpha="character"))
  expect_error(col2hex(col=c("red", "blue", "green"), alpha=c(0.2, 0.6)))
  expect_error(col2hex(col="red", alpha=256))
  
  # example runs
  expect_identical(col2hex(c("red", "blue")), c("#FF0000FF", "#0000FFFF"))
  expect_identical(col2hex(c("red", "blue"), alpha = 0.7), c("#FF0000B2", "#0000FFB2"))
  expect_identical(col2hex(c("red", "blue"), alpha = 0.7), col2hex(c("red", "blue"), alpha = 178))
  expect_identical(col2hex(5), "#00FFFFFF")
   
})
context("symbolic formula for structural break models")


#Model with duplicate regressors in both x and z
test_that("Formula cannot have duplicate z and x regressors", {
  expect_error(check_formula(y~z1|z1+x1),
                 'Duplicate regressor(s): z1', fixed = TRUE)
})

#Model without y
test_that("Formula needed to be 2-sided", {
  expect_error(check_formula(~x1+z1), 'Invalid formula. Need an expression of the form y ~ z | x or y ~ z', fixed = TRUE)
})

#Model without z regressors
test_that("Formula must have at least 1 z regressors", {
  expect_error(check_formula(y~-1|x1),
               'No z regressors. Use lm() instead', fixed = TRUE)
})


#Model with both intercepts
test_that("Formula with intercepts in both z and x regressors",{
  expect_message(check_formula(y~z|x), 
                 'Formula has both intercept in z and x regressors. Intercept in z regressors will take priority', fixed = TRUE)
})
test_that("utils template works properly", {
    source(system.file("templates", "utils_template.R", package = "swissknife"))
    source(system.file("templates", "tests_utils_template.R",
                       package = "swissknife"))
})
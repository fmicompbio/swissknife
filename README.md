## `swissknife` - Handy code shared in the FMI CompBio group

| Platforms |  OS  | R CMD check | Coverage |
|:----------------:|:----------------:|:----------------:|:----------------:|
| GitHub Actions | Linux/Windows/macOS | [![R build status](https://github.com/fmicompbio/swissknife/workflows/R-CMD-check/badge.svg)](https://github.com/fmicompbio/swissknife/actions) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/swissknife/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/swissknife) |


The `swissknife` package is a collection of useful R functions performing various tasks that
might be re-usable and worth sharing.

You can browse the package contents at [https://fmicompbio.github.io/swissknife/](https://fmicompbio.github.io/swissknife/).

### Code guidelines
- Place your function(s) in an R script named in a way that makes it easy to guess what's inside.
- Use [Roxygen](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html) to document your function(s). Some useful tags: `@title`, `@description`, `@details`, `@param`, `@return`, `@author`, `@examples`, `@export`, `@import`/`@importFrom`.
- If possible, add an example in the Roxygen preamble to show how to run each function.
- If possible, add one or more unit tests to check that the function(s) work as intended (using [testthat](https://cran.r-project.org/web/packages/testthat/index.html)).
- When adding code, bump the version number of the package in the DESCRIPTION file, add an entry in the `NEWS.md` file and to the `_pkgdown.yml` file (`reference` section). When you push your edits to github, github actions will automatically run `pkgdown::build_site()` and deploy it to  [https://fmicompbio.github.io/swissknife/](https://fmicompbio.github.io/swissknife/).
- If your function uses functions defined in other packages, import the latter by adding `@importFrom package function` in the Roxygen preamble. Then document the package (e.g. using `devtools::document()`) to regenerate the NAMESPACE file. It's also considered good practice to refer to these external functions with `package::function` in the code, to avoid name clashes. 
- As far as possible, avoid hardcoding options (e.g., point sizes in plots, or other things that someone might want to change). Instead provide them as arguments to your function. 
- Generally, it's good to keep functions small (and doing a single thing), for easier testing and increased modularity. As an example, if a function does both some calculations and plotting, it can often be better to split it up into two functions, each with a clearly defined scope. 

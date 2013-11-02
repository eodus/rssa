library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("Cadzow");

test_that("Cadzow 1d test", {
  env <- new.env();
  load(system.file("extdata", "cadzow.testdata.rda", package = "Rssa"), envir = env);
  names <- c("co2.td", "fr20.td", "fr50.td", "fr20.nz.td", "fr50.nz.td", 
  	"lofr20.nz.td", "lofr50.nz.td");
  for (name in names) {
    test.cadzow.test.data(test.data = env[[name]]);
  }
});

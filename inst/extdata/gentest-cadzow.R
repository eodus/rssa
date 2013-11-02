library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));

all.svd <- c("svd", "eigen", "propack", "nutrlan")
svd.wo.nutrlan <- c("svd", "eigen", "propack")

co2.td <- make.cadzow.test.data(series = co2,
                                Ls = c(17, 234, 235, 300, 400),
                                rank = 4,
                                numiter = 100,
                                kind = "1d-ssa",
                                svd.method = "e",
                                svd.methods = list(svd.wo.nutrlan, all.svd, all.svd, all.svd, all.svd),
                                tolerance = 2e-6,
                                neig = 20);
test.cadzow.test.data(test.data = co2.td);

finite.rank.r5ex1 <- function(N) {
  tt <- 1:N;
  cos(2*pi*(1:N) / 7) + sin(2*pi*(1:N) / 17) * exp(tt / N * 1.5) + exp(-tt / N * 1.2);
}

fr20 <- finite.rank.r5ex1(20);
fr50 <- finite.rank.r5ex1(50);

fr20.td <- make.cadzow.test.data(series = fr20,
                                 Ls = c(8, 10, 15),
                                 rank = 4,
                                 numiter = 50,
                                 kind = "1d-ssa",
                                 svd.method = "e",
                                 svd.methods = list(svd.wo.nutrlan, svd.wo.nutrlan, svd.wo.nutrlan),
                                 neig = 5);
test.cadzow.test.data(test.data = fr20.td);

fr50.td <- make.cadzow.test.data(series = fr50,
                                 Ls = c(17, 25, 40),
                                 rank = 4,
                                 numiter = 50,
                                 kind = "1d-ssa",
                                 svd.method = "e",
                                 svd.methods = list(svd.wo.nutrlan, svd.wo.nutrlan, svd.wo.nutrlan),
                                 neig = 5);
test.cadzow.test.data(test.data = fr50.td);

set.seed(1);

fr20.nz.td <- make.cadzow.test.data(series = fr20 + rnorm(fr20),
                                    Ls = c(8, 10, 15),
                                    rank = 4,
                                    numiter = 50,
                                    kind = "1d-ssa",
                                    svd.method = "e",
                                    svd.methods = list(svd.wo.nutrlan, svd.wo.nutrlan, svd.wo.nutrlan),
                                    neig = 15);
test.cadzow.test.data(test.data = fr20.nz.td);

set.seed(1);

fr50.nz.td <- make.cadzow.test.data(series = fr50 + rnorm(fr50),
                                    Ls = c(17, 25, 40),
                                    rank = 4,
                                    numiter = 50,
                                    kind = "1d-ssa",
                                    svd.method = "e",
                                    svd.methods = list(svd.wo.nutrlan, svd.wo.nutrlan, svd.wo.nutrlan),
                                    neig = 15);
test.cadzow.test.data(test.data = fr50.nz.td);

set.seed(1);

fr_line <- function(x) {finite.rank.r5ex1(x) + rnorm(x)}

lofr20.nz.td <- make.cadzow.test.data(series = sapply(rep(20, 3), fr_line),
                                      Ls = c(8, 10, 15),
                                      rank = 4,
                                      numiter = 50,
                                      kind = "mssa",
                                      svd.method = "e",
                                      svd.methods = list(svd.wo.nutrlan, svd.wo.nutrlan, svd.wo.nutrlan),
                                      neig = 15);
test.cadzow.test.data(test.data = lofr20.nz.td);

set.seed(1);

lofr50.nz.td <- make.cadzow.test.data(series = sapply(rep(50, 3), fr_line),
                                      Ls = c(17, 25, 40),
                                      rank = 4,
                                      numiter = 50,
                                      kind = "mssa",
                                      svd.method = "e",
                                      svd.methods = list(svd.wo.nutrlan, svd.wo.nutrlan, svd.wo.nutrlan),
                                      neig = 15);
test.cadzow.test.data(test.data = lofr50.nz.td);

save(co2.td, fr20.td, fr50.td, fr20.nz.td, fr50.nz.td, lofr20.nz.td, lofr50.nz.td,
     file = system.file("extdata", "cadzow.testdata.rda", package = "Rssa"),
     compress = "xz", compression_level = 9);

#TODO eps check
#TODO correction check
#TODO more mssa tests

#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Anton Korobeynikov <asl@math.spbu.ru>
#
#   This program is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public
#   License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option)
#   any later version.
#
#   This program is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied
#   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#   PURPOSE.  See the GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public
#   License along with this program; if not, write to the
#   Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
#   MA 02139, USA.

mcadzow <- function(x, rank,
                    eps = 1e-6, numiter = 0,
                    ..., cache = TRUE) {
  # Obtain the initial reconstruction of rank r
  r <- reconstruct(x, groups = list(1:rank), ..., cache = cache)
  stopifnot(length(r) == 1)
  F <- r[[1]]
  # Do the actual iterations until the convergence (or stoppping due to number
  # of iterations)
  it <- 0
  repeat {
    s <- clone(x, copy.cache = FALSE, copy.storage = FALSE)
    .set(s, "Fattr", attributes(F))
    .set(s, "Fclass", class(F))
    .set(s, "weights", .get(x, "weights"))
    .set(s, "fmask", .get(x, "fmask"))
    FL <- .to.series.list(r[[1]], na.rm = TRUE)
    .set(s, "F", FL)
    .set(s, "Iattr", lapply(FL, attributes))
    r <- reconstruct(s, groups = list(1:rank), ..., cache = cache)
    stopifnot(length(r) == 1)
    rF <- r[[1]]

    it <- it + 1
    if ((numiter > 0 && it >= numiter) || max((F-rF)^2) < eps)
      break
    F <- rF
  }

  F
}

cadzow.ssa <- function(x, rank,
                       eps = 1e-6, numiter = 0,
                       ..., cache = TRUE) {
  # Obtain the initial reconstruction of rank r
  r <- reconstruct(x, groups = list(1:rank), ..., cache = cache)
  stopifnot(length(r) == 1)
  F <- r[[1]]

  # Do the actual iterations until the convergence (or stoppping due to number
  # of iterations)
  it <- 0
  repeat {
    s <- clone(x, copy.cache = FALSE, copy.storage = FALSE)
    .set(s, "F", F)
    r <- reconstruct(s, groups = list(1:rank), ..., cache = cache)
    stopifnot(length(r) == 1)
    rF <- r[[1]]

    it <- it + 1
    if ((numiter > 0 && it >= numiter) || max((F-rF)^2) < eps)
      break
    F <- rF
  }

  F
}

cadzow.1d.ssa <- function(x, rank,
                          correct = TRUE,
                          eps = 1e-6, numiter = 0,
                          ..., cache = TRUE) {
  # Get the result w/o any correction
  fcall <- match.call(expand.dots = FALSE)
  fcall[[1]] <- cadzow.ssa
  fcall$correct <- NULL
  F <- eval(fcall, parent.frame())

  # Correct the stuff if requested
  if (correct) {
    h1 <- hankel(F, x$window)
    h2 <- hankel(.get(x, "F"), x$window)

    F <- sum(h1 * h2) / sum(h2 * h2) * F
  }

  F
}

cadzow.mssa <- function(x, rank,
                        correct = TRUE,
                        eps = 1e-6, numiter = 0,
                        ..., cache = TRUE) {
  # Get the result w/o any correction
  fcall <- match.call(expand.dots = FALSE)
  fcall[[1]] <- mcadzow
  fcall$correct <- NULL
  F <- eval(fcall, parent.frame())

  # Correct the stuff if requested
  if (correct) {
    h1 <- unlist(.to.series.list(F, na.rm = TRUE))
    h2 <- unlist(.to.series.list(.get(x, "F"), na.rm = TRUE))
    
    F <- sum(h1 * h2 * .hweights(x)) / sum(h2 * h2 * .hweights(x)) * F
  }

  F
}

cadzow <- function(x, ...)
  UseMethod("cadzow")

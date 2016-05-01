#' @import mvtnorm
#' @import purrr
NULL

#' Simulate Parameters from a Model
#'
#' @param x A model object
#' @param ... further arguments passed to or from other methods.
#'
#' @export
postsim <- function(x, ...) {
  UseMethod("postsim")
}

#' Simulate Expected Values from a Model
#'
#' @param x A model object
#' @param ... further arguments passed to or from other methods.
#'
#' @export
postsimev <- function(x, ...) {
  UseMethod("postsimev")
}

#' Simulate Responses from a Model
#'
#' @param x A model object
#' @param ... further arguments passed to or from other methods.
#'
#' @export
postsimy <- function(x, ...) {
  UseMethod("postsimy")
}

#' Posterior Simulation of Partial Effects
#'
#' Calculate the partial (discrete marginal) effects of a model.
#'
#' @details
#'
#' This function calculates:
#' \deqn{(E(y|X_1)) - E(y|X_2)) / \delta, }{((E(y|X1) - E(y|X2)) / delta, }
#' where \eqn{\delta}{delta} is a difference to scale by.
#'
#' @param x Model object
#' @param data1,data2 Data frames
#' @param n Number of iterations
#' @param delta Width of difference. See Details.
#' @param ... further arguments passed to or from other methods.
#' @export
postsim_partialfx <- function(x, data1, data2, n, delta, ...) {
  UseMethod("postsim_partialfx")
}

#' (Average) Partial Effects
#'
#' The partial effects are the effects of discrete change
#' in a variable.
#'
#' @param x A model object
#' @param data1,data2 Model data
#' @param delta Width of difference. See Details.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
partialfx <- function(x, data1, data2, delta, ...) {
  UseMethod("partialfx")
}

#' @rdname partialfx
#' @export
avg_partialfx <- function(x, data1, data2, delta, ...) {
  UseMethod("avg_partialfx")
}

#' Calculate (Average) Marginal Effects of a Model
#'
#' The partial effects are the effects of discrete change
#' in a variable.
#'
#' @param x A model object
#' @param ... further arguments passed to or from other methods.
#'
#' @export
marfx <- function(x, ...) {
  UseMethod("marfx")
}

# Set Numerical Differentiation width
# See http://www.karenkopecky.net/Teaching/eco613614/Notes_NumericalDifferentiation.pdf
numdiff_width <- function(x) {
  # Add and subtract x to remove roundoff error
  x + max(abs(x), 1) * sqrt(.Machine$double.eps) - x
}

mfx_preprocess <- function(x, data, variable, level) {
  data2 <- data
  v <- data[[variable]]
  if (is.null(v)) {
    stop(sprintf("variable %s not found in data", variable))
  }
  if (is.ordered(v)) {
    data2[[variable]] <- ordered(min(as.integer(v) + 1L, nlevels(v)),
                                 levels(v))
    delta <- 1
  } else if (is.factor(v)) {
    data[[variable]] <- factor(levels(v)[1], levels = levels(v))
    if (!is.null(level)) {
      data2[[variable]] <- factor(level, levels = levels(v))
    }
    delta <- 1
  } else if (is.logical(v)) {
    data2[[variable]] <- TRUE
    data[[variable]] <- FALSE
    delta <- 1
  } else if (is.integer(v)) {
    delta <- 1L
    data2[[variable]] <- v + delta
  } else if (is.numeric(v)) {
    delta <- numdiff_width(v)
    data2[[variable]] <- v + delta
  } else {
    stop(sprintf("Variables of class %s not supported", class(v)))
  }
  list(data1 = data, data2 = data2, delta = delta)
}

#' @export
marfx.default <- function(x, variable, level = NULL,
                          data = stats::model.frame(x),  ...) {
  prep <- mfx_preprocess(x, data, variable, level)
  mfx <- partialfx(x, data1 = prep$data1, data2 = prep$data2,
                   delta = prep$delta, ...)
  mfx
}

#' @rdname marfx
#' @export
avg_marfx <- function(x, ...) {
  UseMethod("avg_marfx")
}

#' @export
avg_marfx.default <- function(x, variable, level = NULL,
                              data = stats::model.frame(x), ...) {
  prep <- mfx_preprocess(x, data, variable, level)
  mfx <- avg_partialfx(x, data1 = prep$data1, data2 = prep$data2,
                       delta = prep$delta, ...)
  mfx
}

sim_summary <- function(x, confint = 0.95, estimate = NULL) {
  if (is.numeric(estimate)) {
    point_est <- estimate
  } else {
    FUN <- match.fun(estimate)
    point_est <- apply(x, 1, match.fun(FUN))
  }
  std.error <- apply(x, 1, sd)
  p <- (1 - confint) / 2
  ci <- apply(x, 1, quantile, probs = c(p, 1 - p))
  ret <- cbind(point_est, std.error, ci[1, ], ci[2, ])
  colnames(ret) <- c("estimate", "std.error", "conf.low", "conf.high")
  ret
}


#' Methods for ordered factors
#'
#' Median and quantile methods for ordered factors.
#'
#' @param x An ordered factor
#' @param na.rm logical. If true, any NA and NaN's are removed from \code{x} before computing.
#' @export
#' @importFrom stats median
#' @name ordered
median.ordered <- function(x, na.rm = FALSE) {
  ordered(levels(x)[floor(median(as.integer(x), na.rm = na.rm))], levels(x))
}

#' @param probs numeric vector of probabilities with values in [0,1]. (Values up to 2e-14 outside that range are accepted and moved to the nearby endpoint.)
#' @param left logical. If left, then quantiles are calculated as P(X <= x), else P(X > x).
#' @param names logical; if true, the result has a names attribute. Set to FALSE for speedup with many probs.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @importFrom stats quantile
#' @rdname ordered
quantile.ordered <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
                             left = TRUE, names = TRUE, ...) {
  FUN <- if (left) floor else ceiling
  q <- FUN(quantile(as.integer(x), na.rm = na.rm, type = 2, probs = probs,
           names = names))
  setNames(ordered(levels(x)[q], levels(x)), names(q))
}

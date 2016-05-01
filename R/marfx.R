#' @import mvtnorm
#' @import purrr
NULL

#' Simulate Parameters from a Model
#'
#' @param x A model object
#' @param ... extra arguments
#'
#' @export
postsim <- function(x, ...) {
  UseMethod("postsim")
}

#' Simulate Predicted Values from a Model
#'
#' @param x A model object
#' @param ... extra arguments
#'
#' @export
postsimy <- function(x, ...) {
  UseMethod("postsimy")
}

#' Simulate Expected Values from a Model
#'
#' @param x A model object
#' @param ... extra arguments
#'
#' @export
postsimev <- function(x, ...) {
  UseMethod("postsimev")
}

#' Calculate Expected Values of a Model
#'
#' This method calculates the expected value of common models
#' given known parameter values. Although it is generally easy to calculate expected values
#' of a model after fitting it using the function \code{predict},
#' there is often not an easy way to do this when the parameters
#' of the model are known, rather than estimated.
#'
#' @param param Paramters (coefficients) of the model. Often in the same
#'   format as returned by \code{\link{coef}}.
#' @param ... extra arguments
#'
#' @export
ev <- function(param, ...) {
  UseMethod("ev")
}

#' @export
ev.lm <- function(param, X, ...) {
   as.matrix(X) %*% as.matrix(param)
}

#' Calculate (Average) Partial Effects of a Model
#'
#' The partial effects are the effects of discrete change
#' in a variable.
#'
#' @param x A model object
#' @param data1,data2 Model data
#' @param ... extra arguments
#'
#' @export
partialfx <- function(x, data1, data2, ...) {
  UseMethod("partialfx")
}

sim_partialfx_lm <- function(x, X1, X2, n) {
  obs <- nrow(X1)
  param <- postsim(x, n = n)
  array(as_vector(map(param, function(p, X1, X2) {
    ev.lm(p[["beta"]], X2 - X1)
  }, X1 = X1, X2 = X2), .type = double(obs)), dim = c(obs, n))
}

#' @export
partialfx.lm <- function(x, data1, data2, n = 1000L, confint = 0.95, ...) {
  # Difference
  point_est <- predict(x, newdata = data2) - predict(x, newdata = data1)
  # simulate from posterior to get CI
  mt <- delete.response(terms(x))
  X1 <- model.matrix(mt, data = data1)
  X2 <- model.matrix(mt, data = data2)
  sims <- sim_partialfx_lm(x, X1, X2, n)
  std.error <- apply(sims, 1, sd)
  ret <- cbind(point_est, std.error)
  colnames(ret) <- c("estimate", "std.error")
  ret
}

#' @rdname partialfx
#' @export
avg_partialfx <- function(x, data1, data2, ...) {
  UseMethod("avg_partialfx")
}

#' @export
avg_partialfx.lm <- function(x, data1, data2, n = 1000L, ...) {
  point_est <- mean(predict(x, newdata = data2) - predict(x, newdata = data1))
  # simulate from posterior to get CI
  mt <- delete.response(terms(x))
  X1 <- model.matrix(mt, data = data1)
  X2 <- model.matrix(mt, data = data2)
  sims <- sim_partialfx_lm(x, X1, X2, n)
  std.error <- sd(apply(sims, 2, mean))
  c("estimate" = point_est, "std.error" = std.error)
}

#' Calculate (Average) Marginal Effects of a Model
#'
#' The partial effects are the effects of discrete change
#' in a variable.
#'
#' @param x A model object
#' @param ... extra arguments
#'
#' @export
marfx <- function(x, ...) {
  UseMethod("marfx")
}

#' @export
marfx.lm <- function(x, ...) {
  NULL
}

#' @rdname marfx
#' @export
avg_marfx <- function(x, ...) {
  UseMethod("avg_marfx")
}

#' @export
avg_marfx.lm <- function(x, ...) {
  NULL
}

#' @export
postsim.lm <- function(x, n = 1, ...) {
  x_summary <- summary(x)
  df_residual <- x$df.residual
  beta_hat <- coef(x)
  V_beta_hat <- vcov(x)
  sigma_hat <- x_summary$sigma
  sigma <- sigma_hat / sqrt((df_residual + 1) / rchisq(n, df = df_residual + 1))
  ## TODO parallel process
  lapply(sigma, function(sigma, b, V) {
    list(beta = as.numeric(rmvnorm(1, b, sigma * V)), sigma = sigma)
  }, b = beta_hat, V = V_beta_hat)
}

#' @export
postsimev.lm <- function(x, data = stats::model.frame(x), n = 1000L, ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  map(params, function(p, X) {ev.lm(p[["beta"]], X)}, X = X)
}


#' @export
postsimy.lm <- function(x, n = 1L, data = stats::model.frame(x), ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  map(params, function(p, X) {
    rnorm(nrow(X), ev.lm(p[["beta"]], X), p[["sigma"]])
  }, X = X)
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
#' @param left logical. If left, then quantiles are calculated as P(X <= x),
#'        else P(X > x).
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

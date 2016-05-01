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

#' @export
postsim.lm <- function(x, n = 1L, ...) {
  x_summary <- summary(x)
  df_residual <- x$df.residual
  beta_hat <- coef(x)
  V_beta_hat <- vcov(x)
  sigma_hat <- x_summary$sigma
  #sigma <- sigma_hat / sqrt((df_residual + 1) / rchisq(n, df = df_residual + 1))
  sigma <- rep(sigma_hat, n)
  ## TODO parallel process
  lapply(sigma, function(sigma, b, V) {
    list(beta = as.numeric(rmvnorm(1, b, V)), sigma = sigma)
  }, b = beta_hat, V = V_beta_hat)
}

#' @export
postsim.glm <- function(x, n = 1L, ...) {
  beta_hat <- coef(x)
  V_beta_hat <- vcov(x)
  ## TODO parallel process
  map(array_branch(rmvnorm(n, beta_hat, V_beta_hat), margin = 2),
      function(x) list(beta = x))
}

#' Calculate Expected Values of a Model
#'
#' This method calculates the expected value of common models
#' given known parameter values. Although it is generally easy to calculate expected values
#' of a model after fitting it using the function \code{predict},
#' there is often not an easy way to do this when the parameters
#' of the model are known, rather than estimated.
#'
#' @param param Paramters (coefficients) of the model. Often in the same format as returned by \code{\link{coef}}.
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

#' @export
ev.glm <- function(param, X, family = NULL, ...) {
  mu <- as.matrix(X) %*% as.matrix(param)
  if (!is.null(family)) {
    mu <- family$linkinv(mu)
  }
  mu
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

sim_partialfx_lm <- function(x, data1, data2, n = 1L, delta = 1) {
  mt <- delete.response(terms(x))
  X1 <- model.matrix(mt, data = data1)
  X2 <- model.matrix(mt, data = data2)
  obs <- nrow(X1)
  param <- postsim(x, n = n)
  array(as_vector(map(param, function(p, X1, X2, delta) {
    ev.lm(p[["beta"]], X2 - X1) / delta
  }, X1 = X1, X2 = X2, delta = delta),
  .type = double(obs)), dim = c(obs, n))
}

sim_partialfx_glm <- function(x, data1, data2, n = 1L,
                              response = FALSE, delta = 1) {
  mt <- delete.response(terms(x))
  X1 <- model.matrix(mt, data = data1)
  X2 <- model.matrix(mt, data = data2)
  obs <- nrow(X1)
  param <- postsim(x, n = n)
  array(as_vector(map(param, function(p, X1, X2, response, delta) {
    (ev.lm(p[["beta"]], X2, response = response) -
       ev.lm(p[["beta"]], X1, response = response)) / delta
  }, X1 = X1, X2 = X2, response = response, delta = delta),
  .type = double(obs)), dim = c(obs, n))
}


#' @export
partialfx.lm <- function(x, data1, data2, n = 1000L, confint = 0.95,
                         delta = 1, ...) {
  # Difference
  point_est <- predict(x, newdata = data2) - predict(x, newdata = data1)
  # simulate from posterior to get CI
  sims <- sim_partialfx_lm(x, data1, data2, n, delta)
  std.error <- apply(sims, 1, sd)
  lwr <- point_est + qnorm((1 - confint) / 2) * std.error
  upr <- point_est + qnorm((1 - confint) / 2, lower.tail = FALSE) * std.error
  ret <- cbind(point_est, lwr, upr, std.error)
  colnames(ret) <- c("estimate", "lwr", "upr", "std.error")
  ret
}



#' @export
partialfx.glm <- function(x, data1, data2,
                          n = 1000L,
                          confint = 0.95,
                          response = response,
                          delta = 1, ...) {
  # Difference
  predict_type <- if (response) "response" else "link"
  point_est <- predict(x, newdata = data2, type = predict_type) -
    predict(x, newdata = data1, type = predict_type) / delta
  # simulate from posterior to get CI
  sims <- sim_partialfx_glm(x, data1, data2, n, response, delta)
  std.error <- apply(sims, 1, sd)
  lwr <- point_est + qnorm((1 - confint) / 2) * std.error
  upr <- point_est + qnorm((1 - confint) / 2, lower.tail = FALSE) * std.error
  ret <- cbind(point_est, lwr, upr, std.error)
  colnames(ret) <- c("estimate", "lwr", "upr", "std.error")
  ret
}


#' @rdname partialfx
#' @export
avg_partialfx <- function(x, data1, data2, ...) {
  UseMethod("avg_partialfx")
}


#' @export
avg_partialfx.lm <- function(x, data1, data2, n = 1000L,
                             confint = 0.95,
                             weight = 1,
                             delta = 1,
                             ...) {
  point_est <- weighted.mean(predict(x, newdata = data2) -
                               predict(x, newdata = data1))
  # simulate from posterior to get CI
  sims <- sim_partialfx_lm(x, data1, data2, n, delta)
  std.error <- sd(apply(sims, 2, weighted.mean, weight = weight))
  p <- (1 - confint) / 2
  lwr <- point_est + qnorm(p) * std.error
  upr <- point_est + qnorm(p, lower.tail = FALSE) * std.error
  c("estimate" = point_est, "lwr" = lwr, "upr" = upr, "std.error" = std.error)
}

#' @export
avg_partialfx.glm <- function(x, data1, data2, n = 1000L,
                              confint = 0.95,
                              delta = 1,
                              response = TRUE, ...) {
  predict_type <- if (response) "response" else "link"
  point_est <- predict(x, newdata = data2, type = predict_type) -
    predict(x, newdata = data1, type = predict_type)
  # simulate from posterior to get CI
  sims <- sim_partialfx_glm(x, data1, data2, n, response, delta)
  std.error <- sd(apply(sims, 2, mean))
  p <- (1 - confint) / 2
  lwr <- point_est + qnorm(p) * std.error
  upr <- point_est + qnorm(p, lower.tail = FALSE) * std.error
  c("estimate" = point_est, "lwr" = lwr, "upr" = upr, "std.error" = std.error)
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


mfx_preprocess <- function(x, data, variable, level) {
  data2 <- data
  v <- data[[variable]]
  if (is.null(v)) {
    stop(sprintf("variable %s not found in data", variable))
  }
  if (is.ordered(data[["variable"]])) {
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
    data2[[variable]] <- v + 1L
    delta <- 1
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
                          data = stats::model.frame(x),
                          ...) {
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
                         data = stats::model.frame(x),
                         ...) {
  prep <- mfx_preprocess(x, data, variable, level)
  mfx <- avg_partialfx(x, data1 = prep$data1, data2 = prep$data2,
                       delta = prep$delta, ...)
  mfx
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

#' @export
postsimev.lm <- function(x, data = stats::model.frame(x), n = 1000L, ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  map(params, function(p, X) {ev.lm(p[["beta"]], X)}, X = X)
}

#' @export
postsimev.glm <- function(x, data = stats::model.frame(x), response = TRUE, n = 1000L, ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  family <- if (response) x$family else NULL
  map(params, function(p, X, family) {ev.glm(p[["beta"]], X, family)}, X = X)
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

#' @export
postsimy.lm <- function(x, n = 1L, data = stats::model.frame(x), ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  map(params, function(p, X) {
    rnorm(nrow(X), ev.lm(p[["beta"]], X), p[["sigma"]])
  }, X = X)
}


# @export
# postsimy.glm <- function(x, n = 1L, data = stats::model.frame(x), ...) {
#   X <- model.matrix(delete.response(terms(x)), data = data)
#   params <- postsim.glm(x, n = n, data = data)
#   map(params, function(p, X, object) {
#     object$coef <- p[["beta"]]
#     object$linear.predictors <- ev.glm(p[["beta"]], X, response = FALSE)
#     object$fitted.value <- object$family$linkinv(object$linear.predictors)
#     simulate(object, nsim = 1)
#   }, X = X, object = x)
# }
# Needs to be adapted for all types of glm families
# Can't directly use the simulate method because it requires fitted values,
# which are only in the glm.

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

# Set Numerical Differentiation width
# See http://www.karenkopecky.net/Teaching/eco613614/Notes_NumericalDifferentiation.pdf
numdiff_width <- function(x) {
  # Add and subtract x to remove roundoff error
  x + max(abs(x), 1) * sqrt(.Machine$double.eps) - x
}

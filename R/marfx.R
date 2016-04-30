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
#' @param data Model data
#' @param data2 Model data
#' @param ... extra arguments
#'
#' @export
partialfx <- function(x, data, data2, ...) {
  UseMethod("partialfx")
}

#' @export
partialfx.lm <- function(x, data, data2, n = 1000, ...) {
  # Difference
  point_est <- predict(x, newdata = data2) - predict(x, newdata = data2)
  # If response not deleted, it will be looked for
  f <- delete.response(terms(x))
  X1 <- model.matrix(f, data = data)
  X2 <- model.matrix(f, data = data2)
  obs <- nrow(X1)
  param <- postsim(x, n = n)
  array(as_vector(map(param, function(p, X1, X2) {
    ev.lm(p[["beta"]], X2 - X1)
  }, X1 = X1, X2 = X2), double(obs)), c(double(obs), n))
}

#' @rdname partialfx
#' @export
avg_partialfx <- function(x, data, ...) {
  UseMethod("avg_partialfx")
}

#' @export
avg_partialfx.lm <- function(x, data = stats::model.frame(x), ...) {
  NULL
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
postsimev.lm <- function(x, data = stats::model.frame(x), n = 1000, ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  map(params, function(p, X) {ev.lm(p[["beta"]], X)}, X = X)
}


#' @export
postsimy.lm <- function(x, n = 1, data = stats::model.frame(x), ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  map(params, function(p, X) {
    rnorm(nrow(X), ev.lm(p[["beta"]], X), p[["sigma"]])
  }, X = X)
}

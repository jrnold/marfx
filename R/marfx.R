#' @import mvtnorm
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
#' @param ... extra arguments
#'
#' @export
partialfx <- function(x, ...) {
  UseMethod("partialfx")
}

partialfx.lm <- function(x, data = stats::model.frame(x), data2 = NULL) {
  NULL
}

#' @rdname partialfx
#' @export
avg_partialfx <- function(x, ...) {
  UseMethod("avg_partialfx")
}

avg_partialfx.lm <- function(x, data = stats::model.frame(x), data2 = NULL) {
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

marfx.lm <- function(x, ...) {
  NULL
}

#' @rdname marfx
#' @export
avg_marfx <- function(x, ...) {
  UseMethod("avg_marfx")
}

avg_marfx.lm <- function(x, ...) {
  NULL
}

postsim.lm <- function(x, n = 1, data = stats::model.frame(x)) {
  x_summary <- summary(x)
  df_residual <- x$df.residual
  beta_hat <- coef(x)
  V_beta_hat <- vcov(x)
  sigma_hat <- x_summary$sigma
  sigma <- sigma_hat / sqrt((df_residual + 1) / rchisq(n, df = df_residual + 1))
  ## TODO parallel process
  lapply(sigma, function(sigma, b, V) {
    rmvnorm(1, beta_hat, sigma * V)
  }, b = beta_hat, V = V_beta_hat)
}

postsimev.lm <- function(x, data = stats::model.frame(x)) {
  NULL
}

postsimy.lm <- function(x, data = stats::model.frame(x)) {
  NULL
}

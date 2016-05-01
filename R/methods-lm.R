#' Methods for \code{lm} Models
#'
#' @name lm-methods
NULL

ev_lm <- function(param, X, ...) {
  as.matrix(X) %*% as.matrix(param)
}

#' @param x object
#' @param data1,data2,data Data
#' @param delta Size of the difference
#' @param confint Confidence interval level
#' @param n Number of iterations
#' @param ... further arguments passed to or from other methods.
#' @rdname lm-methods
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

#' @rdname lm-methods
#' @export
postsim_partialfx.lm <- function(x, data1, data2, n = 1L, delta = 1, ...) {
  mt <- delete.response(terms(x))
  X1 <- model.matrix(mt, data = data1)
  X2 <- model.matrix(mt, data = data2)
  obs <- nrow(X1)
  param <- postsim(x, n = n)
  array(as_vector(map(param, function(p, X1, X2, delta) {
    ev_lm(p[["beta"]], X2 - X1) / delta
  }, X1 = X1, X2 = X2, delta = delta),
  .type = double(obs)), dim = c(obs, n))
}

#' @rdname lm-methods
#' @export
partialfx.lm <- function(x, data1, data2, delta = 1, n = 1000L, confint = 0.95,
                         ...) {
  # Difference
  point_est <- (predict(x, newdata = data2) -
                  predict(x, newdata = data1)) / delta
  # simulate from posterior to get CI
  sims <- postsim_partialfx.lm(x, data1, data2, n, delta)
  sim_summary(sims, confint, estimate = point_est)
}

#' @rdname lm-methods
#' @param weights Weights to apply to the data
#' @export
avg_partialfx.lm <- function(x, data1, data2, delta = 1, n = 1000L,
                             confint = 0.95, weights = NULL, ...) {

  # simulate from posterior to get CI
  sims <- postsim_partialfx.lm(x, data1, data2, n, delta)
  if (!is.null(weights)) {
    point_est <- weighted.mean(predict(x, newdata = data2) -
                                 predict(x, newdata = data1), w = weights) /
      delta
    sims_avg <- apply(sims, 2, weighted.mean, w = weights)
  } else {
    point_est <- mean(predict(x, newdata = data2) -
                        predict(x, newdata = data1)) / delta
    sims_avg <- apply(sims, 2, mean)
  }
  sim_summary(array(sims_avg, c(1, n)), confint, estimate = point_est)
}

#' @rdname lm-methods
#' @export
postsimev.lm <- function(x, data = stats::model.frame(x), n = 1000L, ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  map(params, function(p, X) {ev_lm(p[["beta"]], X)}, X = X)
}

#' @rdname lm-methods
#' @export
postsimy.lm <- function(x, n = 1L, data = stats::model.frame(x), ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.lm(x, n = n, data = data)
  map(params, function(p, X) {
    rnorm(nrow(X), ev_lm(p[["beta"]], X), p[["sigma"]])
  }, X = X)
}

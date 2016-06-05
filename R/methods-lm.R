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
#' @param V Variance-covariance matrix for the coefficients. This
#'   allows for the substitution of "robust" covariance matrices.
#' @param ... further arguments passed to or from other methods.
#' @rdname lm-methods
#' @export
simpar.lm <- function(x, n = 1L, V = NULL, ...) {
  summ <- summary(x)
  #.n <- summ$df[1] + summ$df[2]
  #.k <- summ$df[1]
  beta_hat <- coef(x)
  sigma_hat <- summ$sigma
  #sigma <- sigma_hat * sqrt((.n - .k) / rchisq(n, .n - .k))
  sigma <- rep(sigma_hat, n)
  if (is.null(V)) V <- vcov(x)
  map(sigma, function(sigma, b, V) {
    list(beta = setNames(as.numeric(rmvnorm(1, b, V)),
                         names(b)), sigma = sigma)
  }, b = beta_hat, V = V)
}

#' @rdname lm-methods
#' @export
simfdfx.lm <- function(x, data1, data2, n = 1L, delta = 1, ...) {
  mt <- delete.response(terms(x))
  X1 <- model.matrix(mt, data = data1)
  X2 <- model.matrix(mt, data = data2)
  obs <- nrow(X1)
  param <- simpar(x, n = n, ...)
  array(as_vector(map(param, function(p, X1, X2, delta) {
    ev_lm(p[["beta"]], X2 - X1) / delta
  }, X1 = X1, X2 = X2, delta = delta),
  .type = double(obs)), dim = c(obs, n))
}

#' @rdname lm-methods
#' @param weights Weights to apply to the data
#' @export
afdfx.lm <- function(x, data1, data2, delta = 1, n = 1000L,
                             confint = 0.95, weights = NULL, ...) {
  # simulate from posterior to get CI
  sims <- simfdfx(x, data1, data2, n, delta, ...)
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
simev.lm <- function(x, data = stats::model.frame(x), n = 1000L, ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- simpar.lm(x, n = n, data = data, ...)
  map(params, function(p, X) {ev_lm(p[["beta"]], X)}, X = X)
}

#' @rdname lm-methods
#' @export
simy.lm <- function(x, n = 1L,
                        data = stats::model.frame(x), ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- simpar.lm(x, n = n, data = data, ...)
  map(params, function(p, X) {
    rnorm(nrow(X), ev_lm(p[["beta"]], X), p[["sigma"]])
  }, X = X)
}

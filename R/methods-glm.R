#' Methods for \code{glm} Models
#'
#' @name glm-methods
NULL

ev_glm <- function(param, X, family = NULL, ...) {
  mu <- as.matrix(X) %*% as.matrix(param)
  if (!is.null(family)) {
    mu <- family$linkinv(mu)
  }
  mu
}

#' @param x object
#' @param data1,data2,data Data
#' @param delta Size of the difference
#' @param confint Confidence interval level
#' @param n Number of iterations
#' @param response If response is true, then partial effects are on the
#'   the scale of the response variable. If false, then the partial effects
#'   are on the scale of the linear predictors.
#' @param ... further arguments passed to or from other methods.
#' @rdname glm-methods
#' @export
partialfx.glm <- function(x, data1, data2, delta = 1, n = 1000L,
                          confint = 0.95, response = TRUE, ...) {
  # Difference
  predict_type <- if (response) "response" else "link"
  point_est <- (predict(x, newdata = data2, type = predict_type) -
                predict(x, newdata = data1, type = predict_type)) / delta
  # simulate from posterior to get CI
  sims <- postsim_partialfx.glm(x, data1, data2, n, response, delta)
  sim_summary(sims, confint, estimate = point_est)
}

#' @param weights Weights to apply to the data
#' @rdname glm-methods
#' @export
avg_partialfx.glm <- function(x, data1, data2, delta = 1, n = 1000L,
                              confint = 0.95, response = TRUE, weights = NULL,
                              ...) {
  predict_type <- if (response) "response" else "link"
  sims <- postsim_partialfx.glm(x, data1, data2, n, response = response,
                            delta = delta)
  if (!is.null(weights)) {
    point_est <-
      weighted.mean(predict(x, newdata = data2, type = predict_type) -
                      predict(x, newdata = data1, type = predict_type),
                    w = weights)
    sims_avg <- apply(sims, 2, weighted.mean, w = weights)
  } else {
    point_est <- mean(predict(x, newdata = data2, type = predict_type) -
                        predict(x, newdata = data1, type = predict_type))
    sims_avg <- apply(sims, 2, mean)
  }
  sim_summary(array(sims_avg, dim = c(1, n)), confint, estimate = point_est)
}



#' @rdname glm-methods
#' @export
postsimev.glm <- function(x, data = stats::model.frame(x), response = TRUE,
                          n = 1000L, ...) {
  X <- model.matrix(delete.response(terms(x)), data = data)
  params <- postsim.glm(x, n = n, data = data)
  family <- if (response) x$family else NULL
  map(params, function(p, X, family) {ev_glm(p[["beta"]], X, family)}, X = X)
}

# postsimy_glm_binomial <- function(n, eta, family = binomial(), weights = NULL, ...) {
#   rbinom(length(prob), size = weights, prob = family$linkinv(eta)) / size
# }
#
# postsimy_glm_gaussian <- function(n, eta, family = gaussian(), weights = NULL, sigma = 1) {
#   rnorm(length(eta), mean = family$linkinv(eta), sd = sigma / sqrt(weights))
# }
#
#
# postsimy.glm <- function(x, n = 1L, data = stats::model.frame(x), ...) {
#   X <- model.matrix(delete.response(terms(x)), data = data)
#   params <- postsim.glm(x, n = n, data = data, response = FALSE)
# }

#' @rdname glm-methods
#' @export
postsim_partialfx.glm <- function(x, data1, data2, n = 1L,
                                  delta = 1, response = FALSE, ...) {
  mt <- delete.response(terms(x))
  X1 <- model.matrix(mt, data = data1)
  X2 <- model.matrix(mt, data = data2)
  obs <- nrow(X1)
  param <- postsim(x, n = n)
  array(as_vector(map(param, function(p, X1, X2, response, delta) {
    (ev_glm(p[["beta"]], X2, response = response) -
       ev_glm(p[["beta"]], X1, response = response)) / delta
  }, X1 = X1, X2 = X2, response = response, delta = delta),
  .type = double(obs)), dim = c(obs, n))
}

#' @rdname glm-methods
#' @param V The variance-covariance matrix of the coefficients. This arguments
#'   allows for the substitution of "robust" covariance matrices.
#' @export
postsim.glm <- function(x, n = 1L, V = NULL, ...) {
  summ <- summary(x)
  beta_hat <- coef(x)
  sigma <- rep(sqrt(summ$dispersion), n)
  if (is.null(V)) V <- vcov(x)
  ## TODO parallel process
  map(array_branch(rmvnorm(n, beta_hat, V), margin = 2),
      function(x, sigma) list(beta = x, sigma = sigma),
      sigma = sigma)
}


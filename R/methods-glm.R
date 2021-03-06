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
fdfx.glm <- function(x, data1, data2, delta = 1, n = 1000L,
                          confint = 0.95, response = TRUE, ...) {
  # Difference
  predict_type <- if (response) "response" else "link"
  point_est <- (predict(x, newdata = data2, type = predict_type) -
                predict(x, newdata = data1, type = predict_type)) / delta
  # simulate from posterior to get CI
  sims <- simfdfx.glm(x, data1, data2, n = n,
                      response = response, delta = delta)
  sim_summary(sims, confint, estimate = point_est)
}

#' @param weights Weights to apply to the data
#' @rdname glm-methods
#' @export
afdfx.glm <- function(x, data1, data2, delta = 1, n = 1000L,
                      confint = 0.95, response = TRUE, weights = NULL, ...) {
  predict_type <- if (response) "response" else "link"
  sims <- simfdfx.glm(x, data1, data2, n, response = response,
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
simev.glm <- function(x, data = stats::model.frame(x), response = TRUE,
                          n = 1000L, ...) {
  X <- preprocess_data_lm(x, data)
  params <- simpar.glm(x, n = n, data = data)
  family <- if (response) x$family else NULL
  map(params, function(p, X, family) {ev_glm(p[["beta"]], X, family)}, X = X)
}

# simy_glm_binomial <- function(n, eta, family = binomial(), weights = NULL, ...) {
#   rbinom(length(prob), size = weights, prob = family$linkinv(eta)) / size
# }
#
# simy_glm_gaussian <- function(n, eta, family = gaussian(), weights = NULL, sigma = 1) {
#   rnorm(length(eta), mean = family$linkinv(eta), sd = sigma / sqrt(weights))
# }
#
#
# simy.glm <- function(x, n = 1L, data = stats::model.frame(x), ...) {
#   X <- model.matrix(delete.response(terms(x)), data = data)
#   params <- simpar.glm(x, n = n, data = data, response = FALSE)
# }

#' @rdname glm-methods
#' @export
simfdfx.glm <- function(x, data1, data2, n = 1L,
                                  delta = 1, response = FALSE, ...) {
  X1 <- preprocess_data_lm(x, data1)
  X2 <- preprocess_data_lm(x, data2)
  obs <- nrow(X1)
  param <- simpar(x, n = n)
  family <- if (response) {
    x$family
  } else {
    NULL
  }
  array(as_vector(map(param, function(p, X1, X2, family, delta) {
    (ev_glm(p[["beta"]], X2, family) -
       ev_glm(p[["beta"]], X1, family)) / delta
  }, X1 = X1, X2 = X2, family = family, delta = delta),
  .type = double(obs)), dim = c(obs, n))
}

#' @rdname glm-methods
#' @param V The variance-covariance matrix of the coefficients. This arguments
#'   allows for the substitution of "robust" covariance matrices.
#' @export
simpar.glm <- function(x, n = 1L, V = NULL, ...) {
  summ <- summary(x)
  beta_hat <- coef(x)
  sigma <- sqrt(summ$dispersion)
  if (is.null(V)) V <- vcov(x)
  ## TODO parallel process
  rerun(n, list(beta = setNames(as.numeric(rmvnorm(1, beta_hat, V)),
                                names(beta_hat)), sigma = sigma))
}


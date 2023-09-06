#' Simulates a generalized binary Vector Autoregressive (gbVAR) time series.
#'
#' @param n Sample size
#' @param const Constant vector, Phi0
#' @param phi Matrix of gbVAR coefficient probabilities of multinomial
#'   distribution with parameter 1.
#' @param p The order of the gbVAR process.
#' @param theta Matrix of innovation coefficient probabilities of multinomial
#'   distribution with parameter 1.
#' @param innov d-dimensional binary innovation process
#' @param mu Mean vector of innovation.
#' @param n.start The number of initial data to be omitted. Default is 100.
#'
#' @return time series and innovation
#' @import stats
#' @export
#'
#' @examples
#' A = matrix(c(0.49, 0.35, -0.43, -0.39), nrow = 2, byrow = T)
#' B = diag(c(0.16, 0.18), 2)
#' mu_e = c(0.4, 0.8)
#' result = gbVARsim(n = 100, phi = A, p = 1, theta = B, mu = mu_e)
#' plot(result$series[1,], cex = 0)
#' lines(result$series[1,])
#'
gbVARsim <-
  function(n,
           const = NULL,
           phi,
           p,
           theta,
           innov = NULL,
           mu,
           n.start = 100)
  {
    # Check
    if (any(mu > 1 | mu < 0))
      stop("mean of innovation must be from 0 to 1")
    d = length(mu)

    if (!is.null(innov) && ncol(innov) != n)
      stop("length of innovation must be equal to n")
    nT = n + n.start

    if (is.null(innov)) {
      innov = matrix(0, d, nT)
      for (i in 1:d)
        innov[i, ] = rbinom(nT, 1, mu[i])
    }

    if (nrow(phi) != d)
      stop("dimension of phi is wrong")

    if (nrow(theta) != d)
      stop("dimension of theta is wrong")

    if (any(rowSums(abs(phi)) + rowSums(theta) != 1))
      stop("coefficient matrices are wrong")

    ind1 = (phi >= 0) - (phi < 0)
    ind2 = 1 * (phi < 0)

    ist = p + 1
    result = matrix(0, d, nT)
    if (length(const) == 0)
      const = rep(0, d)

    for (it in ist:nT) {
      coeff = matrix(0, ncol(phi) + 1, d)

      for (i in 1:d) {
        coeff[, i] = rmultinom(1, 1, c(abs(phi[i, ]), theta[i, i]))
      }
      coeff = t(coeff)


      result[, it] = const + (coeff[, -(d + 1)] * ind1) %*%
        as.numeric(result[, (it - p):(it - 1)])
      result[, it] = result[, it] + rowSums(coeff[, -(d + 1)] * ind2)
      result[, it] = result[, it] + coeff[, (d + 1)] * innov[, it]

    }

    result = result[, (1 + n.start):nT]
    innov = innov[, (1 + n.start):nT]
    gbVARsim <- list(series = result, noises = innov)
  }

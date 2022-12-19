CalibrateInformationTheory <- function(x, Lambda, scale = TRUE, gamma = 0.5) {
  # Prepare required objects
  n <- nrow(x)
  p <- ncol(x)
  if (scale) {
    S <- stats::cor(x)
  } else {
    S <- stats::cov(x)
  }

  # Loop over lambda grid
  AIC <- BIC <- EBIC <- rep(NA, length(Lambda))
  path <- array(NA, dim = c(nrow(S), ncol(S), length(Lambda)))
  pb <- utils::txtProgressBar(style = 3)
  for (k in 1:length(Lambda)) {
    utils::setTxtProgressBar(pb, k / length(Lambda))

    # Applying the graphical LASSO
    myglasso <- glassoFast::glassoFast(S, rho = Lambda[k])
    omega <- myglasso$wi

    # Computing the log-likelihood of the model
    loglik <- (n / 2) * (log(det(omega)) - sum(diag(omega %*% S)))

    # Computing the number of edges
    df <- 0.5 * (sum(abs(omega) > 0) - sum(abs(diag(omega)) > 0))

    # Computing the AIC and BIC
    AIC[k] <- -2 * loglik + 2 * df
    BIC[k] <- -2 * loglik + df * log(n)
    EBIC[k] <- -2 * loglik + df * (log(n) + 4 * gamma * log(p))
    # BIC=c(BIC, loglik - 0.5 * df * log(n))
    # AIC=c(AIC, loglik - df)

    # Storing the adjacency matrix
    A <- ifelse(myglasso$wi != 0, yes = 1, no = 0)
    A <- A + t(A)
    A <- ifelse(A != 0, yes = 1, no = 0)
    path[, , k] <- A
  }

  return(list(path = path, AIC = AIC, BIC = BIC, EBIC = EBIC))
}


glasso.graphical_model <- function(x, y, q, scale = TRUE, ...) {
  extraargs <- list(...)
  if (scale) {
    empirical.cov <- stats::cor(x)
  } else {
    empirical.cov <- stats::cov(x)
  }
  lams <- extraargs$lams
  if (is.null(lams)) {
    max.cov <- max(abs(empirical.cov[upper.tri(empirical.cov)]))
    lams <- pulsar::getLamPath(max.cov, max.cov * 0.05, len = 40)
  }
  est <- NULL
  est$X <- NULL
  for (k in 1:length(lams)) {
    # est <- QUIC::QUIC(empirical.cov, rho = 1, path = lams, msg = 0)
    myglasso <- glassoFast::glassoFast(empirical.cov, rho = lams[k])
    est$X <- abind::abind(est$X, myglasso$wi, along = 3)
  }
  ut <- upper.tri(empirical.cov)
  qvals <- sapply(1:length(lams), function(idx) {
    m <- est$X[, , idx]
    sum(m[ut] != 0)
  })
  qq <- qvals >= q
  if (!any(qq)) {
    stop("Didn't reach the required number of variables. Try supplying lambda manually")
  }
  lamidx <- which.max(qvals >= q)
  M <- est$X[, , lamidx][ut]
  selected <- (M != 0)
  s <- sapply(1:lamidx, function(idx) {
    m <- est$X[, , idx][ut] != 0
    return(m)
  })
  colnames(s) <- as.character(1:ncol(s))
  return(list(selected = selected, path = s))
}


class(glasso.graphical_model) <- c("function", "graphical_model")


glasso.pulsar <- function(data, lambda, scale = TRUE) {
  x <- data
  # extraargs <- list(...)

  if (scale) {
    empirical.cov <- stats::cor(x)
  } else {
    empirical.cov <- stats::cov(x)
  }
  lams <- lambda
  if (is.null(lams)) {
    max.cov <- max(abs(empirical.cov[upper.tri(empirical.cov)]))
    lams <- pulsar::getLamPath(max.cov, max.cov * 0.05, len = 40)
  }
  path <- list()
  for (k in 1:length(lams)) {
    # est <- QUIC::QUIC(empirical.cov, rho = 1, path = lams, msg = 0)
    myglasso <- glassoFast::glassoFast(empirical.cov, rho = lams[k])
    A <- ifelse(myglasso$wi != 0, yes = 1, no = 0)
    A <- A + t(A)
    A <- ifelse(A != 0, yes = 1, no = 0)
    path <- c(path, list(A))
  }
  return(list(path = path))
}


# ErrorControl <- function(stability, simul, time,
#                          pi_list = seq(0.6, 0.9, by = 0.05)) {
#   # Re-formatting the list of pi values
#   pi_list <- as.character(pi_list)
#
#   # Calculating selection performances for different pi values
#   perf_mb <- NULL
#   for (k in 1:length(pi_list)) {
#     pi_thr <- pi_list[k]
#
#     # Extracting the ID corresponding to pi
#     pi_id <- which(as.character(stability$params$pi_list) == pi_thr)
#
#     # Extracting the ID corresponding to lambda
#     lambda_id <- max(which(stability$PFER_2d[, pi_id] <= PFER_thr))
#
#     # Extracting selection status of calibrated model
#     if (inherits(stability, "variable_selection")) {
#       selected <- SelectedVariables(stability = stability, argmax_id = c(lambda_id, pi_id))
#     } else {
#       selected <- Adjacency(stability = stability, argmax_id = rbind(c(lambda_id, pi_id)))
#     }
#
#     # Computing selection performances
#     perf <- data.frame(c(
#       pi = pi_thr,
#       SelectionPerformance(theta = selected, theta_star = simul$theta),
#       time = time
#     ))
#
#     # Storing selection performances
#     perf_mb <- rbind(perf_mb, perf)
#   }
#   return(perf_mb)
# }


glmnet.lasso_model <- function(x, y, q, lams = NULL, type = c("conservative", "anticonservative"), ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package ", sQuote("glmnet"), " needed but not available")
  }
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"), " is coerced to a model matrix without intercept")
    x <- stats::model.matrix(~ . - 1, x)
  }
  type <- match.arg(type)
  fit <- suppressWarnings(glmnet::glmnet(x, y,
    pmax = q,
    lambda = lams,
    ...
  ))
  selected <- stats::predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  return(list(selected = ret, path = sequence))
}

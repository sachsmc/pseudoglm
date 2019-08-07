#' Simulate competing risks data
#'
#' @export

gen_data <- function(n = 200, alpha0 = .2, alpha1 = .05, eta = .2, gamma = .25) {

    Z <- rbinom(n, 1, .5)

    U1 <- runif(n)
    U2 <- runif(n)

    T1 <- U1 / (alpha0 + alpha1 * Z)
    T2 <- U2 / (eta)

    TT <- pmin(T1, T2)

    lambda <- -log(1 - gamma) / pmin(TT, 1)

    C <- rexp(n, lambda)

    obsT <- pmin(TT, C)
    indi <- ifelse(C < TT, 0,
                   ifelse(T1 < T2, 1, 2))

    data.frame(obsT, trueObs = T1 < T2 & T1 < 1, indi, Z)

}

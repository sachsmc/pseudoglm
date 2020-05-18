#' Generate survival data for analysis
#'
#' Set up the coefficients for generating competing risks data. The coefficients
#' correspond to the covariate vector with distribution (constant, binomial(.5), N(0,1), N(0,1))
#'
#' @param scenario Character constant determining the scenario
#' @param beta.cens Coefficient vector for the association with the censoring distribution
#'
#' @export
#'
#' @return A list with components for each of the coefficient vectors/matrices
#'


switch_scenario <- function(scenario, beta.cens = c(1,0,0,0)) {


    ## NULL scenario
    if(scenario == "0") {
        b1 <- matrix(c(1, 0, 0, 0), nrow = 4)
        b2 <- matrix(c(1, 0, 0, 0), nrow = 4)

        g1 <- 1.5
        g2 <- 1.5
    } else if (scenario == "1") {
        b1 <- matrix(c(1, .5, 0, 0), nrow = 4)
        b2 <- matrix(c(1, 0, 0, 0), nrow = 4)

        g1 <- 1.5
        g2 <- 1.5
    } else if (scenario == "2") {
        b1 <- matrix(c(1, 1, .5, 0), nrow = 4)
        b2 <- matrix(c(1, 0, 0, 0), nrow = 4)

        g1 <- 1.5
        g2 <- 1.5
    } else if (scenario == "3") {
        b1 <- matrix(c(1, 1.5, .5, 0), nrow = 4)
        b2 <- matrix(c(1, .5, 0, 1), nrow = 4)

        g1 <- 1.5
        g2 <- 1.5

    } else if (scenario == "4") {
        b1 <- matrix(c(1, 1.5, 0, 0), nrow = 4)
        b2 <- matrix(c(1, -1.5, 0, 0), nrow = 4)

        g1 <- 1.5
        g2 <- 1.5

    } else if (scenario == "5") {
        b1 <- matrix(c(1, 2.5, .5, .5), nrow = 4)
        b2 <- matrix(c(1, 1, -1, -1), nrow = 4)

        g1 <- 1.5
        g2 <- 1.5


    } else stop("Scenario does not exist")

    list(b1 = b1, b2 = b2, g1 = g1, g2 = g2, beta.cens = beta.cens)

}


#' Cumulative hazard function under the proportional hazards Weibull model
#'
#' @param t Time at which to be evaluated (must be >= 0)
#' @param Z Covariate matrix
#' @param beta Coefficients to be applied to covariate matrix
#' @param gamma Shape parameter
#'
#' @return Cumulative hazards evaluated at t

Haz11<-function(t,Z, beta,gamma){

    lambda <- exp(Z %*% beta)
    (t / lambda)^gamma

}


#' Hazard function under the proportional hazards Weibull model
#'
#' @param t Time at which to be evaluated (must be >= 0)
#' @param Z Covariate matrix
#' @param beta Coefficients to be applied to covariate matrix
#' @param gamma Shape parameter
#'
#' @return Hazards evaluated at t

haz11 <- function(t, Z, beta, gamma) {

    lambda <- exp(Z %*% beta)
    gamma * t^(gamma - 1) * (1/lambda)^gamma

}

#' Overall survival probability
#'
#' Overall survival probability under a competing risks model with two causes
#'
#' @param t Time at which to be evaluated (must be >= 0)
#' @param Z Covariate matrix
#' @param beta1 Coefficients to be applied to covariate matrix for cause 1
#' @param gamma1 Shape parameter for cause 1
#' @param beta2 Coefficients to be applied to covariate matrix for cause 2
#' @param gamma2 Shape parameter for cause 2
#'
#' @return Survival probabilities evaluated at t

St <- function(t, Z,beta1,gamma1,beta2,gamma2){

    exp(-(Haz11(t,Z, beta1,gamma1)+Haz11(t,Z, beta2,gamma2)))

}

#' Numeric inverse
#'
#' @param fn Function to invert
#' @param mix_x Lower limit of the range of fn
#' @param max_x Upper limit of the range of fn
#'
#' @return A function corresponding to the inverse of fn


inverse <- function(fn, min_x, max_x) {

    fn_inv = function(y){
        uniroot((function(x){fn(x) - y}), lower=min_x, upper=max_x)[1]$root
    }
    return(Vectorize(fn_inv))
}


#' Sample event times
#'
#' Overall survival times in the competing risks model using the
#' probability integral transform.
#'
#' @param mat covariate matrix
#' @param b1 Coefficients for cause 1
#' @param g1 Shape for cause 1
#' @param b2 Coefficients for cause 2
#' @param g2 Shape for cause 2
#'
#' @return A vector of overall survival times
#'

sampletimesfunc <- function(mat,b1,g1,b2,g2) {
    n <- dim(mat)[1]
    UU <- runif(n)

    sapply(1:n, function(i){

        inverse(function(t) 1 - St(t, mat[i,], beta1 = b1,gamma1=g1, beta2=b2, gamma2=g2), 0, 1e8)(UU[i])

    })
}



#' Generate a dataset
#'
#' Simulates a dataset under the competing risks model with two causes
#' and censoring.
#'
#' @param n Sample size
#' @param beta1 Coefficients for cause 1
#' @param gamma1 Shape for cause 1
#' @param beta2 Coefficients for cause 2
#' @param gamma2 Shape for cause 2
#' @param beta.cens Coefficients for censoring
#' @param cens.rate Target censoring rate, used for testing
#'
#' @return A data frame with the simulated data
#' @export
#'

genset_func<-function(n, beta1, gamma1, beta2, gamma2, beta.cens = c(0,0,0),
                      cens.rate = .2){

    tr<-rbinom(n,1,0.5)
    X3 <- matrix(rnorm(n * 2, mean = 0, sd = 1), ncol = 2)


    if(!is.matrix(beta1)) {
        beta100<-c(beta1, rep(0,4-length(beta1)))
        beta10<-as.matrix(beta100,nrow=4)

        beta200<-c(beta2, rep(0,4-length(beta2)))
        beta20<-as.matrix(beta200,nrow=4)
    } else {
        beta10 <- beta1
        beta20 <- beta2
    }
    matx<-cbind(1,tr,X3)

    Touts<-sampletimesfunc(matx,beta10,gamma1,beta20,gamma2)

    cen0 <- uniroot(function(c0) {
        lambda <- exp(matx %*% c(c0, beta.cens))
        Cen<- rweibull(n, scale = lambda, shape = gamma1*2)
        cens.rate - mean(Cen < Touts)
    }, c(-1e4, 100))$root

    lambda <- exp(matx %*% c(cen0, beta.cens))
    Cen<- rweibull(n, scale = lambda, shape = gamma1*2)

    type1<-rbinom(n,1,haz11(Touts,matx,beta10,gamma1)/(haz11(Touts,matx,beta10,gamma1)+
                                                           haz11(Touts,matx,beta20,gamma2)))



    type1<-ifelse(type1 == 1, 1, 2)

    ID<-1:length(type1)

    Toutfin<-pmin(Touts, Cen)

    type1b<-ifelse(Cen<Touts, 0, type1)

    dat1<-as.data.frame(cbind(ID, Toutfin, type1b, tr, X3, Touts, type1))

    return(dat1)
}

#' Simulate a dataset and run the analyses
#'
#' Wrapper for the used scenarios
#'
#' @param n Sample size
#' @param scenario Character identifying the scenario
#' @param beta.cens Coefficients for the censoring distribution
#' @param link Link function

simulate_data <- function(n, scenario = "0", beta.cens = c(0,0,0), cens.rate = .2,
                          link = "identity") {


    b <- switch_scenario(scenario, beta.cens)
    data <- genset_func(n, b$b1, b$g1, b$b2, b$g2, b$beta.cens, cens.rate = cens.rate)

    tmax <- quantile(data$Touts, .9)

    results <- list(
        ci.jack = cumincglm(
            formula = Hist(Toutfin, type1b) ~ tr + V5 + V6,
            time = tmax,
            cause = "1",
            link = link,
            data = data),
        ci.infjack = cumincglm.infjack(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6,
                                       time = tmax, cause = "1", link = link, data = data),
        ci.ipcwaalen = cumincglm.ipcw(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6,
                                      time = tmax, cause = "1", link = link, data = data,
                                      model.censoring = "aareg"),
        ci.ipcwcoxph = cumincglm.ipcw(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6,
                                      time = tmax, cause = "1", link = link, data = data,
                                      model.censoring = "coxph"),
        rmean.jack = rmeanglm(
            formula = Hist(Toutfin, type1b) ~ tr + V5 + V6,
            time = tmax,
            cause = "1",
            link = link,
            data = data),
        rmean.infjack = rmeanglm.infjack(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6,
                                         time = tmax, cause = "1", link = link, data = data),
        rmean.ipcwaalen = rmeanglm.ipcw(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6,
                                     time = tmax, cause = "1", link = link, data = data,
                                     model.censoring = "aareg"),
        rmean.ipcwcoxph = rmeanglm.ipcw(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6,
                                     time = tmax, cause = "1", link = link, data = data,
                                     model.censoring = "coxph")
    )


    nmes <- as.data.frame(do.call(rbind, strsplit(names(results), ".", fixed = TRUE)))
    colnames(nmes) <- c("parameter", "method")

    res2 <- do.call(rbind, lapply(results, function(fit) {

        data.frame(estimate = coefficients(fit)[2],
        p.value = summary(fit)$coefficients[2, 4])

    }))

    res3 <- cbind(res2, nmes)

    res3$n <- n
    res3$cens.rate = cens.rate
    res3$scenario = scenario
    res3$link = link
    res3$beta.cens <- paste(beta.cens, collapse = "-")
    res3

}


#' Get true values of the regression coefficients
#'
#' @param scenario Character string identifying the scenario
#' @param tmax Time horizon
#' @param link Link function
#' @param nchunk Sample size for each chunk
#' @param chunks Number of chunks
#'
#' @return A list with the true coefficients for the cumulative incidence and restricted mean
#' @export

true_values <- function(scenario = "0", tmax, link = "identity",
                        nchunk = 500, chunks = 100) {


    b <- switch_scenario(scenario)

    cumincfunc <- function(t, Xi) {

        obj <- function(u) {
            sapply(u, function(uu) St(uu, Xi, b$b1, b$g1, b$b2, b$g2) *
                       haz11(uu, Xi, b$b1, b$g1))
        }
        integrate(obj, 0, t)$value


    }

    bigdat <- NULL
    n <- nchunk

    for(i in 1:chunks){
        tr<-rbinom(n,1,0.5)
        X3 <- matrix(rnorm(n * 2, mean = 0, sd = 1), ncol = 2)
        bigX <- cbind(1, tr, X3)
        cuminc <- sapply(1:n, function(i) {

            cumincfunc(tmax, bigX[i,])

        })

        rmean <- sapply(1:n, function(i) {

            integrate(function(t){
                sapply(t, function(u) cumincfunc(u, bigX[i,]))
                }, 0, tmax)$value

        })
        bigdat <- rbind(bigdat, data.frame(cuminc = cuminc, rmean = rmean, bigX[, -1]))
    }

    ci.fit <- glm(cuminc ~ tr + V2 + V3, data = bigdat, family = quasi(link = link, variance = "constant"))
    rmean.fit <- glm(rmean ~ tr + V2 + V3, data = bigdat, family = quasi(link = link, variance = "constant"))

    list(ci.coef = ci.fit$coefficients, rmean.coef = rmean.fit$coefficients)


}



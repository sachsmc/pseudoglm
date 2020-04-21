#' Generate survival data for analysis
#'
#' 500 subjects, 20 variables, about XX% cumulative incidence. 5 variables truly associated,
#' some linear, some wonky associations, an interaction, and some variables correlated.
#'
#' n Sample size
#' Character designating the simulation scenario:
#'                 0 = NULL, no association,
#'                 A = simple, one variable is linearly associated with IRIS,
#'                 B = spline of 5
#'                 C = wonky, interactions, nonlinearities, etc.
#'                 D = same as C
#' Proportion of missing binary Y values (random censoring)
#'
#'A data frame with X variables, censored survival times (competing risk w death), and true cumulative incidence at Tstar weeks
#'
#'
#'
#'

switch_scenario <- function(scenario) {


    ## NULL scenario
    if(scenario == "0") {
        b1 <- matrix(c(1, 0, 0, 0), nrow = 4)
        b2 <- matrix(c(1, 0, 0, 0), nrow = 4)
        beta.cens <- c(0,0,0,0)
        g1 <- 1.5
        g2 <- 1.5
    } else if (scenario == "1") {
        b1 <- matrix(c(1, 2, 0, 0), nrow = 4)
        b2 <- matrix(c(1, 0, 0, 0), nrow = 4)
        beta.cens <- c(0,0,0,0)
        g1 <- 1.5
        g2 <- 1.5
    } else if (scenario == "2") {
        ## Scenario 2:
        b1 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "3") {

        ## Scenario 3:

        b1 <- matrix(c(1, 0, 0, .5, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "4") {

        b1 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "5") {

        b1 <- matrix(c(1, 0, 0, 0.5, 1, 1, 0, 0), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, 1, 1, 0, 0), nrow = 8)
        g1 <- 3
        g2 <- 3

    } else if (scenario == "6") {

        b1 <- matrix(c(1, 0, 0, 0.5, 1, 1, 0, 0), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, 1, 1, 0, 0), nrow = 8)
        g1 <- 3
        g2 <- 3

    } else if (scenario == "7"){

        b1 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "8"){

        b1 <- matrix(c(1, 0, 0, .5, 0, 0, 0, .5), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else stop("Scenario does not exist")


    list(b1 = b1, b2 = b2, g1 = g1, g2 = g2, beta.cens = beta.cens)

}



Haz11<-function(t,Z, beta,gamma){

    lambda <- exp(Z %*% beta)
    (t / lambda)^gamma

}

haz11 <- function(t, Z, beta, gamma) {

    lambda <- exp(Z %*% beta)
    gamma * t^(gamma - 1) * (1/lambda)^gamma

}



St <- function(t, Z,beta1,gamma1,beta2,gamma2){

    exp(-(Haz11(t,Z, beta1,gamma1)+Haz11(t,Z, beta2,gamma2)))

}



inverse <- function(fn, min_x, max_x,p) {

    fn_inv = function(y){
        uniroot((function(x){fn(x) - y}), lower=min_x, upper=max_x)[1]$root
    }
    return(Vectorize(fn_inv))
}


sampletimesfunc <- function(mat,b1,g1,b2,g2) {
    n <- dim(mat)[1]
    UU <- runif(n)

    sapply(1:n, function(i){

        inverse(function(t) 1 - St(t, mat[i,], beta1 = b1,gamma1=g1, beta2=b2, gamma2=g2), 0, 1e8)(UU[i])

    })
}





genset_func<-function(n, beta1, gamma1, beta2, gamma2, beta.cens = c(0,0,0,0),
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

    lambda <- exp(matx %*% beta.cens)
    Cen<- rweibull(n, scale = lambda, shape = gamma1*2)

    ## adjust censoring rate
    adj.cenrate <-  cens.rate - mean(Cen < Touts)
    if(adj.cenrate < 0) { # need to remove censoring
        adjdex <- sample(which(Cen < Touts), size = floor(abs(adj.cenrate) * n))
        Cen[adjdex] <- Touts[adjdex] + .1
    } else { # need to add censoring
        adjdex <- sample(which(Cen > Touts), size = floor(adj.cenrate * n))
        Touts[adjdex] <- Cen[adjdex] + .1
    }

    type1<-rbinom(n,1,haz11(Touts,matx,beta10,gamma1)/(haz11(Touts,matx,beta10,gamma1)+
                                                           haz11(Touts,matx,beta20,gamma2)))



    type1<-ifelse(type1 == 1, 1, 2)

    ID<-1:length(type1)

    Toutfin<-pmin(Touts, Cen)

    type1b<-ifelse(Cen<Touts, 0, type1)

    dat1<-as.data.frame(cbind(ID, Toutfin, type1b, tr, X3, Touts, type1))

    return(dat1)
}


simulate_data <- function(n, scenario = "0", cens.rate = .2) {


    b <- switch_scenario(scenario)
    data <- genset_func(n, b$b1, b$g1, b$b2, b$g2, b$beta.cens, cens.rate)
    data

}


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



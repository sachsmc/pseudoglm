#' Generalized linear models for cumulative incidence
#'
#' Using the infintesimal jackknife pseudo observations for the cumulative incidence, this function
#' then runs a generalized linear model and estimates the variance correctly according to
#' Overgaard et al (2018). The link function can be "identity" for estimating differences
#' in the cumulative incidence, "log" for estimating ratios, and any of the other link
#' functions supported by \link[stats]{quasi}.
#'
#' @return A pseudoglm object, with its own methods for print, summary, and vcov. It inherits from
#' glm, so predict and other glm methods are supported.
#'
#' @param formula A formula specifying the model. The left hand side must be a \link[survival]{Surv} object. The right hand side is the usual linear combination of covariates.
#' @param time Numeric constant specifying the time at which the cumulative incidence or survival probability effect estimates are desired.
#' @param cause Character constant specifying the cause indicator of interest.
#' @param link Link function for the cumulative incidence regression model.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param ... Other arguments passed to \link[stats]{glm} such as weights, subset, etc.
#'
#' @export
#'
cumincglm.infjack <- function(formula, time, cause = "1", link = "identity", data, ...) {


    marginal.estimate2 <- survival::survfit(update.formula(formula, . ~ 1),
                                            data = data, influence = TRUE)

    tdex <- sapply(time, function(x) max(which(marginal.estimate2$time <= x)))

    pstate <- marginal.estimate2$pstate[tdex, which(marginal.estimate2$states == cause)]


    ## S(t) + (n-1)[S(t) -S_{-i}(t)]
    jackk2 <- matrix(pstate, nrow = marginal.estimate2$n, ncol = length(time), byrow = TRUE) +
        (marginal.estimate2$n - 1) *
        (marginal.estimate2$influence.pstate[, tdex + 1,
                                             which(marginal.estimate2$states == cause)])

    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))
    newdata[["pseudo.vals"]] <-  c(jackk2)
    newdata[["pseudo.time"]] <- rep(time, each = nrow(jackk2))
    newdata[["cluster.id"]] <- rep(1:nrow(data), length(time))

    if(!inherits(marginal.estimate2, "survfitms")) newdata[["pseudo.vals"]] <- 1 - newdata[["pseudo.vals"]]

    newdata[["startmu"]] <- rep(mean(newdata$pseudo.vals), nrow(newdata))

    fit.lin <- stats::glm(update.formula(formula, pseudo.vals ~ .),
                          family = quasi(link = link, variance = "constant"),
                          mustart = startmu,
                          data = newdata, x = TRUE, ...)

    ## update variance estimate

    if(!inherits(marginal.estimate2, "survfitms")) {

        datamat <- cbind(marginal.estimate$model.response[, "time"],
                         marginal.estimate$model.response[, "status"] != 0, ## not censored indicator
                         as.character(marginal.estimate$model.response[, "status"]) == cause,
                         marginal.estimate$model.response[, "status"] == 0 ## censored indicator
        )

    } else {

        dmatframe <- as.matrix(model.frame(update.formula(formula, .~1), data = newdata)[, 1])

        datamat <- cbind(dmatframe[, "time"],
                         dmatframe[, "status"] != 0,
                         dmatframe[, "status"] == cause,
                         dmatframe[, "status"] == 0)

    }


    #    fit.lin$corrected.vcov <- ovg.vcov
    fit.lin$datamat <- datamat
    fit.lin$time <- time
    fit.lin$cause <- cause
    fit.lin$link <- link
    fit.lin$cluster.id <- newdata[["cluster.id"]]
    fit.lin$type <- "cuminc"

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin

}


#' Generalized linear models for restricted mean survival
#'
#' Using the infintesimal jackknife pseudo observations for the restricted mean
#' survival, this function then runs a generalized linear model
#'  The link function can be "identity" for estimating differences,
#'  "log" for estimating ratios, and any of the other link
#' functions supported by \link[stats]{quasi}.
#'
#' @return A pseudoglm object, with its own methods for print, summary, and vcov. It inherits from
#' glm, so predict and other glm methods are supported.
#'
#' @param formula A formula specifying the model. The left hand side must be a \link[survival]{Surv} object. The right hand side is the usual linear combination of covariates.
#' @param time Numeric constant specifying the maximum time at which the mean is restricted.
#' @param cause Character constant specifying the cause indicator of interest.
#' @param link Link function for the regression model.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param ... Other arguments passed to \link[stats]{glm} such as weights, subset, etc.
#'
#' @export
#'
rmeanglm.infjack <- function(formula, time, cause = "1", link = "identity", data, ...) {


    marginal.estimate2 <- survival::survfit(update.formula(formula, . ~ 1),
                                            data = data, influence = TRUE)

    jackk2 <- pseudo_rmst(marginal.estimate2, time)

    if(!is.null(marginal.estimate2$states)) {
        jackk2 <- jackk2[, which(marginal.estimate2$states == cause)]
    }


    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))
    newdata[["pseudo.vals"]] <-  c(jackk2)
    newdata[["pseudo.time"]] <- rep(time, each = length(jackk2))
    newdata[["cluster.id"]] <- rep(1:nrow(data), length(time))

    newdata[["startmu"]] <- rep(mean(newdata$pseudo.vals), nrow(newdata))

    fit.lin <- stats::glm(update.formula(formula, pseudo.vals ~ .),
                          family = quasi(link = link, variance = "constant"),
                          mustart = startmu,
                          data = newdata, x = TRUE, ...)

    ## update variance estimate

    if(!inherits(marginal.estimate2, "survfitms")) {

        datamat <- cbind(marginal.estimate2$model.response[, "time"],
                         marginal.estimate2$model.response[, "status"] != 0, ## not censored indicator
                         as.character(marginal.estimate2$model.response[, "status"]) == cause,
                         marginal.estimate2$model.response[, "status"] == 0 ## censored indicator
        )

    } else {

        dmatframe <- as.matrix(model.frame(update.formula(formula, .~1), data = newdata)[, 1])

        datamat <- cbind(dmatframe[, "time"],
                         dmatframe[, "status"] != 0,
                         dmatframe[, "status"] == cause,
                         dmatframe[, "status"] == 0)

    }


    #    fit.lin$corrected.vcov <- ovg.vcov
    fit.lin$datamat <- datamat
    fit.lin$time <- time
    fit.lin$cause <- cause
    fit.lin$link <- link
    fit.lin$cluster.id <- newdata[["cluster.id"]]
    fit.lin$type <- "rmean"

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin

}

#' Compute restricted mean survival
#'
#' Written by Terry Therneau (2019)
#'
#' @param fit Survfit
#' @param tmax Max time
pseudo_rmst <- function(fit, tmax) {
    # extract the RMST, and the leverage of each subject on the RMST
    if (missing(tmax)) tmax <- max(sfit0$time)
    if (length(tmax) !=1) stop("tmax must be a single value")

    sfit0 <- survfit0(fit)      # add the time=0 point to the result
    rsum <- function(y, x= sfit0$time) {  # sum of rectangles
        keep <- which(x < tmax)
        width <- diff(c(x[keep], tmax))
        sum(width * y[keep])
    }

    if (is.null(sfit0$states)) { # ordinary survival
        rmst <- rsum(sfit0$surv, sfit0$time)
        ijack <- apply(sfit0$influence.surv, 1, rsum)
        rmst + (sfit0$n -1)* ijack
    }
    else {
        rmst <- apply(sfit0$pstate, 2, rsum)
        ijack <- apply(sfit0$influence.pstate, c(1,3), rsum)
        rep(rmst, each= nrow(ijack)) + (nrow(ijack)-1)*ijack
    }
}


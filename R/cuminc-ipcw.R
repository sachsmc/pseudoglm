#' Generalized linear models for cumulative incidence
#'
#' Using pseudo observations for the cumulative incidence, this function then runs a generalized
#' linear model and estimates the variance correctly according to Overgaard et al (2018). The
#' link function can be "identity" for estimating differences in the cumulative incidence, "log"
#' for estimating ratios, and any of the other link functions supported by \link[stats]{quasi}.
#'
#' @return A pseudoglm object, with its own methods for print, summary, and vcov. It inherits from glm, so predict and other glm methods are supported.
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
cumincglm.ipcw <- function(formula, time, cause = "1", link = "identity", data,
                           model.censoring = "aareg", formula.censoring = NULL, ...) {

    outcome <- model.response(model.frame(update.formula(formula, . ~ 1), data = data))
    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))

    newdata$.Ci <- as.numeric(outcome[, "status"] == 0)
    newdata$.Tci <- outcome[, "time"]
    if(is.null(formula.censoring)) {
        cens.formula <- update.formula(formula, survival::Surv(.Tci, .Ci) ~ .)
    } else {
        cens.formula <- formula.censoring
    }


    if(model.censoring == "aareg") {

        predmat <- model.matrix(cens.formula, data = newdata)

        fitcens <- survival::aareg(cens.formula, data = newdata)

        tdex <- sapply(pmin(newdata$.Tci, time), function(t) max(c(1, which(fitcens$times <= t))))
        Gi <- rep(NA, length(tdex))
        for(i in 1:length(tdex)) {

            Gi[i] <- prod(1 - c(fitcens$coefficient[1:tdex[i], ] %*% t(predmat[i, , drop = FALSE])))

        }

    } else if(model.censoring == "coxph") {

        fitcens <- survival::coxph(cens.formula, data = newdata, x = TRUE)
        coxsurv <- survival::survfit(fitcens, newdata = newdata)
        tdex <- sapply(pmin(newdata$.Tci, time), function(t) max(c(1, which(coxsurv$time <= t))))
        Gi <- coxsurv$surv[cbind(tdex,1:ncol(coxsurv$surv))]

    }
    Vi <- as.numeric(outcome[, "time"] < time & outcome[, "status"] == cause)
    Ii <- as.numeric(outcome[, "time"] >= time | outcome[, "status"] != 0)

    nn <- length(Vi)
    theta.n <- mean(Ii * Vi / Gi)

    XXi <- Vi * Ii / Gi
    POi <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))

    newdata[["pseudo.vals"]] <-  POi
    newdata[["pseudo.time"]] <- rep(time, each = length(POi))
    newdata[["cluster.id"]] <- rep(1:nrow(data), length(time))

    newdata[["startmu"]] <- rep(mean(newdata$pseudo.vals), nrow(newdata))

    fit.lin <- stats::glm(update.formula(formula, pseudo.vals ~ .),
                          family = quasi(link = link, variance = "constant"),
                          mustart = startmu,
                          data = newdata, x = TRUE, ...)

    ## update variance estimate
    datamat <- cbind(outcome[, "time"],
                         outcome[, "status"] != 0, ## not censored indicator
                         as.character(outcome[, "status"]) == cause,
                         outcome[, "status"] == 0 ## censored indicator
                    )



    #    fit.lin$corrected.vcov <- ovg.vcov
    fit.lin$datamat <- datamat
    fit.lin$time <- time
    fit.lin$cause <- cause
    fit.lin$link <- link
    fit.lin$cluster.id <- newdata[["cluster.id"]]

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin


}

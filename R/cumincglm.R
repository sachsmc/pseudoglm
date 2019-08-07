#' Generalized linear models for cumulative incidence
#'
#' @param formula A formula specifying the model. The left hand side must be a \link[prodlim]{Hist} object. The right hand side is the usual linear combination of covariates.
#' @param time Numeric constant specifying the time at which the cumulative incidence or survival probability effect estimates are desired.
#' @param cause Character constant specifying the cause indicator of interest.
#' @param link Link function for the cumulative incidence regression model.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param ... Other arguments passed to \link[stats]{glm} such as weights, subset, etc.
#'
#' @export
#'
cumincglm <- function(formula, time, cause = "1", link = "identity", data, ...) {

    marginal.estimate <- prodlim::prodlim(update.formula(formula, . ~ 1), data = data)
    newdata <- data



    newdata[["pseudo.vals"]] <- c(prodlim::jackknife(marginal.estimate, times = time, cause = cause))
    if(marginal.estimate$model == "survival") newdata[["pseudo.vals"]] <- 1 - newdata[["pseudo.vals"]]

    newdata[["startmu"]] <- rep(mean(newdata$pseudo.vals), nrow(newdata))

    fit.lin <- stats::glm(update.formula(formula, pseudo.vals ~ .),
                         family = quasi(link = link, variance = "constant"),
                         mustart = startmu,
                         data = newdata, x = TRUE, ...)

    ## update variance estimate

    if(marginal.estimate$model == "survival") {

        datamat <- cbind(marginal.estimate$model.response[, "time"],
                         marginal.estimate$model.response[, "status"] != 0, ## not censored indicator
                         as.character(marginal.estimate$model.response[, "status"]) == cause,
                         marginal.estimate$model.response[, "status"] == 0 ## censored indicator
        )

    } else {
        datamat <- cbind(marginal.estimate$model.response[, "time"],
                         marginal.estimate$model.response[, "status"] != 0, ## not censored indicator
                         c(attr(marginal.estimate$model.response, "states"), "0")[marginal.estimate$model.response[, "event"]] == cause,
                         marginal.estimate$model.response[, "status"] == 0 ## censored indicator
        )
    }

    noobs <- len <- nrow(datamat)

    datord <- order(datamat[, 1], -datamat[, 2], -datamat[, 3], -datamat[, 4])
    datamat <- datamat[datord, ]

    datain <- datamat[datamat[, 1] <= time, ]

    len2 <- nrow(datain)

    timejump <- datain[1:(len2 - 1), 1] != datain[2:len2, 1]
    times_id <- c(which(timejump == TRUE), len2)

    n_times <- length(times_id)
    times_nr <- cumsum(c(1, timejump))

    Y_lag <- (noobs:1)[times_id]

    N_all <- datain[, 2]
    N_1 <- datain[, 3]
    N_0 <- datain[, 4]

    # Calculating the overall survival function (missing values mean a population die out and should be interpreted as 0)
    H_0 <- cumsum(N_0)/noobs
    H_0_dif = c(H_0[times_id[1]],  ((H_0[times_id])[2:n_times] - (H_0[times_id])[1:(n_times-1)]))


    H_1 = cumsum(N_1)/noobs
    H_1_dif = c(H_1[times_id[1]] , ((H_1[times_id])[2:n_times] - (H_1[times_id])[1:(n_times-1)]))

    H_all = cumsum(N_all)/noobs
    H_all_dif = c(H_all[times_id[1]] , ((H_all[times_id])[2:n_times] - (H_all[times_id])[1:(n_times-1)]))

    H_lag = Y_lag / noobs
    H_lag_inv = 1 / H_lag
    H_lag_inv[is.na(H_lag_inv)] <- 0


    Lambda_0_dif = H_0_dif / H_lag
    Lambda_1_dif = H_1_dif / H_lag
    Lambda_all_dif = H_all_dif / H_lag
    Lambda_0_dif_comp_inv = 1 / (1 - Lambda_0_dif)
    Lambda_0_dif_comp_inv[is.na(Lambda_0_dif_comp_inv)] <- 0

    S = exp(cumsum(log(1 - Lambda_all_dif)))
    S[is.na(S)] <- 0
    S_lag = c(1 , S[1:(n_times-1)])

    G = exp(cumsum(log(1 - Lambda_0_dif)))
    G[is.na(G)] <- 0
    G_lag = c(1 , G[1:(n_times-1)])
    G_lag_inv = 1 / G_lag
    G_lag_inv[is.na(G_lag_inv)] <- 0

    F_1 = cumsum(S_lag * Lambda_1_dif)

    d_phi_1 = c((N_1 * G_lag_inv[times_nr]) , rep(0, noobs - len2))
    d_phi_2 = c(N_0 * (F_1[n_times] - F_1[times_nr]) * Lambda_0_dif_comp_inv[times_nr] * H_lag_inv[times_nr] -
                    cumsum(H_0_dif * (F_1[n_times] - F_1) * Lambda_0_dif_comp_inv * H_lag_inv^2)[times_nr] ,
                rep(-sum(H_0_dif * (F_1[n_times] - F_1) * Lambda_0_dif_comp_inv * H_lag_inv^2), noobs-len2))

    beta = fit.lin$coefficients
    k = length(beta)

    z = fit.lin$x[datord,]  ## what if the model contains po at multiple time points?

    muhat = fit.lin$family$linkinv(z %*% beta)
    Ahat = z * as.vector(fit.lin$family$mu.eta(z %*% beta))
    Ahat_red = Ahat[1:len2,]
    mu_derivhat = z * as.vector(fit.lin$family$mu.eta(z %*% beta))

    a_len <- a1_len <- a2_len <- a3_len <-
        b_1 <- b_2 <- b_3 <- b_4 <- b_5 <- b_6 <-
        matrix(NA, nrow = noobs, ncol = k)
    H_0z <- H_0z_dif <- H_1z <- H_1z_dif <- H_z <- H_z_lag <- matrix(NA, nrow = len2, ncol = k)


    ## compute variance
    for(j in 1:k) {
        a_len[,j] = Ahat[,j] * (d_phi_1-muhat+d_phi_2)
        a1_len[,j] = Ahat[,j]  * (d_phi_1)
        a2_len[,j] = Ahat[,j] * (-muhat)
        a3_len[,j] = Ahat[,j] * (d_phi_2)

        H_0z[,j] = cumsum(Ahat_red[,j] * N_0)/noobs
        H_0z_dif[,j] = c(H_0z[times_id[1],j] , ((H_0z[times_id,j])[2:n_times] - (H_0z[times_id,j])[1:(n_times-1)]))

        H_1z[,j] = cumsum(Ahat_red[,j] * N_1)/noobs
        H_1z_dif[,j] = c(H_1z[times_id[1],j] , ((H_1z[times_id,j])[2:n_times] - (H_1z[times_id,j])[1:(n_times-1)]))

        H_z[,j] = (sum(Ahat[,j]) - cumsum(Ahat_red[,j]))/noobs
        H_z_lag[,j] = c(sum(Ahat[,j]) / noobs , (H_z[1:(len2-1),j])[times_id[1:(n_times-1)]])

        temp1 = cumsum( H_0z_dif[,j] * H_lag_inv * Lambda_0_dif_comp_inv - H_z_lag[,j] * H_lag_inv^2 * Lambda_0_dif_comp_inv * H_0_dif )
        temp1_lag = c(0 , temp1[1:(n_times-1)])
        temp1_dif = temp1 - temp1_lag

        temp2 = cumsum( H_0z_dif[,j] * (F_1[n_times] - F_1) * H_lag_inv * Lambda_0_dif_comp_inv -
                            H_z_lag[,j] * (F_1[n_times] - F_1) * H_lag_inv^2 * Lambda_0_dif_comp_inv  * H_0_dif )
        temp3 = cumsum( H_1z_dif[,j] / G_lag )


        b_1[,j] = c(N_1 * G_lag_inv[times_nr] * temp1_lag[times_nr] , rep(0, len-len2))

        b_2[,j] = c(N_0 * (temp3[n_times] - temp3[times_nr]) * H_lag_inv[times_nr] * Lambda_0_dif_comp_inv[times_nr] -
                        cumsum( H_0_dif * (temp3[n_times] - temp3) * H_lag_inv^2 * Lambda_0_dif_comp_inv )[times_nr] ,
                    rep(-sum( H_0_dif * (temp3[n_times] - temp3) * H_lag_inv^2 * Lambda_0_dif_comp_inv ), len-len2))

        temp4 = cumsum( H_1_dif * G_lag_inv * temp1 )

        b_3[,j] = c(N_0 * (temp4[n_times] - temp4[times_nr]) * H_lag_inv[times_nr] * Lambda_0_dif_comp_inv[times_nr] -
                        cumsum( H_0_dif * (temp4[n_times] - temp4) * H_lag_inv^2 * Lambda_0_dif_comp_inv )[times_nr] ,
                    rep(-sum( H_0_dif * (temp4[n_times] -temp4) / H_lag^2 / (1-Lambda_0_dif) ), len-len2 ))


        b_4[,j] = c(- N_0 * (F_1[n_times] - F_1[times_nr]) * H_z_lag[times_nr,j] * H_lag_inv[times_nr]^2 * Lambda_0_dif_comp_inv[times_nr] +
                        cumsum( H_0_dif * (F_1[n_times] - F_1) * H_z_lag[,j] * H_lag_inv^3 * Lambda_0_dif_comp_inv )[times_nr] ,
                    rep(sum( H_0_dif * (F_1[n_times] - F_1) * H_z_lag[,j] * H_lag_inv^3 * Lambda_0_dif_comp_inv ), len-len2 ))


        b_5[,j] = c(-cumsum( H_0z_dif[,j] * (F_1[n_times] - F_1) * H_lag_inv^2 * Lambda_0_dif_comp_inv -
                                 H_z_lag[,j] * (F_1[n_times] - F_1) * H_lag_inv^3 * Lambda_0_dif_comp_inv * H_0_dif )[times_nr] ,
                    rep(- sum( H_0z_dif[,j] * (F_1[n_times] - F_1) * H_lag_inv^2 * Lambda_0_dif_comp_inv -
                                   H_z_lag[,j] * (F_1[n_times] - F_1) * H_lag_inv^3 * Lambda_0_dif_comp_inv * H_0_dif ), len-len2 ))

        b_6[,j] = c(N_0 * (F_1[n_times] - F_1[times_nr]) * temp1_dif[times_nr] * H_lag_inv[times_nr] * Lambda_0_dif_comp_inv[times_nr] -
                        cumsum( H_0_dif * (F_1[n_times] - F_1) * temp1_dif * H_lag_inv^2 * Lambda_0_dif_comp_inv )[times_nr] ,
                    rep(-sum( H_0_dif * (F_1[n_times] - F_1) * temp1_dif * H_lag_inv^2 * Lambda_0_dif_comp_inv ), len-len2 ))


    }

    b <- b_1 + b_2 + b_3 + b_4 + b_5 + b_6
    Sigma <- t(a_len) %*% a_len / noobs + t(a_len) %*% b / noobs +
        t(b) %*% a_len / noobs + t(b) %*% b / noobs


    Minvhat <- solve(t(Ahat) %*% mu_derivhat / nrow(Ahat))

    ovg.vcov <- Minvhat %*% Sigma %*% t(Minvhat) / noobs

    list(fit.lin, ovg.vcov)

}

## is estimated variance close to true?

nsim <- 1000
estvar.z <- rep(NA, nsim)
estbet.z <- rep(NA, nsim)

for(i in 1:nsim) {
    test <- gen_data(n=200)
    fit <- cumincglm(
        formula = Hist(obsT, indi) ~ Z,
        time = 2,
        cause = "2",
        link = "log",
        data = test)

    estvar.z[i] <- sqrt(diag(vcov(fit)))[2]
    estbet.z[i] <- fit$coefficients[2]

}


sqrt(var(estbet.z))
mean(estvar.z)

## works for standard survival outcome?


nsim <- 1000
estvar.z <- rep(NA, nsim)
estbet.z <- rep(NA, nsim)

for(i in 1:nsim) {
    test <- gen_data(n=200)
    test$indi[test$indi == 2] <- 1

    fit <- cumincglm(formula = Hist(obsT, indi) ~ Z, time = 2, cause = "1", link = "log", data = test)

    estvar.z[i] <- sqrt(diag(vcov(fit)))[2]
    estbet.z[i] <- fit$coefficients[2]

}

sqrt(var(estbet.z))
mean(estvar.z)

## comparison of variance types


nsim <- 1000
estvar <- as.data.frame(matrix(NA, nrow = nsim * 3, ncol = 3))
estbet<- matrix(NA, nrow = nsim, ncol = 2)

j <- 1
for(i in 1:nsim) {
    test <- gen_data(n=200)

    fit <- cumincglm(formula = Hist(obsT, indi) ~ Z,
                     time = 1, cause = "1", link = "log", data = test)

    estvar[j,] <- c(sqrt(diag(vcov(fit, type = "corrected"))), 1)
    estvar[j + 1, ] <- c(sqrt(diag(vcov(fit, type = "robust"))), 2)
    estvar[j + 2, ] <- c(sqrt(diag(vcov(fit, type = "naive"))), 3)

    j <- j + 3

    estbet[i,] <- fit$coefficients

}

colnames(estvar) <- c("intercept", "Z", "type")

sqrt(diag(cov(estbet)))

library(dplyr)
estvar %>% group_by(type) %>% summarize(mint = mean(intercept), mz = mean(Z))


library(ggplot2)
ggplot(estvar, aes(x = factor(type), y = Z)) + geom_boxplot()




nsim <- 1000
estbet.z <- matrix(NA, ncol = 4, nrow = nsim)
se.inf <- matrix(NA, ncol = 4, nrow = nsim)
for(i in 1:nsim) {
    test <- gen_data(n=200)

    fit2 <- cumincglm.infjack(
        formula = Surv(obsT, factor(indi)) ~ Z,
        time = 1,
        cause = "1",
        link = "identity",
        data = test)

    fit <- cumincglm(formula = Hist(obsT, indi) ~ Z,
                     time = 1, cause = "1", link = "identity", data = test)



    estbet.z[i, ] <- c(fit2$coefficients, fit$coefficients)
    se.inf[i, ] <- sqrt(diag(vcov(fit2, type = "cluster")))

}

se.emp <- sqrt(diag(cov(estbet.z)))
par(mfrow = c(2, 2))
for(i in 1:4) {

    hist(se.inf[, i])
    abline(v = se.emp[i], col = "red")

}



nsim <- 1000
estbet.z <- matrix(NA, ncol = 4, nrow = nsim)

for(i in 1:nsim) {
    test <- gen_data(n=200)

    fit2 <- cumincglm.infjack(
        formula = Surv(obsT, factor(indi)) ~ Z,
        time = 1,
        cause = "1",
        link = "identity",
        data = test)

    fit <- cumincglm(formula = Hist(obsT, indi) ~ Z,
                     time = 1, cause = "1", link = "identity", data = test)

    fit3 <- cumincglm.ipcw(formula = Surv(obsT, factor(indi)) ~ Z,
                     time = 1, cause = "1", link = "identity", data = test)


    fit4 <- cumincglm.ipcw(formula = Surv(obsT, factor(indi)) ~ Z,
                           time = 1, cause = "1", link = "identity", data = test,
                           model.censoring = "coxph")

    estbet.z[i, ] <- c(infjack = fit2$coefficients[2],
                       regjack = fit$coefficients[2],
                       ipcw.aareg =  fit3$coefficients[2],
                       ipcw.coxph = fit4$coefficients[2])

}

estbeta <- as.data.frame(estbet.z)
colnames(estbeta) = c("infjack", "regjack", "ipcw.aareg", "ipcw.cox")

library(ggplot2)
ggplot(tidyr::gather(estbeta), aes(x = value, y = key)) + geom_violin()
summary(estbeta)



nsim <- 500
estbet.z <- matrix(NA, ncol = 5, nrow = nsim)

for(i in 1:nsim) {
    test <- generate_data(n=400, scenario = "D")

    fit2 <- cumincglm.infjack(
        formula = Surv(Tout, factor(delta)) ~ X1 + X2,
        time = 26.5,
        cause = "1",
        link = "identity",
        data = test)

    fit <- cumincglm(formula = Hist(Tout, delta) ~ X1 + X2,
                     time = 26.5, cause = "1", link = "identity", data = test)

    fit3 <- cumincglm.ipcw(formula = Surv(Tout, factor(delta)) ~ X1 + X2,
                           time = 26.5, cause = "1", link = "identity", data = test)


    fit4 <- cumincglm.ipcw(formula = Surv(Tout, factor(delta)) ~ X1 + X2,
                           time = 26.5, cause = "1", link = "identity", data = test,
                           model.censoring = "coxph")

    fittrue <- lm(trueT ~ X1 + X2, data = test)

    estbet.z[i, ] <- c(infjack = fit2$coefficients[2],
                       regjack = fit$coefficients[2],
                       ipcw.aareg =  fit3$coefficients[2],
                       ipcw.coxph = fit4$coefficients[2],
                       true = fittrue$coefficients[2])

}

estbeta <- as.data.frame(estbet.z[, 1:4] - estbet.z[, 5])
colnames(estbeta) = c("infjack", "regjack", "ipcw.aareg", "ipcw.cox")

library(ggplot2)
ggplot(tidyr::gather(estbeta), aes(x = value, y = key))  +
    geom_vline(xintercept = 0, col = "red") +
    geom_boxplot()
summary(estbeta)

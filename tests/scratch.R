
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

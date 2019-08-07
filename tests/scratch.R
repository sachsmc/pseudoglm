
## is estimated variance close to true?

nsim <- 1000
estvar.z <- rep(NA, nsim)
estbet.z <- rep(NA, nsim)

for(i in 1:nsim) {
    test <- gen_data(n=200)
    fit <- cumincglm(formula = Hist(obsT, indi) ~ Z, time = 2, cause = "1", link = "log", data = test)

    estvar.z[i] <- sqrt(diag(fit[[2]]))[2]
    estbet.z[i] <- fit[[1]]$coefficients[2]

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

    estvar.z[i] <- sqrt(diag(fit[[2]]))[2]
    estbet.z[i] <- fit[[1]]$coefficients[2]

}

sqrt(var(estbet.z))
mean(estvar.z)


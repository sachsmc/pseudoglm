
estbetas1 <- matrix(NA, nrow = 100, ncol = 4)
estbetas2 <- matrix(NA, nrow = 100, ncol = 4)
for(i in 1:100) {


    test <- simulate_data(200, scenario = "1")
    estbetas1[i, ] <- cumincglm.infjack(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6, data = test,
                      cause = "1", time = 3)$coefficients
    estbetas2[i, ] <- rmeanglm.infjack(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6, data = test,
                      cause = "1", time = 3)$coefficients


}

true.values <- true_values("1", 3, nchunk = 400, chunks = 5)


hist(estbetas1[, 2])
abline(v = true.values$ci.coef[2])

hist(estbetas2[, 2])
abline(v = true.values$rmean.coef[2])

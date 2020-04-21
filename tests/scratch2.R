
estbetas1 <- matrix(NA, nrow = 100, ncol = 4)
estbetas2 <- matrix(NA, nrow = 100, ncol = 4)
for(i in 1:100) {


    test <- simulate_data(200, scenario = "1")

    estbetas1[i, ] <- cumincglm.infjack(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6, data = test,
                      cause = "1", time = 6, link = "log")$coefficients
    estbetas2[i, ] <- rmeanglm.infjack(Surv(Toutfin, factor(type1b)) ~ tr + V5 + V6, data = test,
                      cause = "1", time = 6, link = "log")$coefficients

    if(estbetas2[i, 2] < -1000 | estbetas1[i , 2] < -1000) break

}

true.values <- true_values("1", 6, nchunk = 400, chunks = 5, link = "log")


hist(estbetas1[-96, 2])
abline(v = true.values$ci.coef[2])

hist(estbetas2[-96, 2])
abline(v = true.values$rmean.coef[2])

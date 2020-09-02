devtools::load_all(".")
library(parallel)

setup <- expand.grid(nn = c(500, 1000),
                     cens.rate = c(0.3, .5, .8),
                     b2.cens = c(0, 1),
                     b3.cens = c(0, .5),
                     b4.cens = c(0, .5),
                     gamma.cens = c(1,3),
                     link = c("identity"),
                     scenario = c("2", "3", "4"),
                     stringsAsFactors = FALSE)


library(parallel)
cl <- makeCluster(8)

init <- clusterEvalQ(cl, {
    devtools::load_all(".")

    setup <- expand.grid(nn = c(500, 1000),
                         cens.rate = c(0.3, .5, .8),
                         b2.cens = c(0, 1),
                         b3.cens = c(0, .5),
                         b4.cens = c(0, .5),
                         gamma.cens = c(1,3),
                         link = c("identity"),
                         scenario = c("2", "3", "4"),
                         stringsAsFactors = FALSE)

})



for(j in 1:nrow(setup)) {

    clusterExport(cl, "j")
    res <- do.call(rbind, clusterApplyLB(cl, 101:1000, function(i){

        inres <- tryCatch(simulate_data(setup$nn[j], scenario = setup$scenario[j],
                                        beta.cens = c(setup$b2.cens[j],setup$b3.cens[j],setup$b4.cens[j]),
                                        cens.rate = setup$cens.rate[j], gamma.cens = setup$gamma.cens[j],
                                        link = setup$link[j]),
                          error = function(e) {
                              data.frame(estimate = NA, std.err.rob = NA, std.err.cor = NA,
                                         std.err.nai = NA,
                                         parameter = NA, method = "failed",
                                         n = setup$nn[j],
                                         cens.rate = setup$cens.rate[j], scenario = setup$scenario[j],
                                         link = setup$link[j],
                                         beta.cens = paste(c(setup$b2.cens[j],setup$b3.cens[j],setup$b4.cens[j]),
                                                                               collapse = "-"),
                                         gamma.cens = setup$gamma.cens[j])
                          })

        inres$replicate <- i
        inres

    }))

    if(FALSE) {
        write.csv(res, "sims-type-2.csv", row.names = FALSE)
    } else {
        write.table(res, "sims-type-2.csv", row.names = FALSE, append = TRUE,
                    col.names = FALSE, sep = ",")
    }
}

stopCluster(cl)


true.values <- NULL
for(li in c("identity")) {
    for(sc in c( "2", "3", "4")){

        tv <- true_values(scenario = sc, link = li,
                          nchunk = 500, chunks = 100)

        true.values <- rbind(true.values,
                             data.frame(scenario = sc, link = li,
                                        tv.cuminc = tv$ci.coef[2],
                                        tv.rmean = tv$rmean.coef[2]))
    }}


write.csv(true.values, "true-coefficients.csv")

### analysis in scratch2


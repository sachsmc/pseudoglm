devtools::load_all(".")
library(parallel)

setup <- expand.grid(nn = c(500, 1000),
                     cens.rate = c(0.3, .5, .8),
                     b2.cens = c(0, 1),
                     b3.cens = c(0, .5),
                     b4.cens = c(0, .5),
                     gamma.cens = c(1,3),
                     link = c("identity"), stringsAsFactors = FALSE)


library(parallel)
cl <- makeCluster(8)

init <- clusterEvalQ(cl, {
    devtools::load_all(".")

    setup <- expand.grid(nn = c(500, 1000),
                         cens.rate = c(0.3, .5, .8),
                         b2.cens = c(0, 1),
                         b3.cens = c(0, .5),
                         b4.cens = c(0, .5),
                         gamma.cens = c(1, 3),
                         link = c("identity"), stringsAsFactors = FALSE)
})



for(j in 1:nrow(setup)) {

    clusterExport(cl, "j")
    res <- do.call(rbind, clusterApplyLB(cl, 1:100, function(i){

        inres <- tryCatch(simulate_data(setup$nn[j], scenario = "0",
                                        beta.cens = c(setup$b2.cens[j],setup$b3.cens[j],setup$b4.cens[j]),
                                        cens.rate = setup$cens.rate[j], gamma.cens = setup$gamma.cens[j],
                                        link = setup$link[j]),
                          error = function(e) {
                              data.frame(estimate = NA, std.err.rob = NA, std.err.cor = NA,
                                         std.err.nai = NA,
                                         parameter = NA, method = "failed",
                                         n = setup$nn[j],
                                         cens.rate = setup$cens.rate[j], scenario = "0",
                                         link = setup$link[j],
                                         beta.cens = paste(c(setup$b2.cens[j],setup$b3.cens[j],setup$b4.cens[j]),
                                                                                 collapse = "-"),
                                         gamma.cens = setup$gamma.cens[j])
                          })

        inres$replicate <- i
        inres

    }))

    if(j == 1) {
        write.csv(res, "tests/sims-type-1.csv", row.names = FALSE)
    } else {
        write.table(res, "tests/sims-type-1.csv", row.names = FALSE, append = TRUE,
                    col.names = FALSE, sep = ",")
    }
}

stopCluster(cl)
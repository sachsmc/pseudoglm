devtools::load_all(".")
library(parallel)

setup <- expand.grid(nn = c(500, 1000, 2500),
                     cens.rate = c(0.1, .2, .3),
                     b2.cens = c(0, .5, 1),
                     link = c("identity", "log"), stringsAsFactors = FALSE)


library(parallel)
cl <- makeCluster(8)

init <- clusterEvalQ(cl, {
    devtools::load_all(".")

    setup <- expand.grid(nn = c(500, 1000, 2500),
                         cens.rate = c(0.1, 0.2, 0.3),
                         b2.cens = c(0, .5, 1),
                         link = c("identity", "log"), stringsAsFactors = FALSE)
})



for(j in 1:nrow(setup)) {

    clusterExport(cl, "j")
    res <- do.call(rbind, clusterApplyLB(cl, 1:100, function(i){

        inres <- tryCatch(simulate_data(setup$nn[j], scenario = "0", beta.cens = c(0,setup$b2.cens[j],0),
                                        cens.rate = setup$cens.rate[j], link = setup$link[j]),
                          error = function(e) {
                              data.frame(estimate = NA, p.value = NA, parameter = NA, method = "failed",
                                         n = setup$nn[j],
                                         cens.rate = setup$cens.rate[j], scenario = "0",
                                         link = setup$link[j], beta.cens = paste(c(0, setup$b2.cens[j], 0),
                                                                                 collapse = "-"))
                          })

        inres$replicate <- i
        inres

    }))

    if(j == 1) {
        write.csv(res, "sims-type-1.csv", row.names = FALSE)
    } else {
        write.table(res, "sims-type-1.csv", row.names = FALSE, append = TRUE,
                    col.names = FALSE, sep = ",")
    }
}

stopCluster(cl)

library(data.table)
library(ggplot2)
library(xtable)

sim1 <- data.table(read.csv("tests/sims-type-1.csv"))

res1 <- sim1[!is.na(estimate), .(mest = median(estimate, na.rm = TRUE),
                 lower = quantile(estimate, .05) ,
                 upper = quantile(estimate, .95)),
     by = .(scenario, n, parameter, cens.rate, beta.cens, gamma.cens, link, method)]

res1[, method := factor(method, levels = rev(c("jack", "strat", "ipcwcoxph", "ipcwaalen")),
                        ordered = TRUE)]


sim1[, method := factor(method, levels = rev(c("jack", "strat", "ipcwcoxph", "ipcwaalen")),
                        ordered = TRUE)]


type1 <- sim1[!is.na(estimate), .(type1err = c(
  mean(2 * pnorm(-abs(estimate / std.err.rob)) < .05, na.rm = TRUE),
  mean(2 * pnorm(-abs(estimate / std.err.cor)) < .05, na.rm = TRUE),
  mean(2 * pnorm(-abs(estimate / std.err.nai)) < .05, na.rm = TRUE)),
  se.method = c("robust", "corrected", "naive")
  ), .(scenario, n, parameter, cens.rate, beta.cens, gamma.cens, link, method)]



tabout <- dcast(type1[parameter == "ci"],
                cens.rate + beta.cens + method + se.method ~  link + n, value.var = "type1err")


tabout[, method := factor(method, levels = (c("jack", "infjack", "ipcwcoxph", "ipcwaalen")),
                        ordered = TRUE)]
setkey(tabout, cens.rate, beta.cens, method, se.method)

knitr::kable(tabout[, -c(1:2)],
      format = "latex", digits = 3)


tabout <- dcast(type1[parameter == "rmean"],
                cens.rate + beta.cens + method + se.method ~  link + n, value.var = "type1err")


tabout[, method := factor(method, levels = (c("jack", "infjack", "ipcwcoxph", "ipcwaalen")),
                          ordered = TRUE)]
setkey(tabout, cens.rate, beta.cens, method, se.method)

knitr::kable(tabout[, -c(1:2)],
             format = "latex", digits = 3)



library(data.table)
library(ggplot2)

sim2 <- data.table(read.csv("tests/sims-type-2.csv"))
truev <- data.table(read.csv("tests/true-coefficients.csv"))
truev <- rbind(truev, data.table(X = "tr3", scenario = 0, link = "identity",
                                 tv.cuminc = 0, tv.rmean = 0))


sim2 <- merge(sim2,
              rbind(truev[, .(scenario, link, trueval = tv.cuminc, parameter = "ci")],
                    truev[, .(scenario, link, trueval = tv.rmean, parameter = "rmean")]),
              by = c("parameter", "link", "scenario"))

sim2[, method := factor(method, levels = rev(c("jack", "strat", "ipcwcoxph", "ipcwaalen")),
                        ordered = TRUE)]


## bias

sim2[, bias := (estimate - trueval) / trueval]

res2 <- sim2[!is.na(estimate), .(mest = median(bias),
                                 lower = quantile(bias, .05) ,
                                 upper = quantile(bias, .95)),
             by = .(scenario, n, parameter, cens.rate, beta.cens, gamma.cens, link, method)]



sim2[, emp.sd := sd(estimate, na.rm = TRUE),
     by = .(parameter, method, n, cens.rate, scenario, link, beta.cens, gamma.cens)]

empsd2 <- sim2[, .(emp.sd = sd(estimate, na.rm = TRUE)),
     by = .(parameter, method, n, cens.rate, scenario, link, beta.cens, gamma.cens)]


empsd2$method <- factor(empsd2$method, levels = (c("jack", "strat", "ipcwcoxph", "ipcwaalen")),
                          ordered = TRUE)

#empsd2$beta.cens <- substr(empsd2$beta.cens, 3, 3)
empsd2$scen <- with(empsd2,paste0(scenario, ";", beta.cens, ";", cens.rate, ";", gamma.cens))

empsdres <- dcast(empsd2, link + parameter + scen  ~ n + method , value.var = "emp.sd")


knitr::kable(empsdres[, -c(1:2)], format = "latex", digits = 3)



## coverage
covertrue <- function(est, se, true) {

  (true >= est - 1.96 * se) &
    (true <= est + 1.96 * se)

}

emperr <- sim2[, .(estvar = c(median((std.err.rob - emp.sd) / emp.sd, na.rm = TRUE),
         median((std.err.cor - emp.sd) / emp.sd),
         median((std.err.nai - emp.sd) / emp.sd)),
         estvar.low = c(quantile((std.err.rob - emp.sd) / emp.sd, .05, na.rm = TRUE),
                        quantile((std.err.cor - emp.sd) / emp.sd, .05),
                        quantile((std.err.nai - emp.sd) / emp.sd, .05)),
         estvar.high = c(quantile((std.err.rob - emp.sd) / emp.sd, .95, na.rm = TRUE),
                         quantile((std.err.cor - emp.sd) / emp.sd, .95),
                         quantile((std.err.nai - emp.sd) / emp.sd, .95)),
         coverage = c(mean(covertrue(estimate, std.err.rob, trueval), na.rm = TRUE),
                      mean(covertrue(estimate, std.err.cor, trueval)),
                      mean(covertrue(estimate, std.err.nai, trueval))),
         varmethod = c("robust", "corrected", "naive")),
         by = .(parameter, method, n, cens.rate, scenario, link, beta.cens, gamma.cens)]



table(res2$beta.cens)

bias.emp.sd <- merge(
rbind(dcast(res2[scenario == 4 & cens.rate == .5 &
       n == 500 & gamma.cens == 1][beta.cens %in% c("0-0-0", "1-0-0", "1-0.5-0", "1-0.5-0.5")],
      parameter + beta.cens + gamma.cens ~ method, value.var = "mest"),
dcast(res2[scenario == 4 & cens.rate == .5 &
             n == 500 & gamma.cens == 3][beta.cens %in% c("1-0.5-0", "1-0.5-0.5", "0-0.5-0.5")],
      parameter + beta.cens + gamma.cens ~ method, value.var = "mest")),
rbind(dcast(empsd2[scenario == 4 & cens.rate == .5 &
                   n == 500 & gamma.cens == 1][beta.cens %in% c("0-0-0", "1-0-0", "1-0.5-0", "1-0.5-0.5")],
            parameter + beta.cens + gamma.cens ~ method, value.var = "emp.sd"),
      dcast(empsd2[scenario == 4 & cens.rate == .5 &
                   n == 500 & gamma.cens == 3][beta.cens %in% c("1-0.5-0", "1-0.5-0.5", "0-0.5-0.5")],
            parameter + beta.cens + gamma.cens ~ method, value.var = "emp.sd")),
by = c("parameter", "beta.cens", "gamma.cens"), sort = FALSE)


mary <- lapply(levels(res2$method), function(j) {

  sprintf("%.3f (%.3f)", bias.emp.sd[[paste0(j, ".x")]], bias.emp.sd[[paste0(j, ".y")]])

})
names(mary) <- levels(res2$method)
btab <- cbind(bias.emp.sd[, .(parameter, beta.cens, gamma.cens)], do.call(cbind, mary))

btab <- btab[order(btab$parameter),]
bcen.look <- c("indep", "strat", "ignorable", "missp1", "midssp2")
names(bcen.look) <- c("0-0-0", "1-0-0", "1-0.5-0", "1-0.5-0.5", "0-0.5-0.5")
btab$beta.cens <-bcen.look[btab$beta.cens]
btab$gamma.cens <- ifelse(btab$gamma.cens == 1, "PH", "non-PH")
print(xtable(btab), include.rownames = FALSE)

cover <- rbind(dcast(emperr[scenario == 4 & cens.rate == .5 &
                     n == 500 & gamma.cens == 1][beta.cens %in% c("0-0-0", "1-0-0", "1-0.5-0", "1-0.5-0.5")],
            parameter + beta.cens + gamma.cens + varmethod ~ method , value.var = "coverage"),
      dcast(emperr[scenario == 4 & cens.rate == .5 &
                     n == 500 & gamma.cens == 3][beta.cens %in% c("1-0.5-0", "1-0.5-0.5", "0-0.5-0.5")],
            parameter + beta.cens + gamma.cens + varmethod ~ method , value.var = "coverage"))
cover$beta.cens <- bcen.look[cover$beta.cens]
cover$gamma.cens <- ifelse(cover$gamma.cens == 1, "PH", "non-PH")
print(xtable(cover[parameter == "ci", -1]), include.rownames = FALSE)

bias.type1 <- rbind(dcast(res1[cens.rate == .5 &
             n == 500 & gamma.cens == 1][beta.cens %in% c("0-0-0", "1-0-0", "1-0.5-0", "1-0.5-0.5")],
      parameter + beta.cens + gamma.cens ~ method, value.var = "mest"),
      dcast(res1[cens.rate == .5 &
                   n == 500 & gamma.cens == 3][beta.cens %in% c("0-0-0", "1-0-0", "1-0.5-0", "1-0.5-0.5")],
            parameter + beta.cens + gamma.cens ~ method, value.var = "mest"))


t1err <- rbind(dcast(type1[cens.rate == .5 &
                            n == 500 & gamma.cens == 1][beta.cens %in% c("0-0-0", "1-0-0", "1-0.5-0", "1-0.5-0.5")],
                     parameter + beta.cens + gamma.cens + se.method ~ method, value.var = "type1err"),
               dcast(type1[cens.rate == .5 &
                            n == 500 & gamma.cens == 3][beta.cens %in% c("0-0-0", "1-0-0", "1-0.5-0", "1-0.5-0.5")],
                     parameter + beta.cens + gamma.cens + se.method ~ method, value.var = "type1err"))


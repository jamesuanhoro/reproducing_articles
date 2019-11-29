source("sls_mods.R")
library(dplyr)
library(sandwich)
library(lmtest)
library(stargazer)

dat <- read.csv("eclsk_subset.csv")

coef(summary(fit.logit <- glm(profread ~ ., binomial, dat)))
coeftest(fit.logit <- glm(profread ~ ., binomial, dat))
coeftest(fit.ols <- lm(profread ~ ., dat), vcov. = vcovHC, type = "HC4")
fit.sls <- s.ols(fit.ols, print = T)
coeftest(fit.sls$sls, vcov. = vcovHC, type = "HC4")
coeftest(fit.pois <- glm(profread ~ ., poisson, dat), vcov. = vcovHC, type = "HC4")
fit.spois <- s.pois(fit.pois, print = T)
coeftest(fit.spois$sls, vcov. = vcovHC, type = "HC4")

ecls.mods <- list(fit.logit, fit.ols, fit.sls$sls, fit.pois, fit.spois$sls)
ecls.mods.inf <- list(
  se = list(
    coeftest(fit.logit <- glm(profread ~ ., binomial, dat))[, 2],
    coeftest(fit.ols <- lm(profread ~ ., dat), vcov. = vcovHC, type = "HC4")[, 2],
    coeftest(fit.sls$sls, vcov. = vcovHC, type = "HC4")[, 2],
    coeftest(fit.pois <- glm(profread ~ ., poisson, dat), vcov. = vcovHC, type = "HC4")[, 2],
    coeftest(fit.spois$sls, vcov. = vcovHC, type = "HC4")[, 2]),
  p = list(
    coeftest(fit.logit <- glm(profread ~ ., binomial, dat))[, 4],
    coeftest(fit.ols <- lm(profread ~ ., dat), vcov. = vcovHC, type = "HC4")[, 4],
    coeftest(fit.sls$sls, vcov. = vcovHC, type = "HC4")[, 4],
    coeftest(fit.pois <- glm(profread ~ ., poisson, dat), vcov. = vcovHC, type = "HC4")[, 4],
    coeftest(fit.spois$sls, vcov. = vcovHC, type = "HC4")[, 4])
)

# Table 4 in paper
stargazer(
  ecls.mods,
  digits = 2, title = "Analysis of dichotomized reading proficiency", header = FALSE,
  dep.var.caption = "Outcome variable: Reading proficiency", dep.var.labels.include = FALSE,
  model.names = FALSE, model.numbers = FALSE, report = "vcs*",
  se = ecls.mods.inf$se, p = ecls.mods.inf$p,
  column.labels = c("Logistic", "OLS", "SLS", "Poisson", "Seq-Poisson"),
  no.space = TRUE, column.sep.width = "2.5pt",
  omit.stat = c("ll", "aic", "adj.rsq", "rsq", "f", "ser"),
  star.cutoffs = c(.05, .01, .001), out = "ecls_reg_res_0.tex", font.size = "small",
  notes = c("Standard errors for the linear and Poisson models are HC4",
            "robust standard errors."), notes.align = "l", type = "latex")

library(ggplot2)

med.dat <- as.data.frame(t(apply(dat[, -1], 2, median)))
med.dat <- med.dat[rep(1, nrow(dat)), ]
med.dat$wksesl <- dat$wksesl # ; med.dat$p1ageent_66 <- mean(dat$p1ageent_66)
med.dat <- cbind(1, as.matrix(med.dat))
dat.new <- data.frame(
  SES = dat$wksesl, logit0 = plogis(med.dat %*% coef(fit.logit)),
  logit1 = plogis(med.dat %*% coef(fit.logit)),
  ols = punif(med.dat %*% coef(fit.ols)), sls = punif(med.dat %*% coef(fit.sls$sls)),
  pois = punif(exp(med.dat %*% coef(fit.pois))), spois = punif(exp(med.dat %*% coef(fit.spois$sls))))
dat.new <- tidyr::gather(dat.new, method, fitted, logit0:spois)
dat.new$RD <- ifelse(dat.new$method %in% c("logit0", "ols", "sls"), "Risk difference", "Risk ratio")
dat.new$approach <- ifelse(grepl("logit", dat.new$method), "Logit", ifelse(
  dat.new$method %in% c("ols", "pois"), "Standard", "Sequential"))
dat.new$approach <- factor(
  dat.new$approach, levels = c("Logit", "Standard", "Sequential"), ordered = TRUE)
saveRDS(dat.new, "models_ecls_plot.RDS")

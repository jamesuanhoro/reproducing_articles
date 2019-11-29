library(SimDesign)
library(sandwich)
library(lmtest)

gen_unif_bern <- function(n, a = -1, b = 1, p = .5, r = .3) {
  r <- r / (dnorm(qnorm(p)) / sqrt(p * (1 - p)))
  U <- pnorm(t(chol(matrix(c(1, r, r, 1), 2))) %*%
               matrix(rnorm(n * 2), nrow = 2))
  cbind(qunif(U[1, ], a, b), qbinom(U[2, ], 1, p))
}

s.ols <- function(fit.ols) {
  count <- 1
  dat.ols <- model.frame(fit.ols)
  n.org <- nrow(dat.ols)
  fitted <- fit.ols$fitted.values
  while (any(fitted > 1 | fitted < 0)) {
    count <- count + 1
    dat.ols <- dat.ols[!(fitted > 1 | fitted < 0), ]
    m.ols <- lm(paste(names(dat.ols)[1], "~ ."), dat.ols)
    fitted <- m.ols$fitted.values
  }
  # dat.ols$var_inv <- as.numeric(1 / (fitted * (1 - fitted)))
  m.ols <- lm(paste(names(dat.ols)[1], "~ ."), dat.ols)
  # m.ols <- lm(y ~ xb + xc, dat.ols, weights = var_inv)
  print(paste("Number of iterations:", count))
  return(list(sls = m.ols, count = count,
              n.org = n.org, n.fin = nrow(dat.ols)))
}

# condition = list(n = 500, eff.b = 0, eff.c = .5)

design_7 <- expand.grid(n = 500, eff.b = c(0, .1), eff.c = c(0, .25, .5, .75))
saveRDS(design_7, "design_7.rds")

# 0-0 should be fine through till
# 0-.5 should be borderline but fine
# x-.5 should be problematic if |x| > 0

generate_7 <- function(condition, fixed_objects = NULL) {
  dat <- with(condition, {
    X <- gen_unif_bern(n)
    data.frame(
      xc = X[, 1], xb = X[, 2],
      y = ((0.5 + eff.b * X[, 2] + eff.c * X[, 1] + runif(n, -.5, .5)) > 0.5) + 0)
  })
}

analyze_7 <- function(condition, dat, fixed_objects = NULL) {
  ols.0 <- lm(y ~ xb + xc, dat)
  # logit.0 <- glm(y ~ xb + xc, binomial, dat)
  ols.0.c <- unname(coef(ols.0))
  ols.0.s <- coeftest(ols.0, vcov. = vcovHC, type = "HC4")
  ols.0.ci <- coefci(ols.0, vcov. = vcovHC, type = "HC4")
  sigma.0 <- sigma(ols.0)

  ols.0s <- s.ols(ols.0)$sls
  ols.0s.c <- unname(coef(ols.0s))
  ols.0s.s <- coeftest(ols.0s, vcov. = vcovHC, type = "HC4")
  ols.0s.ci <- coefci(ols.0s, vcov. = vcovHC, type = "HC4")
  sigma.0s <- sigma(ols.0s)
  # as.vector(tcrossprod(solve(crossprod(model.matrix(ols.0))),
  #                      model.matrix(ols.0)) %*% dat$y)

  # Misspecified sigma affects statistical inference
  # sigma.0 * sqrt(diag(solve(crossprod(model.matrix(ols.0)))))

  fitted.0 <- unname(fitted(ols.0))
  perc.exceed <- 1 - mean(fitted.0 < 1 & fitted.0 > 0)

  ecr.0s.i <- ECR(ols.0s.ci[1, ], .5)
  ecr.0s.b <- ECR(ols.0s.ci[2, ], condition$eff.b)
  ecr.0s.c <- ECR(ols.0s.ci[3, ], condition$eff.c)
  edr.0s.i <- EDR(ols.0s.s[1, 4])
  edr.0s.b <- EDR(ols.0s.s[2, 4])
  edr.0s.c <- EDR(ols.0s.s[3, 4])
  bias.0s.i <- bias(ols.0s.c[1], .5)
  bias.0s.b <- bias(ols.0s.c[2], condition$eff.b)
  bias.0s.c <- bias(ols.0s.c[3], condition$eff.c)

  ecr.0.i <- ECR(ols.0.ci[1, ], .5)
  ecr.0.b <- ECR(ols.0.ci[2, ], condition$eff.b)
  ecr.0.c <- ECR(ols.0.ci[3, ], condition$eff.c)

  edr.0.i <- EDR(ols.0.s[1, 4])
  edr.0.b <- EDR(ols.0.s[2, 4])
  edr.0.c <- EDR(ols.0.s[3, 4])

  bias.0.i <- bias(ols.0.c[1], .5)
  bias.0.b <- bias(ols.0.c[2], condition$eff.b)
  bias.0.c <- bias(ols.0.c[3], condition$eff.c)

  stat <- c(
    exceed = perc.exceed,
    ecr.0.i = ecr.0.i, ecr.0.b = ecr.0.b, ecr.0.c = ecr.0.c,
    ecr.0s.i = ecr.0s.i, ecr.0s.b = ecr.0s.b, ecr.0s.c = ecr.0s.c,
    edr.0.i = edr.0.i, edr.0.b = edr.0.b, edr.0.c = edr.0.c,
    edr.0s.i = edr.0s.i, edr.0s.b = edr.0s.b, edr.0s.c = edr.0s.c,
    bias.0.i = bias.0.i, bias.0.b = bias.0.b, bias.0.c = bias.0.c,
    bias.0s.i = bias.0s.i, bias.0s.b = bias.0s.b, bias.0s.c = bias.0s.c,
    sigma.0 = sigma.0, sigma.0s = sigma.0s
  )
}

summaries_7 <- function(condition, results, fixed_objects = NULL) {
  ret <- colMeans(results)
  ret
}

results_7 <- runSimulation(
  design = design_7, replications = 4999, parallel = TRUE,
  generate = generate_7, analyse = analyze_7, summarise = summaries_7, edit = "none",
  save = TRUE, save_results = TRUE, save_generate_data = FALSE, progress = TRUE,
  filename = "./simdata_7", packages = c("sandwich", "lmtest"))

library(dplyr)

(results_7 <- readRDS("simdata_7.rds"))

{
  ecr_7 <- select(results_7, contains("eff"), exceed, contains("ecr"))
  ecr_7 <- tidyr::gather(ecr_7, variable, ecr, contains("ecr"))
  ecr_7$model <- (unlist(lapply(strsplit(
    ecr_7$variable, ".", TRUE), "[[", 2)))
  ecr_7$variable <- unlist(lapply(strsplit(
    ecr_7$variable, ".", TRUE), "[[", 3))
  ecr_7$variable <- ifelse(ecr_7$variable == "i", "intercept", ifelse(
    ecr_7$variable == "b", "binary", "continuous"))
  ecr_7 <- filter(ecr_7, variable != "intercept")
  ecr_7$exceed <- round(ecr_7$exceed * 100, 2)
  names(ecr_7)[1:2] <- paste0("coef.", c("bin", "contin"))
  ecr_7$c75 <- "Does not greatly exceed unit interval"
  ecr_7$c75[ecr_7$coef.contin == .75] <- "Greatly exceeds unit interval"
  ecr_7$exp <- .95
  ecr_7$exp[ecr_7$coef.contin == .75 & ecr_7$model == 0] <- NA
}
ecr_7

{
  bias_7 <- select(results_7, contains("eff"), exceed, contains("bias"))
  bias_7 <- tidyr::gather(bias_7, variable, bias, contains("bias"))
  bias_7$model <- unlist(lapply(strsplit(
    bias_7$variable, ".", TRUE), "[[", 2))
  bias_7$model[bias_7$model == 0] <- "OLS"
  bias_7$model[bias_7$model == "0s"] <- "Seq. OLS"
  bias_7$variable <- unlist(lapply(strsplit(
    bias_7$variable, ".", TRUE), "[[", 3))
  bias_7$variable <- ifelse(bias_7$variable == "i", "intercept", ifelse(
    bias_7$variable == "b", "binary", "continuous"))
  bias_7 <- filter(bias_7, variable != "intercept")
  bias_7$exp <- 0
  bias_7$exceed <- round(bias_7$exceed * 100, 2)
  names(bias_7)[1:2] <- paste0("coef.", c("bin", "contin"))
}
bias_7

saveRDS(list(bias_7 = bias_7, ecr_7 = ecr_7), "res_7.RDS")

{
  edr_7 <- select(results_7, contains("eff"), exceed, contains("edr.0"))
  edr_7 <- tidyr::gather(edr_7, variable, edr, contains("edr.0"))
  edr_7$model <- unlist(lapply(strsplit(edr_7$variable, ".", TRUE), "[[", 2))
  edr_7$model[edr_7$model == 0] <- "OLS"
  edr_7$model[edr_7$model == "0s"] <- "Seq. OLS"
  edr_7$variable <- gsub(".", "", gsub("edr.0.", "", edr_7$variable), fixed = TRUE)
  edr_7$variable <- ifelse(edr_7$variable == "i", "intercept", ifelse(
    edr_7$variable == "b", "binary", "continuous"))
  edr_7 <- filter(edr_7, variable != "intercept")
  edr_7$exp <- NA
  edr_7$exp[edr_7$variable == "binary" & edr_7$eff.b == 0] <- .05
  edr_7$exp[edr_7$variable == "continuous" & edr_7$eff.c == 0] <- .05
  edr_7$hyp <- ifelse(is.na(edr_7$exp), "Alt", "Null")
  edr_7$hyp <- paste("Hypothesis:", edr_7$hyp)
  edr_7$hyp <- factor(edr_7$hyp, levels = paste("Hypothesis:", c("Null", "Alt")))
  edr_7$perf <- ifelse(edr_7$eff.c == .75, "poor", ifelse(
    edr_7$eff.c == .5 & edr_7$eff.b == .1, "poor", ifelse(
      edr_7$eff.c == .5 & edr_7$eff.b == 0, "borderline", "good"
    )))
  edr_7$perf <- factor(edr_7$perf, levels = c("good", "borderline", "poor"))
  names(edr_7)[1:2] <- paste0("coef.", c("bin", "contin"))
}
edr_7

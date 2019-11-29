library(SimDesign)
library(sandwich)
library(lmtest)

gen_unif_bern <- function(n, a = -1, b = 1, p = .5, r = .3) {
  r <- r / (dnorm(qnorm(p)) / sqrt(p * (1 - p)))
  U <- pnorm(t(chol(matrix(c(1, r, r, 1), 2))) %*%
               matrix(rnorm(n * 2), nrow = 2))
  cbind(qunif(U[1, ], a, b), qbinom(U[2, ], 1, p))
}

s.pois <- function(fit.pois) {
  count <- 1
  dat.pois <- model.frame(fit.pois)
  n.org <- nrow(dat.pois)
  fitted <- fit.pois$fitted.values
  while (any(fitted > 1)) {
    count <- count + 1
    dat.pois <- dat.pois[!(fitted > 1 | fitted < 0), ]
    m.pois <- glm(paste(names(dat.pois)[1], "~ ."), poisson, dat.pois)
    fitted <- m.pois$fitted.values
  }
  m.pois <- glm(paste(names(dat.pois)[1], "~ ."), poisson, dat.pois)
  print(paste("Number of iterations:", count))
  return(list(sls = m.pois, count = count,
              n.org = n.org, n.fin = nrow(dat.pois)))
}

# condition <- list(n = 5000, eff.b = .5, eff.c = .25)

design_8 <- expand.grid(n = 500, eff.b = c(0, .5), eff.c = c(0, .5, 1))
saveRDS(design_8, "design_8.rds")

# 0-0 should be fine through till
# 0-.5 should be borderline but fine
# x-.5 should be problematic if |x| > 0

generate_8 <- function(condition, fixed_objects = NULL) {
  dat <- with(condition, {
    X <- gen_unif_bern(n)
    eps <- rexp(n, 1) - 1
    data.frame(
      xc = X[, 1], xb = X[, 2],
      y = ((-1 + eff.b * X[, 2] + eff.c * X[, 1] + eps) > -1) + 0)
  })
}

analyze_8 <- function(condition, dat, fixed_objects = NULL) {
  pois.0 <- glm(y ~ xb + xc, poisson, dat)
  # logit.0 <- glm(y ~ xb + xc, binomial, dat)
  pois.0.c <- unname(coef(pois.0))
  pois.0.s <- coeftest(pois.0, vcov. = vcovHC, type = "HC4")
  pois.0.ci <- coefci(pois.0, vcov. = vcovHC, type = "HC4")

  pois.0s <- s.pois(pois.0)$sls
  pois.0s.c <- unname(coef(pois.0s))
  pois.0s.s <- coeftest(pois.0s, vcov. = vcovHC, type = "HC4")
  pois.0s.ci <- coefci(pois.0s, vcov. = vcovHC, type = "HC4")
  # as.vector(tcrossprod(solve(crossprod(model.matrix(pois.0))),
  #                      model.matrix(pois.0)) %*% dat$y)

  # Misspecified sigma affects statistical inference
  # sigma.0 * sqrt(diag(solve(crossprod(model.matrix(pois.0)))))

  fitted.0 <- unname(fitted(pois.0))
  perc.exceed <- 1 - mean(fitted.0 < 1 & fitted.0 > 0)

  ecr.0s.i <- ECR(pois.0s.ci[1, ], -1)
  ecr.0s.b <- ECR(pois.0s.ci[2, ], condition$eff.b)
  ecr.0s.c <- ECR(pois.0s.ci[3, ], condition$eff.c)
  edr.0s.i <- EDR(pois.0s.s[1, 4])
  edr.0s.b <- EDR(pois.0s.s[2, 4])
  edr.0s.c <- EDR(pois.0s.s[3, 4])
  bias.0s.i <- bias(pois.0s.c[1], -1)
  bias.0s.b <- bias(pois.0s.c[2], condition$eff.b)
  bias.0s.c <- bias(pois.0s.c[3], condition$eff.c)

  ecr.0.i <- ECR(pois.0.ci[1, ], -1)
  ecr.0.b <- ECR(pois.0.ci[2, ], condition$eff.b)
  ecr.0.c <- ECR(pois.0.ci[3, ], condition$eff.c)

  edr.0.i <- EDR(pois.0.s[1, 4])
  edr.0.b <- EDR(pois.0.s[2, 4])
  edr.0.c <- EDR(pois.0.s[3, 4])

  bias.0.i <- bias(pois.0.c[1], -1)
  bias.0.b <- bias(pois.0.c[2], condition$eff.b)
  bias.0.c <- bias(pois.0.c[3], condition$eff.c)

  stat <- c(
    exceed = perc.exceed,
    ecr.0.i = ecr.0.i, ecr.0.b = ecr.0.b, ecr.0.c = ecr.0.c,
    ecr.0s.i = ecr.0s.i, ecr.0s.b = ecr.0s.b, ecr.0s.c = ecr.0s.c,
    edr.0.i = edr.0.i, edr.0.b = edr.0.b, edr.0.c = edr.0.c,
    edr.0s.i = edr.0s.i, edr.0s.b = edr.0s.b, edr.0s.c = edr.0s.c,
    bias.0.i = bias.0.i, bias.0.b = bias.0.b, bias.0.c = bias.0.c,
    bias.0s.i = bias.0s.i, bias.0s.b = bias.0s.b, bias.0s.c = bias.0s.c
  )
}

summaries_8 <- function(condition, results, fixed_objects = NULL) {
  ret <- colMeans(results)
  ret
}

results_8 <- runSimulation(
  design = design_8, replications = 4999, parallel = TRUE,
  generate = generate_8, analyse = analyze_8, summarise = summaries_8, edit = "none",
  save = TRUE, save_results = TRUE, save_generate_data = FALSE, progress = TRUE,
  filename = "./simdata_8", packages = c("sandwich", "lmtest"))

library(dplyr)

(results_8 <- readRDS("simdata_8.rds"))

{
  ecr_8 <- select(results_8, contains("eff"), exceed, contains("ecr"))
  ecr_8 <- tidyr::gather(ecr_8, variable, ecr, contains("ecr"))
  ecr_8$model <- (unlist(lapply(strsplit(
    ecr_8$variable, ".", TRUE), "[[", 2)))
  ecr_8$model[ecr_8$model == 0] <- "Poisson"
  ecr_8$model[ecr_8$model == "0s"] <- "Seq. Poisson"
  ecr_8$variable <- unlist(lapply(strsplit(
    ecr_8$variable, ".", TRUE), "[[", 3))
  ecr_8$variable <- ifelse(ecr_8$variable == "i", "intercept", ifelse(
    ecr_8$variable == "b", "binary", "continuous"))
  ecr_8 <- filter(ecr_8, variable != "intercept")
  ecr_8$exceed <- round(ecr_8$exceed * 100, 2)
  names(ecr_8)[1:2] <- paste0("coef.", c("bin", "contin"))
  ecr_8$c75 <- "Does not greatly exceed unit interval"
  ecr_8$c75[ecr_8$coef.contin == .75] <- "Greatly exceeds unit interval"
  ecr_8$exp <- .95
  ecr_8$exp[ecr_8$coef.contin == .75 & ecr_8$model == 0] <- NA
}
ecr_8

{
  bias_8 <- select(results_8, contains("eff"), exceed, contains("bias"))
  bias_8 <- tidyr::gather(bias_8, variable, bias, contains("bias"))
  bias_8$model <- unlist(lapply(strsplit(
    bias_8$variable, ".", TRUE), "[[", 2))
  bias_8$model[bias_8$model == 0] <- "Poisson"
  bias_8$model[bias_8$model == "0s"] <- "Seq. Poisson"
  bias_8$variable <- unlist(lapply(strsplit(
    bias_8$variable, ".", TRUE), "[[", 3))
  bias_8$variable <- ifelse(bias_8$variable == "i", "intercept", ifelse(
    bias_8$variable == "b", "binary", "continuous"))
  bias_8 <- filter(bias_8, variable != "intercept")
  bias_8$exp <- 0
  bias_8$exceed <- round(bias_8$exceed * 100, 2)
  names(bias_8)[1:2] <- paste0("coef.", c("bin", "contin"))
}
bias_8

saveRDS(list(bias_8 = bias_8, ecr_8 = ecr_8), "res_8.RDS")

{
  edr_8 <- select(results_8, contains("eff"), exceed, contains("edr.0"))
  edr_8 <- tidyr::gather(edr_8, variable, edr, contains("edr.0"))
  edr_8$model <- unlist(lapply(strsplit(edr_8$variable, ".", TRUE), "[[", 2))
  edr_8$model[edr_8$model == 0] <- "Poisson"
  edr_8$model[edr_8$model == "0s"] <- "Seq. Poisson"
  edr_8$variable <- gsub(".", "", gsub("edr.0.", "", edr_8$variable), fixed = TRUE)
  edr_8$variable <- ifelse(edr_8$variable == "i", "intercept", ifelse(
    edr_8$variable == "b", "binary", "continuous"))
  edr_8 <- filter(edr_8, variable != "intercept")
  edr_8$exp <- NA
  edr_8$exp[edr_8$variable == "binary" & edr_8$eff.b == 0] <- .05
  edr_8$exp[edr_8$variable == "continuous" & edr_8$eff.c == 0] <- .05
  edr_8$hyp <- ifelse(is.na(edr_8$exp), "Alt", "Null")
  edr_8$hyp <- paste("Hypothesis:", edr_8$hyp)
  edr_8$hyp <- factor(edr_8$hyp, levels = paste("Hypothesis:", c("Null", "Alt")))
  edr_8$perf <- ifelse(edr_8$eff.c == .75, "poor", ifelse(
    edr_8$eff.c == .5 & edr_8$eff.b == .1, "poor", ifelse(
      edr_8$eff.c == .5 & edr_8$eff.b == 0, "borderline", "good"
    )))
  edr_8$perf <- factor(edr_8$perf, levels = c("good", "borderline", "poor"))
  names(edr_8)[1:2] <- paste0("coef.", c("bin", "contin"))
}
ecr_8

library(sandwich); library(lmtest)

s.pois <- function(fit.pois, print = FALSE) {
  count <- 1
  dat.pois <- model.frame(fit.pois)
  n.org <- nrow(dat.pois)
  fitted <- fit.pois$fitted.values
  form <- formula(fit.pois)
  while (any(fitted > 1)) {
    count <- count + 1
    dat.pois <- dat.pois[!(fitted > 1 | fitted < 0), ]
    m.pois <- glm(form, poisson, dat.pois)
    fitted <- m.pois$fitted.values
  }
  m.pois <- glm(form, poisson, dat.pois)
  # m.pois$fitted.values <- exp(as.numeric(model.matrix(fit.pois) %*% coef(m.pois)))
  # m.pois$fitted.values[m.pois$fitted.values > 1] <- 1
  deleted <- n.org - nrow(dat.pois)
  if (print) print(paste("Number of deleted cases:", deleted))
  return(list(sls = m.pois, count = count, dat = dat.pois,
              n.del = deleted, n.fin = nrow(dat.pois)))
}

s.pois.1 <- function(fit.pois, print = FALSE) {
  count <- 1
  dat.pois <- model.matrix(fit.pois)[, -1]
  y <- fit.pois$model[, 1]
  n.org <- nrow(dat.pois)
  fitted <- fit.pois$fitted.values
  while (any(fitted > 1)) {
    count <- count + 1
    dat.pois <- dat.pois[!(fitted > 1 | fitted < 0), ]
    y <- y[!(fitted > 1 | fitted < 0)]
    m.pois <- glm(y ~ dat.pois, poisson)
    fitted <- m.pois$fitted.values
  }
  m.pois <- glm(y ~ dat.pois, poisson)
  deleted <- n.org - nrow(dat.pois)
  if (print) print(paste("Number of deleted cases:", deleted))
  return(list(sls = m.pois, count = count, dat = dat.pois,
              n.del = deleted, n.fin = nrow(dat.pois), y = y))
}

s.ols <- function(fit.ols, print = FALSE) {
  count <- 1
  dat.ols <- model.frame(fit.ols)
  n.org <- nrow(dat.ols)
  fitted <- fit.ols$fitted.values
  form <- formula(fit.ols)
  while (any(fitted > 1 | fitted < 0)) {
    count <- count + 1
    dat.ols <- dat.ols[!(fitted > 1 | fitted < 0), ]
    m.ols <- lm(form, dat.ols)
    fitted <- m.ols$fitted.values
  }
  # dat.ols$var_inv <- as.numeric(1 / (fitted * (1 - fitted)))
  # m.wls <- lm(form, dat.ols, weights = var_inv)
  # m.wls$fitted.values <- punif(as.numeric(model.matrix(fit.ols) %*% coef(m.wls)))
  m.ols <- lm(form, dat.ols)
  # m.ols$fitted.values <- punif(as.numeric(model.matrix(fit.ols) %*% coef(m.ols)))
  deleted <- n.org - nrow(dat.ols)
  if (print) print(paste("Number of deleted cases:", deleted))
  return(list(sls = m.ols, count = count, dat = dat.ols,
              n.del = deleted, n.fin = nrow(dat.ols)))
}

s.ols.1 <- function(fit.ols, print = FALSE) {
  count <- 1
  dat.ols <- model.matrix(fit.ols)[, -1]
  y <- fit.ols$model[, 1]
  n.org <- nrow(dat.ols)
  fitted <- fit.ols$fitted.values
  while (any(fitted > 1 | fitted < 0)) {
    count <- count + 1
    dat.ols <- dat.ols[!(fitted > 1 | fitted < 0), ]
    y <- y[!(fitted > 1 | fitted < 0)]
    m.ols <- lm(y ~ dat.ols)
    fitted <- m.ols$fitted.values
  }
  # m.ols <- lm(y ~ dat.ols)
  deleted <- n.org - nrow(dat.ols)
  if (print) print(paste("Number of deleted cases:", deleted))
  return(list(sls = m.ols, count = count, dat = dat.ols,
              n.del = deleted, n.fin = nrow(dat.ols), y = y))
}

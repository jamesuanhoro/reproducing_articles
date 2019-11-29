library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)

# Figure 1

cunif <- function (xs) sapply(xs, function (x) ifelse(x < 0, 0, ifelse(x > 1, 1, x)))
ggplot(data.frame(x = seq(-.5, 1.5, .0001)), aes(x, cunif(x))) +
  geom_line() + theme_bw() + scale_y_continuous(labels = scales::percent) +
  labs(x = "Fitted values from regression (Xbeta)", y = "P(y = 1 | X)")

# Figure 2

cpwdisc <- function (xs) sapply(xs, function (x) ifelse(x > 0, 0, x))
cexp <- function (xs) sapply(xs, function (x) ifelse(x > 0, 1, exp(x)))
dat <- data.frame(x = seq(-6, 1.5, .01))
dat$logp <- cpwdisc(dat$x); dat$p <- cexp(dat$x)
dat <- tidyr::gather(dat, outcome, val, logp, p)
dat$outcome <- ifelse(dat$outcome == "p", "P(y = 1 | X)", "log(P(y = 1 | X))")
ggplot(dat, aes(x, val)) + geom_line() + theme_bw() + scale_y_continuous() +
  facet_wrap(~ outcome, scales = "free_y", strip.position = "left") +
  labs(x = "Fitted values from regression", y = "")

# Retrieve simulation 1 and 2 results and format data for plotting

## Simulation 1
des_7 <- readRDS("design_7.rds") %>%
  rename(coef.bin = eff.b, coef.contin = eff.c) %>%
  mutate(Condition = 1:nrow(.))
res_2 <- readRDS("res_7.RDS")
bias_2 <- merge(res_2$bias_7, des_7)
bias_2$exceed.f <- ifelse(
  bias_2$exceed == 0, "a", ifelse(bias_2$exceed > 20, "c", "b"))
bias_2$exceed.f <- factor(
  bias_2$exceed.f, letters[1:3], labels = c("None", "2% - 8%", "> 20%"),
  ordered = TRUE)
bias_2$model[grepl("Seq", bias_2$model)] <- "SLS"
ecr_2 <- merge(res_2$ecr_7, des_7)
ecr_2$exceed.f <- ifelse(
  ecr_2$exceed == 0, "a", ifelse(ecr_2$exceed > 20, "c", "b"))
ecr_2$exceed.f <- factor(
  ecr_2$exceed.f, letters[1:3], labels = c("None", "2% - 8%", "> 20%"),
  ordered = TRUE)
ecr_2$model <- ifelse(grepl("0s", ecr_2$model), "SLS", "OLS")

## Simulation 2
des_8 <- readRDS("design_8.rds") %>%
  rename(coef.bin = eff.b, coef.contin = eff.c) %>%
  mutate(Condition = 1:nrow(.))
res_3 <- readRDS("res_8.RDS")
bias_3 <- merge(res_3$bias_8, des_8)
bias_3$exceed.f <- ifelse(
  bias_3$exceed == 0, "a", ifelse(bias_3$exceed > 8, "c", "b"))
bias_3$exceed.f <- factor(
  bias_3$exceed.f, letters[1:3], labels = c("None", "1% - 2%%", "11.5%"),
  ordered = TRUE)
ecr_3 <- merge(res_3$ecr_8, des_8)
ecr_3$exceed.f <- ifelse(
  ecr_3$exceed == 0, "a", ifelse(ecr_3$exceed > 8, "c", "b"))
ecr_3$exceed.f <- factor(
  ecr_3$exceed.f, letters[1:3], labels = c("None", "1% - 2%%", "11.5%"),
  ordered = TRUE)

## Figure 3

bias_2$exp.1 <- -.05
bias_2$exp.2 <- .05
bias_2_sub <- bias_2 %>% filter(exceed != 0)
bias_2_sub$rel_bias <- ifelse(
  grepl("bin", bias_2_sub$variable), bias_2_sub$bias / as.numeric(as.character(bias_2_sub$coef.bin)),
  bias_2_sub$bias / as.numeric(as.character(bias_2_sub$coef.contin)))
bias_2_sub <- bias_2_sub[is.finite(bias_2_sub$rel_bias), ]
lay <- rbind(c(1, 1, 1, 2, 2, 2, 2), c(3, 3, 3, 4, 4, 4, 4))
p1 <- ggplot(filter(bias_2_sub, variable == "binary"),
             aes(factor(Condition), rel_bias, fill = model)) +
  geom_point(aes(size = exceed.f, shape = model), position = position_dodge(.75)) +
  geom_point(size = .01, position = position_dodge(.75)) +
  facet_wrap(~ paste("Variable:", variable), ncol = 2, scales = "free_x") +
  theme_classic() + geom_hline(aes(yintercept = exp), linetype = 3) +
  geom_hline(aes(yintercept = exp.1), linetype = 2) +
  geom_hline(aes(yintercept = exp.2), linetype = 2) +
  scale_shape_manual(values = c(1, 3)) + coord_cartesian(ylim = c(-.35, .05)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, seq(-.55, .5, .1))) +
  labs(y = "Relative bias (%)", fill = "Model", shape = "Model",
       size = "% fitted values outside unit interval", x = "", tag = "A") +
  theme(legend.position = "top", legend.spacing.x = ) +
  guides(fill = FALSE, size = FALSE)
p2 <- ggplot(filter(bias_2_sub, variable == "continuous"),
             aes(factor(Condition), rel_bias, fill = model)) +
  geom_point(aes(size = exceed.f, shape = model), position = position_dodge(.75)) +
  geom_point(size = .01, position = position_dodge(.75)) +
  facet_wrap(~ paste("Variable:", variable), ncol = 2, scales = "free_x") +
  theme_classic() + geom_hline(aes(yintercept = exp), linetype = 3) +
  geom_hline(aes(yintercept = exp.1), linetype = 2) +
  geom_hline(aes(yintercept = exp.2), linetype = 2) +
  scale_shape_manual(values = c(1, 3)) + coord_cartesian(ylim = c(-.35, .05)) +
  scale_y_continuous(labels = scales::percent, breaks = c(),
                     sec.axis = sec_axis(~ ., labels = scales::percent, breaks = c(0, seq(-.55, .5, .1)))) +
  labs(y = "", fill = "Model", shape = "Model",
       size = "% fitted outside 0 - 1",
       x = "Condition number", tag = "") +
  theme(legend.position = "top", axis.title.x = element_text(hjust = -.2)) +
  guides(fill = FALSE, shape = FALSE)
ecr_2$exp.1 <- .925
ecr_2$exp.2 <- .975
ecr_2_sub <- filter(ecr_2, exceed != 0)
p3 <- ggplot(filter(ecr_2_sub, variable == "binary"),
             aes(factor(Condition), ecr, fill = model)) +
  geom_point(aes(size = exceed.f, shape = model), position = position_dodge(.75)) +
  geom_point(size = .01, position = position_dodge(.75)) +
  facet_wrap(~ paste("Variable:", variable), ncol = 2, scales = "free_x") +
  theme_classic() + geom_hline(aes(yintercept = exp), linetype = 3) +
  geom_hline(aes(yintercept = exp.1), linetype = 2) +
  geom_hline(aes(yintercept = exp.2), linetype = 2) +
  scale_shape_manual(values = c(1, 3)) +
  scale_y_continuous(labels = percent, breaks = c(.95, seq(.025, .975, .05))) +
  labs(y = "Empirical coverage rate (95% CI)", fill = "Model", shape = "Model",
       size = "% fitted outside 0 - 1", x = "", tag = "B") +
  coord_cartesian(ylim = c(.80, .975)) +
  theme(legend.position = "bottom") +
  guides(fill = FALSE, size = FALSE)
p4 <- ggplot(filter(ecr_2_sub, variable == "continuous"),
             aes(factor(Condition), ecr, fill = model)) +
  geom_point(aes(size = exceed.f, shape = model), position = position_dodge(.75)) +
  geom_point(size = .01, position = position_dodge(.75)) +
  facet_wrap(~ paste("Variable:", variable), ncol = 2, scales = "free_x") +
  theme_classic() + geom_hline(aes(yintercept = exp), linetype = 3) +
  geom_hline(aes(yintercept = exp.1), linetype = 2) +
  geom_hline(aes(yintercept = exp.2), linetype = 2) +
  scale_shape_manual(values = c(1, 3)) +
  scale_y_continuous(labels = percent, breaks = c(), sec.axis = sec_axis(
    ~ ., labels = percent, breaks = c(.95, seq(.025, .975, .05)))) +
  coord_cartesian(ylim = c(.80, .975)) +
  labs(y = "", fill = "Model", shape = "Model", tag = "",
       size = "% fitted outside 0 - 1", x = "Condition number") +
  theme(legend.position = "bottom", axis.title.x = element_text(hjust = -.2)) +
  guides(fill = FALSE, shape = FALSE)
grid.arrange(p1, p2, p3, p4, layout_matrix = lay)

## Figure 4

bias_3$exp.1 <- -.05
bias_3$exp.2 <- .05
bias_3_sub <- bias_3 %>% filter(exceed != 0)
bias_3_sub$rel_bias <- ifelse(
  grepl("bin", bias_3_sub$variable), bias_3_sub$bias / as.numeric(as.character(bias_3_sub$coef.bin)),
  bias_3_sub$bias / as.numeric(as.character(bias_3_sub$coef.contin)))
bias_3_sub <- bias_3_sub[is.finite(bias_3_sub$rel_bias), ]
lay <- rbind(c(1, 1, 1, 1, 2, 2, 2, 2, 2), c(3, 3, 3, 3, 4, 4, 4, 4, 4))
p5 <- ggplot(filter(bias_3_sub, variable == "binary"),
             aes(factor(Condition), rel_bias, fill = model)) +
  geom_point(aes(size = exceed.f, shape = model), position = position_dodge(.75)) +
  geom_point(size = .01, position = position_dodge(.75)) +
  facet_wrap(~ paste("Variable:", variable), ncol = 2, scales = "free_x") +
  theme_classic() + geom_hline(aes(yintercept = exp), linetype = 3) +
  geom_hline(aes(yintercept = exp.1), linetype = 2) +
  geom_hline(aes(yintercept = exp.2), linetype = 2) +
  coord_cartesian(ylim = c(-.2, .05)) +
  scale_shape_manual(values = c(1, 3)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, seq(-.55, .05, .05))) +
  labs(y = "Relative bias (%)", fill = "Model", shape = "Model",
       size = "% fitted outside 0 - 1", x = "", tag = "A") +
  theme(legend.position = "top") +
  guides(fill = FALSE, size = FALSE, shape = guide_legend(keywidth = .05, default.unit = "inch"))
p6 <- ggplot(filter(bias_3_sub, variable == "continuous"),
             aes(factor(Condition), rel_bias, fill = model)) +
  geom_point(aes(size = exceed.f, shape = model), position = position_dodge(.75)) +
  geom_point(size = .01, position = position_dodge(.75)) +
  facet_wrap(~ paste("Variable:", variable), ncol = 2, scales = "free_x") +
  theme_classic() + geom_hline(aes(yintercept = exp), linetype = 3) +
  geom_hline(aes(yintercept = exp.1), linetype = 2) +
  geom_hline(aes(yintercept = exp.2), linetype = 2) +
  coord_cartesian(ylim = c(-.2, .05)) +
  scale_shape_manual(values = c(1, 3)) +
  scale_y_continuous(labels = scales::percent, breaks = c(), sec.axis = sec_axis(
    ~ ., labels = scales::percent, breaks = c(0, seq(-.55, .05, .05)))) +
  labs(y = "", fill = "Model", shape = "Model", tag = "",
       size = "% fitted outside 0 - 1", x = "Condition number") +
  theme(legend.position = "top", axis.title.x = element_text(hjust = -.2)) +
  guides(fill = FALSE, shape = FALSE)
ecr_3$exp.1 <- .925
ecr_3$exp.2 <- .975
ecr_3_sub <- filter(ecr_3, exceed != 0)
p7 <- ggplot(filter(ecr_3_sub, variable == "binary"),
             aes(factor(Condition), ecr, fill = model)) +
  geom_point(aes(size = exceed.f, shape = model), position = position_dodge(.75)) +
  geom_point(size = .01, position = position_dodge(.75)) +
  facet_wrap(~ paste("Variable:", variable), ncol = 2, scales = "free_x") +
  theme_classic() + geom_hline(aes(yintercept = exp), linetype = 3) +
  geom_hline(aes(yintercept = exp.1), linetype = 2) +
  geom_hline(aes(yintercept = exp.2), linetype = 2) +
  scale_shape_manual(values = c(1, 3)) +
  scale_y_continuous(labels = percent, breaks = c(.975, .925, seq(.05, .85, .1))) +
  coord_cartesian(ylim = c(.45, .975)) +
  labs(y = "Empirical coverage rate (95% CI)", fill = "Model", shape = "Model",
       size = "% fitted outside 0 - 1", x = "", tag = "B") +
  theme(legend.position = "bottom") +
  guides(fill = FALSE, size = FALSE, shape = guide_legend(keywidth = .05, default.unit = "inch"))
p8 <- ggplot(filter(ecr_3_sub, variable == "continuous"),
             aes(factor(Condition), ecr, fill = model)) +
  geom_point(aes(size = exceed.f, shape = model), position = position_dodge(.75)) +
  geom_point(size = .01, position = position_dodge(.75)) +
  facet_wrap(~ paste("Variable:", variable), ncol = 2, scales = "free_x") +
  theme_classic() + geom_hline(aes(yintercept = exp), linetype = 3) +
  geom_hline(aes(yintercept = exp.1), linetype = 2) +
  geom_hline(aes(yintercept = exp.2), linetype = 2) +
  scale_shape_manual(values = c(1, 3)) +
  scale_y_continuous(labels = percent, breaks = c(), sec.axis = sec_axis(
    ~ ., labels = percent, breaks = c(.975, .925, seq(.05, .85, .1)))) +
  coord_cartesian(ylim = c(.45, .975)) +
  labs(y = "", fill = "Model", shape = "Model",
       size = "% fitted outside 0 - 1", tag = "",
       x = "Condition number") +
  theme(legend.position = "bottom", axis.title.x = element_text(hjust = -.2)) +
  guides(fill = FALSE, shape = FALSE)
grid.arrange(p5, p6, p7, p8, layout_matrix = lay)

# Figure 5

dat.new <- readRDS("models_ecls_plot.RDS")
ggplot(dat.new, aes(SES, fitted, linetype = approach)) + geom_line() +
  theme_classic() + scale_y_continuous(labels = scales::percent) + facet_wrap(~ RD) +
  coord_cartesian(ylim = 0:1) + scale_linetype_manual(values = c(1, 3, 2)) +
  labs(x = "Family SES prior to kindergarten", y = "Predicted probabilities",
       linetype = "Approach") + theme(legend.position = "top")

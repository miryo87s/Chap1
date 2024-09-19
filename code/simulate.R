librarian::shelf(boot, data.table, MASS, ggplot2, sandwich, ggthemes, MetBrewer, dplyr, tidyr, hdm, randomForest, quiet = T)

#### DGP for CATE ####

df_sim1 <- function(n, tau = 1) {
  x <- matrix(rnorm(n * 2), ncol = 2)
  
  w <- rbinom(n, 1, 0.5)
  
  beta <- rep(1, 2)
  error <- rnorm(n)
  y <- 1 + tau * (x[, 1] + x[, 2]) * w + as.matrix(x) %*% beta + error
  
  as.data.frame(cbind(y, w, x)) %>% setNames(c("y", "w", paste0("x", 1:2)))
}

df_sim2 <- function(n, tau = 1) {
  x <- matrix(rnorm(n * 2), ncol = 2)

  k <- - 1 - 0.4 * x[, 1] + 0.2 * x[, 2]
  pr <- 1 / (1 + exp(-k))
  w <- rbinom(n, 1, pr)
  
  beta <- rep(1, 2)
  error <- rnorm(n)
  y <- 1 + tau * (x[, 1] + x[, 2]) * w + x %*% beta + error
  
  as.data.frame(cbind(y, w, x)) %>% setNames(c("y", "w", paste0("x", 1:2)))
}

df_sim3 <- function(n, tau = 1) {
  x <- matrix(rnorm(n * 2), ncol = 2)
  
  k <- - 1 - 0.4 * x[, 1] + 0.2 * x[, 2]
  pr <- 1 / (1 + exp(-k))
  w <- rbinom(n, 1, pr)
  
  beta <- rep(1, 2)
  error <- rnorm(n)
  y <- 1 + tau * (x[, 1] + x[, 2] + x[, 1] * x[, 2]) * w + x %*% beta + error
  
  as.data.frame(cbind(y, w, x)) %>% setNames(c("y", "w", paste0("x", 1:2)))
}

df_sim4 <- function(n, tau = 1) {
  x <- matrix(rnorm(n * 2), ncol = 2)
  
  k <- - 1 - 0.4 * x[, 1] + 0.2 * x[, 2]
  pr <- 1 / (1 + exp(-k))
  w <- rbinom(n, 1, pr)
  
  beta <- rep(1, 2)
  error <- rnorm(n)
  y <- 1 + tau * (ifelse(x[, 1] <= 0, -1, 1)) * w + x %*% beta + error
  
  as.data.frame(cbind(y, w, x)) %>% setNames(c("y", "w", paste0("x", 1:2)))
}


# DGP: Linear projection fails, add qaudratic term
df_sim5 <- function(n, tau = 1) {
  x <- matrix(rnorm(n * 2), ncol = 2)
  
  k <- - 1 - 0.4 * x[, 1] + 0.2 * x[, 2]
  pr <- 1 / (1 + exp(-k))
  w <- rbinom(n, 1, pr)
  
  beta <- rep(1, 2)
  error <- rnorm(n)
  y <- 1 + tau * (ifelse(x[, 1] < -0.5 | x[, 1] > 0.5, -1, 1)) * w + x %*% beta + error
  
  as.data.frame(cbind(y, w, x)) %>% setNames(c("y", "w", paste0("x", 1:2)))
}

df_sim6 <- function(n, tau = 1) {
  x <- matrix(rnorm(n * 2), ncol = 2)
  
  k <- - 1 - 0.4 * x[, 1] + 0.2 * x[, 2]
  pr <- 1 / (1 + exp(-k))
  w <- rbinom(n, 1, pr)
  
  beta <- rep(1, 2)
  error <- rnorm(n)
  y <- 1 + tau * (cos(x[, 1]) + cos(x[, 2])) * w + x %*% beta + error
  
  as.data.frame(cbind(y, w, x)) %>% setNames(c("y", "w", paste0("x", 1:2)))
}

#### DGP for CLATE ####

df_sim7 <- function(n, tau = 1) {
  x <- rnorm(n)
  k <- 1 - 0.3 * x
  pr.z <- 1 / (1 + exp(-k))
  z <- rbinom(n, 1, pr.z)
  
  pz <- ifelse(z == 1, 0.5, -0.5)
  uv <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.5, 0.5, 1), ncol = 2))
  w <- ifelse(uv[, 2] <= pz, 1, 0)
  
  y <- 1 + x + tau * x * w + uv[, 1]
  
  as.data.frame(cbind(y, w, z, x)) %>% setNames(c("y", "w", "z", "x"))
}

#### Test Functions ####

dr_test <- function(df, degree = 1, p = 2) {
  mod.pr <- glm(w ~ . - y, data = df, family = binomial(link = "logit"))
  mod.or <- lm(y ~ . * w, data = df)
  propensity <- mod.pr$fitted.values
  cond.mean <- mod.or$fitted.values
  res.pr <- mod.pr$residuals
  res.or <- mod.or$residuals
  
  df.treatall <- df.treatnone <- df
  df.treatall$w <- 1
  df.treatnone$w <- 0
  m1 <- predict(mod.or, df.treatall)
  m0 <- predict(mod.or, df.treatnone)
  
  df.dr <- df %>% 
    mutate(dr = (w - propensity) * (y - cond.mean) / (propensity * (1 - propensity)) + m1 - m0, 
           cond.mean = cond.mean, propensity = propensity)
  
  if (p == 1) {
    mod.cate <- lm(dr ~ poly(x1, degree, raw = TRUE), data = df.dr)
  } else {
    mod.cate <- lm(dr ~ poly(x1, degree, raw = TRUE) + poly(x2, degree, raw = TRUE), data = df.dr)
  }
  beta <- coef(mod.cate)[-1]
  Sig <- vcov(mod.cate)[-1, -1]
  Tstats <- t(beta) %*% solve(Sig) %*% beta
  list(pvalue = 1 - pchisq(Tstats, df = length(beta)), Tstats = Tstats, df.dr = df.dr)
}

dr_iv_test <- function(df) {
  mod.pr.z <- glm(z ~ x, data = df, family = binomial(link = "logit"))
  mod.or.w <- lm(w ~ z, data = df)
  mod.or.y <- lm(y ~ z * x, data = df)
  propensity <- mod.pr.z$fitted.values
  cond.mean.w <- mod.or.w$fitted.values
  cond.mean.y <- mod.or.y$fitted.values
  res.pr <- mod.pr.z$residuals
  res.or.w <- mod.or.w$residuals
  res.or.y <- mod.or.y$residuals
  
  df.treatall <- df.treatnone <- df
  df.treatall$z <- 1
  df.treatnone$z <- 0
  m1.w <- predict(mod.or.w, df.treatall)
  m0.w <- predict(mod.or.w, df.treatnone)
  m1.y <- predict(mod.or.y, df.treatall)
  m0.y <- predict(mod.or.y, df.treatnone)
  
  df.dr <- df %>% 
    mutate(dr.w = (z - propensity) * (w - cond.mean.w) / (propensity * (1 - propensity)) + m1.w - m0.w, 
           dr.y = (z - propensity) * (y - cond.mean.y) / (propensity * (1 - propensity)) + m1.y - m0.y,
           cond.mean.w = cond.mean.w,
           cond.mean.y = cond.mean.y,
           propensity = propensity)
  
  n <- nrow(df.dr)
  design.mat <- rbind(cbind(rep(1, n), rep(0, n), as.matrix(df.dr$x), matrix(rep(0, n), n, 1)), 
                      cbind(rep(0, n), rep(1, n), matrix(rep(0, n), n, 1), as.matrix(df.dr$x)))
  outcome <- with(df.dr, c(dr.y, dr.w))
  dr.reg <- as.data.frame(cbind(outcome, design.mat))
  mod.clate <- lm(outcome ~ . - 1, data = dr.reg)
  
  est <- coef(mod.clate)
  cov <- vcovHC(mod.clate, type = "HC0", cluster = c(rep(1, nrowdf.dr), rep(0, nrwo(df.dr))))
  
  h <- est[3] - est[1] * est[4] / est[2]
  h.prime <- c(-est[4]/est[2], est[1]*est[4]/est[2]^2, 1, -est[1]/est[2])
  vcov.h <- t(h.prime) %*% cov %*% h.prime
  
  Tstats <- h / sqrt(vcov.h)
  pvalue <- 2 * (1 - pnorm(abs(Tstats)))
  
  list(pvalue = pvalue, Tstats = Tstats, df.dr = df.dr)
}

#### Power Analysis for CATE ####

power_function <- function(effect, N, nsim, dgp, degree = 1) {
  result <- replicate(nsim, dr_test(dgp(N, tau = effect), degree)$pvalue <= 0.05)
  mean(result)
}

effect.size <- seq(-1, 1, 0.1)
nsim <- 10000

set.seed(1)

power.sim1.100 <- sapply(effect.size, \(i) power_function(i, N = 100, nsim = nsim, df_sim1))
power.sim1.300 <- sapply(effect.size, \(i) power_function(i, N = 300, nsim = nsim, df_sim1))
power.sim1.1000 <- sapply(effect.size, \(i) power_function(i, N = 1000, nsim = nsim, df_sim1))

power.sim2.100 <- sapply(effect.size, \(i) power_function(i, N = 100, nsim = nsim, df_sim2))
power.sim2.300 <- sapply(effect.size, \(i) power_function(i, N = 300, nsim = nsim, df_sim2))
power.sim2.1000 <- sapply(effect.size, \(i) power_function(i, N = 1000, nsim = nsim, df_sim2))

power.sim3.100 <- sapply(effect.size, \(i) power_function(i, N = 100, nsim = nsim, df_sim3))
power.sim3.300 <- sapply(effect.size, \(i) power_function(i, N = 300, nsim = nsim, df_sim3))
power.sim3.1000 <- sapply(effect.size, \(i) power_function(i, N = 1000, nsim = nsim, df_sim3))

power.sim4.100 <- sapply(effect.size, \(i) power_function(i, N = 100, nsim = nsim, df_sim4))
power.sim4.300 <- sapply(effect.size, \(i) power_function(i, N = 300, nsim = nsim, df_sim4))
power.sim4.1000 <- sapply(effect.size, \(i) power_function(i, N = 1000, nsim = nsim, df_sim4))

power.sim1 <- data.frame(power.sim1.100, power.sim1.300, power.sim1.1000) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:3, names_to = "N", values_to = "Power")

power.sim2 <- data.frame(power.sim2.100, power.sim2.300, power.sim2.1000) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:3, names_to = "N", values_to = "Power")

power.sim3 <- data.frame(power.sim3.100, power.sim3.300, power.sim3.1000) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:3, names_to = "N", values_to = "Power")

power.sim4 <- data.frame(power.sim4.100, power.sim4.300, power.sim4.1000) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:3, names_to = "N", values_to = "Power")

power.sim <- cbind(power.sim1, power.sim2$Power, power.sim3$Power, power.sim4$Power) %>%
  setNames(c("effect", "N", paste0("design", 1:4))) %>%
  pivot_longer(cols = 3:6, names_to = "Design", values_to = "Power") %>%
  mutate(N = recode(N, "power.sim1.100" = 100, "power.sim1.300" = 300, "power.sim1.1000" = 1000))

# p1 <- ggplot(power.sim, aes(x = effect, y = Power, group = as.factor(N))) +
#   geom_line(aes(linetype = as.factor(N)), linewidth = 0.6) + 
#   facet_wrap(~Design, nrow = 1, ncol = 4) + 
#   ylab("Simulated Rejection Rate") +
#   scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 1.00)) +
#   scale_linetype_manual(values=c("twodash", "dotted", "solid"), name = "Sample Size") + 
#   theme_bw() + guides(linetype = guide_legend(reverse = T)) +
#   theme(strip.background = element_blank(), 
#         legend.title = element_text(size = 7), 
#         legend.text = element_text(size = 6))

#### Power Analysis for CLATE

power_function_iv <- function(effect, N, nsim, dgp = df_sim7) {
  result <- replicate(nsim, dr_iv_test(dgp(N, tau = effect))$pvalue <= 0.05)
  mean(result)
}

power.sim7.100 <- sapply(effect.size, \(i) power_function_iv(i, N = 400, nsim = nsim, df_sim7))
power.sim7.300 <- sapply(effect.size, \(i) power_function_iv(i, N = 1200, nsim = nsim, df_sim7))
power.sim7.1000 <- sapply(effect.size, \(i) power_function_iv(i, N = 4000, nsim = nsim, df_sim7))

power.sim7 <- data.frame(power.sim7.100, power.sim7.300, power.sim7.1000) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:3, names_to = "N", values_to = "Power")

power.sim.alt <- cbind(power.sim1, power.sim2$Power, power.sim3$Power, power.sim4$Power, power.sim7$Power) %>%
  setNames(c("effect", "N", paste0("design", 1:4), "IV")) %>%
  pivot_longer(cols = 3:7, names_to = "Design", values_to = "Power") %>%
  mutate(N = recode(N, "power.sim1.100" = 100, "power.sim1.300" = 300, "power.sim1.1000" = 1000))

p1 <- ggplot(power.sim.alt, aes(x = effect, y = Power, group = as.factor(N))) +
  geom_line(aes(linetype = as.factor(N)), linewidth = 0.6) + 
  facet_wrap(~Design, nrow = 1, ncol = 5) + 
  ylab("Simulated Rejection Rate") +
  scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 1.00)) +
  scale_linetype_manual(values=c("twodash", "dotted", "solid"), name = "Sample Size") + 
  theme_gdocs() +
  theme(strip.background = element_blank(),
        text = element_text(size = 10))
p1

#### Size Analysis ####

set.seed(1)

size_function <- function(crit, N, nsim, dgp) {
  result <- replicate(nsim, dr_test(dgp(N, tau = 0))$pvalue <= crit)
  mean(result)
}

size_function_iv <- function(crit, N, nsim, dgp) {
  result <- replicate(nsim, dr_iv_test(dgp(N, tau = 0))$pvalue <= crit)
  mean(result)
}

pcrit <- c(0.1, 0.05, 0.01)

size1.100 <- sapply(pcrit, \(p) size_function(p, N = 100, nsim = nsim, df_sim1))
size1.300 <- sapply(pcrit, \(p) size_function(p, N = 300, nsim = nsim, df_sim1))
size1.1000 <- sapply(pcrit, \(p) size_function(p, N = 1000, nsim = nsim, df_sim1))

size2.100 <- sapply(pcrit, \(p) size_function(p, N = 100, nsim = nsim, df_sim2))
size2.300 <- sapply(pcrit, \(p) size_function(p, N = 300, nsim = nsim, df_sim2))
size2.1000 <- sapply(pcrit, \(p) size_function(p, N = 1000, nsim = nsim, df_sim2))

size3.100 <- sapply(pcrit, \(p) size_function(p, N = 100, nsim = nsim, df_sim3))
size3.300 <- sapply(pcrit, \(p) size_function(p, N = 300, nsim = nsim, df_sim3))
size3.1000 <- sapply(pcrit, \(p) size_function(p, N = 1000, nsim = nsim, df_sim3))

size4.100 <- sapply(pcrit, \(p) size_function(p, N = 100, nsim = nsim, df_sim4))
size4.300 <- sapply(pcrit, \(p) size_function(p, N = 300, nsim = nsim, df_sim4))
size4.1000 <- sapply(pcrit, \(p) size_function(p, N = 1000, nsim = nsim, df_sim4))

size7.400 <- sapply(pcrit, \(p) size_function_iv(p, N = 400, nsim = nsim, df_sim7))
size7.1200 <- sapply(pcrit, \(p) size_function_iv(p, N = 1200, nsim = nsim, df_sim7))
size7.4000 <- sapply(pcrit, \(p) size_function_iv(p, N = 4000, nsim = nsim, df_sim7))

size.sim <- rbind(size1.100, size1.300, size1.1000, 
                  size2.100, size2.300, size2.1000, 
                  size3.100, size3.300, size3.1000, 
                  size4.100, size4.300, size4.1000, 
                  size7.400, size7.1200, size7.4000)






#### Add Quadratic terms ####

set.seed(0)

power1.sim5 <- sapply(effect.size, \(i) power_function(i, N = 1000, nsim = nsim, df_sim5))
power1.sim6 <- sapply(effect.size, \(i) power_function(i, N = 1000, nsim = nsim, df_sim6))
power2.sim5 <- sapply(effect.size, \(i) power_function(i, N = 1000, nsim = nsim, df_sim5, degree = 2))
power2.sim6 <- sapply(effect.size, \(i) power_function(i, N = 1000, nsim = nsim, df_sim6, degree = 2))

power.sim5 <- data.frame(power1.sim5, power2.sim5) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:2, names_to = "model", values_to = "Power")

power.sim6 <- data.frame(power1.sim6, power2.sim6) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:2, names_to = "model", values_to = "Power")

power2.sim <- cbind(power.sim5, power.sim6$Power) %>%
  setNames(c("effect", "model", paste0("design", 5:6))) %>%
  pivot_longer(cols = 3:4, names_to = "Design", values_to = "Power") %>%
  mutate(model = recode(model, "power1.sim5" = 0, "power2.sim5" = 1))

p3 <- ggplot(power2.sim, aes(x = effect, y = Power, color = as.factor(model))) +
  geom_line(linewidth = 0.6) + 
  facet_wrap(~Design, nrow = 1, ncol = 2) + 
  scale_color_discrete(name = "Quadratic") +
  scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 1.00)) +
  theme_bw() + guides(linetype = guide_legend(reverse = T)) +
  theme(strip.background = element_blank(), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6))

####  Comparison with Crump test ####

set.seed(0)

effect.size <- seq(-1, 1, 0.1)
nsim <- 10000

crump_test <- function(df, degree = 1, p = 2) {
  mod0 <- lm(y ~ . - w, data = df %>% filter(w == 0))
  mod1 <- lm(y ~ . - w, data = df %>% filter(w == 1))
  
  coef0 <- coef(mod0)[-1]
  coef1 <- coef(mod1)[-1]
  cov.mat0 <- vcov(mod0)[-1, -1]
  cov.mat1 <- vcov(mod1)[-1, -1]
  
  Tstats <- t(matrix(coef0 - coef1)) %*% solve(cov.mat0 + cov.mat1) %*% matrix(coef0 - coef1)
  pvalue <- 1 - pchisq(Tstats, df = length(variable.names(mod1)))
  list(pvalue = 1 - pchisq(Tstats, df = length(beta)), Tstats = Tstats)
}

power_crump <- function(effect, N, nsim, dgp) {
  result <- replicate(nsim, crump_test(dgp(N, tau = effect))$pvalue <= 0.05)
  mean(result)
}

crump.sim1.300 <- sapply(effect.size, \(i) power_crump(i, 300, nsim, df_sim1))
crump.sim2.300 <- sapply(effect.size, \(i) power_crump(i, 300, nsim, df_sim2))
crump.sim3.300 <- sapply(effect.size, \(i) power_crump(i, 300, nsim, df_sim3))
crump.sim4.300 <- sapply(effect.size, \(i) power_crump(i, 300, nsim, df_sim4))

# adjustment_test <- function(df, degree = 1, p = 2) {
#   df <- df %>% mutate_at(names(.)[- c(1, 2)], list(d = ~ w * (. - mean(.))))
# 
#   mod <- lm(y ~ ., data = df)
#   beta <- tail(coef(mod), p)
#   Sig <- vcov(mod)[- c(1:(p + 2)), - c(1:(p + 2))]
# 
#   Tstats <- t(beta) %*% solve(Sig) %*% beta
#   list(pvalue = 1 - pchisq(Tstats, df = length(beta)), Tstats = Tstats)
# }
# 
# power_adjustment <- function(effect, N, nsim, dgp) {
#   result <- replicate(nsim, adjustment_test(dgp(N, tau = effect))$pvalue <= 0.05)
#   mean(result)
# }
# 
# adjustment.sim1.300 <- sapply(effect.size, \(i) power_adjustment(i, 300, nsim, df_sim1))
# adjustment.sim2.300 <- sapply(effect.size, \(i) power_adjustment(i, 300, nsim, df_sim2))
# adjustment.sim3.300 <- sapply(effect.size, \(i) power_adjustment(i, 300, nsim, df_sim3))
# adjustment.sim4.300 <- sapply(effect.size, \(i) power_adjustment(i, 300, nsim, df_sim4))

# power.comparison1 <- data.frame(adjustment.sim1.300, crump.sim1.300, power.sim1.300) %>% 
#   mutate(effect = effect.size) %>%
#   pivot_longer(cols = 1:3, names_to = "method", values_to = "Power")
# 
# power.comparison2 <- data.frame(adjustment.sim2.300, crump.sim2.300, power.sim2.300) %>% 
#   mutate(effect = effect.size) %>%
#   pivot_longer(cols = 1:3, names_to = "method", values_to = "Power")
# 
# power.comparison3 <- data.frame(adjustment.sim3.300, crump.sim3.300, power.sim3.300) %>% 
#   mutate(effect = effect.size) %>%
#   pivot_longer(cols = 1:3, names_to = "method", values_to = "Power")
# 
# power.comparison4 <- data.frame(adjustment.sim4.300, crump.sim4.300, power.sim4.300) %>% 
#   mutate(effect = effect.size) %>%
#   pivot_longer(cols = 1:3, names_to = "method", values_to = "Power")

power.comparison1 <- data.frame(crump.sim1.300, power.sim1.300) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:2, names_to = "method", values_to = "Power")

power.comparison2 <- data.frame(crump.sim2.300, power.sim2.300) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:2, names_to = "method", values_to = "Power")

power.comparison3 <- data.frame(crump.sim3.300, power.sim3.300) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:2, names_to = "method", values_to = "Power")

power.comparison4 <- data.frame(crump.sim4.300, power.sim4.300) %>% 
  mutate(effect = effect.size) %>%
  pivot_longer(cols = 1:2, names_to = "method", values_to = "Power")

power.comparison <- cbind(power.comparison1, power.comparison2$Power, power.comparison3$Power, power.comparison4$Power) %>%
  setNames(c("effect", "method", paste0("design", 1:4))) %>%
  pivot_longer(cols = 3:6, names_to = "Design", values_to = "Power") %>%
  mutate(method = recode(method, "crump.sim1.300" = "Crump et al.", "power.sim1.300" = "AIPW"))

p.comparison <- ggplot(power.comparison, aes(x = effect, y = Power, group = method)) +
  geom_line(aes(linetype = method), linewidth = 0.6) + 
  facet_wrap(~Design, nrow = 1, ncol = 4) + 
  ylab("Simulated Rejection Rate") +
  scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 1.00)) + 
  theme_gdocs() +
  theme(strip.background = element_blank(), 
        text = element_text(size = 10))
p.comparison

#### Tests with DoubleML ####

df_sim_ml1 <- function(n, tau = 1) {
  x <- matrix(rnorm(n * 100), ncol = 100)
  beta_w <- c(rep(0.1, 10), rep(0, 90))
  
  k <- 1 + x %*% beta_w
  pr <- 1 / (1 + exp(-k))
  w <- rbinom(n, 1, pr)
  
  beta <- c(rep(1, 10), rep(0, 90))
  error <- rnorm(n)
  y <- 1 + tau * x %*% beta * w + x %*% beta + error
  
  as.data.frame(cbind(y, w, x)) %>% setNames(c("y", "w", paste0("x", 1:100)))
}

df_sim_ml2 <- function(n, tau = 1) {
  x <- matrix(runif(n * 10), ncol = 10)
  beta_w <- rep(0.1, 10)
  
  k <- 1 + x %*% beta_w
  pr <- 1 / (1 + exp(-k))
  w <- rbinom(n, 1, pr)
  
  beta <- rep(1, 10)
  error <- rnorm(n)
  y <- 1 + tau * x %*% beta * w + sin(x %*% beta) + error
  
  as.data.frame(cbind(y, w, x)) %>% setNames(c("y", "w", paste0("x", 1:10)))
}

dml_test1 <- function(df, wreg = rlasso, yreg = rlasso, nfold = 5) {
  nobs <- nrow(df) #number of observations
  foldid <- rep.int(1:nfold, times = ceiling(nobs / nfold))[sample.int(nobs)] #define folds indices
  I <- split(1:nobs, foldid)  #split observation indices into folds
  
  x <- as.matrix(df[, -c(1:2)])
  w <- df$w
  y <- df$y
  
  propensity <- m0 <- m1 <- rep(NA, nobs)
  for(b in 1:length(I)){
    mod.pr <- wreg(x[-I[[b]],], w[-I[[b]]])
    mod.or1 <- yreg(subset(x[-I[[b]],], w[-I[[b]]] == 1), subset(y[-I[[b]]], w[-I[[b]]] == 1))
    mod.or0 <- yreg(subset(x[-I[[b]],], w[-I[[b]]] == 0), subset(y[-I[[b]]], w[-I[[b]]] == 0))
    
    propensity[I[[b]]] <- predict(mod.pr, x[I[[b]],], type = "response") 
    m1[I[[b]]] <- predict(mod.or1, x[I[[b]],], type = "response")
    m0[I[[b]]] <- predict(mod.or0, x[I[[b]],], type = "response")
  }
  df.dr <- df %>% 
    mutate(dr = (y - m1) * w / propensity - (y - m0) * (1 - w)/(1 - propensity) + m1 - m0, 
           m1 = m1, m0 = m0, propensity = propensity)
  
  mod.cate <- lm(dr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, data = df.dr)
  beta <- coef(mod.cate)[-1]
  Sig <- vcov(mod.cate)[-1, -1]
  Tstats <- t(beta) %*% solve(Sig) %*% beta
  list(pvalue = 1 - pchisq(Tstats, df = length(beta)), Tstats = Tstats, df.dr = df.dr)
}

# dml_test2 <- function(df, wreg = randomForest, yreg = randomForest, nfold = 5) {
#   nobs <- nrow(df) #number of observations
#   foldid <- rep.int(1:nfold, times = ceiling(nobs / nfold))[sample.int(nobs)] #define folds indices
#   I <- split(1:nobs, foldid)  #split observation indices into folds
#   
#   x <- as.matrix(df[, -c(1:2)])
#   w <- as.factor(df$w)
#   y <- df$y
#   
#   propensity <- m0 <- m1 <- rep(NA, nobs)
#   for(b in 1:length(I)){
#     mod.pr <- wreg(x[-I[[b]],], w[-I[[b]]])
#     mod.or1 <- yreg(subset(x[-I[[b]],], w[-I[[b]]] == 1), subset(y[-I[[b]]], w[-I[[b]]] == 1))
#     mod.or0 <- yreg(subset(x[-I[[b]],], w[-I[[b]]] == 0), subset(y[-I[[b]]], w[-I[[b]]] == 0))
#     
#     propensity[I[[b]]] <- predict(mod.pr, x[I[[b]],], type = "prob")[, 2]
#     m1[I[[b]]] <- predict(mod.or1, x[I[[b]],], type = "response")
#     m0[I[[b]]] <- predict(mod.or0, x[I[[b]],], type = "response")
#   }
#   df.dr <- df %>% 
#     mutate(dr = (y - m1) * w / propensity - (y - m0) * (1 - w)/(1 - propensity) + m1 - m0, 
#            m1 = m1, m0 = m0, propensity = propensity)
#   
#   mod.cate <- lm(dr ~ . , data = df.dr)
#   beta <- coef(mod.cate)[-1]
#   Sig <- vcov(mod.cate)[-1, -1]
#   Tstats <- t(beta) %*% solve(Sig) %*% beta
#   list(pvalue = 1 - pchisq(Tstats, df = length(beta)), Tstats = Tstats, df.dr = df.dr)
# }

dml_test2 <- function(df, wreg = xgboost, yreg = xgboost, nfold = 5) {
  nobs <- nrow(df) #number of observations
  foldid <- rep.int(1:nfold, times = ceiling(nobs / nfold))[sample.int(nobs)] #define folds indices
  I <- split(1:nobs, foldid)  #split observation indices into folds
  
  x <- as.matrix(df[, -c(1:2)])
  w <- df$w
  y <- df$y
  
  propensity <- m0 <- m1 <- rep(NA, nobs)
  for(b in 1:length(I)){
    mod.pr <- wreg(x[-I[[b]],], w[-I[[b]]], max.depth = 2, eta = 1, nthread = 2, nrounds = 5, objective = "binary:logistic")
    mod.or1 <- yreg(subset(x[-I[[b]],], w[-I[[b]]] == 1), subset(y[-I[[b]]], w[-I[[b]]] == 1), max.depth = 2, eta = 0.5, nthread = 2, nrounds = 5)
    mod.or0 <- yreg(subset(x[-I[[b]],], w[-I[[b]]] == 0), subset(y[-I[[b]]], w[-I[[b]]] == 0), max.depth = 2, eta = 0.5, nthread = 2, nrounds = 5)
    
    propensity[I[[b]]] <- predict(mod.pr, x[I[[b]],])
    m1[I[[b]]] <- predict(mod.or1, x[I[[b]],])
    m0[I[[b]]] <- predict(mod.or0, x[I[[b]],])
  }
  df.dr <- df %>% 
    mutate(dr = (y - m1) * w / propensity - (y - m0) * (1 - w)/(1 - propensity) + m1 - m0, 
           m1 = m1, m0 = m0, propensity = propensity)
  
  mod.cate <- lm(dr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, data = df.dr)
  beta <- coef(mod.cate)[-1]
  Sig <- vcov(mod.cate)[-1, -1]
  Tstats <- t(beta) %*% solve(Sig) %*% beta
  list(pvalue = 1 - pchisq(Tstats, df = length(beta)), Tstats = Tstats, df.dr = df.dr)
}


dml_size_function <- function(dml_test, crit, N, nsim, dgp) {
  result <- replicate(nsim, dml_test(dgp(N, tau = 0))$pvalue <= crit)
  mean(result)
}

pcrit <- c(0.1, 0.05, 0.01)
nsim <- 10000

set.seed(0)

dml1.size.500 <- sapply(pcrit, \(p) dml_size_function(dml_test1, p, 500, nsim, df_sim_ml1))
dml1.size.1000 <- sapply(pcrit, \(p) dml_size_function(dml_test1, p, 1000, nsim, df_sim_ml1))
dml1.size.2000 <- sapply(pcrit, \(p) dml_size_function(dml_test1, p, 2000, nsim, df_sim_ml1))
dml1.size.5000 <- sapply(pcrit, \(p) dml_size_function(dml_test1, p, 5000, nsim, df_sim_ml1))

dml2.size.500 <- sapply(pcrit, \(p) dml_size_function(dml_test2, p, 500, nsim, df_sim_ml2))
dml2.size.1000 <- sapply(pcrit, \(p) dml_size_function(dml_test2, p, 1000, nsim, df_sim_ml2))
dml2.size.2000 <- sapply(pcrit, \(p) dml_size_function(dml_test2, p, 2000, nsim, df_sim_ml2))
dml2.size.5000 <- sapply(pcrit, \(p) dml_size_function(dml_test2, p, 5000, nsim, df_sim_ml2))





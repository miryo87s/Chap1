librarian::shelf(MatchIt, hdm, glmnet, purrr, ggplot2, sandwich, car, dplyr, quiet = T)

#### 401k ####

data(pension)
df <- pension %>% mutate(net_tfa = net_tfa * 1e-3, inc = inc * 1e-3) %>%
  select(net_tfa, p401, e401, age, male, inc, educ, fsize, marr, twoearn, db, hown, pira)

fm.full.y <- net_tfa ~ e401 * (age + male + inc + educ + fsize + marr + twoearn + db + hown + pira)
fm.full.w <- p401 ~ e401 * (age + male + inc + educ + fsize + marr + twoearn + db + hown + pira)
fm.z <- e401 ~ age + male + inc + educ + fsize + marr + twoearn + db + hown + pira
fm.cate.y <- dr.y ~ age + male + inc + educ + fsize + marr + twoearn + db + hown + pira
fm.cate.w <- dr.w ~ age + male + inc + educ + fsize + marr + twoearn + db + hown + pira
mod.or.y <- lm(fm.full.y, data = df)
mod.or.w <- lm(fm.full.w, data = df)
mod.pr <- glm(fm.z, data = df, family = binomial(link = "logit"))

propensity <- mod.pr$fitted.values
cond.mean.y <- mod.or.y$fitted.values
cond.mean.w <- mod.or.w$fitted.values

df.treatall <- df.treatnone <- df
df.treatall$e401 <- 1
df.treatnone$e401 <- 0
m1y <- predict(mod.or.y, df.treatall)
m0y <- predict(mod.or.y, df.treatnone)
m1w <- predict(mod.or.w, df.treatall) %>% ifelse(. > 1, 1, .)
m0w <- 0

# LATE test
df.dr <- df %>% 
  mutate(dr.y = (e401 - propensity) * (net_tfa - cond.mean.y) / (propensity * (1 - propensity)) + m1y - m0y, 
         dr.w = (e401 - propensity) * (p401 - cond.mean.w) / (propensity * (1 - propensity)) + m1w - m0w) %>%
  filter(dr.y > -4000)

mod.cate.y <- lm(fm.cate.y, data = df.dr)
mod.cate.w <- lm(fm.cate.w, data = df.dr)

n <- nrow(df.dr)
p <- 10
X <- df.dr %>% select(age, male, inc, educ, fsize, marr, twoearn, db, hown, pira)
Xu <- rbind(cbind(rep(1, n), rep(0, n), as.matrix(X), matrix(rep(0, n * 10), n, 10)), 
            cbind(rep(0, n), rep(1, n), matrix(rep(0, n * 10), n, 10), as.matrix(X)))
Y <- c(df.dr$dr.y, df.dr$dr.w)
ivdf <- as.data.frame(cbind(Y, Xu))
ivBLP <- lm(Y ~ . - 1, data = ivdf)
R <- cbind(rep(0, p), rep(0, p), diag(1, p), diag(-1, p))
test <- linearHypothesis(ivBLP, R, rep(0, p), 
                         vcov = vcovHC(ivBLP, type = "HC0", cluster = c(rep(1, nrow(df.dr)), rep(0, nrow(df.dr)))), 
                         test = "Chisq")

# High-Dim test

xbasic <- pension %>% select(age, male, inc, educ, fsize, marr, twoearn, db, hown, pira)
x <- pension %>% mutate(inc = inc * 1e-3) %>%
  select(age, male, inc, educ, fsize, marr, twoearn, db, hown, pira, i1, i2, i3, i4, i5, i6, i7) %>%
  mutate(agesq = age^2, fsizesq = fsize^2, incsq = inc^2) %>%
  model.matrix( ~ .^2, data = .) 
y <- pension$net_tfa * 1e-3
w <- pension$p401

DML_AIPW <- function(x, w, y, wreg, yreg, nfold = 2, bx) {
  w0.ind <- which(w == 0)
  w1.ind <- which(w == 1)
  
  nobs <- nrow(x)
  foldid <- rep.int(1:nfold, times = ceiling(nobs / nfold))[sample.int(nobs)]
  I <- split(1:nobs, foldid)
  dr <- rep(NA, nobs)
  cat("fold: ")
  for(b in 1:length(I)) {
    wfit <- wreg(x[-I[[b]], ], w[-I[[b]]])
    y0fit <- yreg(x[-union(I[[b]], w1.ind), ], y[-union(I[[b]], w1.ind)])
    y1fit <- yreg(x[-union(I[[b]], w0.ind), ], y[-union(I[[b]], w0.ind)])
    
    what <- predict(wfit, x[I[[b]], ], type = "response")
    y1hat <- predict(y1fit, x[I[[b]], ], type = "response")
    y0hat <- predict(y0fit, x[I[[b]], ], type = "response")
    
    dr[I[[b]]] <- w[I[[b]]] * (y[I[[b]]] - y1hat) / what - (1 - w[I[[b]]]) * (y[I[[b]]] - y0hat) / (1 - what) + y1hat - y0hat
    cat(b, " ")
  }
  df.dr <- data.frame(cbind(dr, bx))
  blp <- lm(dr ~ ., data = df.dr)
  return(blp)
}

wreg <- function(x, w) {rlasso(x, w, post = T)}
yreg <- function(x, y) {rlasso(x, y, post = T)}
set.seed(0)
dml.blp <- DML_AIPW(x, w, y, wreg, yreg, nfold = 5, bx = xbasic)

coefs <- names(coef(dml.blp))
test <- linearHypothesis(dml.blp, coefs[-1], test = "Chisq")



#### One Child Policy & Mental Health ####

CFPS_2010 <- readr::read_csv("data/CFPS_2010.csv")
cfps <- CFPS_2010 %>% filter(fincome < 10^6, fage - age > 14, mage - age > 14)
df <- cfps %>% filter(confidence >= 1, anxiety >= 1, desperation >= 1) %>% 
  mutate(qk601 = qk601 / 100000, fincome = fincome / 100000)

covariates <- c("age", "urban", "gender", "fincome", "fage", "mage", "feduy", "meduy", "divorce", "remarriage",
                "han", "qe1", "qk601")
outcomes <- c("confidence", "anxiety", "desperation")
treatment <- "onechild"

fm.w <- paste(treatment, "~", paste(covariates, collapse = " + "))
fm.full.y <- lapply(outcomes, \(i) paste(i, " ~ onechild * (", paste(covariates, collapse = " + "), ")"))
fm.blp <- paste("dr ~ ", paste(covariates, collapse = " + "))

cate.blp <- function(fm.full.y, fm.w, fm.blp, df) {
  mod.or <- lm(fm.full.y, data = df)
  mod.pr <- glm(fm.w, data = df, family = binomial(link = "logit"))
  
  propensity <- mod.pr$fitted.values
  cond.mean <- mod.or$fitted.values
  res.pr <- mod.pr$residuals
  res.or <- mod.or$residuals
  
  df.treatall <- df.treatnone <- df
  df.treatall$onechild <- 1
  df.treatnone$onechild <- 0
  m1 <- predict(mod.or, df.treatall)
  m0 <- predict(mod.or, df.treatnone)
  
  df.dr <- df %>% 
    mutate(dr = (onechild - propensity) * (anxiety - cond.mean) / (propensity * (1 - propensity)) + m1 - m0) %>%
    filter(propensity > 0.05)
  mod.blp <- lm(fm.blp, data = df.dr)
  return(list(mod.blp, mod.pr))
}

blp.results <- lapply(fm.full.y, \(i) cate.blp(i, fm.w, fm.blp, df))

mod.blp <- blp.results[[2]][[1]]
coefs <- names(coef(mod.blp))
test0 <- linearHypothesis(mod.blp, coefs, test = "Chisq")
test1 <- linearHypothesis(mod.blp, coefs[-1], test = "Chisq")

# Propensity density
propden <- ggplot(data.frame(propensity = blp.results[[1]][[2]]$fitted.values), aes(x = propensity)) +
  geom_density() + 
  xlim(0, 1) +
  theme_bw()


# Crump test

crump <- function(fm, df, test = "cate") {
  mod1 <- lm(fm, data = df %>% filter(onechild == 1))
  mod0 <- lm(fm, data = df %>% filter(onechild == 0))
  
  if (test == "cate") {
    coef0 <- coef(mod0)
    coef1 <- coef(mod1)
    cov.mat0 <- vcov(mod0)
    cov.mat1 <- vcov(mod1)
  } else if (test == "hte") {
    coef0 <- coef(mod0)[-1]
    coef1 <- coef(mod1)[-1]
    cov.mat0 <- vcov(mod0)[-1, -1]
    cov.mat1 <- vcov(mod1)[-1, -1]
  }

  Tstats <- t(matrix(coef0 - coef1)) %*% solve(cov.mat0 + cov.mat1) %*% matrix(coef0 - coef1)
  pvalue <- 1 - pchisq(Tstats, df = length(variable.names(mod1)))
  return(c(Tstats, pvalue))
}

fm.y <- lapply(outcomes, \(i) paste(i, "~", paste(covariates, collapse = " + ")))
fm.y.quad <- lapply(fm.y, \(i) paste(i, " + I(age^2) + I(fincome^2) + I(fage^2) + I(mage^2) + I(feduy^2) + I(meduy^2) + I(qk601^2)"))

crump.test0 <- map(fm.y, crump, df = df, test = "cate")
crump.quad.test0 <- map(fm.y.quad, crump, df = df, test = "cate")

crump.test <- map(fm.y, crump, df = df, test = "hte")
crump.quad.test <- map(fm.y.quad, crump, df = df, test = "hte")





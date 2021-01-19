##### R Code Syntax for Ratcliff et al. paper titled "An Examination of the Longitudinal Stability of Psychological Measures Contained within the U.S. Army's Global Assessment Tool (GAT)"
## NOTE: for questions about the data or using the PDE, contact Nathaniel Ratcliff, email = nr3xe@virginia.edu

### Load GAT 1.0 Data

## Libraries used -----------------------------------------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(psych)
library(pastecs)
library(lavaan)

### Functions used -----------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to run through calculating alpha and omega reliabilities using variable name in output
bindr <- function(s1, s2, s3, s4, s5)
{
  a <- get(x = s1)
  b <- get(x = s2)
  c <- get(x = s3)
  d <- get(x = s4)
  e <- get(x = s5)
  r1.1 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(a), type = "alpha", interval.type = "parallel", conf.level = .95))) %>%
    mutate(n = sum(!is.na(data$adapt.1.scale))) %>% mutate(var = s1) 
  r1.2 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(a), type = "omega", interval.type = "mlr", conf.level = .95)))  %>%
    mutate(n = sum(!is.na(data$adapt.1.scale)))  %>% mutate(var =  s1)   
  r2.1 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(b), type = "alpha", interval.type = "parallel", conf.level = .95))) %>%
    mutate(n = sum(!is.na(data$adapt.2.scale)))  %>% mutate(var =  s2)  
  r2.2 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(b), type = "omega", interval.type = "mlr", conf.level = .95)))  %>%
    mutate(n = sum(!is.na(data$adapt.2.scale)))  %>% mutate(var =  s2)  
  r3.1 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(c), type = "alpha", interval.type = "parallel", conf.level = .95))) %>%
    mutate(n = sum(!is.na(data$adapt.3.scale)))  %>% mutate(var =  s3)  
  r3.2 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(c), type = "omega", interval.type = "mlr", conf.level = .95)))  %>%
    mutate(n = sum(!is.na(data$adapt.3.scale)))  %>% mutate(var =  s3) 
  r4.1 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(d), type = "alpha", interval.type = "parallel", conf.level = .95))) %>%
    mutate(n = sum(!is.na(data$adapt.4.scale)))  %>% mutate(var =  s4)  
  r4.2 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(d), type = "omega", interval.type = "mlr", conf.level = .95)))  %>%
    mutate(n = sum(!is.na(data$adapt.4.scale)))  %>% mutate(var =  s4)  
  r5.1 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(e), type = "alpha", interval.type = "parallel", conf.level = .95))) %>%
    mutate(n = sum(!is.na(data$adapt.5.scale)))  %>% mutate(var =  s5)  
  r5.2 <- suppressWarnings(as.data.frame(MBESS::ci.reliability(tidyr::drop_na(e), type = "omega", interval.type = "mlr", conf.level = .95)))  %>%
    mutate(n = sum(!is.na(data$adapt.5.scale)))  %>% mutate(var =  s5)  
  bindrel <- rbind(r1.1, r1.2, r2.1, r2.2, r3.1, r3.2, r4.1, r4.2, r5.1, r5.2)
  bindrel %>% dplyr::select(var, type, est, se, ci.lower, ci.upper, conf.level, n)
}

# Function for one-way, 5 timepoint repeated measures ANOVA with generalized eta squared + CIs
rANOVA.sum <- function(x1, x2, x3, x4, x5) {
  repmodel1 <- lm(cbind(x1, x2, x3, x4, x5) ~ 1)
  trial1 <- factor(c('Time 1', 'Time 2', 'Time 3', ' Time 4', 'Time 5'), ordered = FALSE)
  rANOVA1 <- car::Anova(repmodel1, idata = data.frame(trial1), idesign = ~ trial1, type = 'III', digits = 4)
  rANOVA1 <- summary(rANOVA1, multivariate = FALSE)
  a1 <- rANOVA1$univariate.tests %>% as.matrix() %>% unclass() %>% as.data.frame()
  a2 <- rANOVA1$pval.adjustments %>% as.matrix() %>% unclass() %>% as.data.frame()
  a3 <- rANOVA1$sphericity.tests %>% as.matrix() %>% unclass() %>% as.data.frame()
  a4 <- cbind(a1, a2, a3)
  a4 <- a4 %>% dplyr::rename(SS_n = `Sum Sq`, DF_n = `num Df`, SS_d = `Error SS`, DF_d = `den Df`, F.val = `F value`, p.val = `Pr(>F)`, GGe = `GG eps`, GG.p.val = `Pr(>F[GG])`,
                             HFe = `HF eps`, HFe.p.val = `Pr(>F[HF])`, spherie = `Test statistic`, spherie.p.val = `p-value`)
  ges.ci <- MOTE::ges.partial.SS.mix(dfm = a4[2,2], dfe = a4[2,4], ssm = a4[2,1], sss = a4[1,3], sse = a4[2,3], Fvalue = a4[2,5], a = 0.10)
  a4 <- a4 %>% dplyr::mutate(ges = ges.ci$ges, ci.lower = ges.ci$geslow, ci.upper = ges.ci$geshigh)
  print(a4, digits = 4)
}

# Create p-value matrix for correlation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}  

# Full model syntax generator for within SEM with compound symmetry variances
syntax.full.CorComp <- function(data, t1, t2, t3, t4, t5) {
  full <- paste0(
    '# full model;',
    t1, ' ~ i1*1;',
    t2, ' ~ i2*1;',
    t3, ' ~ i3*1;',
    t4, ' ~ i4*1;',
    t5, ' ~ i5*1;',
    
    t1, ' ~~ v1*',t1,";",
    t2, ' ~~ v1*',t2,";",
    t3, ' ~~ v1*',t3,";",
    t4, ' ~~ v1*',t4,";",
    t5, ' ~~ v1*',t5,";",
    
    t1, ' ~~ c*',t2,' + c*', t3, ' + c*', t4, ' + c*', t5, ';',
    t2, ' ~~ c*',t3,' + c*', t4, ' + c*', t5, ';',
    t3, ' ~~ c*',t4,' + c*', t5, ';',
    t4, ' ~~ c*',t5, ';')
  return(as.character(full))
}

# Null model syntax generator for within SEM with compound symmetry variances; equal intercepts/means
syntax.null.CorComp <-  function(data, t1, t2, t3, t4, t5) {
  null <- paste0(
    '# null model;',
    t1, ' ~ i1*1;',
    t2, ' ~ i1*1;',
    t3, ' ~ i1*1;',
    t4, ' ~ i1*1;',
    t5, ' ~ i1*1;',
    
    t1, ' ~~ v1*',t1,";",
    t2, ' ~~ v1*',t2,";",
    t3, ' ~~ v1*',t3,";",
    t4, ' ~~ v1*',t4,";",
    t5, ' ~~ v1*',t5,";",
    
    t1, ' ~~ c*',t2,' + c*', t3, ' + c*', t4, ' + c*', t5, ';',
    t2, ' ~~ c*',t3,' + c*', t4, ' + c*', t5, ';',
    t3, ' ~~ c*',t4,' + c*', t5, ';',
    t4, ' ~~ c*',t5, ';')
  return(as.character(null))
}

# Function to calculate degrees of freedom error within based on number of valid observations and levels of fixed factor
dfe.fun <- function(t1, t2, t3, t4, t5, n, k) {
  # Number of non-NA observation count
  s.t1 <- sum(!is.na(t1))
  s.t2 <- sum(!is.na(t2))
  s.t3 <- sum(!is.na(t3))
  s.t4 <- sum(!is.na(t4))
  s.t5 <- sum(!is.na(t5))
  total <- s.t1 + s.t2 + s.t3 + s.t4 + s.t5
  # Calc of df error
  dfe <- total - (n - 1) - (k - 1)
  return(dfe)
}

# Function to get F, eta2, and 90% CI
mod.fit.Fun <- function(x, dfn, dfe) {
  modc <- x %>% dplyr::mutate(F.val = `Chisq diff` / `Df diff`)
  F.val <- modc[2,8]
  modc <- modc %>% dplyr::mutate(eta2g = (F.val*dfn)/((F.val*dfn) + dfe))
  lims <- MBESS::conf.limits.ncf(F.value = F.val, conf.level = .90, df.1 <- dfn, df.2 <- dfe)
  Lower.lim <- lims$Lower.Limit/(lims$Lower.Limit + df.1 + df.2 + 1)
  Upper.lim <- lims$Upper.Limit/(lims$Upper.Limit + df.1 + df.2 + 1)
  modc <- modc %>% dplyr::mutate(lower.CI = Lower.lim) %>% dplyr::mutate(upper.CI = Upper.lim)
  return(modc)
}

# Function to calculate Cohen's w
cohenw.Fun <- function(x, n) {
  modc <- x %>% dplyr::mutate(w = sqrt(`Chisq diff` / (n*`Df diff`)))
  return(modc)
}


# Function to calculate degrees of freedom error within based on number of valid observations and levels of fixed factor
# data = data frame of all individual items 
dfe.item.fun <- function(data, n, k) {
  # Number of non-NA observation count
  obs <- sum(!is.na(data))
  # Calc of df error
  dfe <- obs - (n - 1) - (k - 1)
  return(dfe)
}

# Function to get F, eta2, and 90% CI (multiple comparison rows)
mod.fit.Fun2 <- function(x, dfn, dfe) {
  modc <- x %>% dplyr::mutate(F.val = `Chisq diff` / `Df diff`)
  modc <- modc %>% dplyr::mutate(eta2g = (F.val*dfn)/((F.val*dfn) + dfe))
  F.val.1 <- modc[2,8]
  F.val.2 <- modc[3,8]
  lims.1 <- MBESS::conf.limits.ncf(F.value = F.val.1, conf.level = .90, df.1 <- dfn, df.2 <- dfe)
  Lower.lim.1 <- lims.1$Lower.Limit/(lims.1$Lower.Limit + dfn + dfe + 1)
  Upper.lim.1 <- lims.1$Upper.Limit/(lims.1$Upper.Limit + dfn + dfe + 1)
  lims.2 <- MBESS::conf.limits.ncf(F.value = F.val.2, conf.level = .90, df.1 <- dfn, df.2 <- dfe)
  Lower.lim.2 <- lims.2$Lower.Limit/(lims.2$Lower.Limit + dfn + dfe + 1)
  Upper.lim.2 <- lims.2$Upper.Limit/(lims.2$Upper.Limit + dfn + dfe + 1)
  modc <- modc %>% dplyr::mutate(lower.CI = NA) %>% dplyr::mutate(upper.CI = NA)
  modc$lower.CI[2] <- Lower.lim.1
  modc$upper.CI[2] <- Upper.lim.1
  modc$lower.CI[3] <- Lower.lim.2
  modc$upper.CI[3] <- Upper.lim.2
  return(modc)
}

# Function to get eta2, and 90% CI from RM-MLM output
mlm.fit.Fun <- function(x) {
  dfn <- x[2,1]
  dfe <- x[2,2]
  F.val <- x[2,3]
  modc <- x %>% dplyr::mutate(eta2g = (F.val*dfn)/((F.val*dfn) + dfe))
  lims <- MBESS::conf.limits.ncf(F.value = F.val, conf.level = .90, df.1 <- dfn, df.2 <- dfe)
  Lower.lim <- lims$Lower.Limit/(lims$Lower.Limit + df.1 + df.2 + 1)
  Upper.lim <- lims$Upper.Limit/(lims$Upper.Limit + df.1 + df.2 + 1)
  modc <- modc %>% dplyr::mutate(lower.CI = Lower.lim) %>% dplyr::mutate(upper.CI = Upper.lim)
  return(modc)
}

### Recode Variables -----------------------------------------------------------------------------------------------------------------------------------------------------------------

# Recoding schemes
rcs1 <- "23=1; 24=2; 25=3; 26=4; 27=5; else=NA" # First recoding scheme
rcs2 <- "67=0; 8=1; 9=2; 10=3; 11=4; 12=5; 13=6; 14=7; 15=8; 16=9; 17=10; else=NA" # Second recoding scheme
rcs3 <- "7=0; 6=1; else=NA" # Third recoding scheme
rcs4 <- "33=1; 34=2; 25=3; 26=4; 27=5; else=NA" # Fourth recoding scheme
rcs5 <- "18=1; 19=2; 20=3; 21=4; 22=5; else=NA" # Fifth recoding scheme
rcs6 <- "28=1; 29=2; 30=3; 31=4; 32=5; else=NA" # Sixth recoding scheme
rcs7 <- "35=1; 36=2; 37=3; 38=4; 39=5; else=NA" # Seventh recoding scheme
rcs8 <- "11=1; 12=1; 13=1; 14=1; 21=2; 22=2; 23=2; 24=1; 25=2; 26=2; 27=2; 28=9; 31=2; 32=1; 41=3; 42=3; 43=3; 44=4;
              45=4; 46=3; 51=5; 52=6; 61=7; 62=7; 63=8; 64=8; 65=8; 99=9; else=NA" # Eighth coding scheme
rcs9 <- "None=0; Sometimes=1; Alot=3; else=NA" # Ninth coding scheme
rcs10 <- "Never=0; Monthly or Less=1; 2-4 times a month=2; 2-3 times a week=3; 4 or more times a week=4; else=NA" # Tenth coding scheme
rcs11 <- "Never=0; Less than monthly=1; Monthly=2; Weekly=3; else=NA" # Eleventh coding scheme
rcs12 <- "1-2=1; 3-4=2; 5-6=3; 7-9=4; 10 or more=5; else=NA" # Twelfth coding scheme

# Recode adaptability items across time points
data <- ryouready::recode2(data, vars = c("Q64.1", "Q64.2", "Q64.3", "Q64.4", "Q64.5",
                                          "Q66.1", "Q66.2", "Q66.3", "Q66.4", "Q66.5",
                                          "Q67.1", "Q67.2", "Q67.3", "Q67.4", "Q67.5"), recodes = rcs1)

# Reverse score adaptability items
data[c("Q66.1", "Q66.2", "Q66.3", "Q66.4", "Q66.5")] <- 
  lapply(data[c("Q66.1", "Q66.2", "Q66.3", "Q66.4", "Q66.5")], function(x) 6 - x)

# Recode bad coping items across time points
data <- ryouready::recode2(data, vars = c("Q71.1", "Q71.2", "Q71.3", "Q71.4", "Q71.5",
                                          "Q76.1", "Q76.2", "Q76.3", "Q76.4", "Q76.5",
                                          "Q79.1", "Q79.2", "Q79.3", "Q79.4", "Q79.5"), recodes = rcs1)

# Reverse score passive coping items
data[c("Q71.1", "Q71.2", "Q71.3", "Q71.4", "Q71.5",
       "Q76.1", "Q76.2", "Q76.3", "Q76.4", "Q76.5",
       "Q79.1", "Q79.2", "Q79.3", "Q79.4", "Q79.5")] <- 
  lapply(data[c("Q71.1", "Q71.2", "Q71.3", "Q71.4", "Q71.5",
                "Q76.1", "Q76.2", "Q76.3", "Q76.4", "Q76.5",
                "Q79.1", "Q79.2", "Q79.3", "Q79.4", "Q79.5")], function(x) 6 - x)

# Recode catastrophizing items across time points
data <- ryouready::recode2(data, vars = c("Q175.1", "Q175.2", "Q175.3", "Q175.4", "Q175.5",
                                          "Q176.1", "Q176.2", "Q176.3", "Q176.4", "Q176.5",
                                          "Q54.1", "Q54.2", "Q54.3", "Q54.4", "Q54.5",
                                          "Q55.1", "Q55.2", "Q55.3", "Q55.4", "Q55.5",
                                          "Q56.1", "Q56.2", "Q56.3", "Q56.4", "Q56.5",
                                          "Q57.1", "Q57.2", "Q57.3", "Q57.4", "Q57.5",
                                          "Q58.1", "Q58.2", "Q58.3", "Q58.4", "Q58.5"), recodes = rcs1)

# Recode character items across time points
data <- ryouready::recode2(data, vars = c("Q30.1", "Q30.2", "Q30.3", "Q30.4", "Q30.5",
                                          "Q31.1", "Q31.2", "Q31.3", "Q31.4", "Q31.5",
                                          "Q32.1", "Q32.2", "Q32.3", "Q32.4", "Q32.5",
                                          "Q33.1", "Q33.2", "Q33.3", "Q33.4", "Q33.5",
                                          "Q34.1", "Q34.2", "Q34.3", "Q34.4", "Q34.5",
                                          "Q35.1", "Q35.2", "Q35.3", "Q35.4", "Q35.5",
                                          "Q36.1", "Q36.2", "Q36.3", "Q36.4", "Q36.5",
                                          "Q37.1", "Q37.2", "Q37.3", "Q37.4", "Q37.5",
                                          "Q38.1", "Q38.2", "Q38.3", "Q38.4", "Q38.5",
                                          "Q39.1", "Q39.2", "Q39.3", "Q39.4", "Q39.5",
                                          "Q40.1", "Q40.2", "Q40.3", "Q40.4", "Q40.5",
                                          "Q41.1", "Q41.2", "Q41.3", "Q41.4", "Q41.5",
                                          "Q42.1", "Q42.2", "Q42.3", "Q42.4", "Q42.5",
                                          "Q43.1", "Q43.2", "Q43.3", "Q43.4", "Q43.5",
                                          "Q44.1", "Q44.2", "Q44.3", "Q44.4", "Q44.5",
                                          "Q45.1", "Q45.2", "Q45.3", "Q45.4", "Q45.5",
                                          "Q46.1", "Q46.2", "Q46.3", "Q46.4", "Q46.5",
                                          "Q47.1", "Q47.2", "Q47.3", "Q47.4", "Q47.5",
                                          "Q48.1", "Q48.2", "Q48.3", "Q48.4", "Q48.5",
                                          "Q49.1", "Q49.2", "Q49.3", "Q49.4", "Q49.5",
                                          "Q50.1", "Q50.2", "Q50.3", "Q50.4", "Q50.5",
                                          "Q51.1", "Q51.2", "Q51.3", "Q51.4", "Q51.5",
                                          "Q52.1", "Q52.2", "Q52.3", "Q52.4", "Q52.5",
                                          "Q53.1", "Q53.2", "Q53.3", "Q53.4", "Q53.5"), recodes = rcs2)

# Recode active coping items across time points
data <- ryouready::recode2(data, vars = c("Q69.1", "Q69.2", "Q69.3", "Q69.4", "Q69.5",
                                          "Q70.1", "Q70.2", "Q70.3", "Q70.4", "Q70.5",
                                          "Q72.1", "Q72.2", "Q72.3", "Q72.4", "Q72.5",
                                          "Q74.1", "Q74.2", "Q74.3", "Q74.4", "Q74.5",
                                          "Q78.1", "Q78.2", "Q78.3", "Q78.4", "Q78.5"), recodes = rcs1)

# Recode life meaning items across time points
data <- ryouready::recode2(data, vars = c("Q82.1", "Q82.2", "Q82.3", "Q82.4", "Q82.5",
                                          "Q84.1", "Q84.2", "Q84.3", "Q84.4", "Q84.5",
                                          "Q86.1", "Q86.2", "Q86.3", "Q86.4", "Q86.5",
                                          "Q90.1", "Q90.2", "Q90.3", "Q90.4", "Q90.5",
                                          "Q92.1", "Q92.2", "Q92.3", "Q92.4", "Q92.5"), recodes = rcs4)

# Recode loneliness items across time points
data <- ryouready::recode2(data, vars = c("Q181.1", "Q181.2", "Q181.3", "Q181.4", "Q181.5",
                                          "Q185.1", "Q185.2", "Q185.3", "Q185.4", "Q185.5",
                                          "Q187.1", "Q187.2", "Q187.3", "Q187.4", "Q187.5"), recodes = rcs5)

# Reverse score loneliness items
data[c("Q185.1", "Q185.2", "Q185.3", "Q185.4", "Q185.5",
       "Q187.1", "Q187.2", "Q187.3", "Q187.4", "Q187.5")] <- 
  lapply(data[c("Q185.1", "Q185.2", "Q185.3", "Q185.4", "Q185.5",
                "Q187.1", "Q187.2", "Q187.3", "Q187.4", "Q187.5")], function(x) 6 - x)

# Recode negative affect items across time points
data <- ryouready::recode2(data, vars = c("Q156.1", "Q156.2", "Q156.3", "Q156.4", "Q156.5",
                                          "Q157.1", "Q157.2", "Q157.3", "Q157.4", "Q157.5",
                                          "Q160.1", "Q160.2", "Q160.3", "Q160.4", "Q160.5",
                                          "Q161.1", "Q161.2", "Q161.3", "Q161.4", "Q161.5",
                                          "Q164.1", "Q164.2", "Q164.3", "Q164.4", "Q164.5",
                                          "Q165.1", "Q165.2", "Q165.3", "Q165.4", "Q165.5",
                                          "Q167.1", "Q167.2", "Q167.3", "Q167.4", "Q167.5",
                                          "Q168.1", "Q168.2", "Q168.3", "Q168.4", "Q168.5",
                                          "Q169.1", "Q169.2", "Q169.3", "Q169.4", "Q169.5",
                                          "Q174.1", "Q174.2", "Q174.3", "Q174.4", "Q174.5",
                                          "Q177.1", "Q177.2", "Q177.3", "Q177.4", "Q177.5"), recodes = rcs5)

# Recode optimism items across time points
data <- ryouready::recode2(data, vars = c("Q93.1", "Q93.2", "Q93.3", "Q93.4", "Q93.5",
                                          "Q94.1", "Q94.2", "Q94.3", "Q94.4", "Q94.5",
                                          "Q97.1", "Q97.2", "Q97.3", "Q97.4", "Q97.5",
                                          "Q98.1", "Q98.2", "Q98.3", "Q98.4", "Q98.5"), recodes = rcs6)

# Reverse score optimism items
data[c("Q94.1", "Q94.2", "Q94.3", "Q94.4", "Q94.5",
       "Q97.1", "Q97.2", "Q97.3", "Q97.4", "Q97.5")] <- 
  lapply(data[c("Q94.1", "Q94.2", "Q94.3", "Q94.4", "Q94.5",
                "Q97.1", "Q97.2", "Q97.3", "Q97.4", "Q97.5")], function(x) 6 - x)

# Recode org trust items across time points
data <- ryouready::recode2(data, vars = c("Q113.1", "Q113.2", "Q113.3", "Q113.4", "Q113.5",
                                          "Q115.1", "Q115.2", "Q115.3", "Q115.4", "Q115.5",
                                          "Q117.1", "Q117.2", "Q117.3", "Q117.4", "Q117.5",
                                          "Q119.1", "Q119.2", "Q119.3", "Q119.4", "Q119.5",
                                          "Q124.1", "Q124.2", "Q124.3", "Q124.4", "Q124.5"), recodes = rcs7)

# Recode positive affect items across time points
data <- ryouready::recode2(data, vars = c("Q155.1", "Q155.2", "Q155.3", "Q155.4", "Q155.5",
                                          "Q158.1", "Q158.2", "Q158.3", "Q158.4", "Q158.5",
                                          "Q159.1", "Q159.2", "Q159.3", "Q159.4", "Q159.5",
                                          "Q162.1", "Q162.2", "Q162.3", "Q162.4", "Q162.5",
                                          "Q163.1", "Q163.2", "Q163.3", "Q163.4", "Q163.5",
                                          "Q166.1", "Q166.2", "Q166.3", "Q166.4", "Q166.5",
                                          "Q170.1", "Q170.2", "Q170.3", "Q170.4", "Q170.5",
                                          "Q171.1", "Q171.2", "Q171.3", "Q171.4", "Q171.5",
                                          "Q172.1", "Q172.2", "Q172.3", "Q172.4", "Q172.5",
                                          "Q173.1", "Q173.2", "Q173.3", "Q173.4", "Q173.5"), recodes = rcs5)

# Recode work engagement items across time points
data <- ryouready::recode2(data, vars = c("Q100.1", "Q100.2", "Q100.3", "Q100.4", "Q100.5",
                                          "Q103.1", "Q103.2", "Q103.3", "Q103.4", "Q103.5",
                                          "Q104.1", "Q104.2", "Q104.3", "Q104.4", "Q104.5",
                                          "Q106.1", "Q106.2", "Q106.3", "Q106.4", "Q106.5"), recodes = rcs1)

### Create scale variables by time occasion -------------------------------------------------------------------------------------------------------------------------------------------- 

# GAT: Adaptability; 3 items
data <- data %>% mutate(adapt.1.scale = rowMeans(data[c("Q64.1", "Q66.1", "Q67.1")], na.rm = TRUE))
data <- data %>% mutate(adapt.2.scale = rowMeans(data[c("Q64.2", "Q66.2", "Q67.2")], na.rm = TRUE))
data <- data %>% mutate(adapt.3.scale = rowMeans(data[c("Q64.3", "Q66.3", "Q67.3")], na.rm = TRUE))
data <- data %>% mutate(adapt.4.scale = rowMeans(data[c("Q64.4", "Q66.4", "Q67.4")], na.rm = TRUE))
data <- data %>% mutate(adapt.5.scale = rowMeans(data[c("Q64.5", "Q66.5", "Q67.5")], na.rm = TRUE))

# GAT: Passive Coping; 3 items
data <- data %>% mutate(pcope.1.scale = rowMeans(data[c("Q71.1", "Q76.1", "Q79.1")], na.rm = TRUE))
data <- data %>% mutate(pcope.2.scale = rowMeans(data[c("Q71.2", "Q76.2", "Q79.2")], na.rm = TRUE))
data <- data %>% mutate(pcope.3.scale = rowMeans(data[c("Q71.3", "Q76.3", "Q79.3")], na.rm = TRUE))
data <- data %>% mutate(pcope.4.scale = rowMeans(data[c("Q71.4", "Q76.4", "Q79.4")], na.rm = TRUE))
data <- data %>% mutate(pcope.5.scale = rowMeans(data[c("Q71.5", "Q76.5", "Q79.5")], na.rm = TRUE))

# GAT: Catastrophizing; 7 items
data <- data %>% mutate(catastro.1.scale = rowMeans(data[c("Q175.1", "Q176.1", "Q54.1", "Q55.1", "Q56.1", "Q57.1", "Q58.1")], na.rm = TRUE))
data <- data %>% mutate(catastro.2.scale = rowMeans(data[c("Q175.2", "Q176.2", "Q54.2", "Q55.2", "Q56.2", "Q57.2", "Q58.2")], na.rm = TRUE))
data <- data %>% mutate(catastro.3.scale = rowMeans(data[c("Q175.3", "Q176.3", "Q54.3", "Q55.3", "Q56.3", "Q57.3", "Q58.3")], na.rm = TRUE))
data <- data %>% mutate(catastro.4.scale = rowMeans(data[c("Q175.4", "Q176.4", "Q54.4", "Q55.4", "Q56.4", "Q57.4", "Q58.4")], na.rm = TRUE))
data <- data %>% mutate(catastro.5.scale = rowMeans(data[c("Q175.5", "Q176.5", "Q54.5", "Q55.5", "Q56.5", "Q57.5", "Q58.5")], na.rm = TRUE))

# GAT: Character; 24 items
data <- data %>% mutate(chr.1.scale = rowMeans(data[c("Q30.1", "Q31.1", "Q32.1", "Q33.1", "Q34.1", "Q35.1", "Q36.1", "Q37.1", "Q38.1", "Q39.1", 
                                                      "Q40.1", "Q41.1", "Q42.1", "Q43.1", "Q44.1", "Q45.1", "Q46.1", "Q47.1", "Q48.1", "Q49.1", 
                                                      "Q50.1", "Q51.1", "Q52.1", "Q53.1")], na.rm = TRUE))
data <- data %>% mutate(chr.2.scale = rowMeans(data[c("Q30.2", "Q31.2", "Q32.2", "Q33.2", "Q34.2", "Q35.2", "Q36.2", "Q37.2", "Q38.2", "Q39.2", 
                                                      "Q40.2", "Q41.2", "Q42.2", "Q43.2", "Q44.2", "Q45.2", "Q46.2", "Q47.2", "Q48.2", "Q49.2", 
                                                      "Q50.2", "Q51.2", "Q52.2", "Q53.2")], na.rm = TRUE))
data <- data %>% mutate(chr.3.scale = rowMeans(data[c("Q30.3", "Q31.3", "Q32.3", "Q33.3", "Q34.3", "Q35.3", "Q36.3", "Q37.3", "Q38.3", "Q39.3", 
                                                      "Q40.3", "Q41.3", "Q42.3", "Q43.3", "Q44.3", "Q45.3", "Q46.3", "Q47.3", "Q48.3", "Q49.3", 
                                                      "Q50.3", "Q51.3", "Q52.3", "Q53.3")], na.rm = TRUE))
data <- data %>% mutate(chr.4.scale = rowMeans(data[c("Q30.4", "Q31.4", "Q32.4", "Q33.4", "Q34.4", "Q35.4", "Q36.4", "Q37.4", "Q38.4", "Q39.4", 
                                                      "Q40.4", "Q41.4", "Q42.4", "Q43.4", "Q44.4", "Q45.4", "Q46.4", "Q47.4", "Q48.4", "Q49.4", 
                                                      "Q50.4", "Q51.4", "Q52.4", "Q53.4")], na.rm = TRUE))
data <- data %>% mutate(chr.5.scale = rowMeans(data[c("Q30.5", "Q31.5", "Q32.5", "Q33.5", "Q34.5", "Q35.5", "Q36.5", "Q37.5", "Q38.5", "Q39.5", 
                                                      "Q40.5", "Q41.5", "Q42.5", "Q43.5", "Q44.5", "Q45.5", "Q46.5", "Q47.5", "Q48.5", "Q49.5", 
                                                      "Q50.5", "Q51.5", "Q52.5", "Q53.5")], na.rm = TRUE))

# GAT: Depression; 10 items
data <- data %>% mutate(depress.1.scale = rowMeans(data[c("Q142.1", "Q143.1", "Q144.1", "Q145.1", "Q146.1", "Q147.1", "Q148.1", "Q149.1", "Q150.1", "Q151.1")], na.rm = TRUE))
data <- data %>% mutate(depress.2.scale = rowMeans(data[c("Q142.2", "Q143.2", "Q144.2", "Q145.2", "Q146.2", "Q147.2", "Q148.2", "Q149.2", "Q150.2", "Q151.2")], na.rm = TRUE))
data <- data %>% mutate(depress.3.scale = rowMeans(data[c("Q142.3", "Q143.3", "Q144.3", "Q145.3", "Q146.3", "Q147.3", "Q148.3", "Q149.3", "Q150.3", "Q151.3")], na.rm = TRUE))
data <- data %>% mutate(depress.4.scale = rowMeans(data[c("Q142.4", "Q143.4", "Q144.4", "Q145.4", "Q146.4", "Q147.4", "Q148.4", "Q149.4", "Q150.4", "Q151.4")], na.rm = TRUE))
data <- data %>% mutate(depress.5.scale = rowMeans(data[c("Q142.5", "Q143.5", "Q144.5", "Q145.5", "Q146.5", "Q147.5", "Q148.5", "Q149.5", "Q150.5", "Q151.5")], na.rm = TRUE))

# GAT: Active Coping; 5 items
data <- data %>% mutate(acope.1.scale = rowMeans(data[c("Q69.1", "Q70.1", "Q72.1", "Q74.1", "Q78.1")], na.rm = TRUE))
data <- data %>% mutate(acope.2.scale = rowMeans(data[c("Q69.2", "Q70.2", "Q72.2", "Q74.2", "Q78.2")], na.rm = TRUE))
data <- data %>% mutate(acope.3.scale = rowMeans(data[c("Q69.3", "Q70.3", "Q72.3", "Q74.3", "Q78.3")], na.rm = TRUE))
data <- data %>% mutate(acope.4.scale = rowMeans(data[c("Q69.4", "Q70.4", "Q72.4", "Q74.4", "Q78.4")], na.rm = TRUE))
data <- data %>% mutate(acope.5.scale = rowMeans(data[c("Q69.5", "Q70.5", "Q72.5", "Q74.5", "Q78.5")], na.rm = TRUE))

# GAT: Life Meaning; 5 items
data <- data %>% mutate(lifemean.1.scale = rowMeans(data[c("Q82.1", "Q84.1", "Q86.1", "Q90.1", "Q92.1")], na.rm = TRUE))
data <- data %>% mutate(lifemean.2.scale = rowMeans(data[c("Q82.2", "Q84.2", "Q86.2", "Q90.2", "Q92.2")], na.rm = TRUE))
data <- data %>% mutate(lifemean.3.scale = rowMeans(data[c("Q82.3", "Q84.3", "Q86.3", "Q90.3", "Q92.3")], na.rm = TRUE))
data <- data %>% mutate(lifemean.4.scale = rowMeans(data[c("Q82.4", "Q84.4", "Q86.4", "Q90.4", "Q92.4")], na.rm = TRUE))
data <- data %>% mutate(lifemean.5.scale = rowMeans(data[c("Q82.5", "Q84.5", "Q86.5", "Q90.5", "Q92.5")], na.rm = TRUE))

# GAT: Loneliness; 3 items
data <- data %>% mutate(lone.1.scale = rowMeans(data[c("Q181.1", "Q185.1", "Q187.1")], na.rm = TRUE))
data <- data %>% mutate(lone.2.scale = rowMeans(data[c("Q181.2", "Q185.2", "Q187.2")], na.rm = TRUE))
data <- data %>% mutate(lone.3.scale = rowMeans(data[c("Q181.3", "Q185.3", "Q187.3")], na.rm = TRUE))
data <- data %>% mutate(lone.4.scale = rowMeans(data[c("Q181.4", "Q185.4", "Q187.4")], na.rm = TRUE))
data <- data %>% mutate(lone.5.scale = rowMeans(data[c("Q181.5", "Q185.5", "Q187.5")], na.rm = TRUE))

# GAT: Negative Affect; 11 items
data <- data %>% mutate(negaffect.1.scale = rowMeans(data[c("Q156.1", "Q157.1", "Q160.1", "Q161.1", "Q164.1", "Q165.1", "Q167.1", "Q168.1", "Q169.1", "Q174.1", "Q177.1")], na.rm = TRUE))
data <- data %>% mutate(negaffect.2.scale = rowMeans(data[c("Q156.2", "Q157.2", "Q160.2", "Q161.2", "Q164.2", "Q165.2", "Q167.2", "Q168.2", "Q169.2", "Q174.2", "Q177.2")], na.rm = TRUE))
data <- data %>% mutate(negaffect.3.scale = rowMeans(data[c("Q156.3", "Q157.3", "Q160.3", "Q161.3", "Q164.3", "Q165.3", "Q167.3", "Q168.3", "Q169.3", "Q174.3", "Q177.3")], na.rm = TRUE))
data <- data %>% mutate(negaffect.4.scale = rowMeans(data[c("Q156.4", "Q157.4", "Q160.4", "Q161.4", "Q164.4", "Q165.4", "Q167.4", "Q168.4", "Q169.4", "Q174.4", "Q177.4")], na.rm = TRUE))
data <- data %>% mutate(negaffect.5.scale = rowMeans(data[c("Q156.5", "Q157.5", "Q160.5", "Q161.5", "Q164.5", "Q165.5", "Q167.5", "Q168.5", "Q169.5", "Q174.5", "Q177.5")], na.rm = TRUE))

# GAT: Optimism; 4 items
data <- data %>% mutate(optimism.1.scale = rowMeans(data[c("Q93.1", "Q94.1", "Q97.1", "Q98.1")], na.rm = TRUE))
data <- data %>% mutate(optimism.2.scale = rowMeans(data[c("Q93.2", "Q94.2", "Q97.2", "Q98.2")], na.rm = TRUE))
data <- data %>% mutate(optimism.3.scale = rowMeans(data[c("Q93.3", "Q94.3", "Q97.3", "Q98.3")], na.rm = TRUE))
data <- data %>% mutate(optimism.4.scale = rowMeans(data[c("Q93.4", "Q94.4", "Q97.4", "Q98.4")], na.rm = TRUE))
data <- data %>% mutate(optimism.5.scale = rowMeans(data[c("Q93.5", "Q94.5", "Q97.5", "Q98.5")], na.rm = TRUE))

# GAT: Organizational Trust; 5 items
data <- data %>% mutate(orgtrust.1.scale = rowMeans(data[c("Q113.1", "Q115.1", "Q117.1", "Q119.1", "Q124.1")], na.rm = TRUE))
data <- data %>% mutate(orgtrust.2.scale = rowMeans(data[c("Q113.2", "Q115.2", "Q117.2", "Q119.2", "Q124.2")], na.rm = TRUE))
data <- data %>% mutate(orgtrust.3.scale = rowMeans(data[c("Q113.3", "Q115.3", "Q117.3", "Q119.3", "Q124.3")], na.rm = TRUE))
data <- data %>% mutate(orgtrust.4.scale = rowMeans(data[c("Q113.4", "Q115.4", "Q117.4", "Q119.4", "Q124.4")], na.rm = TRUE))
data <- data %>% mutate(orgtrust.5.scale = rowMeans(data[c("Q113.5", "Q115.5", "Q117.5", "Q119.5", "Q124.5")], na.rm = TRUE))

# GAT: Positive Affect; 10 items
data <- data %>% mutate(posaffect.1.scale = rowMeans(data[c("Q155.1", "Q158.1", "Q159.1", "Q162.1", "Q163.1", "Q166.1", "Q170.1", "Q171.1", "Q172.1", "Q173.1")], na.rm = TRUE))
data <- data %>% mutate(posaffect.2.scale = rowMeans(data[c("Q155.2", "Q158.2", "Q159.2", "Q162.2", "Q163.2", "Q166.2", "Q170.2", "Q171.2", "Q172.2", "Q173.2")], na.rm = TRUE))
data <- data %>% mutate(posaffect.3.scale = rowMeans(data[c("Q155.3", "Q158.3", "Q159.3", "Q162.3", "Q163.3", "Q166.3", "Q170.3", "Q171.3", "Q172.3", "Q173.3")], na.rm = TRUE))
data <- data %>% mutate(posaffect.4.scale = rowMeans(data[c("Q155.4", "Q158.4", "Q159.4", "Q162.4", "Q163.4", "Q166.4", "Q170.4", "Q171.4", "Q172.4", "Q173.4")], na.rm = TRUE))
data <- data %>% mutate(posaffect.5.scale = rowMeans(data[c("Q155.5", "Q158.5", "Q159.5", "Q162.5", "Q163.5", "Q166.5", "Q170.5", "Q171.5", "Q172.5", "Q173.5")], na.rm = TRUE))

# GAT: Work Engagement; 4 items
data <- data %>% mutate(wkengage.1.scale = rowMeans(data[c("Q100.1", "Q103.1", "Q104.1", "Q106.1")], na.rm = TRUE))
data <- data %>% mutate(wkengage.2.scale = rowMeans(data[c("Q100.2", "Q103.2", "Q104.2", "Q106.2")], na.rm = TRUE))
data <- data %>% mutate(wkengage.3.scale = rowMeans(data[c("Q100.3", "Q103.3", "Q104.3", "Q106.3")], na.rm = TRUE))
data <- data %>% mutate(wkengage.4.scale = rowMeans(data[c("Q100.4", "Q103.4", "Q104.4", "Q106.4")], na.rm = TRUE))
data <- data %>% mutate(wkengage.5.scale = rowMeans(data[c("Q100.5", "Q103.5", "Q104.5", "Q106.5")], na.rm = TRUE))

### Scale Reliability Analysis --------------------------------------------------------------------------------------------------------------------------------------------

## Create scale data groupings

# GAT: Adaptability; 3 items
adapt_rel.1 <- dplyr::select(data, Q64.1, Q66.1, Q67.1)
adapt_rel.2 <- dplyr::select(data, Q64.2, Q66.2, Q67.2)
adapt_rel.3 <- dplyr::select(data, Q64.3, Q66.3, Q67.3)
adapt_rel.4 <- dplyr::select(data, Q64.4, Q66.4, Q67.4)
adapt_rel.5 <- dplyr::select(data, Q64.5, Q66.5, Q67.5)
adapt_rel.x <- dplyr::select(data, Q64.1, Q66.1, Q67.1, Q64.2, Q66.2, Q67.2, Q64.3, Q66.3, Q67.3, Q64.4, Q66.4, Q67.4, Q64.5, Q66.5, Q67.5)

# GAT: Passive Coping; 3 items
pcope_rel.1 <- dplyr::select(data, Q71.1, Q76.1, Q79.1)
pcope_rel.2 <- dplyr::select(data, Q71.2, Q76.2, Q79.2)
pcope_rel.3 <- dplyr::select(data, Q71.3, Q76.3, Q79.3)
pcope_rel.4 <- dplyr::select(data, Q71.4, Q76.4, Q79.4)
pcope_rel.5 <- dplyr::select(data, Q71.5, Q76.5, Q79.5)
pcope_rel.x <- dplyr::select(data, Q71.1, Q76.1, Q79.1, Q71.2, Q76.2, Q79.2, Q71.3, Q76.3, Q79.3, Q71.4, Q76.4, Q79.4, Q71.5, Q76.5, Q79.5)

# GAT: Catastrophizing; 7 items
catastro_rel.1 <- dplyr::select(data, Q175.1, Q176.1, Q54.1, Q55.1, Q56.1, Q57.1, Q58.1)
catastro_rel.2 <- dplyr::select(data, Q175.2, Q176.2, Q54.2, Q55.2, Q56.2, Q57.2, Q58.2)
catastro_rel.3 <- dplyr::select(data, Q175.3, Q176.3, Q54.3, Q55.3, Q56.3, Q57.3, Q58.3)
catastro_rel.4 <- dplyr::select(data, Q175.4, Q176.4, Q54.4, Q55.4, Q56.4, Q57.4, Q58.4)
catastro_rel.5 <- dplyr::select(data, Q175.5, Q176.5, Q54.5, Q55.5, Q56.5, Q57.5, Q58.5)
catastro_rel.x <- dplyr::select(data, Q175.1, Q176.1, Q54.1, Q55.1, Q56.1, Q57.1, Q58.1, Q175.2, Q176.2, Q54.2, Q55.2, Q56.2, Q57.2, Q58.2, Q175.3, Q176.3, Q54.3, Q55.3, Q56.3, Q57.3, Q58.3, 
                                Q175.4, Q176.4, Q54.4, Q55.4, Q56.4, Q57.4, Q58.4, Q175.5, Q176.5, Q54.5, Q55.5, Q56.5, Q57.5, Q58.5)

# GAT: Character; 24 items
chr_rel.1 <- dplyr::select(data, Q30.1, Q31.1, Q32.1, Q33.1, Q34.1, Q35.1, Q36.1, Q37.1, Q38.1, Q39.1, Q40.1, Q41.1, 
                           Q42.1, Q43.1, Q44.1, Q45.1, Q46.1, Q47.1, Q48.1, Q49.1, Q50.1, Q51.1, Q52.1, Q53.1)
chr_rel.2 <- dplyr::select(data, Q30.2, Q31.2, Q32.2, Q33.2, Q34.2, Q35.2, Q36.2, Q37.2, Q38.2, Q39.2, Q40.2, Q41.2, 
                           Q42.2, Q43.2, Q44.2, Q45.2, Q46.2, Q47.2, Q48.2, Q49.2, Q50.2, Q51.2, Q52.2, Q53.2)
chr_rel.3 <- dplyr::select(data, Q30.3, Q31.3, Q32.3, Q33.3, Q34.3, Q35.3, Q36.3, Q37.3, Q38.3, Q39.3, Q40.3, Q41.3, 
                           Q42.3, Q43.3, Q44.3, Q45.3, Q46.3, Q47.3, Q48.3, Q49.3, Q50.3, Q51.3, Q52.3, Q53.3)
chr_rel.4 <- dplyr::select(data, Q30.4, Q31.4, Q32.4, Q33.4, Q34.4, Q35.4, Q36.4, Q37.4, Q38.4, Q39.4, Q40.4, Q41.4, 
                           Q42.4, Q43.4, Q44.4, Q45.4, Q46.4, Q47.4, Q48.4, Q49.4, Q50.4, Q51.4, Q52.4, Q53.4)
chr_rel.5 <- dplyr::select(data, Q30.5, Q31.5, Q32.5, Q33.5, Q34.5, Q35.5, Q36.5, Q37.5, Q38.5, Q39.5, Q40.5, Q41.5, 
                           Q42.5, Q43.5, Q44.5, Q45.5, Q46.5, Q47.5, Q48.5, Q49.5, Q50.5, Q51.5, Q52.5, Q53.5)
chr_rel.x <- dplyr::select(data, Q30.1, Q31.1, Q32.1, Q33.1, Q34.1, Q35.1, Q36.1, Q37.1, Q38.1, Q39.1, Q40.1, Q41.1, 
                           Q42.1, Q43.1, Q44.1, Q45.1, Q46.1, Q47.1, Q48.1, Q49.1, Q50.1, Q51.1, Q52.1, Q53.1, Q30.2, Q31.2, Q32.2, Q33.2, Q34.2, Q35.2, Q36.2, Q37.2, Q38.2, Q39.2, Q40.2, Q41.2, 
                           Q42.2, Q43.2, Q44.2, Q45.2, Q46.2, Q47.2, Q48.2, Q49.2, Q50.2, Q51.2, Q52.2, Q53.2, Q30.3, Q31.3, Q32.3, Q33.3, Q34.3, Q35.3, Q36.3, Q37.3, Q38.3, Q39.3, Q40.3, Q41.3, 
                           Q42.3, Q43.3, Q44.3, Q45.3, Q46.3, Q47.3, Q48.3, Q49.3, Q50.3, Q51.3, Q52.3, Q53.3, Q30.4, Q31.4, Q32.4, Q33.4, Q34.4, Q35.4, Q36.4, Q37.4, Q38.4, Q39.4, Q40.4, Q41.4, 
                           Q42.4, Q43.4, Q44.4, Q45.4, Q46.4, Q47.4, Q48.4, Q49.4, Q50.4, Q51.4, Q52.4, Q53.4, Q30.5, Q31.5, Q32.5, Q33.5, Q34.5, Q35.5, Q36.5, Q37.5, Q38.5, Q39.5, Q40.5, Q41.5, 
                           Q42.5, Q43.5, Q44.5, Q45.5, Q46.5, Q47.5, Q48.5, Q49.5, Q50.5, Q51.5, Q52.5, Q53.5)

# GAT: Depression; 10 items
depress_rel.1 <- dplyr::select(data, Q142.1, Q143.1, Q144.1, Q145.1, Q146.1, Q147.1, Q148.1, Q149.1, Q150.1, Q151.1)
depress_rel.2 <- dplyr::select(data, Q142.2, Q143.2, Q144.2, Q145.2, Q146.2, Q147.2, Q148.2, Q149.2, Q150.2, Q151.2)
depress_rel.3 <- dplyr::select(data, Q142.3, Q143.3, Q144.3, Q145.3, Q146.3, Q147.3, Q148.3, Q149.3, Q150.3, Q151.3)
depress_rel.4 <- dplyr::select(data, Q142.4, Q143.4, Q144.4, Q145.4, Q146.4, Q147.4, Q148.4, Q149.4, Q150.4, Q151.4)
depress_rel.5 <- dplyr::select(data, Q142.5, Q143.5, Q144.5, Q145.5, Q146.5, Q147.5, Q148.5, Q149.5, Q150.5, Q151.5)
depress_rel.x <- dplyr::select(data, Q142.1, Q143.1, Q144.1, Q145.1, Q146.1, Q147.1, Q148.1, Q149.1, Q150.1, Q151.1, Q142.2, Q143.2, Q144.2, Q145.2, Q146.2, Q147.2, Q148.2, Q149.2, Q150.2, Q151.2,
                               Q142.3, Q143.3, Q144.3, Q145.3, Q146.3, Q147.3, Q148.3, Q149.3, Q150.3, Q151.3, Q142.4, Q143.4, Q144.4, Q145.4, Q146.4, Q147.4, Q148.4, Q149.4, Q150.4, Q151.4,
                               Q142.5, Q143.5, Q144.5, Q145.5, Q146.5, Q147.5, Q148.5, Q149.5, Q150.5, Q151.5)

# GAT: Active Coping; 5 items
acope_rel.1 <- dplyr::select(data, Q69.1, Q70.1, Q72.1, Q74.1, Q78.1)
acope_rel.2 <- dplyr::select(data, Q69.2, Q70.2, Q72.2, Q74.2, Q78.2)
acope_rel.3 <- dplyr::select(data, Q69.3, Q70.3, Q72.3, Q74.3, Q78.3)
acope_rel.4 <- dplyr::select(data, Q69.4, Q70.4, Q72.4, Q74.4, Q78.4)
acope_rel.5 <- dplyr::select(data, Q69.5, Q70.5, Q72.5, Q74.5, Q78.5)
acope_rel.x <- dplyr::select(data, Q69.1, Q70.1, Q72.1, Q74.1, Q78.1, Q69.2, Q70.2, Q72.2, Q74.2, Q78.2, Q69.3, Q70.3, Q72.3, Q74.3, Q78.3, Q69.4, Q70.4, Q72.4, Q74.4, Q78.4, 
                             Q69.5, Q70.5, Q72.5, Q74.5, Q78.5)

# GAT: Life Meaning; 5 items
lifemean_rel.1 <- dplyr::select(data, Q82.1, Q84.1, Q86.1, Q90.1, Q92.1)
lifemean_rel.2 <- dplyr::select(data, Q82.2, Q84.2, Q86.2, Q90.2, Q92.2)
lifemean_rel.3 <- dplyr::select(data, Q82.3, Q84.3, Q86.3, Q90.3, Q92.3)
lifemean_rel.4 <- dplyr::select(data, Q82.4, Q84.4, Q86.4, Q90.4, Q92.4)
lifemean_rel.5 <- dplyr::select(data, Q82.5, Q84.5, Q86.5, Q90.5, Q92.5)
lifemean_rel.x <- dplyr::select(data, Q82.1, Q84.1, Q86.1, Q90.1, Q92.1, Q82.2, Q84.2, Q86.2, Q90.2, Q92.2,  Q82.3, Q84.3, Q86.3, Q90.3, Q92.3, Q82.4, Q84.4, Q86.4, Q90.4, Q92.4,
                                Q82.5, Q84.5, Q86.5, Q90.5, Q92.5)

# GAT: Loneliness; 3 items
lone_rel.1 <- dplyr::select(data, Q181.1, Q185.1, Q187.1)
lone_rel.2 <- dplyr::select(data, Q181.2, Q185.2, Q187.2)
lone_rel.3 <- dplyr::select(data, Q181.3, Q185.3, Q187.3)
lone_rel.4 <- dplyr::select(data, Q181.4, Q185.4, Q187.4)
lone_rel.5 <- dplyr::select(data, Q181.5, Q185.5, Q187.5)
lone_rel.x <- dplyr::select(data, Q181.1, Q185.1, Q187.1, Q181.2, Q185.2, Q187.2, Q181.3, Q185.3, Q187.3, Q181.4, Q185.4, Q187.4, Q181.5, Q185.5, Q187.5)

# GAT: Negative Affect; 11 items
negaffect_rel.1 <- dplyr::select(data, Q156.1, Q157.1, Q160.1, Q161.1, Q164.1, Q165.1, Q167.1, Q168.1, Q169.1, Q174.1, Q177.1)
negaffect_rel.2 <- dplyr::select(data, Q156.2, Q157.2, Q160.2, Q161.2, Q164.2, Q165.2, Q167.2, Q168.2, Q169.2, Q174.2, Q177.2)
negaffect_rel.3 <- dplyr::select(data, Q156.3, Q157.3, Q160.3, Q161.3, Q164.3, Q165.3, Q167.3, Q168.3, Q169.3, Q174.3, Q177.3)
negaffect_rel.4 <- dplyr::select(data, Q156.4, Q157.4, Q160.4, Q161.4, Q164.4, Q165.4, Q167.4, Q168.4, Q169.4, Q174.4, Q177.4)
negaffect_rel.5 <- dplyr::select(data, Q156.5, Q157.5, Q160.5, Q161.5, Q164.5, Q165.5, Q167.5, Q168.5, Q169.5, Q174.5, Q177.5)
negaffect_rel.x <- dplyr::select(data, Q156.1, Q157.1, Q160.1, Q161.1, Q164.1, Q165.1, Q167.1, Q168.1, Q169.1, Q174.1, Q177.1, Q156.2, Q157.2, Q160.2, Q161.2, Q164.2, Q165.2, Q167.2, Q168.2, Q169.2, Q174.2, Q177.2,
                                 Q156.3, Q157.3, Q160.3, Q161.3, Q164.3, Q165.3, Q167.3, Q168.3, Q169.3, Q174.3, Q177.3, Q156.4, Q157.4, Q160.4, Q161.4, Q164.4, Q165.4, Q167.4, Q168.4, Q169.4, Q174.4, Q177.4,
                                 Q156.5, Q157.5, Q160.5, Q161.5, Q164.5, Q165.5, Q167.5, Q168.5, Q169.5, Q174.5, Q177.5)

# GAT: Optimism; 4 items
optimism_rel.1 <- dplyr::select(data, Q93.1, Q94.1, Q97.1, Q98.1)
optimism_rel.2 <- dplyr::select(data, Q93.2, Q94.2, Q97.2, Q98.2)
optimism_rel.3 <- dplyr::select(data, Q93.3, Q94.3, Q97.3, Q98.3)
optimism_rel.4 <- dplyr::select(data, Q93.4, Q94.4, Q97.4, Q98.4)
optimism_rel.5 <- dplyr::select(data, Q93.5, Q94.5, Q97.5, Q98.5)
optimism_rel.x <- dplyr::select(data,  Q93.1, Q94.1, Q97.1, Q98.1, Q93.2, Q94.2, Q97.2, Q98.2, Q93.3, Q94.3, Q97.3, Q98.3, Q93.4, Q94.4, Q97.4, Q98.4, Q93.5, Q94.5, Q97.5, Q98.5)

# GAT: Organizational Trust; 5 items
orgtrust_rel.1 <- dplyr::select(data, Q113.1, Q115.1, Q117.1, Q119.1, Q124.1)
orgtrust_rel.2 <- dplyr::select(data, Q113.2, Q115.2, Q117.2, Q119.2, Q124.2)
orgtrust_rel.3 <- dplyr::select(data, Q113.3, Q115.3, Q117.3, Q119.3, Q124.3)
orgtrust_rel.4 <- dplyr::select(data, Q113.4, Q115.4, Q117.4, Q119.4, Q124.4)
orgtrust_rel.5 <- dplyr::select(data, Q113.5, Q115.5, Q117.5, Q119.5, Q124.5)
orgtrust_rel.x <- dplyr::select(data, Q113.1, Q115.1, Q117.1, Q119.1, Q124.1, Q113.2, Q115.2, Q117.2, Q119.2, Q124.2, Q113.3, Q115.3, Q117.3, Q119.3, Q124.3, Q113.4, Q115.4, Q117.4, Q119.4, Q124.4,
                                Q113.5, Q115.5, Q117.5, Q119.5, Q124.5)

# GAT: Positive Affect; 10 items
posaffect_rel.1 <- dplyr::select(data, Q155.1, Q158.1, Q159.1, Q162.1, Q163.1, Q166.1, Q170.1, Q171.1, Q172.1, Q173.1)
posaffect_rel.2 <- dplyr::select(data, Q155.2, Q158.2, Q159.2, Q162.2, Q163.2, Q166.2, Q170.2, Q171.2, Q172.2, Q173.2)
posaffect_rel.3 <- dplyr::select(data, Q155.3, Q158.3, Q159.3, Q162.3, Q163.3, Q166.3, Q170.3, Q171.3, Q172.3, Q173.3)
posaffect_rel.4 <- dplyr::select(data, Q155.4, Q158.4, Q159.4, Q162.4, Q163.4, Q166.4, Q170.4, Q171.4, Q172.4, Q173.4)
posaffect_rel.5 <- dplyr::select(data, Q155.5, Q158.5, Q159.5, Q162.5, Q163.5, Q166.5, Q170.5, Q171.5, Q172.5, Q173.5)
posaffect_rel.x <- dplyr::select(data, Q155.1, Q158.1, Q159.1, Q162.1, Q163.1, Q166.1, Q170.1, Q171.1, Q172.1, Q173.1, Q155.2, Q158.2, Q159.2, Q162.2, Q163.2, Q166.2, Q170.2, Q171.2, Q172.2, Q173.2,
                                 Q155.3, Q158.3, Q159.3, Q162.3, Q163.3, Q166.3, Q170.3, Q171.3, Q172.3, Q173.3, Q155.4, Q158.4, Q159.4, Q162.4, Q163.4, Q166.4, Q170.4, Q171.4, Q172.4, Q173.4,
                                 Q155.5, Q158.5, Q159.5, Q162.5, Q163.5, Q166.5, Q170.5, Q171.5, Q172.5, Q173.5)

# GAT: Work Enagement; 4 items
wkengage_rel.1 <- dplyr::select(data, Q100.1, Q103.1, Q104.1, Q106.1)
wkengage_rel.2 <- dplyr::select(data, Q100.2, Q103.2, Q104.2, Q106.2)
wkengage_rel.3 <- dplyr::select(data, Q100.3, Q103.3, Q104.3, Q106.3)
wkengage_rel.4 <- dplyr::select(data, Q100.4, Q103.4, Q104.4, Q106.4)
wkengage_rel.5 <- dplyr::select(data, Q100.5, Q103.5, Q104.5, Q106.5)
wkengage_rel.x <- dplyr::select(data, Q100.1, Q103.1, Q104.1, Q106.1, Q100.2, Q103.2, Q104.2, Q106.2, Q100.3, Q103.3, Q104.3, Q106.3, Q100.4, Q103.4, Q104.4, Q106.4,
                                Q100.5, Q103.5, Q104.5, Q106.5)

## Calculate Cronbach's Alpha and Omega Total (four examples) 

# Table for Adaptability
r1 <- bindr(s1 = "adapt_rel.1", s2 = "adapt_rel.2", s3 = "adapt_rel.3", s4 = "adapt_rel.4", s5 = "adapt_rel.5")

# Table for Passive Coping
r2 <- bindr(s1 = "pcope_rel.1", s2 = "pcope_rel.2", s3 = "pcope_rel.3", s4 = "pcope_rel.4", s5 = "pcope_rel.5")

# Table for Catastrophizing
r3 <- bindr(s1 = "catastro_rel.1", s2 = "catastro_rel.2", s3 = "catastro_rel.3", s4 = "catastro_rel.4", s5 = "catastro_rel.5")

# Table for Character
r4 <- bindr(s1 = "chr_rel.1", s2 = "chr_rel.2", s3 = "chr_rel.3", s4 = "chr_rel.4", s5 = "chr_rel.5")

### Correlation Table --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

data.cor.0 <- data %>% dplyr::select(adapt.1.scale, adapt.2.scale, adapt.3.scale, adapt.4.scale, adapt.5.scale,
                                     pcope.1.scale, pcope.2.scale, pcope.3.scale, pcope.4.scale, pcope.5.scale, 
                                     catastro.1.scale, catastro.2.scale, catastro.3.scale, catastro.4.scale, catastro.5.scale, 
                                     chr.1.scale, chr.2.scale, chr.3.scale, chr.4.scale, chr.5.scale, depress.1.scale, depress.2.scale, 
                                     depress.3.scale, depress.4.scale, depress.5.scale, acope.1.scale, acope.2.scale, acope.3.scale, 
                                     acope.4.scale, acope.5.scale, lifemean.1.scale, lifemean.2.scale,  lifemean.3.scale, lifemean.4.scale, 
                                     lifemean.5.scale, lone.1.scale, lone.2.scale, lone.3.scale, lone.4.scale, lone.5.scale, negaffect.1.scale, 
                                     negaffect.2.scale, negaffect.3.scale, negaffect.4.scale, negaffect.5.scale, optimism.1.scale, optimism.2.scale,
                                     optimism.3.scale, optimism.4.scale, optimism.5.scale, orgtrust.1.scale, orgtrust.2.scale, orgtrust.3.scale, orgtrust.4.scale, 
                                     orgtrust.5.scale, posaffect.1.scale, posaffect.2.scale, posaffect.3.scale, posaffect.4.scale, posaffect.5.scale, 
                                     wkengage.1.scale, wkengage.2.scale, wkengage.3.scale, wkengage.4.scale, wkengage.5.scale)

# Create correlation matrix
t0.scale <- cor(data.cor.0, use = "pairwise.complete.obs", method = "pearson")                          
t0.scale.p <- cor.mtest(data.cor.0, conf.level = .95)  # Correlation p-value matrix


### Descriptive Statistics Analysis -------------------------------------------------------------------------------------------------------------------------------------------------------

data.desc <- data %>% subset(select = c(adapt.1.scale, adapt.2.scale, adapt.3.scale, adapt.4.scale, adapt.5.scale,
                                        pcope.1.scale, pcope.2.scale, pcope.3.scale, pcope.4.scale, pcope.5.scale, 
                                        catastro.1.scale, catastro.2.scale, catastro.3.scale, catastro.4.scale, catastro.5.scale, 
                                        chr.1.scale, chr.2.scale, chr.3.scale, chr.4.scale, chr.5.scale, depress.1.scale, depress.2.scale, 
                                        depress.3.scale, depress.4.scale, depress.5.scale, acope.1.scale, acope.2.scale, acope.3.scale, 
                                        acope.4.scale, acope.5.scale, lifemean.1.scale, lifemean.2.scale,  lifemean.3.scale, lifemean.4.scale, 
                                        lifemean.5.scale, lone.1.scale, lone.2.scale, lone.3.scale, lone.4.scale, lone.5.scale, negaffect.1.scale, 
                                        negaffect.2.scale, negaffect.3.scale, negaffect.4.scale, negaffect.5.scale, optimism.1.scale, optimism.2.scale,
                                        optimism.3.scale, optimism.4.scale, optimism.5.scale, orgtrust.1.scale, orgtrust.2.scale, orgtrust.3.scale, orgtrust.4.scale, 
                                        orgtrust.5.scale, posaffect.1.scale, posaffect.2.scale, posaffect.3.scale, posaffect.4.scale, posaffect.5.scale, 
                                        wkengage.1.scale, wkengage.2.scale, wkengage.3.scale, wkengage.4.scale, wkengage.5.scale))

# Descriptive stats table for scale variables in data frame
dat <- pastecs::stat.desc(data.desc, norm = FALSE, p = .95)

### Repeated Measures ANOVA (RM-ANOVA) Analysis (four examples) ------------------------------------------------------------------------------------------------------

# Adaptability
rA1 <- rANOVA.sum(x1 = data$adapt.1.scale, x2 = data$adapt.2.scale, x3 = data$adapt.3.scale, x4 = data$adapt.4.scale, x5 = data$adapt.5.scale)
rA1 <- rA1 %>% dplyr::mutate(var = "Adaptability") %>% dplyr::mutate(model = c("intercept", "trial")) %>% dplyr::select(var, model, dplyr::everything())

# Passive Coping
rA2 <- rANOVA.sum(x1 = data$pcope.1.scale, x2 = data$pcope.2.scale, x3 = data$pcope.3.scale, x4 = data$pcope.4.scale, x5 = data$pcope.5.scale)
rA2 <- rA2 %>% dplyr::mutate(var = "Passive Coping") %>% dplyr::mutate(model = c("intercept", "trial")) %>% dplyr::select(var, model, dplyr::everything())

# Catastrophizing
rA3 <- rANOVA.sum(x1 = data$catastro.1.scale, x2 = data$catastro.2.scale, x3 = data$catastro.3.scale, x4 = data$catastro.4.scale, x5 = data$catastro.5.scale)
rA3 <- rA3 %>% dplyr::mutate(var = "Catastrophizing") %>% dplyr::mutate(model = c("intercept", "trial")) %>% dplyr::select(var, model, dplyr::everything())

# Character
rA4 <- rANOVA.sum(x1 = data$chr.1.scale, x2 = data$chr.2.scale, x3 = data$chr.3.scale, x4 = data$chr.4.scale, x5 = data$chr.5.scale)
rA4 <- rA4 %>% dplyr::mutate(var = "Character") %>% dplyr::mutate(model = c("intercept", "trial")) %>% dplyr::select(var, model, dplyr::everything())


### Repeated Measures SEM (RM-SEM) Analysis (two examples) ---------------------------------------------------------------------------------------------------------------

## Adaptability
sem1.n <- syntax.null.CorComp(data = data, t1 = 'adapt.1.scale', t2 = 'adapt.2.scale', t3 = 'adapt.3.scale', t4 = 'adapt.4.scale', t5 = 'adapt.5.scale')
sem1.f <- syntax.full.CorComp(data = data, t1 = 'adapt.1.scale', t2 = 'adapt.2.scale', t3 = 'adapt.3.scale', t4 = 'adapt.4.scale', t5 = 'adapt.5.scale')
fit1.n <- lavaan::sem(sem1.n, data = data, meanstructure = FALSE, estimator = "ML", missing = "fiml")
fit1.f <- lavaan::sem(sem1.f, data = data, meanstructure = FALSE, estimator = "ML", missing = "fiml")
comp1 <- round(cbind(full = inspect(fit1.f, 'fit.measures'), null = inspect(fit1.n, 'fit.measures')), 3) %>% data.frame() %>% tibble::rownames_to_column() %>% dplyr::mutate(var = "Adaptability")
# Compare models 
modc1 <- as.data.frame(lavaan::lavTestLRT(fit1.f, fit1.n))
listn.fun(t1 = data$adapt.1.scale, t2 = data$adapt.2.scale, t3 = data$adapt.3.scale, t4 = data$adapt.4.scale, t5 = data$adapt.5.scale) #  n = 95277
dfe.fun(t1 = data$adapt.1.scale, t2 = data$adapt.2.scale, t3 = data$adapt.3.scale, t4 = data$adapt.4.scale, t5 = data$adapt.5.scale, n = 95277, k = 5) # 177378
mod1 <- tibble::as_tibble(mod.fit.Fun(x = modc1, dfn = 4, dfe = 177378)) %>% dplyr::mutate(name = "dfe_Num_Obs") %>% dplyr::mutate(name = "Adaptability") # df observations eta2
mod1b <- tibble::as_tibble(mod.fit.Fun(x = modc1, dfn = 4, dfe = 95277))  %>% dplyr::mutate(name = "dfe_Num_n") %>% dplyr::mutate(name = "Adaptability") # df sample size eta2
mod1c <- tibble::as_tibble(cohenw.Fun(x = modc1, n = 95277)) %>% dplyr::mutate(name = "Adaptability") # cohen's w

## Passive Coping
sem2.n <- syntax.null.CorComp(data = data, t1 = 'pcope.1.scale', t2 = 'pcope.2.scale', t3 = 'pcope.3.scale', t4 = 'pcope.4.scale', t5 = 'pcope.5.scale')
sem2.f <- syntax.full.CorComp(data = data, t1 = 'pcope.1.scale', t2 = 'pcope.2.scale', t3 = 'pcope.3.scale', t4 = 'pcope.4.scale', t5 = 'pcope.5.scale')
fit2.n <- lavaan::sem(sem2.n, data = data, meanstructure = FALSE, estimator = "ML", missing = "fiml")
fit2.f <- lavaan::sem(sem2.f, data = data, meanstructure = FALSE, estimator = "ML", missing = "fiml")
comp2 <- round(cbind(full = inspect(fit2.f, 'fit.measures'), null = inspect(fit2.n, 'fit.measures')), 3) %>% data.frame() %>% tibble::rownames_to_column() %>% dplyr::mutate(var = "Passive Coping")
# Compare models 
modc2 <- as.data.frame(lavaan::lavTestLRT(fit2.f, fit2.n))
listn.fun(t1 = data$pcope.1.scale, t2 = data$pcope.2.scale, t3 = data$pcope.3.scale, t4 = data$pcope.4.scale, t5 = data$pcope.5.scale) # n = 95277
dfe.fun(t1 = data$pcope.1.scale, t2 = data$pcope.2.scale, t3 = data$pcope.3.scale, t4 = data$pcope.4.scale, t5 = data$pcope.5.scale, n = 95277, k = 5) #177376
mod2 <- tibble::as_tibble(mod.fit.Fun(x = modc2, dfn = 4, dfe = 177376)) %>% dplyr::mutate(name = "dfe_Num_Obs") %>% dplyr::mutate(name = "Passive Coping")
mod2b <- tibble::as_tibble(mod.fit.Fun(x = modc2, dfn = 4, dfe = 95277))  %>% dplyr::mutate(name = "dfe_Num_n") %>% dplyr::mutate(name = "Passive Coping")
mod2c <- tibble::as_tibble(cohenw.Fun(x = modc2, n = 95277)) %>% dplyr::mutate(name = "Passive Coping")

### Measurement Invariance (MI) Analysis (one example) --------------------------------------------------------------------------------------------------------------------------

## Adaptability
# Configural invariance: Equality of factor structure (baseline model)
config1 <- '
# Creating first-order factors (no constraints)
  t1 =~ 1*Q64.1 + Q66.1 + Q67.1;
  t2 =~ 1*Q64.2 + Q66.2 + Q67.2;
  t3 =~ 1*Q64.3 + Q66.3 + Q67.3;
  t4 =~ 1*Q64.4 + Q66.4 + Q67.4;
  t5 =~ 1*Q64.5 + Q66.5 + Q67.5;
# Item intercepts (no constraints)
  Q64.1~1; Q66.1~1; Q67.1~1;
  Q64.2~1; Q66.2~1; Q67.2~1;
  Q64.3~1; Q66.3~1; Q67.3~1;
  Q64.4~1; Q66.4~1; Q67.4~1;
  Q64.5~1; Q66.5~1; Q67.5~1;
# Item residual covariances (no constraints)
  Q64.1 ~~ Q64.2 + Q64.3 + Q64.4 + Q64.5;
  Q66.1 ~~ Q66.2 + Q66.3 + Q66.4 + Q66.5;
  Q67.1 ~~ Q67.2 + Q67.3 + Q67.4 + Q67.5;
# First-order factor means (all fixed to zero)
  t1~0;
  t2~0;
  t3~0;
  t4~0;
  t5~0; 
# First-order factor covariance (no constraints)
  t1~~t2; t1~~t3; t1~~t4; t1~~t5;
  t2~~t3; t2~~t4; t2~~t5;
  t3~~t4; t3~~t5;
  t4~~t5;
'
# Weak (metric) invariance: Equality of factor loadings
weak1 <- '
# Creating first-order factors (equal)
  t1 =~ 1*Q64.1 + l1*Q66.1 + l2*Q67.1;
  t2 =~ 1*Q64.2 + l1*Q66.2 + l2*Q67.2;
  t3 =~ 1*Q64.3 + l1*Q66.3 + l2*Q67.3;
  t4 =~ 1*Q64.4 + l1*Q66.4 + l2*Q67.4;
  t5 =~ 1*Q64.5 + l1*Q66.5 + l2*Q67.5;
# Item intercepts (no constraints)
  Q64.1~1; Q66.1~1; Q67.1~1;
  Q64.2~1; Q66.2~1; Q67.2~1;
  Q64.3~1; Q66.3~1; Q67.3~1;
  Q64.4~1; Q66.4~1; Q67.4~1;
  Q64.5~1; Q66.5~1; Q67.5~1;
# Item residual covariances (no constraints)
  Q64.1 ~~ Q64.2 + Q64.3 + Q64.4 + Q64.5;
  Q66.1 ~~ Q66.2 + Q66.3 + Q66.4 + Q66.5;
  Q67.1 ~~ Q67.2 + Q67.3 + Q67.4 + Q67.5;
# First-order factor means (all fixed to zero)
  t1~0;
  t2~0;
  t3~0;
  t4~0;
  t5~0; 
# First-order factor covariance (no constraints)
  t1~~t2; t1~~t3; t1~~t4; t1~~t5;
  t2~~t3; t2~~t4; t2~~t5;
  t3~~t4; t3~~t5;
  t4~~t5;
'
# Strong (scalar) invariance: Equality of indicator intercepts
strong1 <- '
# Creating first-order factors (equal)
  t1 =~ 1*Q64.1 + l1*Q66.1 + l2*Q67.1;
  t2 =~ 1*Q64.2 + l1*Q66.2 + l2*Q67.2;
  t3 =~ 1*Q64.3 + l1*Q66.3 + l2*Q67.3;
  t4 =~ 1*Q64.4 + l1*Q66.4 + l2*Q67.4;
  t5 =~ 1*Q64.5 + l1*Q66.5 + l2*Q67.5;
# Item intercepts (equal)
  Q64.1~i1*1; Q66.1~i2*1; Q67.1~i3*1;
  Q64.2~i1*1; Q66.2~i2*1; Q67.2~i3*1;
  Q64.3~i1*1; Q66.3~i2*1; Q67.3~i3*1;
  Q64.4~i1*1; Q66.4~i2*1; Q67.4~i3*1;
  Q64.5~i1*1; Q66.5~i2*1; Q67.5~i3*1;
# Item residual covariances (no constraints)
  Q64.1 ~~ Q64.2 + Q64.3 + Q64.4 + Q64.5;
  Q66.1 ~~ Q66.2 + Q66.3 + Q66.4 + Q66.5;
  Q67.1 ~~ Q67.2 + Q67.3 + Q67.4 + Q67.5;
# First-order factor means (all fixed to zero)
  t1~0*1;
  t2~1;
  t3~1;
  t4~1;
  t5~1; 
# First-order factor covariance (no constraints)
  t1~~t2; t1~~t3; t1~~t4; t1~~t5;
  t2~~t3; t2~~t4; t2~~t5;
  t3~~t4; t3~~t5;
  t4~~t5;
'

# Fit the models
mi.fit1a <- lavaan::cfa(config1, data = data, meanstructure = TRUE, estimator = "ML", missing = "fiml", mimic = "mplus")
mi.fit1b <- lavaan::cfa(weak1, data = data, meanstructure = TRUE, estimator = "ML", missing = "fiml", mimic = "mplus")
mi.fit1c <- lavaan::cfa(strong1, data = data, meanstructure = TRUE, estimator = "ML", missing = "fiml", mimic = "mplus")

# Compare the models
mi.modc1 <- as.data.frame(lavaan::lavTestLRT(mi.fit1a, mi.fit1b, mi.fit1c))
dfe.item.fun(adapt_rel.x, n = 95277, k = 5) # 722694
mi.mod1 <- tibble::as_tibble(mod.fit.Fun2(x = mi.modc1, dfn = 8, dfe = 722694)) %>% dplyr::mutate(name = "Adaptability")
mi.mod1b <- tibble::as_tibble(mod.fit.Fun2(x = mi.modc1, dfn = 8, dfe = 95277)) %>% dplyr::mutate(name = "Adaptability")
mi.mod1c <- tibble::as_tibble(cohenw.Fun(x = mi.modc1, n = 95277)) %>% dplyr::mutate(name = "Adaptability")
mi.comp1 <- round(cbind(configural = inspect(mi.fit1a, 'fit.measures'), weak = inspect(mi.fit1b, 'fit.measures'), strong = inspect(mi.fit1c, 'fit.measures')), 3)

### Repeated Measures CFA (rCFA) Analysis (one example) ----------------------------------------------------------------------------------------------------------------------

## Adaptability
# Model where factor means are allowed to be freely estimated
fac.mean.full1 <- '
# Creating first-order factors (constrained to be equal)
  t1 =~ 1*Q64.1 + l1*Q66.1 + l2*Q67.1;
  t2 =~ 1*Q64.2 + l1*Q66.2 + l2*Q67.2;
  t3 =~ 1*Q64.3 + l1*Q66.3 + l2*Q67.3;
  t4 =~ 1*Q64.4 + l1*Q66.4 + l2*Q67.4;
  t5 =~ 1*Q64.5 + l1*Q66.5 + l2*Q67.5;
# Item intercepts (constrained to be equal)
  Q64.1~i1*1; Q66.1~i2*1; Q67.1~i3*1;
  Q64.2~i1*1; Q66.2~i2*1; Q67.2~i3*1;
  Q64.3~i1*1; Q66.3~i2*1; Q67.3~i3*1;
  Q64.4~i1*1; Q66.4~i2*1; Q67.4~i3*1;
  Q64.5~i1*1; Q66.5~i2*1; Q67.5~i3*1;
# Item residual covariances (no constraints)
  Q64.1 ~~ Q64.2 + Q64.3 + Q64.4 + Q64.5;
  Q66.1 ~~ Q66.2 + Q66.3 + Q66.4 + Q66.5;
  Q67.1 ~~ Q67.2 + Q67.3 + Q67.4 + Q67.5;
# First-order factor means (freely estimated)
  t1~0*1;
  t2~1;
  t3~1;
  t4~1;
  t5~1; 
# First-order factor variance (constrained to be equal)
  t1 ~~ v1*t1;
  t2 ~~ v1*t2;
  t3 ~~ v1*t3;
  t4 ~~ v1*t4;
  t5 ~~ v1*t5;
# First-order factor covariance (constrained to be equal)
  t1 ~~ c*t2 + c*t3 + c*t4 + c*t5;
  t2 ~~ c*t3 + c*t4 + c*t5;
  t3 ~~ c*t4 + c*t5;
  t4 ~~ c*t5;
'
# Model where factor means are constrained to be equal
fac.mean.null1 <- '
# Creating first-order factors (constrained to be equal)
  t1 =~ 1*Q64.1 + l1*Q66.1 + l2*Q67.1;
  t2 =~ 1*Q64.2 + l1*Q66.2 + l2*Q67.2;
  t3 =~ 1*Q64.3 + l1*Q66.3 + l2*Q67.3;
  t4 =~ 1*Q64.4 + l1*Q66.4 + l2*Q67.4;
  t5 =~ 1*Q64.5 + l1*Q66.5 + l2*Q67.5;
# Item intercepts (constrained to be equal)
  Q64.1~i1*1; Q66.1~i2*1; Q67.1~i3*1;
  Q64.2~i1*1; Q66.2~i2*1; Q67.2~i3*1;
  Q64.3~i1*1; Q66.3~i2*1; Q67.3~i3*1;
  Q64.4~i1*1; Q66.4~i2*1; Q67.4~i3*1;
  Q64.5~i1*1; Q66.5~i2*1; Q67.5~i3*1;
# Item residual covariances (no constraints)
  Q64.1 ~~ Q64.2 + Q64.3 + Q64.4 + Q64.5;
  Q66.1 ~~ Q66.2 + Q66.3 + Q66.4 + Q66.5;
  Q67.1 ~~ Q67.2 + Q67.3 + Q67.4 + Q67.5;
# First-order factor means (constrained to be equal)
  t1~0*1;
  t2~0*1;
  t3~0*1;
  t4~0*1;
  t5~0*1; 
# First-order factor variance (constrained to be equal)
  t1 ~~ v1*t1;
  t2 ~~ v1*t2;
  t3 ~~ v1*t3;
  t4 ~~ v1*t4;
  t5 ~~ v1*t5;
# First-order factor covariance (constrained to be equal)
  t1 ~~ c*t2 + c*t3 + c*t4 + c*t5;
  t2 ~~ c*t3 + c*t4 + c*t5;
  t3 ~~ c*t4 + c*t5;
  t4 ~~ c*t5;
'

rcfa.fit1a <- lavaan::cfa(fac.mean.full1, data = data, meanstructure = TRUE, estimator = "ML", missing = "fiml", mimic = "mplus")
sum1a <- lavaan::summary(rcfa.fit1a, standardized = TRUE, fit.measures = TRUE)

rcfa.fit1b <- lavaan::cfa(fac.mean.null1, data = data, meanstructure = TRUE, estimator = "ML", missing = "fiml", mimic = "mplus")
sum1b <- lavaan::summary(rcfa.fit1b, standardized = TRUE, fit.measures = TRUE)

rcfa.modc1 <- as.data.frame(lavaan::lavTestLRT(rcfa.fit1a, rcfa.fit1b))
dfe.item.fun(adapt_rel.x, n = 95277, k = 5) # 722694
rcfa.mod1 <- tibble::as_tibble(mod.fit.Fun(x = rcfa.modc1, dfn = 4, dfe = 722694)) %>% dplyr::mutate(name = "Adaptability")
rcfa.mod1b <- tibble::as_tibble(mod.fit.Fun(x = rcfa.modc1, dfn = 4, dfe = 95277)) %>% dplyr::mutate(name = "Adaptability")
rcfa.mod1c <- tibble::as_tibble(cohenw.Fun(x = rcfa.modc1, n = 95277)) %>% dplyr::mutate(name = "Adaptability")
rcfa.comp1 <- round(cbind(full = inspect(rcfa.fit1a, 'fit.measures'), null = inspect(rcfa.fit1b, 'fit.measures')), 3)

### Repeated Measures Multi-Level Model (RM-MLM) Treating Time as an Individually-Varying Continuous Metric (one example) ----------------------------

## Set all completion dates to be dates
data$COMPLETEDDATE.1 <- as.Date(data$COMPLETEDDATE.1, format = "%Y-%m-%d") 
data$COMPLETEDDATE.2 <- as.Date(data$COMPLETEDDATE.2, format = "%Y-%m-%d") 
data$COMPLETEDDATE.3 <- as.Date(data$COMPLETEDDATE.3, format = "%Y-%m-%d") 
data$COMPLETEDDATE.4 <- as.Date(data$COMPLETEDDATE.4, format = "%Y-%m-%d") 
data$COMPLETEDDATE.5 <- as.Date(data$COMPLETEDDATE.5, format = "%Y-%m-%d") 

## Calculate time for each measurement since first date for each participant
data$GAT.contutime1 <- lubridate::time_length(lubridate::as.period(lubridate::interval(lubridate::as_date(data$COMPLETEDDATE.1), lubridate::as_date(data$COMPLETEDDATE.1)), unit = "year"), unit = "year")
data$GAT.contutime2 <- lubridate::time_length(lubridate::as.period(lubridate::interval(lubridate::as_date(data$COMPLETEDDATE.1), lubridate::as_date(data$COMPLETEDDATE.2)), unit = "year"), unit = "year")
data$GAT.contutime3 <- lubridate::time_length(lubridate::as.period(lubridate::interval(lubridate::as_date(data$COMPLETEDDATE.1), lubridate::as_date(data$COMPLETEDDATE.3)), unit = "year"), unit = "year")
data$GAT.contutime4 <- lubridate::time_length(lubridate::as.period(lubridate::interval(lubridate::as_date(data$COMPLETEDDATE.1), lubridate::as_date(data$COMPLETEDDATE.4)), unit = "year"), unit = "year")
data$GAT.contutime5 <- lubridate::time_length(lubridate::as.period(lubridate::interval(lubridate::as_date(data$COMPLETEDDATE.1), lubridate::as_date(data$COMPLETEDDATE.5)), unit = "year"), unit = "year")

## Adaptability
# Select right dv and time vars
mlm01 <- data %>% dplyr::select(PID_PDE, adapt.1.scale, adapt.2.scale, adapt.3.scale, adapt.4.scale, adapt.5.scale, GAT.contutime1, GAT.contutime2, GAT.contutime3, GAT.contutime4, GAT.contutime5) %>% na.omit

# Transform data to long format with time-varying measures: time, GAT measure
mlm.L01 <- stats::reshape(mlm01, direction = 'long', 
                          varying = c("adapt.1.scale", "GAT.contutime1", "adapt.2.scale", "GAT.contutime2", "adapt.3.scale", "GAT.contutime3", "adapt.4.scale", "GAT.contutime4", "adapt.5.scale", "GAT.contutime5"), 
                          timevar = 'occasion',
                          times = c('1', '2', '3', '4', '5'),
                          v.names = c("y", "time"),
                          idvar = "PID_PDE") %>% tibble::remove_rownames() %>% dplyr::rename(time = y, y = time)

# Modeling
# random intercept only model
mlm.mod.ri01 <- nlme::lme(fixed = y ~ time, random = ~1|PID_PDE, method = "ML", cor = nlme::corCompSymm(form = ~time|PID_PDE), data = mlm.L01, na.action = na.omit, control = nlme::lmeControl(opt = "optim"))
# random intercept plus random slope model
mlm.mod.rs01 <- nlme::lme(fixed = y ~ time, random = ~1 + time|PID_PDE, method = "ML", cor = nlme::corCompSymm(form = ~time|PID_PDE), data = mlm.L01, na.action = na.omit, control = nlme::lmeControl(opt = "optim"))
summary(mlm.mod.ri01) # get model summary intercept only
summary(mlm.mod.rs01) # get model summary intercept plus random slope

mlm.fit.ri01 <- anova(mlm.mod.ri01) %>% as.data.frame() # get F-test for fixed effect
mlm.fit.rs01 <- anova(mlm.mod.rs01) %>% as.data.frame() # get F-test for fixed effect
mlm.fit.Fun(mlm.fit.ri01) # get generalized eta-squared effect size for fixed effect with 90% CIs
mlm.fit.Fun(mlm.fit.rs01) # get generalized eta-squared effect size for fixed effect with 90% CIs

# Calculate pseudo r-squared
# Calculates pseudo R2 effect size to test magnitude of change for addition of random slope to model
rcompanion::nagelkerke(fit = mlm.mod.rs01, null = mlm.mod.ri01) # pseudo.R.squared

### Random Intercept and Random Slope Analysis using SEM (one example) ------------------------------------------------------------------------------------------------

## Intercept model with slope fixed to zero
mlmsem.int <- '
int =~ 1*adapt.1.scale + 1*adapt.2.scale + 1*adapt.3.scale + 1*adapt.4.scale + 1*adapt.5.scale;
slope =~ 0*adapt.1.scale + 1* adapt.2.scale + 2* adapt.3.scale + 3* adapt.4.scale + 4*adapt.5.scale;
# intercepts (fixed effects);
int ~ 1;
slope ~ 1;
# random intercept;
int ~~ int;
# fixed slope;
slope ~~ 0*slope # no variance;
int ~~ 0*slope # no covariance;
# force same variance for all (compound symmetry);
adapt.1.scale ~~ v1* adapt.1.scale;
adapt.2.scale ~~ v1* adapt.2.scale;
adapt.3.scale ~~ v1* adapt.3.scale;
adapt.4.scale ~~ v1* adapt.4.scale;
adapt.5.scale ~~ v1* adapt.5.scale;
'
fit1.int <- lavaan::lavaan(mlmsem.int, data = data, meanstructure = FALSE, estimator = "ML", missing = "fiml")
lavaan::summary(fit1.int)

## Intercept model with random slope
mlmsem.slp <- '
int =~ 1*adapt.1.scale + 1*adapt.2.scale + 1*adapt.3.scale + 1*adapt.4.scale + 1*adapt.5.scale;
slope =~ 0*adapt.1.scale + 1* adapt.2.scale + 2* adapt.3.scale + 3* adapt.4.scale + 4*adapt.5.scale;
# intercepts (fixed effects);
int ~ 1;
slope ~ 1;
# random intercept with random slope;
int ~~ int;
slope ~~ slope;
int ~~ slope;
# force same variance for all (compound symmetry);
adapt.1.scale ~~ v1* adapt.1.scale;
adapt.2.scale ~~ v1* adapt.2.scale;
adapt.3.scale ~~ v1* adapt.3.scale;
adapt.4.scale ~~ v1* adapt.4.scale;
adapt.5.scale ~~ v1* adapt.5.scale;
'
fit1.slp <- lavaan::lavaan(mlmsem.slp, data = data, meanstructure = FALSE, estimator = "ML", missing = "fiml")
lavaan::summary(fit1.slp)

# compare nested models
lavaan::lavTestLRT(fit1.int, fit1.slp)

### Figure Visualizations  ------------------------------------------------------------------------------------------------

## Figure 1: Mean Plot by Time Occasion for both studies
# Load in data
GATc.mp <- read.csv("Mean_Plot_data.csv", header = TRUE)
# Re-order variables
GATc.mp$var <- factor(GATc.mp$var, levels = c("Adaptability", "Active Coping", "Passive Coping", "Character", "Catastrophizing", "Depression", "Loneliness", "Life Meaning",
                                              "Optimism", "Positive Affect", "Negative Affect", "Org. Trust", "Work Engagement"))
# Create dummy data set to modify y-axis range using the mean statement below
blank_data3 <- data.frame(var = c("Adaptability", "Adaptability", "Active Coping", "Active Coping", "Passive Coping", "Passive Coping", "Character", "Character", 
                                  "Catastrophizing", "Catastrophizing", "Depression", "Depression", "Loneliness", "Loneliness", "Life Meaning", "Life Meaning", 
                                  "Optimism", "Optimism", "Positive Affect", "Positive Affect", "Negative Affect", "Negative Affect", "Org. Trust", "Org. Trust", "Work Engagement", "Work Engagement"), 
                          vers = "Model", time = 1, mean = c(3.5, 4.5, 3.5, 4.5, 2, 3, 7.5, 8.5, 1.5, 2.5, 1, 2, 1.75, 2.75, 3.5, 4.5, 3.25, 4.25, 3.5, 4.5, 1.5, 2.5, 3.5, 4.5, 3.5, 4.5))
override.linetype <- c("solid", "dashed")
#dodge <- position_dodge(width = 0.2)
#png(width = 1900, height = 1709, res = 150)
ggplot2::ggplot(GATc.mp, aes(x = time, y = mean)) +
  geom_line(aes(color = vers, linetype = vers), size = .75) + 
  geom_point(aes(color = vers), size = 2.5) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper, color = vers), alpha = .75, width = 0.07, cex = 1) +
  geom_blank(data = blank_data3, aes(x = time, y = mean)) +
  labs(linetype = "", x = "Time Occassion", y = "Mean Value on Measure") +
  facet_wrap(~ var, strip.position = "top", nrow = 6, scales = "free_y") + # or "fixed" shows too much variation
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, color = "black", size = 10, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "black", size = 8, face = "italic"),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", colour = "black", size = 20, vjust = .1),
        axis.title.y = element_text(face = "bold", colour = "black", size = 20, vjust = 2),
        strip.text.x = element_text(hjust = 0.5, vjust = 0, angle = 0, face = "bold", color = "black", size = 14),
        strip.background = element_rect(fill = "#f2f1f1"),
        panel.grid.major.y = element_blank(), panel.background = element_rect(fill = "#fbfbfb"), panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  scale_colour_manual(values = cbbPalette3) +
  guides(colour = FALSE) +
  guides(linetype = guide_legend(override.aes = list(color = cbbPalette3, size = 1.4))) 
#dev.off()

## Figure 2: Forest-Facet Effect Size Plot Study 1 (GAT 1.0)

cbbPalette <- c("#E69F00", "#57b4e9", "#d55e00")
override.shape_c <- c(17, 15, 16)
#png(width = 1900, height = 1709, res = 150)
ggplot2::ggplot(data = GAT1.es, aes(x = group, y = ges, ymin = ci.lower, ymax = ci.upper, shape = group)) +
  geom_pointrange(aes(color = group, shape = group, stroke = 2.5)) +
  geom_hline(aes(yintercept = 0), color = "black", linetype ="dotted", size = .75) +
  geom_hline(aes(yintercept = 0.02), color = "black", linetype ="dashed", size = .75) +
  labs(color = "Model Type", x = "GAT Measure", y = expression(Effect ~ Size ~ `(`* eta[G]^2 *`)`)) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper, color = group), alpha = 1, width = 1.2, cex = 1) + 
  scale_shape_manual(values = c(17, 15, 16)) +
  facet_wrap(~var, strip.position = "left", nrow = 13, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.002, 0.075), breaks = scales::pretty_breaks(n = 7)) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, color = "black", size = 10, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "black", size = 8, face = "italic"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 20),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(face = "bold", colour = "black", size = 26, vjust = .1),
        axis.title.y = element_text(face = "bold", colour = "black", size = 26),
        strip.text.y = element_text(hjust = 0, vjust = .5, angle = 180, face = "bold", color = "black", size = 20),
        strip.background = element_rect(fill = "#f2f1f1"),
        panel.grid.major.y = element_blank(), panel.background = element_rect(fill = "#fbfbfb"),
        legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  scale_colour_manual(values = cbbPalette) +
  guides(shape = FALSE) +
  guides(colour = guide_legend(override.aes = list(shape = override.shape_c, size = 2))) +
  coord_flip()
#dev.off()

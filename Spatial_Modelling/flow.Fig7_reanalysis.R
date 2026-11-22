################################################################################
# flow.Fig7_reanalysis.R
# Sensitivity analysis for Figure 7 (transect land-cover vs. salinity),
# addressing reviewer concern about pseudoreplication / non-independence.
#
# Replaces: per-transect Pearson correlations -> one-way ANOVA on the coefficients
# With:     a hierarchical, temporally-autocorrelated model fit to the raw
#           transect x year panel, testing whether the land-cover:salinity
#           relationship differs among subregions (salinity x subregion term).
#
# Inputs (already created by Spatial_Work_Flow_cleanup.R Section 6):
#   summary_data_df : transect land-cover proportions (ID, time,
#                     proportion_Int1 = VZ, proportion_Int2 = DZ)
#   salinity_data   : site, year, annual_avg  (salinity_annual_avg_full.csv)
#   west / central / east : the salinity-site vectors defining each subregion
# Author: reanalysis for revision, 2026
################################################################################
library(nlme)

## ---- 0. Region assignment (matches compute_site_correlations_by_region) -----
##   transect ID > 61  -> west ; 37-61 -> central ; < 37 -> east
subof <- function(id) if (id > 61) "west" else if (id >= 37) "central" else "east"
yrs <- 1991:2024

## ---- 1. Subregional annual-mean salinity (one series per subregion) ---------
reg_mean <- function(sites) {
  s <- salinity_data[salinity_data$site %in% sites & salinity_data$year %in% yrs, ]
  aggregate(annual_avg ~ year, s, mean, na.rm = TRUE)
}
salreg <- do.call(rbind, list(
  cbind(reg_mean(west),    subregion = "west"),
  cbind(reg_mean(central), subregion = "central"),
  cbind(reg_mean(east),    subregion = "east")))
names(salreg)[2] <- "salinity"

## ---- 2. Build the transect x year panel ------------------------------------
panel <- do.call(rbind, lapply(1:100, function(id) {
  sr <- subof(id)
  td <- summary_data_df[summary_data_df$ID == id & summary_data_df$time %in% yrs,
                        c("time", "proportion_Int2", "proportion_Int1")]
  names(td) <- c("year", "DZ", "VZ")
  s  <- salreg[salreg$subregion == sr, c("year", "salinity")]
  m  <- merge(td, s, by = "year"); if (nrow(m) == 0) return(NULL)
  m$ID <- factor(id); m$subregion <- sr; m
}))
panel$subregion <- factor(panel$subregion, levels = c("central", "east", "west"))
panel <- panel[is.finite(panel$DZ) & is.finite(panel$VZ) & is.finite(panel$salinity), ]
# salinity standardized within subregion; response standardized within transect,
# so fixed slopes are interpretable on a correlation-like (partial-r) scale.
zsc <- function(x) if (sd(x, na.rm = TRUE) > 0) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) else x * 0
panel$sal_zs <- ave(panel$salinity, panel$subregion, FUN = zsc)

ctrl <- lmeControl(opt = "optim", maxIter = 300, msMaxIter = 300, returnObject = TRUE)

## ---- 3. Fit models & test the salinity x subregion interaction --------------
reanalyze <- function(resp) {
  d <- panel; d$y <- ave(d[[resp]], d$ID, FUN = zsc)
  # (A) mixed: transect random intercept + continuous AR1 within transect
  full  <- lme(y ~ sal_zs * subregion, random = ~1|ID,
               correlation = corCAR1(0.5, form = ~year|ID), data = d, method = "ML", control = ctrl)
  noint <- lme(y ~ sal_zs + subregion, random = ~1|ID,
               correlation = corCAR1(0.5, form = ~year|ID), data = d, method = "ML", control = ctrl)
  lrt <- anova(noint, full)
  # (B) random-slope check (transect-specific salinity slopes)
  full_s  <- try(update(full,  random = ~sal_zs|ID), silent = TRUE)
  noint_s <- try(update(noint, random = ~sal_zs|ID), silent = TRUE)
  lrt_s <- if (!inherits(full_s,"try-error")) anova(noint_s, full_s) else NULL
  # (C) GLS + AR1 (no random effect) cross-check
  gf <- gls(y ~ sal_zs * subregion, correlation = corCAR1(0.5, form = ~year|ID), data = d, method = "ML", control = ctrl)
  g0 <- gls(y ~ sal_zs + subregion, correlation = corCAR1(0.5, form = ~year|ID), data = d, method = "ML", control = ctrl)
  lrt_g <- anova(g0, gf)
  # per-subregion standardized association (~ autocorrelation-adjusted partial r)
  mm <- lme(y ~ subregion + subregion:sal_zs, random = ~1|ID,
            correlation = corCAR1(0.5, form = ~year|ID), data = d, method = "REML", control = ctrl)
  b <- fixef(mm); V <- vcov(mm); idx <- grep(":sal_zs$", names(b))
  assoc <- data.frame(subregion = sub("subregion","",sub(":sal_zs","",names(b)[idx])),
                      assoc = round(b[idx],3), se = round(sqrt(diag(V)[idx]),3),
                      p = signif(2*pnorm(-abs(b[idx]/sqrt(diag(V)[idx]))),3), row.names = NULL)
  cat("\n=========================== ", resp, " ===========================\n")
  cat("[A] Mixed (1|transect)+CAR1  interaction LRT: L=",round(lrt$L.Ratio[2],2),
      " df=",lrt$df[2]-lrt$df[1]," p=",signif(lrt$`p-value`[2],3),"\n",sep="")
  if(!is.null(lrt_s)) cat("[B] Random-slope           interaction LRT: L=",round(lrt_s$L.Ratio[2],2),
      " p=",signif(lrt_s$`p-value`[2],3),"\n",sep="")
  cat("[C] GLS+CAR1               interaction LRT: L=",round(lrt_g$L.Ratio[2],2),
      " p=",signif(lrt_g$`p-value`[2],3),"\n",sep="")
  cat("Autocorrelation-adjusted per-subregion association (partial-r scale):\n"); print(assoc)
  invisible(list(lrt=lrt, assoc=assoc, model=mm))
}
res_DZ <- reanalyze("DZ")
res_VZ <- reanalyze("VZ")

## ---- 4. Naive baseline for comparison (the criticized test) -----------------
naive <- function(propcol){
  rr <- do.call(rbind, lapply(1:100, function(id){
    sr<-subof(id); rm<-salreg[salreg$subregion==sr,]
    td<-summary_data_df[summary_data_df$ID==id & summary_data_df$time %in% rm$year,c("time",propcol)]
    names(td)<-c("year","v"); m<-merge(rm,td,by="year")
    if(nrow(m)>1 && sd(m$v)>0) data.frame(ID=id,group=sr,r=cor(m$salinity,m$v)) else NULL}))
  rr}
nd<-naive("proportion_Int2"); nv<-naive("proportion_Int1")
cat("\n--- Naive per-transect ANOVA-on-correlations (for reference) ---\n")
cat("DZ: F =",round(summary(aov(r~group,nd))[[1]]$`F value`[1],1),
    " p =",signif(summary(aov(r~group,nd))[[1]]$`Pr(>F)`[1],3),
    "| VZ: F =",round(summary(aov(r~group,nv))[[1]]$`F value`[1],1),
    " p =",signif(summary(aov(r~group,nv))[[1]]$`Pr(>F)`[1],3),"\n")

################################################################################
## PART B — Non-stationarity / spurious-correlation robustness (reviewer #2)
## Both salinity and land cover trend over time; correlation between two
## non-stationary series can be spurious. We (1) quantify trends, (2) re-test
## with a linear time term (detrended), and (3) first-difference the series.
################################################################################

## ---- B1. Trend diagnostics -------------------------------------------------
cat("\n--- Subregional salinity linear trend (PSU/yr) ---\n")
for (sr in c("central","east","west")) {
  s <- salreg[salreg$subregion == sr, ]; s <- s[order(s$year), ]
  f <- lm(salinity ~ year, s)
  cat(sprintf("  %-8s slope=%+.4f  p=%.3g\n", sr, coef(f)[2], summary(f)$coefficients[2,4]))
}
trend_share <- function(col) {
  pv <- sapply(split(panel, panel$ID), function(d) {
    d <- d[order(d$year), ]; if (nrow(d) < 5 || sd(d[[col]]) == 0) return(NA)
    summary(lm(d[[col]] ~ d$year))$coefficients[2,4] })
  mean(pv < 0.05, na.rm = TRUE) }
cat(sprintf("  Transects with significant trend: DZ %.0f%%, VZ %.0f%%\n",
            100*trend_share("DZ"), 100*trend_share("VZ")))

## ---- B2. Detrended model: salinity effect net of the common time trend -----
detrended <- function(resp) {
  d <- panel; d$y <- ave(d[[resp]], d$ID, FUN = zsc); d$yr_c <- d$year - mean(d$year)
  full  <- lme(y ~ sal_zs * subregion + yr_c, random = ~1|ID,
               correlation = corCAR1(0.5, form = ~year|ID), data = d, method = "ML", control = ctrl)
  noint <- lme(y ~ sal_zs + subregion + yr_c, random = ~1|ID,
               correlation = corCAR1(0.5, form = ~year|ID), data = d, method = "ML", control = ctrl)
  a <- anova(noint, full)
  cat(sprintf("  %s detrended: salinity x subregion LRT L=%.2f p=%.3g\n", resp, a$L.Ratio[2], a$`p-value`[2]))
}
cat("\n--- Detrended (year covariate) interaction tests ---\n"); detrended("DZ"); detrended("VZ")

## ---- B3. First differences: relate annual CHANGES (trend-free) --------------
dfp <- panel[order(panel$ID, panel$year), ]
for (v in c("DZ","VZ","salinity")) dfp[[paste0("d",v)]] <- ave(dfp[[v]], dfp$ID, FUN = function(x) c(NA, diff(x)))
dfp$dyr <- ave(dfp$year, dfp$ID, FUN = function(x) c(NA, diff(x)))
dfp <- dfp[is.finite(dfp$dsalinity) & dfp$dyr == 1, ]         # consecutive years only
dfp$dSal_z <- ave(dfp$dsalinity, dfp$subregion, FUN = zsc)
firstdiff <- function(resp) {
  d <- dfp; d$dy <- d[[paste0("d",resp)]]
  full  <- lme(dy ~ dSal_z * subregion, random = ~1|ID, data = d, method = "ML", control = ctrl)
  noint <- lme(dy ~ dSal_z + subregion, random = ~1|ID, data = d, method = "ML", control = ctrl)
  a <- anova(noint, full)
  m <- lme(dy ~ subregion + subregion:dSal_z, random = ~1|ID, data = d, method = "REML", control = ctrl)
  b <- fixef(m); V <- vcov(m); idx <- grep(":dSal_z$", names(b))
  as <- data.frame(subregion = sub("subregion","",sub(":dSal_z","",names(b)[idx])),
                   assoc = round(b[idx],3), p = signif(2*pnorm(-abs(b[idx]/sqrt(diag(V)[idx]))),3), row.names=NULL)
  cat(sprintf("\n  %s first-difference: interaction LRT L=%.2f p=%.3g\n", resp, a$L.Ratio[2], a$`p-value`[2]))
  print(as)
}
cat("\n--- First-difference (trend-free) interaction tests & per-subregion associations ---\n")
firstdiff("DZ"); firstdiff("VZ")
cat(sprintf("(first-difference n = %d transect-year changes)\n", nrow(dfp)))
## EOF Part B

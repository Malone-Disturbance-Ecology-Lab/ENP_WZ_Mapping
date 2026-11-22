# Region-adjusted VZ–salinity threshold analysis

This R Markdown workflow tests whether the VZ–salinity breakpoint remains after accounting for regional baseline differences and region-specific temporal trends.

The analysis assumes that `mixed_data` has already been created in the main workflow and contains:

- `proportion_Int1`
- `salinity`
- `region`
- `year`
- `ID`
- `RegionYear`
- `year_f`

---

## 1. Prepare data

```{r prepare-data}
vz_region <- mixed_data %>%
  filter(
    !is.na(proportion_Int1),
    !is.na(salinity),
    !is.na(region),
    !is.na(year),
    !is.na(ID)
  ) %>%
  mutate(
    ID = factor(ID),
    region = factor(region, levels = c("West", "Central", "East")),
    year_f = factor(year, levels = sort(unique(year))),
    year_c = year - mean(year, na.rm = TRUE),
    RegionYear = interaction(region, year, drop = TRUE)
  )

candidate_bp <- seq(
  quantile(vz_region$salinity, 0.10, na.rm = TRUE),
  quantile(vz_region$salinity, 0.90, na.rm = TRUE),
  by = 0.25
)

candidate_bp
```

---

## 2. Fit region-adjusted threshold models

The threshold model includes:

- a common piecewise salinity response,
- regional baseline differences,
- region-specific temporal trends,
- transect random intercepts,
- region-year random intercepts,
- AR(1) temporal autocorrelation within transects.

```{r fit-threshold-models}
models_region <- vector("list", length(candidate_bp))

aic_region <- data.frame(
  threshold = candidate_bp,
  AIC = NA_real_
)

for (i in seq_along(candidate_bp)) {

  bp <- candidate_bp[i]

  dat <- vz_region %>%
    mutate(
      salinity_bp = salinity - bp,
      post_threshold = pmax(0, salinity - bp)
    )

  models_region[[i]] <- glmmTMB(
    proportion_Int1 ~
      salinity_bp +
      post_threshold +
      region +
      year_c * region +
      (1 | ID) +
      (1 | RegionYear) +
      ar1(year_f + 0 | ID),
    data = dat,
    family = gaussian()
  )

  aic_region$AIC[i] <- AIC(models_region[[i]])

  cat(
    i, "/", length(candidate_bp),
    " BP =", round(bp, 2),
    " AIC =", round(aic_region$AIC[i], 2),
    "\n"
  )
}
```

---

## 3. Select the best breakpoint

```{r select-best-breakpoint}
aic_region <- aic_region %>%
  mutate(
    delta_AIC = AIC - min(AIC)
  ) %>%
  arrange(AIC)

best_bp <- aic_region$threshold[1]

best_i <- which.min(
  abs(candidate_bp - best_bp)
)

best_region_model <- models_region[[best_i]]

best_bp
aic_region
summary(best_region_model)
```

---

## 4. Fit the comparable linear model

This model has the same regional, temporal, random-effect, and AR(1) structure as the threshold model, but assumes a single linear salinity effect.

```{r linear-model}
linear_region_model <- glmmTMB(
  proportion_Int1 ~
    salinity +
    region +
    year_c * region +
    (1 | ID) +
    (1 | RegionYear) +
    ar1(year_f + 0 | ID),
  data = vz_region,
  family = gaussian()
)

AIC(
  linear_region_model,
  best_region_model
)

AIC(linear_region_model) -
  AIC(best_region_model)
```

---

## 5. Calculate below- and above-threshold slopes

Below the breakpoint:

\[
\text{slope}_{below} = \beta_{salinity}
\]

Above the breakpoint:

\[
\text{slope}_{above} =
\beta_{salinity} + \beta_{post-threshold}
\]

The `post_threshold` coefficient tests the change in slope at the breakpoint.

```{r calculate-slopes}
b <- fixef(best_region_model)$cond
V <- vcov(best_region_model)$cond

slope_test <- function(terms) {

  L <- setNames(
    rep(0, length(b)),
    names(b)
  )

  L[terms] <- 1

  est <- sum(L * b)

  se <- sqrt(
    as.numeric(
      t(L) %*% V %*% L
    )
  )

  z <- est / se

  data.frame(
    estimate = est * 100,
    SE = se * 100,
    lower95 = (est - 1.96 * se) * 100,
    upper95 = (est + 1.96 * se) * 100,
    z = z,
    p = 2 * pnorm(
      abs(z),
      lower.tail = FALSE
    )
  )
}

below <- slope_test(
  "salinity_bp"
)

above <- slope_test(
  c(
    "salinity_bp",
    "post_threshold"
  )
)

change <- slope_test(
  "post_threshold"
)

below
above
change
```

---

## 6. Save results

This writes the main results to one compact text file.

```{r save-results}
output <- paste0(
  "/Users/bz294/Documents/WhiteZone/Figure/Backup/",
  "VZ_threshold_region_adjusted_results.txt"
)

capture.output({

  cat(
    "REGION-ADJUSTED VZ-SALINITY THRESHOLD ANALYSIS\n\n"
  )

  cat(
    "Best breakpoint:",
    best_bp,
    "PSU\n\n"
  )

  cat("BELOW THRESHOLD\n")
  print(below)

  cat("\nABOVE THRESHOLD\n")
  print(above)

  cat("\nCHANGE IN SLOPE\n")
  print(change)

  cat("\nAIC COMPARISON\n")

  print(
    AIC(
      linear_region_model,
      best_region_model
    )
  )

  cat(
    "\nAIC improvement:",
    AIC(linear_region_model) -
      AIC(best_region_model),
    "\n"
  )

  cat("\nBEST MODEL SUMMARY\n")
  print(summary(best_region_model))

  cat("\nBREAKPOINT AIC PROFILE\n")
  print(aic_region)

}, file = output)

output
```

---

## 7. Optional: save the AIC profile as CSV

```{r save-aic-profile}
write.csv(
  aic_region,
  "/Users/bz294/Documents/WhiteZone/Figure/Backup/VZ_threshold_region_AIC_profile.csv",
  row.names = FALSE
)
```

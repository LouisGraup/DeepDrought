library(tidyverse)
library(lubridate)
library(mgcv)

df_Pfyn = read_csv("../../Data/Pfyn/PFY_daily_all_swp.csv")

# filter for treatments and variables
df_Pfyn = df_Pfyn |> filter(yr > 2015, yr < 2024, para=="SWP", treatment != "irrigation") |> 
  mutate(SWP = value) |> select(date, node, depth, SWP, treatment)

# monthly means
df_monthly <- df_Pfyn %>% mutate(month_date = floor_date(date, "month")) |> 
  group_by(node, treatment, depth, month_date) %>%
  summarise(SWP = mean(SWP, na.rm = TRUE), .groups = "drop") %>%
  mutate(year = year(month_date) + (month(month_date) - 0.5) / 12,
         month = month(month_date), node = factor(node), 
         treatment = factor(treatment, levels = c("control", "irrigation_stop")))

# fit GAM
fit_swp <- function(dat, soil_depth) {
  
  dat_depth <- dat %>%
    filter(depth == soil_depth) %>%
    droplevels()
  
  gam(SWP ~ treatment +
      s(year, by = treatment, k = 7) +
      s(month, bs = "cc", k = 8) +
      s(node, bs = "re"),
    data = dat_depth,
    method = "REML",
    knots = list(month = c(0.5, 12.5))
  )
}

mod_10 <- fit_swp(df_monthly, 10)
mod_80 <- fit_swp(df_monthly, 80)

summary(mod_10)
summary(mod_80)

# difference between control and irrigation-stop
predict_difference <- function(model, dat, soil_depth) {
  
  # Use the observed temporal range
  years <- seq(
    min(dat$year[dat$depth == soil_depth]),
    max(dat$year[dat$depth == soil_depth]),
    length.out = 500
  )
  
  # Use the same reference month and plot
  # for both treatment predictions
  newdata <- expand.grid(
    year = years,
    treatment = c("control", "irrigation_stop"),
    month = 7
  )
  
  newdata$node <- factor(
    levels(dat$node)[1],
    levels = levels(dat$node)
  )
  
  # Model matrix
  X <- predict(model, newdata = newdata, type = "lpmatrix")
  
  n <- length(years)
  
  # Stop minus Control
  X_diff <- X[(n + 1):(2 * n), , drop = FALSE] - X[1:n, , drop = FALSE]
  
  # Estimated contrast
  difference <- as.vector(X_diff %*% coef(model))
  
  # Standard errors of the contrast
  V <- vcov(model, unconditional = TRUE)
  
  se <- sqrt(rowSums((X_diff %*% V) * X_diff))
  
  tibble(depth = soil_depth, year = years, difference = difference,
    lower = difference - 1.96 * se, upper = difference + 1.96 * se)
}

df_diff <- bind_rows(
  predict_difference(mod_10, df_monthly, 10),
  predict_difference(mod_80, df_monthly, 80)
)

ggplot(df_diff, aes(x = year, y = difference)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue4", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ depth, ncol = 1) +
  labs(x = "", y = "SWP difference, Stop - Control (kPa)") +
  theme_bw()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
                   strip.text=element_text(size=12, face="bold"))


find_crossover <- function(dat) {
  
  dat <- dat %>% arrange(year)
  
  idx <- which(dat$difference[-nrow(dat)] * dat$difference[-1] < 0)
  
  if (length(idx) == 0) {
    return(numeric(0))
  }
  
  map_dbl(idx, function(i) {
    
    x1 <- dat$year[i]
    x2 <- dat$year[i + 1]
    
    y1 <- dat$difference[i]
    y2 <- dat$difference[i + 1]
    
    x1 - y1 * (x2 - x1) / (y2 - y1)
  })
}

df_diff %>% group_by(depth) %>%
  group_modify(~ tibble(crossover_year = find_crossover(.x)))
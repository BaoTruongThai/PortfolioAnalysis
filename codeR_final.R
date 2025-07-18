
# Section 1. Introduction
# Section 2. Stock Selection and Data Retrieval

library(quantmod)
library(ggplot2)
library(dplyr)
library(reshape2)
library(readrba)
library(knitr)
library(lubridate)   
library(Rsolnp) 
library(quadprog)

stocks <- c("CBA.AX", "BHP.AX", "CSL.AX", "WES.AX", "WOW.AX", 
            "TLS.AX", "STO.AX", "MQG.AX", "GMG.AX", "APA.AX")

start_date <- "2019-04-18"
end_date <- "2025-04-18"

stock_monthly_list <- list()
stock_return_list <- list()


# Step 1: Get data for 10 stocks
for (stock in stocks) {
  stock_data <- getSymbols(stock, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
  
  stock_monthly <- to.monthly(stock_data, indexAt = "lastof", drop.time = TRUE)
  stock_monthly_list[[stock]] <- stock_monthly
  
  stock_return <- periodReturn(stock_data, period = "monthly", type = "arithmetic")
  stock_return_list[[stock]] <- stock_return
}


# Step 2: Get data for All Ordinaries
aord_data <- getSymbols("^AORD", src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)

aord_monthly <- to.monthly(aord_data, indexAt = "lastof", drop.time = TRUE)

aord_monthly_returns <- periodReturn(aord_data, period = "monthly", type = "arithmetic")



# Step 3: 3-month BABs/NCDs
Rf <- read_rba_seriesid("FIRMMBAB90") %>% filter(date>= start_date)
Rf$Month <- format(Rf$date, "%m")
Rf$Year <- format(Rf$date, "%Y")
Rf$value <- Rf$value/1200
Rf <- Rf[c("date", "Year", "Month", "value")]


# Step 4: Combine data
library(xts)
library(dplyr)

# Convert Rf$date to Date type
Rf$date <- as.Date(Rf$date)

# Create end-of-month dates to match stock data
Rf_xts <- xts(Rf$value, order.by = Rf$date)

# If needed, aggregate to monthly (e.g., take the last observation each month)
Rf_monthly <- apply.monthly(Rf_xts, last)
colnames(Rf_monthly) <- "Risk_free"

all_returns <- do.call(merge, c(stock_return_list, list(Market = aord_monthly_returns, Risk_free = Rf_monthly)))


colnames(all_returns) <- c(stocks, "AORD", "Risk-free")


# Section 3. Calculating statistics
# a. avarage return
avg_return <-as.data.frame(colMeans(all_returns, na.rm = TRUE))
colnames(avg_return) <- "Average Return"

# b. standard deviation
std_dev <- as.data.frame(apply(all_returns, 2, sd, na.rm = TRUE))
colnames(std_dev) <- "Standard Deviation"

# c. covariance
cov_matrix <- as.data.frame(cov(all_returns, use = "pairwise.complete.obs"))
colnames(cov_matrix) <- paste0("Covariance with ", colnames(cov_matrix))

# c. merge avg_return and std_dev
stats <- merge(avg_return, std_dev, by = "row.names")


# Section 4. Plotting to see the performance (return and risk) of 10 stocks compared to All Ordinaries
# First ensure Row.names is treated as a factor
stats$Row.names <- as.factor(stats$Row.names)
# 4.1	Plot the average returns against the standard deviations on a graph
#windows(width = 10, height = 6)

ggplot(stats, aes(x = `Standard Deviation`, y = `Average Return`)) +
  geom_point(aes(color = Row.names), size = 4, alpha = 0.7) +
  scale_color_viridis_d() +
  
  geom_text(aes(label = Row.names), 
            vjust = -1.2, 
            size = 2) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  
  labs(title = "Risk-Return Profile of Australian Stocks",
       x = "Standard Deviation (Risk)",
       y = "Average Monthly Return") +
  
  theme_minimal() +
  theme(legend.position = "none")


# Section 5. Running the regression model for estimating Betas of 10 stocks
# Section 5: Estimate Beta - chosen stock by regression
# 5.1. Using CAPM to predicts the relationship between expected return of each stock and the expected market return
# 5.1 Calculate excess returns for stocks and market

#
excess_returns <- all_returns

# Remove Risk-free from excess_returns 
rf_values <- all_returns[, "Risk-free"]
excess_returns <- excess_returns[, !colnames(excess_returns) %in% "Risk-free"]

# Calculate excess returns for each stock and All ordinaries
for (stock in stocks) {
  excess_returns[, stock] <- all_returns[, stock] - rf_values
}
excess_returns[, "AORD"] <- all_returns[, "AORD"] - rf_values

# 5.2 Create a dataframe to store regression results
regression_results <- data.frame(
  Stock = stocks,
  Alpha = numeric(length(stocks)),
  Beta = numeric(length(stocks)),
  Alpha_p_value = numeric(length(stocks)),
  Beta_p_value = numeric(length(stocks)), # Cnowers
  R_squared = numeric(length(stocks)),
  stringsAsFactors = FALSE
)

# 5.3 Run regression for each stock
for (i in 1:length(stocks)) {
  stock <- stocks[i]
  
  # Ta:!o dataframe riC*ng cho ma;i ca; phia:?u Da; cha:!y ha;i quy
  reg_data <- na.omit(data.frame(
    y = as.numeric(excess_returns[, stock]),
    x = as.numeric(excess_returns[, "AORD"])
  ))
  
  # Da;i tC*n bia:?n Da; da; hia;u
  names(reg_data) <- c("Stock_Excess", "Market_Excess")
  
  # Cha:!y ha;i quy
  model <- lm(Stock_Excess ~ Market_Excess, data = reg_data)
  
  # LF0u ka:?t qua:#
  regression_results$Alpha[i] <- coef(model)[1]  # Intercept (alpha)
  regression_results$Beta[i] <- coef(model)[2]   # Slope (beta)
  regression_results$Alpha_p_value[i] <- summary(model)$coefficients[1, 4]  # p-value for alpha
  regression_results$R_squared[i] <- summary(model)$r.squared
}

# 5.4 Create table with results
regression_table <- regression_results %>%
  mutate(
    Alpha = round(Alpha, 4),
    Beta = round(Beta, 4),
    Alpha_p_value = round(Alpha_p_value, 4),
    Beta_p_value = round(Beta_p_value, 4), # C Nowers
    R_squared = round(R_squared, 4)*100,
    Null_Hypothesis = ifelse(Alpha_p_value < 0.05, "Reject", "Do not reject"),
    H0_Beta = ifelse(Beta_p_value < 0.05, "Reject", "Do not reject")
  )

# 5.5 Create a nice formatted table for the report
kable(regression_table, 
      caption = "CAPM Regression Results for 10 ASX Stocks",
      col.names = c("Stock", "Alpha", "Beta", "Alpha p-value", "Beta p-value", "R-squared", 
                    "H_0: alpha = 0", "H_0: beta = 0")) 

# 5.6 Visualization of beta values
ggplot(regression_table, aes(x = reorder(Stock, Beta), y = Beta)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Beta), vjust = -0.5, size = 3.5) +
  labs(title = "Beta Coefficients for ASX Stocks",
       subtitle = "Based on monthly returns (2019-2025)",
       x = "Stock",
       y = "Beta") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 5.7 Diagnostic plots for one example stock
example_stock <- "CBA.AX"
example_reg_data <- na.omit(data.frame(
  y = as.numeric(excess_returns[, example_stock]),
  x = as.numeric(excess_returns[, "AORD"])
))
names(example_reg_data) <- c("Stock_Excess", "Market_Excess")
example_model <- lm(Stock_Excess ~ Market_Excess, data = example_reg_data)

# Create diagnostic plots
par(mfrow = c(2, 2))
plot(example_model)


# Section 6. Compare with Betas from Financial Review

published_betas_5y<- c(CBA.AX = 0.85, BHP.AX = 0.75, CSL.AX=0.33, WES.AX =0.78, WOW.AX=0.28, 
                       TLS.AX=0.27, STO.AX=0.43, MQG.AX = 1.12, GMG.AX = 0.88, APA.AX = 0.24)

# Create a dataframe with published betas
beta_comparison <- data.frame(
  symbol = names(published_betas_5y),
  published_beta = published_betas_5y,
  estimated_beta = regression_results$Beta,
  residual = regression_results$Beta - published_betas_5y,
  stringsAsFactors = FALSE
)
mse <- mean(beta_comparison$residual^2)
# 6.2 comment on the difference

# 6.3 Plotting the comparison 
beta_comparison <- beta_comparison %>%
  arrange(estimated_beta)



# Dumbbell chart
ggplot(beta_comparison, aes(y = symbol)) +
  # 
  geom_segment(aes(x = published_beta, xend = estimated_beta,
                   y = symbol,       yend = symbol),
               color = "gray80", size = 2) +
  # Betas Published
  geom_point(aes(x = published_beta), color = "steelblue", size = 5) +
  # Betas Estimated
  geom_point(aes(x = estimated_beta), color = "firebrick", size = 5) +
  labs(
    title = "Beta comparison: Published vs Estimated",
    x     = "Beta",
    y     = "Stock",
    caption = "Blue: Published Beta b"
      ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    plot.title  = element_text(face = "bold", hjust = 5)
  )


# Section 7. Porfolio Analysis
# Section 7.1:	Calculate the average return and the covariance of the 10 stocks using the data of first 5 years. 
# 1. Define 5 years training window
train_start <- as.Date(start_date)                     # "2019-04-18"
train_end   <- train_start + years(5)                  # "2024-04-18"

# 2. Subset all_returns xts for that window and only the 10 stocks
#    all_returns has columns: stocks, "AORD", "Risk-free"
train_returns <- all_returns[paste0(train_start, "/", train_end), stocks]

# 3. Calculate average monthly return for each stock over first 5 years
#    (arithmetic mean of monthly returns)
avg_return_5yr <- colMeans(train_returns, na.rm = TRUE)

# 4. Calculate covariance matrix of those returns
cov_matrix_5yr <- cov(train_returns, use = "pairwise.complete.obs")

# 5. Prepare a summary table: Average Return & Standard Deviation
stats_5yr <- data.frame(
  Stock      = names(avg_return_5yr),
  AvgReturn  = round(avg_return_5yr, 4),
  StdDev     = round(sqrt(diag(cov_matrix_5yr)), 4),
  row.names  = NULL,
  stringsAsFactors = FALSE
)


# 6. (Optional) Print the summary nicely
kable(
  stats_5yr,
  caption = "Section 7.1 — Average Monthly Return & Std. Dev. (2019-04-18 to 2024-04-18)",
  col.names = c("Stock", "Avg Return", "Std Dev")
)

# 7. (Optional) Make the covariance matrix available for report
#    You can convert it to a data.frame or keep as matrix for optimization
cov_matrix_5yr_df <- as.data.frame(cov_matrix_5yr)

# 8. Visualization
# Giả sử bạn đã có stats_5yr:
#   Stock | AvgReturn | StdDev


ggplot(stats_5yr, aes(x = StdDev, y = AvgReturn)) +
  # điểm màu xanh dương, kích cỡ lớn để dễ nhìn
  geom_point(size = 5, color = "dodgerblue") +
  # ghi nhãn mã cổ phiếu, font size vừa phải, màu đỏ đậm
  geom_text(aes(label = Stock),
            size = 4, color = "darkred",
            vjust = -1.2, hjust = 0.5) +
  # trục phần trăm
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Risk–Return Profile (5-Year Training Period)",
    subtitle = "Average vs Std Dev (2019-04-18 – 2024-04-18)",
    x     = "Standard Deviation (Monthly Risk)",
    y     = "Average Monthly Return",
    caption = "Source: Yahoo Finance & RBA"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    plot.subtitle= element_text(size = 12),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12)
  )


# Assume cov_melt was created as:
    cov_melt <- melt(cov_matrix_5yr)
    colnames(cov_melt) <- c("Stock1", "Stock2", "Covariance")


ggplot(cov_melt, aes(x = Stock1, y = Stock2, fill = Covariance)) +
  # tiles with light grey borders
  geom_tile(color = "grey90") +
  # divergent palette: blues for low, whites at zero, reds for high
  scale_fill_gradient2(
    low      = "#4575b4",
mid      = "#ffffbf",
high     = "#d73027",
midpoint = 0,
name     = "Covariance"
  ) +
  # add numeric labels inside each cell
  geom_text(aes(label = sprintf("%.4f", Covariance)), size = 3, color = "black") +
  # titles and labels
  labs(
    title    = "Covariance Matrix Heatmap (5-Year Period)",
    subtitle = "Monthly Returns 2019-04-18 to 2024-04-18",
    x        = NULL,
    y        = NULL,
    caption  = "Source: Yahoo Finance & RBA"
  ) +
  # minimal theme with larger base font
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y      = element_text(size = 12),
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 10),
    panel.grid.major = element_blank()
  )

# assume avg_return_5yr (R_exp), cov_matrix_5yr (Σ), rf_raw, all_returns, stocks are already in memory
train_start <- as.Date("2019-04-18")
train_end   <- train_start + years(5)
test_start  <- train_end + 1
test_end    <- as.Date("2025-04-18")

# risk-free mean over training window
Rf_monthly <- Rf %>%
  filter(date >= train_start & date <= train_end) %>%
  pull(value) %>%
  mean(na.rm = TRUE)

# ——————————————
# 2. Max-Sharpe weights
# ——————————————
sharpe_fn <- function(w) {
  ret <- sum(w * avg_return_5yr)
  sd  <- sqrt(t(w) %*% cov_matrix_5yr %*% w)
  - (ret - Rf_monthly) / sd
}
n   <- length(avg_return_5yr)
w0  <- rep(1/n, n)
opt_sh <- solnp(
  pars  = w0,
  fun   = sharpe_fn,
  eqfun = function(w) sum(w) - 1, eqB = 0,
  LB    = rep(0, n), UB = rep(1, n)
)
w_sh <- opt_sh$pars
names(w_sh) <- stocks

# ——————————————
# 3. Min-Variance weights
# ——————————————
Dmat <- 2 * cov_matrix_5yr
dvec <- rep(0, n)
Amat <- cbind(rep(1, n), diag(n))
bvec <- c(1, rep(0, n))
sol_mv <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_mv  <- sol_mv$solution
names(w_mv) <- stocks

# ——————————————
# 4. Out-of-sample returns (6th year)
# ——————————————
test_ret <- all_returns[paste0(test_start, "/", test_end), stocks]
# portfolio returns
p_sh <- as.numeric(test_ret %*% w_sh)
p_mv <- as.numeric(test_ret %*% w_mv)

# build comparison table
comparison <- tibble(
  Asset     = c(stocks, "Max-Sharpe", "Min-Variance"),
  AvgReturn = c(
    colMeans(test_ret, na.rm = TRUE),
    mean(p_sh, na.rm = TRUE),
    mean(p_mv, na.rm = TRUE)
  ),
  StdDev    = c(
    apply(test_ret, 2, sd, na.rm = TRUE),
    sd(p_sh, na.rm = TRUE),
    sd(p_mv, na.rm = TRUE)
  )
)


# 1. Prepare data frames
weights_sh_df <- tibble(
  Stock  = names(w_sh),
  Weight = w_sh
)

weights_mv_df <- tibble(
  Stock  = names(w_mv),
  Weight = w_mv
)

# 2. Max-Sharpe weights chart
ggplot(weights_sh_df, aes(y = reorder(Stock, Weight), x = Weight)) +
  geom_col(fill = "firebrick", width = 0.6) +
  geom_text(aes(label = scales::percent(Weight, accuracy = 0.1)),
            hjust = -0.1, size = 3) +
  scale_x_continuous(labels = scales::percent_format(1),
                     expand = expansion(add = 0.1)) +
  labs(
    title = "Max-Sharpe Portfolio Weights",
    x     = "Weight (%)",
    y     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text  = element_text(size = 12)
  )

# 3. Min-Variance weights chart
ggplot(weights_mv_df, aes(y = reorder(Stock, Weight), x = Weight)) +
  geom_col(fill = "darkgreen", width = 0.6) +
  geom_text(aes(label = scales::percent(Weight, accuracy = 0.1)),
            hjust = -0.1, size = 3) +
  scale_x_continuous(labels = scales::percent_format(1),
                     expand = expansion(add = 0.1)) +
  labs(
    title = "Min-Variance Portfolio Weights",
    x     = "Weight (%)",
    y     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text  = element_text(size = 12)
  )


# 1. Define test window
test_start <- train_end + 1        # 2024-04-19
test_end   <- as.Date(end_date)    # 2025-04-18

# 2. Extract monthly returns for stocks
test_ret <- all_returns[paste0(test_start, "/", test_end), stocks]

# 3. Compute Avg & Std Dev for each stock
stock_perf <- tibble(
  Asset     = stocks,
  AvgReturn = colMeans(test_ret, na.rm = TRUE),
  StdDev    = apply(test_ret, 2, sd, na.rm = TRUE)
)

# 4. Compute portfolio returns & performance
portf_ret <- function(weights) {
  as.numeric(test_ret %*% weights)
}

perf_row <- function(name, weights) {
  r <- portf_ret(weights)
  tibble(
    Asset     = name,
    AvgReturn = mean(r, na.rm = TRUE),
    StdDev    = sd(r,   na.rm = TRUE)
  )
}

sharpe_perf <- perf_row("Max-Sharpe", w_sh)
mv_perf     <- perf_row("Min-Variance", w_mv)

# 5. Combine into one table
comparison <- bind_rows(stock_perf, sharpe_perf, mv_perf) %>%
  mutate(
    Type = case_when(
      Asset == "Max-Sharpe"    ~ "Max-Sharpe",
      Asset == "Min-Variance"  ~ "Min-Variance",
      TRUE                     ~ "Stock"
    )
  )

# 6. Print table
comparison %>%
  arrange(Type, desc(AvgReturn)) %>%
  kable(
    caption = "Out-of-Sample Returns & Risks (6th Year)",
    col.names = c("Asset", "Avg Monthly Return", "Std Dev", "Type")
  )

# 7. Scatter plot Risk–Return comparison
ggplot(comparison, aes(x = StdDev, y = AvgReturn, color = Type)) +
  geom_point(size = 4) +
  geom_text(aes(label = Asset), nudge_y = 0.002, size = 3) +
  scale_color_manual(
    values = c(
      "Stock"        = "steelblue",
      "Max-Sharpe"   = "firebrick",
      "Min-Variance" = "darkgreen"
    )
  ) +
  scale_x_continuous(labels = scales::percent_format(1)) +
  scale_y_continuous(labels = scales::percent_format(1)) +
  labs(
    title    = "6th Year Risk–Return Comparison",
    subtitle = "Stocks vs. Max-Sharpe vs. Min-Variance",
    x        = "Std Dev (Monthly Risk)",
    y        = "Avg Monthly Return"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    legend.title  = element_blank()
  )

print(Rf_monthly)



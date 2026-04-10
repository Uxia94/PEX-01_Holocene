# =========================
# Statistical assessment of δ2Hterr – δ¹³Cterr covariation
# =========================
# Libraries
library(readxl)
library(dplyr)

# 1. Load data

data <- read_excel("Data_figure_5.xlsx", sheet = 1)
#colnames(data)

# 2. Remove missing values

clean_data <- na.omit(data[, c("dD", "d13C")])

# 3. Pearson correlation (raw data)

cor(clean_data$dD,
    clean_data$d13C,
    method = "pearson")

cor.test(clean_data$dD,
         clean_data$d13C,
         method = "pearson")

# 4. Scatt plot + regression line

plot(clean_data$dD,
     clean_data$d13C,
     pch = 19,
     xlab = expression(delta*D),
     ylab = expression(delta^{13}*C))

abline(lm(d13C ~ dD, data = clean_data),
       col = "red",
       lwd = 2)

# 5. Linear detrending

time <- 1:nrow(clean_data) #create time index

# remove linear trends
model_D <- lm(dD ~ time, data = clean_data)
residuals_D <- resid(model_D)

model_C <- lm(d13C ~ time, data = clean_data)
residuals_C <- resid(model_C)

# correlation of detrended series
cor.test(residuals_D, residuals_C)

# 6. First-differencing

diff_D <- diff(clean_data$dD)
diff_C <- diff(clean_data$d13C)

cor.test(diff_D, diff_C)

# 7. Autocorrelation check

acf(clean_data$dD, main = "ACF dD")
acf(clean_data$d13C, main = "ACF d13C")


















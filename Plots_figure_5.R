# =========================
# δ2Hterr z-score
# =========================
# Libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Load data
z <- read_excel("Data_figure_5.xlsx", sheet = 1)

# 2. Convert variables to numeric
z <- z %>%
  mutate(
    Edad_calBP = as.numeric(`Age (calBP)`),
    dD = as.numeric(dD),
  )

# 3. Calculate z-score

z <- z %>%
  mutate(
    dD_z = scale(dD)[,1]
  )

# 4. Reshape data (optional structure)

z_long <- z %>%
  select(`Age (calBP)`, dD_z) %>%
  pivot_longer(
    cols = c(dD_z),
    names_to = "variable",
    values_to = "z"
  )

# 5. Define colors

green <- rgb(0.2, 0.6, 0.3, 0.6) #positive anomalies
purple <- rgb(0.5, 0.3, 0.7, 0.6) #negative anomalies

# 6. Define area-filling function

fill_area <- function(x,y,condition,col){
  r<-rle(condition)
  ends<-cumsum(r$lengths)
  starts<-c(1,head(ends+1,-1))
  
  for(i in seq_along(r$values)){
    if(r$values[i]){
      idx<-starts[i]:ends[i]
      polygon(c(x[idx],rev(x[idx])),
              c(y[idx],rep(0,length(idx))),
              col=col,
              border=NA)}}}

# 7. Define plotting variables

x <- z$`Age (calBP)`
y <- z$dD_z

# 8. Plot

plot(x, y, 
     type = "n", 
     xlab = "Years BP", ylab = "dD (z-score)", xaxt = "n") #initialize plot

axis(1, at = seq(-1000, 12000, by = 1000)) # custom x-axis

fill_area(x, y, y >= 0, green)
fill_area(x, y, y < 0, purple) # fill positive and negative areas

lines(x, y, col = "black", lwd = 1) # draw the line
abline(h = 0, col = "black", lwd = 0.8) # add reference line

# =========================
# NAO-like index
# =========================
# Libraries
library(readxl)
library(dplyr)

# 1. Load data

df <- read_excel("Data_figure_5.xlsx", sheet = 2)

# 2. Define bin size

bin_size <- 500

# 3. Create time bins

df$bin_start <- floor(df$Age..decade_BP. / bin_size) * bin_size
df$bin_end   <- df$bin_start + bin_size

# 4. Compute median per bin

df_bin <- df %>%
  dplyr::group_by(bin_start, bin_end) %>%
  dplyr::summarise(
    nao_med = median(NAO_like, na.rm = TRUE),
    .groups = "drop"
  )

# 5. Plot

plot(NULL,
     xlim = c(-1000, 12000),
     ylim = range(df_bin$nao_med, na.rm = TRUE),
     xlab = "Years BP",
     ylab = "NAO-like index",
     xaxt = "n") # empty plot

axis(1, at = seq(-1000, 12000, by = 1000)) # add custom x-axis

for(i in 1:nrow(df_bin)) {
  
  x0 <- df_bin$bin_start[i]
  x1 <- df_bin$bin_end[i]
  y  <- df_bin$nao_med[i] # draw bars manually
  
  col_fill <- ifelse(y >= 0,
                     rgb(1,0,0,0.6), # red: NAO +
                     rgb(0,0,1,0.6)) # blue: NAO -
  
  polygon(c(x0, x1, x1, x0),
          c(0, 0, y, y),
          col = col_fill,
          border = NA) # draw each bar
}





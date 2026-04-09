# =========================
# XRF_clr
# =========================

# Libraries
library(readxl)
library(compositions)
library(tidyverse)

# 1. Load data

xrf <- read_excel("XRF_Peixao_clean.xlsx")

meta <- xrf %>%
  select(`Depth (cm)`, `Age BP`) 

# 2. Select compositional data

xrf_data <- xrf %>%
  select(-`Depth (cm)`, -`Age BP`)

# 3. clar transformation

xrf_clr <- clr(xrf_data)

# 4. Recombine data

xrf_clr <- cbind(meta, xrf_clr)

# =========================
# PCA & PCA plot
# =========================

# 1. PCA

PCA <- prcomp(xrf_clr[, -(1:2)], center = TRUE, scale. = TRUE)

# 1. Plot
library(ellipse)

scores   <- as.data.frame(PCA$x) 
loadings <- as.data.frame(PCA$rotation)
imp      <- summary(PCA)$importance

pch.group_2025 <- c(rep(21, times=502), rep(22, times=423), rep(24, times=204))
col.group_2025 <- c(rep("green", times=502), rep("gold", times=423), rep("red", times=204))

#dim(PCA$x)[1]
PCA$x[,1:2]

summary_PCA <- summary(PCA)
par(mfrow = c(1, 1))
plot1 <- plot(PCA$x[,1], 
              PCA$x[,2], 
              xlab=paste("PCA 1 (", round(summary_PCA$importance[2]*100, 1), "%)", sep = ""), 
              ylab=paste("PCA 2 (", round(summary_PCA$importance[5]*100, 1), "%)", sep = ""), pch=pch.group_2025, 
              col="black", bg=col.group_2025, cex=1, las=1, asp=1)


#text(PCA$x[,1], PCA$x[,2], labels=1:nrow(PCA$x), pos=3, cex=0.7)

abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")


fac <- 0.8 * min(
  diff(range(scores[,1])) / diff(range(loadings[,1])),
  diff(range(scores[,2])) / diff(range(loadings[,2]))
)


l.x <- loadings[,1] * fac
l.y <- loadings[,2] * fac
arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red2", length=0.15, lwd=0.8)

l.pos <- l.y
lo <- which(l.y < 0)
hi <- which(l.y > 0)

l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")

text(l.x, l.y, labels=row.names(PCA$rotation), col="red", pos=l.pos, cex=0.7)
par(xpd = TRUE)
legend("topleft", legend=c("LH", "MH", "EH"), 
       col="black", pt.bg=c("green", "gold", "red"), 
       pch=c(21, 22, 24), pt.cex=1, cex = 0.5)


library(ellipse)
tab <- matrix(c(PCA$x[,1], PCA$x[,2]), ncol=2)
c1 <- cor(tab[1:502,])
c2 <- cor(tab[503:925,])
c3 <- cor(tab[926:1129,])


polygon(ellipse(c1*(max(abs(PCA$rotation))*1), 
                centre=colMeans(tab[1:502,]), level=0.95), 
        col=adjustcolor("green", alpha.f=0.25), border="green4")
polygon(ellipse(c2*(max(abs(PCA$rotation))*1), 
                centre=colMeans(tab[502:925,]), level=0.95), 
        col=adjustcolor("gold", alpha.f=0.25), border="gold2")
polygon(ellipse(c3*(max(abs(PCA$rotation))*1), 
                centre=colMeans(tab[926:1129,]), level=0.95), 
        col=adjustcolor("red", alpha.f=0.25), border="red4")

# =========================
# Broken stick
# =========================
library(vegan)

bstick(PCA)
screeplot(PCA, bstick = TRUE, type = "lines")

# =========================
# ANOSIM
# =========================
library(dplyr)

# 1. Extract PCA scores

scores <- as.data.frame(PCA$x)

# 2. Add age

scores$depth <- xrf_clr$`Depth (cm)`

# 3. Create groups

scores <- scores %>%
  mutate(
    Group = case_when(
      depth >= 612 ~ "Early",
      depth < 612 & age >= 332 ~ "Mid",
      depth < 332 ~ "Late"
    )
  )

scores$Group <- as.factor(scores$Group)

# 4. Select PC1 and PC2

scores_sel <- scores[, c("PC1", "PC2")]

# 5. Distance

dist_matrix <- dist(scores_sel, method = "euclidean")

# 6. ANOSIM

anosim_result <- anosim(dist_matrix, grouping = scores$Group)

# 6. Results

summary(anosim_result)
plot(anosim_result)

# =========================
# PC1 time series plot
# =========================
library(ggplot2)
library(dplyr)
library(zoo)

# 1. Prepare data
# Extract PC1 scores and corresponding ages

df_PC1 <- data.frame(
  Age_calBP = xrf_clr$`Age BP`,     
  PC1 = PCA$x[, 1]             
)

df_PC1$Age_interp <- na.approx(df_PC1$Age_calBP, na.rm = FALSE)

# 2. Raw PC1 plot

ggplot(df_PC1, aes(x = Age_interp, y = PC1)) +
  
  geom_line(color = "black", size = 1.2) +
  
  scale_x_reverse(
    breaks = seq(11000, -1000, by = -1000)
  ) +
  
  coord_cartesian(xlim = c(11000, -1000)) +
  
  labs(
    x = "Age (cal BP)",
    y = "PC1"
  ) +
  
  theme_classic(base_size = 14)

# 3. 200-year smoothing function

smooth_200yr <- function(x, Age_calBP, window = 200) {
  sapply(seq_along(x), function(i) {
    idx <- which(abs(Age_calBP - Age_calBP[i]) <= window/2)
    mean(x[idx], na.rm = TRUE)
  })
}

df_PC1$PC1_smooth <- smooth_200yr(df_PC1$PC1, df_PC1$Age_interp, 200) # apply smoothing

# 4. Plot raw + smoothed PC1

ggplot(df_PC1, aes(x = Age_interp)) +
  
  geom_line(aes(y = PC1),
            color = "black",
            alpha = 0.4,
            linewidth = 0.8) +
  
  geom_line(aes(y = PC1_smooth),
            color = "blue",
            linewidth = 1.3) +
  
  scale_x_reverse(
    breaks = seq(11000, -1000, by = -1000)
  ) +
  
  coord_cartesian(xlim = c(11000, -1000)) +
  
  labs(
    x = "Age (cal BP)",
    y = "PC1"
  ) +
  
  theme_classic(base_size = 14)

# =========================
# kclr time series plot
# =========================

# 1. Prepare data

df_K <- data.frame(
  Age_calBP = xrf_clr$`Age BP`,
  K_clr = xrf_clr$K
)

df_K$Age_interp <- na.approx(df_K$Age_calBP, na.rm = FALSE)

# 2. Smoothing function

smooth_200yr <- function(x, age, window = 200) {
  sapply(seq_along(x), function(i) {
    idx <- which(abs(age - age[i]) <= window / 2)
    mean(x[idx], na.rm = TRUE)
  })
}

df_K$K_smooth <- smooth_200yr(df_K$K_clr, df_K$Age_interp, 200)

# 3. Plot

ggplot(df_K, aes(x = Age_interp)) +
  
  geom_line(aes(y = K_clr),
            color = "black",
            alpha = 0.4,
            linewidth = 0.8) +
  
  geom_line(aes(y = K_smooth),
            color = "blue",
            linewidth = 1.3) +
  
  scale_x_reverse(
    breaks = seq(11000, -1000, by = -1000)
  ) +
  
  coord_cartesian(xlim = c(11000, -1000)) +
  
  labs(
    x = "Age (cal BP)",
    y = "K (clr)"
  ) +
  
  theme_classic(base_size = 14)








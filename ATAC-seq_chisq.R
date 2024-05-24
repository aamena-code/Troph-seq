### ATAC-seq: basic stats to test significance of distributions of DARs ###

library("gplots")
library ("ggplot2")
library ("ggpubr")
library("vcd")

nDARs <- read.delim (file = "Chi-sq_ATAC.txt", header = T, row.names = 1)

nDARs_subset <- nDARs[1:(nrow(nDARs)-2), ]
chisq <- chisq.test(nDARs_subset)
chisq
library(corrplot)
corrplot(chisq$residuals, is.cor = FALSE, cl.ratio = 0.6, cl.align.text = "l")
residual_df <- as.data.frame(chisq$residuals)

# Create a data frame with residuals and comparisons
residual_df <- data.frame(Comparison = rownames(residuals), residuals)

# Reshape the data for plotting
residual_df <- reshape2::melt(residual_df, id.vars = "Comparison")

ggplot(residual_df, aes(x = Comparison, y = value, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, guide = "legend") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", fill = "Residuals")

# make 2x2 contingency tables
first <- chisq.test(nDARs[1:2,])
sec <- chisq.test(nDARs[3:4, ])
third <- chisq.test(nDARs[5:6, ])

# chi sq for DisttoTSS tables
downDARs_Dist <- read.delim (file = "/chi-sq_downDARs_disttoTSS.txt", header = T, row.names = 1)
upDARs_Dist <- read.delim (file = "/chi-sq_upDARs_disttoTSS.txt", header = T, row.names = 1)
allDARs_Dist <- read.delim (file = "/chi-sq_allDARs_disttoTSS.txt", header = T, row.names = 1)

chisq.test(downDARs_Dist)
chisq.test(upDARs_Dist)
chisq.test(allDARs_Dist)


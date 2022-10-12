# ovCancer-signature - Reproducibility code

#  Copyright (C) 2022 Fernando Palluzzi
#  e-mail: <fernando.palluzzi@gmail.com>
#  Bioinformatics facility 
#  Gemelli Science and Technological Park (GSTeP)
#  Fondazione Policlinico Universitario Agostino Gemelli IRCCS,
#  Largo Agostino Gemelli 8, 00168 Roma, Italy

#  ovCancer-signature is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  ovCancer-signature is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.


##### Ten-genes signature evaluation

library(randomForest)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(stringr)
library(dplyr)


performance <- function(M, tp = "topleft") {
	if (tp == "topleft") {
		TP <- M[1, 1]
		FP <- M[1, 2]
		FN <- M[2, 1]
		TN <- M[2, 2]
	} else if (tp == "bottomright") {
		TP <- M[2, 2]
		FP <- M[2, 1]
		FN <- M[1, 2]
		TN <- M[1, 1]		
	}
	Se <- TP/(TP + FN)
	Sp <- TN/(TN + FP)
	PPV <- TP/(TP + FP)
	F1 <- 2*PPV*Se/(PPV + Se)
	A <- sum(diag(M))/sum(M)
	return(list(Se = Se, Sp = Sp, PPV = PPV, F1 = F1, A = A))
}


# Importing RFC dataset

load("~/ovCancer-signature/data/Provac_RFC10.RData")
data <- rfc10.provac$data


##### Defining training sets

T1 <- data[data$K == "A" | data$K == "B" | data$K == "C",]
T2 <- data[data$K == "A" | data$K == "B" | data$K == "D",]
T3 <- data[data$K == "A" | data$K == "C" | data$K == "D",]
T4 <- data[data$K == "B" | data$K == "C" | data$K == "A",]


# Trees number
n <- 5000


##### Model training

if (FALSE) {
# Top-42 DEGs

model <- formula(paste0(c("as.factor(y) ~ ABDHD5 + ACVR1B + ALOX5AP + ",
                          "C3AR1 + CD4 + CKB + CPZ + CTNNBL1 + ",
                          "CXCL16 + DHX35 + DSN1 + GNG11 + HERC5 + ",
                          "IGFBP7 + IGSF9 + INPP5D + LDLR + LINC01816 + ",
                          "MILR1 + MYO5A + NEXN + NOL4L + PLCB1 + PLCG2 + ",
                          "PODN + RNF24 + RPRD1B + SELPLG + SERPINF1 + ",
                          "SLC15A3 + STOM + TACC1 + TARBP2 + TMEM140 + ",
                          "TNFSF13B + TPP1 + TRPM2 + TSPAN31 + TTI1 + ",
                          "UQCC1 + WDPCP + ZNF738"), collapse = ""))
}

# Best predictors

model <- formula(paste0(c("as.factor(y) ~ UQCC1 + CKB + TSPAN31 + ",
                          "GNG11 + PLCG2 + IGFBP7 + SLC15A3 + CTNNBL1 + ",
                          "RNF24 + TTI1"), collapse = ""))

F1 <- randomForest(model, data = T1, ntree = n, importance = TRUE, mtry = 3)
F2 <- randomForest(model, data = T2, ntree = n, importance = TRUE, mtry = 3)
F3 <- randomForest(model, data = T3, ntree = n, importance = TRUE, mtry = 3)
F4 <- randomForest(model, data = T4, ntree = n, importance = TRUE, mtry = 3)


##### Defining validation sets

V1 <- data[data$K == "D",]
V2 <- data[data$K == "C",]
V3 <- data[data$K == "B",]
V4 <- data[data$K == "A",]


##### Defining ranking by importance

R1 <- data.frame(importance(F1))[, 3:4]
R2 <- data.frame(importance(F2))[, 3:4]
R3 <- data.frame(importance(F3))[, 3:4]
R4 <- data.frame(importance(F4))[, 3:4]

R <- cbind(R1, R2, R3, R4)
R$MDA <- apply(R[, c(1, 3, 5, 7)], 1, mean)
R$MDG <- apply(R[, c(2, 4, 6, 8)], 1, mean)
R <- R[, c(9, 10)]

# Minmax normalization
R$fA <- 100*(R$MDA - min(R$MDA))/(max(R$MDA) - min(R$MDA))
R$fG <- 100*(R$MDG - min(R$MDG))/(max(R$MDG) - min(R$MDG))

# Final score (f)
R$f <- apply(R[, c(3, 4)], 1, mean)
R <- R[order(R$f, decreasing = TRUE),]
R


##### Prediction and performances

C1 <- predict(F1, V1)
C2 <- predict(F2, V2)
C3 <- predict(F3, V3)
C4 <- predict(F4, V4)

TP <- rbind(V1[C1 == V1$y & V1$y == 1,], V2[C2 == V2$y & V2$y == 1,], V3[C3 == V3$y & V3$y == 1,], V4[C4 == V4$y & V4$y == 1,])
TN <- rbind(V1[C1 == V1$y & V1$y == 0,], V2[C2 == V2$y & V2$y == 0,], V3[C3 == V3$y & V3$y == 0,], V4[C4 == V4$y & V4$y == 0,])
FP <- rbind(V1[C1 == 1 & V1$y == 0,], V2[C2 == 1 & V2$y == 0,], V3[C3 == 1 & V3$y == 0,], V4[C4 == 1 & V4$y == 0,])
FN <- rbind(V1[C1 == 0 & V1$y == 1,], V2[C2 == 0 & V2$y == 1,], V3[C3 == 0 & V3$y == 1,], V4[C4 == 0 & V4$y == 1,])

confusion <- matrix(c(nrow(TP), nrow(FP), nrow(FN), nrow(TN)), nrow = 2, byrow = TRUE)

rownames(confusion) <- c("Pred1", "Pred0")
colnames(confusion) <- c("Real1", "Real0")
confusion

P <- performance(confusion)
P


##### RNA-seq DEGs Heatmap (top-42 markers)

x <- log2(exprs[, 2:15] + 1)
x <- x[exprs$symbol %in% top,]
x <- t(apply(x, 1, scale))
rownames(x) <- exprs$symbol[exprs$symbol %in% top]
colnames(x) <- colnames(exprs[, c(2:15)])
x <- x[, c(1:3, 7, 4:6, 8:10, 14, 11:13)]
colnames(x) <- c("1", "2", "3", "4", "5", "6", "7",
                 "8", "9", "10", "11", "12", "13", "14")

pdf("Provac_RNAseq_Top42.pdf", width = 15, height = 10)
colors <- list(Phenotype = c("Sensitive" = "lightblue", "Resistant" = "gold"))
hann <- HeatmapAnnotation(Phenotype = c(rep("Sensitive", 7), rep("Resistant", 7)),
                          col = colors,
                          show_legend = c(FALSE, FALSE))
col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("blue", "blue", "white", "red", "red"))
Heatmap(as.matrix(x), name = "Standardized\nlog2(counts)",
        top_annotation = hann,
        row_names_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 24),
        column_names_rot = 0,
        heatmap_legend_param = list(title_gp = gpar(col = "black", fontsize = 16),
                                    labels_gp = gpar(col = "black", fontsize = 16),
                                    direction = "horizontal"))
pheno <- Legend(labels = c("Sensitive", "Resistant"),
                legend_gp = gpar(fill = c("lightblue", "gold")),
                title = "Phenotype",
                title_gp = gpar(col = "black", fontsize = 16),
                labels_gp = gpar(col = "black", fontsize = 16),
                direction = "horizontal")
draw(pheno, x = unit(36, "cm"), y = unit(15, "cm"))
dev.off()


##### Volcano plots (DEGs)

# RNA-seq DEGs

exprs <- read.delim("~/ovCancer-signature/data/Provac_ResistantVsSensitive_raw.txt", stringsAsFactors = FALSE)

exprs$Regulation <- "NR"
exprs$Regulation[exprs$log2FoldChange > 1 & exprs$padj < 0.05] <- "UP"
exprs$Regulation[exprs$log2FoldChange < -1 & exprs$padj < 0.05] <- "DOWN"
cols <- c("UP" = "red2", "DOWN" = "blue", "NR" = "grey60")

top <- exprs$symbol[exprs$Tier == "Top42"]
top.degs <- exprs[exprs$Tier == "Top42" & exprs$padj < 0.001,]
up <- exprs[exprs$log2FoldChange > 1 & exprs$padj < 0.001,]
down <- exprs[exprs$log2FoldChange < -1 & exprs$padj < 0.001,]

pdf("Provac_RNAseq_volcano.pdf", width = 20, height = 10)
ggplot(data = exprs,
       aes(x = log2FoldChange,
           y = log10padj)) +
  theme_bw() +
  theme(panel.border = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 26),
    axis.title = element_text(size = 28, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 28)) +
  geom_point(aes(colour = Regulation),
             alpha = 0.8,
             shape = 16,
             size = 3) +
  scale_colour_manual(values = cols) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  annotate("text", x = -3.8, y = 1.42, size = 8,
           label = "P.adj = 0.05",
           parse = FALSE) +
  labs(x = "log2(Resistant/Sensitive)", y = "-log10(BH-adjusted P-value)") +
  scale_x_continuous(breaks = c(seq(-10, 10, 1)),     
                     limits = c(-4, 4))
dev.off()

# HT RT-qPCR DEGs

W <- read.delim("~/ovCancer-signature/data/Provac_RNAseq_RTqPCR_summary.txt", stringsAsFactors = FALSE)

W$Regulation <- "NR"
W$Regulation[W$wPvalue <= 0.05 & W$shift > 0] <- "UP"
W$Regulation[W$wPvalue <= 0.05 & W$shift < 0] <- "DOWN"

cols <- c("UP" = "red2", "DOWN" = "blue", "NR" = "grey60")

pdf("Provac_RTqPCR_volcano.png", width = 20, height = 10)
lbl <- ifelse(W$symbol %in% c("GNG11", "RNF24", "UQCC1", "CTNNBL1",
                              "PLCG2", "TTI1", "IGFBP7", "SLC15A3",
                              "TSPAN31", "CKB"), W$symbol, "")
ggplot(data = W,
       aes(x = shift,
           y = -log10(wPvalue))) +
  theme_bw() +
  theme(panel.border = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 26),
    axis.title = element_text(size = 22, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26)) +
  geom_point(aes(colour = Regulation),
             alpha = 0.8,
             shape = 16,
             size = 5) +
  scale_colour_manual(values = cols) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed") +
  annotate("text", x = -1.2, y = 1.35, size = 7,
           label = "P-value = 0.05",
           parse = FALSE) +
  geom_text(label = lbl, cex = 6, nudge_x = -0.04, nudge_y = 0.06, show.legend = FALSE) +
  labs(x = "Wilcoxon estimated shift", y = "-log10(P-value)") +
  scale_x_continuous(breaks = c(seq(-2, 2, 0.25)),
                     limits = c(-1.25, 1.25))
dev.off()


##### Polar bar plot of the top-10 markers

x <- data.frame(Gene = c("GNG11", "RNF24", "UQCC1", "CTNNBL1", "PLCG2",
                         "TTI1", "IGFBP7", "SLC15A3", "TSPAN31", "CKB"),
                Estimate = c(-1.09, 0.59, 0.54, 0.53, -0.91, 0.59, -0.42,
                             -0.98, 0.53, 0.94),
                P = c(0.00475, 0.00498, 0.01013, 0.01781, 0.02592,
                      0.02751, 0.02837, 0.02926, 0.04038, 0.04040),
                CI95 = c("-1.71, -0.34", "0.20, 1.12", "0.15, 0.93",
                         "0.05, 1.03", "-1.70, -0.13", "0.06, 1.21",
                         "-0.89, -0.06", "-1.75, -0.17", "0.03, 0.98",
                         "0.05, 1.79"),
                Marker = c("W", "WP", "W", "WP", "W", "WP", "W", "W", "W", "W"),
                Human = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                Primary = c(0, 1, 0, 1, 0, 1, 0, 0, 0, 0),
                Stabilized = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0))
x <- x[x$Human == 1,]
x$N <- apply(x[, 6:7], 1, sum)

pdf("Provac_polar_barplot.pdf", width = 20, height = 18)
plt <- ggplot(x) +
  theme_bw() +
  theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = c(29, 29, 29, 29, 29, 29, 29, 26, 29, 29)),
  axis.text.y = element_blank(),
  axis.title = element_text(size = 36, face = "bold"),
  legend.key.size = unit(1, "cm"),
  legend.text = element_text(size = 26),
  legend.title = element_text(size = 30)) +
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0, 1, 1.3, 2)),
    color = c("grey", "grey", "red3", "grey")
  ) + 
  geom_col(
    aes(
      x = reorder(str_wrap(Gene, 5), N),
      y = N,
      fill = Estimate
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = 0.9
  ) +
  geom_point(
    aes(
      x = reorder(str_wrap(Gene, 7), -log10(P)),
      y = -log10(P)
    ),
    size = 4,
    color = "gray12"
  ) +
  geom_segment(
    aes(
      x = reorder(str_wrap(Gene, 5), -log10(P)),
      y = -1,
      xend = reorder(str_wrap(Gene, 5), -log10(P)),
      yend = 2
    ),
    linetype = "dashed",
    color = "gray12"
  ) + 
  scale_fill_gradientn(
    "Shift",
     colours = c("deepskyblue", "white", "darkorange"),
     limits = c(-1.1, 1.1)
  ) +
  labs(x = "Selected markers", y = "1: Patients,   2: Patients and OV.GEM") +
  coord_polar()
plt
dev.off()

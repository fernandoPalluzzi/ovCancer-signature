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


##### Volcano plots (DEGs)

# 

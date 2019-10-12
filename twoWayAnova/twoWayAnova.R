#Two way ANOVA
#fictional experiment:
#survival time (in hours) after three different toxins (factor A)
#and four different treatments (factor B)
#----------------------------------------------------------------------------------

# Prepare the Data
rm(list = ls());

treat1 <- read.csv("treat1.csv");
treat1 <- as.matrix(treat1[, -1]);

treat2 <- read.csv("treat2.csv");
treat2 <- as.matrix(treat2[, -1]);

treat3 <- read.csv("treat3.csv");
treat3 <- as.matrix(treat3[, -1]);

treat4 <- read.csv("treat4.csv");
treat4 <- as.matrix(treat4[, -1]);


treat <- cbind(treat1, treat2, treat3, treat4);
dim(treat) <- c(nrow(treat1), ncol(treat1), 4);

treatDim <- dim(treat);

#Two way Anova begin

I <- treatDim[1];
J <- treatDim[2];
K <- treatDim[3];
N <- I*J*K;

dfA <- (I - 1);
dfB <- (J - 1);
dfAxB <- ((I - 1) * (J - 1));
dfE <- (I * J * (K - 1));
dfTot <- (I * J * K - 1);

grandMean <- mean(treat);

SQA <- 0; # Square sums of Factor A from grand mean (not the variable A)

for (i in c(1:treatDim[1])) {
  yi_u <- mean(treat[i, , ]);
  SQA <- SQA +(yi_u - grandMean)^2;
}

SQA <- (J * K * SQA);
MQA <- (SQA / dfA);

SQB <- 0;  # Square sums of Factor B from grand mean (not the variable B)

for (j in c(1:treatDim[2])) {
  yj_u <- mean(treat[, j, ]);
  SQB <- SQB + (yj_u - grandMean)^2;
}

SQB <- (I * K * SQB);
MQB <- (SQB / dfB);

SQAxB <- 0; # Square sums of Interaction between factor A and B 

for (i in c(1:treatDim[1])) {
  yi_u <- mean(treat[i, , ]);
  for (j in c(1:treatDim[2])) {
    yj_u <- mean(treat[, j, ]);
    yij_u <- mean(treat[i, j, ]);
    
    SQAxB <- SQAxB +(yij_u - yi_u - yj_u + grandMean)^2;
  }
}

SQAxB <- (SQAxB * K);
MQAxB <- (SQAxB / dfAxB);


SQE <- 0; # Square sums within the groups / residuals

for (i in c(1:treatDim[1])) {
  for (j in c(1:treatDim[2])) {
    yij_u <- mean(treat[i, j, ]);
    for (k in c(1:treatDim[3])) {
      yijk <- treat[i, j, k];
      SQE <- SQE + (yijk - yij_u)^2;
    }
  }
}

MQE <- (SQE / dfE);
SQTot <- (SQA + SQB + SQAxB + SQE);

#Calculate the F Statistics

alpha <- 0.05;

FA <- (MQA / MQE);
FB <- (MQB / MQE);
FAxB <- (MQAxB / MQE);

FTestResA <- FTest(FA, dfA, dfE, alpha);
FTestResB <- FTest(FB, dfB, dfE, alpha);
FTestResAxB <- FTest(FAxB, dfAxB, dfE, alpha);

#control the results

treatControl <- read.csv("BehandlungGiftTabelle.csv");
boxplot(treatControl$hours ~ treatControl$poison + treatControl$treatment);
summary(aov(treatControl$hours ~ treatControl$poison + treatControl$treatment +
              treatControl$poison:treatControl$treatment));



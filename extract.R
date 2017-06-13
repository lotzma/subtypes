## R functions for preparing the data set for analysis

library(deconstructSigs)
library(robust)

# Extract exposure of mutation catalog M to signatures P using non-negative least squares

# Function to extract signatures from mutation matrix

extract.sigs <- function(M, P, norm) {
  # Attempt to use off-the-shelf software
  E <- data.frame(row.names = rownames(M))
  #rel.sigs <- c('Sig1','Sig5','Sig7','Sig11','Sig17')
  rel.sigs <- c('Signature.1','Signature.5','Signature.7','Signature.11','Signature.17')
  E[,rel.sigs]=0
  for (x in rownames(M)) {
    sample <- whichSignatures(tumor.ref = M,
                              signatures.ref = P,
                              sample.id = x,
                              contexts.needed = TRUE,
                              tri.counts.method = norm,
                              associated = rel.sigs,
                              signature.cutoff = 0)
    E[x,rel.sigs] <- sample$weights[,rel.sigs] 
  }
  return(E)
}

ratio2count <- function(M, E, P, method) {
  # Get number of mutations caused by every signatures in P
  E.scaled <- E
  ratio <- tri.counts.genome/tri.counts.exome
  total.muts <- rowSums(M)
  for (x in rownames(E)) { # Samples
    for (y in colnames(E)) { # Signatures
      if (method == 'default') { 
        E.scaled[x,y] <- total.muts[x]*E[x,y]
      }
      else if (method == 'exome') {
        E.scaled[x,y] <- total.muts[x]*E[x,y]
      }
      else if (method == 'genome') {
        E.scaled[x,y] <- total.muts[x]*E[x,y]
      }
      else if (method == 'exome2genome') {
        up.scaled <- sum(norm.it(M[x,], ratio))
        prob.scaled <- sum(norm.it(P[y,], 1/ratio))
        E.scaled[x,y] <- E[x,y]*up.scaled*prob.scaled
      }
    }
  }
  return(E.scaled)
}

# Load data
load('data/tcga_samples.Rdata')
load('data/icgc_samples.Rdata')

mutations.tcga <- data.tcga[,-c(1:6)]
mutations.icgc <- data.icgc[,-c(1:4)]
sigdata.tcga <- data.tcga[,c(1:6)]
sigdata.icgc <- data.icgc[,c(1:4)]

# Create signature matrices
P.genome <- signatures.cosmic
#P.exome <- signatures.cosmic
# Renormalize for exome
#ratio <- tri.counts.genome/tri.counts.exome
#for (x in colnames(P.exome)) {
#  P.exome[x] <- norm.it(P.genome[x], ratio)
#}
#P.exome = P.exome/rowSums(P.exome)
  
# Extract signatures
# Function to extract signatures from mutation matrix
sigs.tcga <- extract.sigs(mutations.tcga, P.genome, 'default')
sigs.icgc <- extract.sigs(mutations.icgc, P.genome, 'genome')

# Rescale to create total mutation counts

#genome.length <- 3088286401
#exome.length <- 44000000
  
sigdata.tcga$TotalSNV <- rowSums(mutations.tcga)
sigdata.tcga$Sig1Rel <- sigs.tcga[,'Signature.1']
sigdata.tcga$Sig1Total <- sigdata.tcga$Sig1Rel*as.vector(sigdata.tcga$TotalSNV)

sigdata.icgc$TotalSNV <- rowSums(mutations.icgc)
sigdata.icgc$Sig1Rel <- sigs.icgc[,'Signature.1']
sigdata.icgc$Sig1Total <- sigdata.icgc$Sig1Rel*as.vector(sigdata.icgc$TotalSNV)

# Remove outliers
Y <- sqrt(2*sigdata.tcga$Sig1Total/sigdata.tcga$Age)
m <- median(Y)
sigdata.tcga <- sigdata.tcga[sigdata.tcga$Sig1Total/sigdata.tcga$Age<m+3,]

#Y <- sqrt(2*sigdata.icgc$Sig1Total/sigdata.icgc$Age)
#m <- median(Y)
#sigdata.icgc <- sigdata.icgc[sigdata.icgc$Sig1Total/sigdata.icgc$Age<m+3,]

# Compute averages in separate data frame
ages <- sort(unique(sigdata.tcga$Age))
averages <- data.frame(Age=ages)
averages$Sig1Mean <- NA
averages$Sig1Med <- NA
for (a in ages) {
    i <- match(a,averages$Age)
    averages[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a),"Sig1Total"])
    averages[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a),"Sig1Total"])
    if (nrow(subset(sigdata.tcga, Age==a))>0) {
      averages[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a), family=poisson))["(Intercept)"])
    }
}
   
# Compute averages in separate data frame (BRAF)
ages <- sort(unique(subset(sigdata.tcga, Cohort=="BRAF")$Age))
averages.braf <- data.frame(Age=ages)
averages.braf$Sig1Mean <- NA
averages.braf$Sig1Med <- NA
averages.braf$Sig1RobMean <- NA

for (a in ages) {
    i <- match(a,averages.braf$Age)
    averages.braf[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="BRAF"),"Sig1Total"])
    averages.braf[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="BRAF"),"Sig1Total"])
    if (nrow(subset(sigdata.tcga, Age==a & Cohort=="BRAF"))>0) {
      averages.braf[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="BRAF"), family=poisson))["(Intercept)"])
    }
}

# Compute averages in separate data frame (NRAS)
ages <- sort(unique(subset(sigdata.tcga, Cohort=="NRAS")$Age))
averages.nras <- data.frame(Age=ages)
averages.nras$Sig1Mean <- NA
averages.nras$Sig1Med <- NA
averages.nras$Sig1RobMean <- NA
for (a in ages) {
    i <- match(a,averages.nras$Age)
    averages.nras[i,c("Sig1Mean")] <- mean(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NRAS"),"Sig1Total"])
    averages.nras[i,c("Sig1Med")] <- median(sigdata.tcga[which(sigdata.tcga$Age==a & sigdata.tcga$Cohort=="NRAS"),"Sig1Total"])
    if (nrow(subset(sigdata.tcga, Age==a & Cohort=="NRAS"))>0) {
      averages.nras[i,c("Sig1RobMean")] <- exp(coef(glmRob(as.integer(Sig1Total)~1, data=subset(sigdata.tcga, Age==a & Cohort=="NRAS"), family=poisson))["(Intercept)"])
    }
}

averages.braf$Cohort <- "BRAF"
averages.nras$Cohort <- "NRAS"
averages.brafnras <- rbind(averages.braf, averages.nras)
    
# Save into file
#write.csv(sigdata.tcga, 'data/tcga_4class.csv', sep=',')
save(sigdata.tcga, file = 'data/tcga_4class.Rdata')
write.table(averages, 'data/averages.csv', sep=',')
write.table(averages.braf, 'data/averages_braf.csv', sep=',')
write.table(averages.nras, 'data/averages_nras.csv', sep=',')
write.table(averages.brafnras, 'data/averages_brafnras.csv', sep=',')

# Save ICGC data
#write.csv(sigdata.icgc, 'data/icgc_4class.csv', sep=',')
save(sigdata.icgc, file = 'data/icgc_4class.Rdata')
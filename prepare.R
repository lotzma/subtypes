## Prepare a data frame with mutation counts and additional information / no signatures

## TCGA Data

# All samples with age and gender from cBioPortal
data.tcga       <- read.table(file = 'data/skcm_tcga_clinical_data.tsv', sep = '\t', header = TRUE)
data.tcga       <- data.tcga[,c('Patient.ID','Diagnosis.Age','Person.Gender', 'Tumor.Site')]
colnames(data.tcga) <- c('ID','Age','Gender', 'Site')

# Remove duplicate samples
data.tcga       <- data.tcga[!duplicated(data.tcga$ID),]

# Remove samples with no age information
data.tcga       <- data.tcga[!is.na(data.tcga$Age),]

# BRAF/NRAS/NF1 information from Simon
data.tcga$Cohort <- rep('W3',nrow(data.tcga))
braf <- read.csv('data/tcga.346.braf.mutant.csv', header=FALSE)
nras <- read.csv('data/tcga.346.nras.mutant.csv', header=FALSE)
nf1 <- read.csv('data/tcga.346.nf1.mutant.csv', header=FALSE)
braf$Cohort <- rep('BRAF', nrow(braf))
nras$Cohort <- rep('NRAS', nrow(nras))
nf1$Cohort <- rep('NF1', nrow(nf1))
allcohort.tcga <- rbind(rbind(braf, nras), nf1)
allcohort.tcga <- allcohort.tcga[,c('V1', 'Cohort')]
# Remove multiple entries
brafnras <- merge(braf,nras,by.x="V1", by.y="V1")[,1]
allcohort.tcga <- allcohort.tcga[! allcohort.tcga$V1 %in% brafnras, ]
data.tcga <- data.tcga[! data.tcga$ID %in% brafnras,]
# Remove remaining duplicates (these are the ones with BRAF/NRAS and NF1)
dups <- duplicated(allcohort.tcga$V1)
allcohort.tcga <- allcohort.tcga[!dups,]
# Enter information into data frame
indices <- match(allcohort.tcga$V1, data.tcga$ID)
indices <- indices[!is.na(indices)]
indices2 <- match(data.tcga[indices,c('ID')], allcohort.tcga$V1)
data.tcga[indices,c('Cohort')] <- allcohort.tcga[indices2,c('Cohort')]

# mc1r information from MC1R paper
mc1r <- read.csv('data/ncomms12064-s3.csv', header=TRUE, skip=1)
data.tcga$MC1R <- rep(NA, nrow(data.tcga))
indices <- match(mc1r$BARCODE, data.tcga$ID)
indices <- indices[!is.na(indices)]
indices2 <- match(data.tcga[indices,'ID'], mc1r$BARCODE)
data.tcga[indices,'MC1R'] <- as.character(mc1r[indices2,c('rgeno')])

# Mutation count from TCGA
load("data/tcga.344.sigs.input.RData")
mutations.tcga  <- sigs.input
ids.tcga        <- row.names(mutations.tcga)

indices          <- match(ids.tcga, data.tcga$ID)
indices          <- indices[!is.na(indices)]
data.tcga        <- data.tcga[indices,]

# Extract rows containing mutations
indices <- match(data.tcga$ID,ids.tcga)
indices <- indices[!is.na(indices)]
mutations.tcga <- mutations.tcga[ids.tcga[indices],]

row.names(data.tcga) <- seq(1,nrow(data.tcga))
colnames(data.tcga) <- c('ID','Age','Gender','Site','Cohort','MC1R')
# Add mutations and principal components
data.tcga <- cbind(data.tcga, mutations.tcga)

## ICGC data

# All samples with age, gender and cohort from ICGC
data.icgc           <- read.csv("data/icgc.au.cutaneous.124.clinical.data.csv")
data.icgc           <- data.icgc[,c('Sample','Age','Sex','Mutant')]
colnames(data.icgc) <- c("ID", "Age", "Gender", "Cohort")

# Mutation count from ICGC
load("data/icgc.au.mel.cutaneous.124.sigs.input.RData")
mutations.icgc  <- sigs.input
ids.icgc        <- row.names(mutations.icgc)
indices <- match(data.icgc$ID,ids.icgc)
indices <- indices[!is.na(indices)]
mutations.icgc <- mutations.icgc[ids.icgc[indices],]
data.icgc <- cbind(data.icgc, mutations.icgc)

# Change names
levels(data.icgc$Cohort) <- c(levels(data.icgc$Cohort), 'W3')
data.icgc$Cohort[data.icgc$Cohort=='triple_wild_type'] <- 'W3'
levels(data.icgc$Gender) <- c(levels(data.icgc$Gender), 'MALE', 'FEMALE')
data.icgc$Gender[data.icgc$Gender=='male'] <- 'MALE'
data.icgc$Gender[data.icgc$Gender=='female'] <- 'FEMALE'
data.icgc <- data.icgc[data.icgc$Gender=='FEMALE' | data.icgc$Gender=='MALE',]
data.icgc <- data.icgc[!(data.icgc$Cohort=='HRAS'),]

# Save into file
#write.csv(data.tcga, 'data/tcga_samples.csv')
#write.csv(data.icgc, 'data/icgc_samples.csv')
save(data.tcga, file = 'data/tcga_samples.Rdata')
save(data.icgc, file = 'data/icgc_samples.Rdata')

# Prepare data for hdp package
rownames(data.tcga) <- data.tcga$ID
tcga.braf <- subset(data.tcga, Cohort=='BRAF')[,-c(1:6)]
tcga.nras <- subset(data.tcga, Cohort=='NRAS')[,-c(1:6)]
tcga.nf1 <- subset(data.tcga, Cohort=='NF1')[,-c(1:6)]
tcga.w3 <- subset(data.tcga, Cohort=='W3')[,-c(1:6)]

tcga.all <- do.call("rbind", list(tcga.braf, tcga.nras, tcga.nf1, tcga.w3))
#write.csv(tcga.all, 'data/tcga_all_sorted.csv')
save(tcga.all, file = 'data/tcga_all_sorted.Rdata')

# Prepare data for hdp package
rownames(data.icgc) <- data.icgc$ID
icgc.braf <- subset(data.icgc, Cohort=='BRAF')[,-c(1:4)]
icgc.nras <- subset(data.icgc, Cohort=='NRAS')[,-c(1:4)]
icgc.nf1 <- subset(data.icgc, Cohort=='NF1')[,-c(1:4)]
icgc.w3 <- subset(data.icgc, Cohort=='W3')[,-c(1:4)]

icgc.all <- do.call("rbind", list(icgc.braf, icgc.nras, icgc.nf1, icgc.w3))
#write.csv(icgc.all, 'data/icgc_all_sorted.csv')
save(icgc.all, file = 'data/icgc_all_sorted.Rdata')
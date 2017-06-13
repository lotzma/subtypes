library(hdp)
library(Matrix)

#Load data set(s)
mut_data <- read.csv("data/icgc_all_sorted.csv", row.names=1)
mut_data <- data.frame(sapply(mut_data, as.numeric))

group_size <- c(56, 42, 12)

# initialise HDP, assign data, activate
hdp <- hdp_init(ppindex = c(0, 1, 1, 1, rep(2, group_size[1]), rep(3, group_size[2]), rep(4, group_size[3])),
                cpindex = c(1, 2, 2, 2, rep(3, group_size[1]), rep(4, group_size[2]), rep(5, group_size[3])),
                hh = rep(1, 96),
                alphaa = rep(1, 5),
                alphab = rep(1, 5))

hdp <- hdp_setdata(hdp, 5:numdp(hdp), mut_data)
hdp <- dp_activate(hdp, 1:numdp(hdp), initcc=4, seed = 740613)

# run four posterior sampling chains
hdplist <- list(hdp, hdp, hdp, hdp)
seeds <- c(3106053, 6418374, 667621, 7699454)
for (i in 1:4){
    hdplist[[i]] <- hdp_posterior(hdplist[[i]], 10000, 500, 50, 3, seed=seeds[i])
}

hdp_multi <- hdp_multi_chain(hdplist)

# extract components
hdp_multi <- hdp_extract_components(hdp_multi)
save(hdp_multi, file="hdp_multi_comp.RData")

# assess posterior sampling chains
par(mfrow=c(2,2))
lapply(chains(hdp_multi), plot_lik, bty='L')
lapply(chains(hdp_multi), plot_numcluster, bty='L')
lapply(chains(hdp_multi), plot_data_assigned, bty='L')
par(mfrow=c(1,1))

# signatures plot
bases <- c("A", "C", "G", "T")
trinuc_context <- paste0(rep(rep(bases, times=6), each=4),
                         rep(c("C", "T"), each=48),
                         rep(bases, times=24))
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

png(filename="sigs.png", width=1800, height=900)
par(mfrow=c(3,3), mar=c(2,3,1,1), cex=1.35)
plot_comp_distn(hdp_multi, comp=0:5,
                grouping=group_factor, col=RColorBrewer::brewer.pal(8, "Set2")[-c(6,7)],
                col_nonsig="grey80", show_group_labels=TRUE, yaxp=c(0, 0.3, 3),
                plot_title = paste("Signature", 0:5), cex.cat=0.9)
#par(mar=c(3,3,1,1))
#plot_comp_distn(hdp_multi, comp=6:8, cat_names=trinuc_context,
#                grouping=group_factor, col=RColorBrewer::brewer.pal(8, "Set2")[-c(6,7)],
#                col_nonsig="grey80", show_group_labels=TRUE, yaxp=c(0, 0.3, 3),
#                plot_title = paste("Signature", 6:8), cex.cat=0.9)
dev.off()

# exposures plots
png(filename="exp_BRAF.png", width=900, height=600)
plot_dp_comp_exposure(hdp_multi, 5 + 1:group_size[1],
                      RColorBrewer::brewer.pal(12, "Set3"),
                      incl_numdata_plot=TRUE, incl_nonsig=TRUE,
                      main_text="BRAF", ylab_numdata="SNV count",
                      ylab_exp="Signature exposure",
                      leg.title="Signature", cex.axis = 1)
dev.off()

png(filename="exp_NRAS.png", width=900, height=600)
plot_dp_comp_exposure(hdp_multi, 5+group_size[1] + 1:group_size[2],
                      RColorBrewer::brewer.pal(12, "Set3"),
                      incl_numdata_plot=TRUE, incl_nonsig=TRUE,
                      main_text="NRAS", ylab_numdata="SNV count",
                      ylab_exp="Signature exposure",
                      leg.title="Signature", cex.axis = 1)
dev.off()

png(filename="exp_NF1.png", width=900, height=600)
plot_dp_comp_exposure(hdp_multi, 5+group_size[1]+group_size[2] + 1:group_size[3],
                      RColorBrewer::brewer.pal(12, "Set3"),
                      incl_numdata_plot=TRUE, incl_nonsig=TRUE,
                      main_text="NF1", ylab_numdata="SNV count",
                      ylab_exp="Signature exposure",
                      leg.title="Signature", cex.axis = 1)
dev.off()

png(filename="exp_group.png", width=600, height=450)
par(mfrow=c(1,1), mar=c(8,4,3,2), cex=1)
plot_dp_comp_exposure(hdp_multi, 2:4,
                      RColorBrewer::brewer.pal(12, "Set3"),
                      incl_numdata_plot=FALSE,
                      main_text="Group signature distribution",
                      ylab_exp="Signature exposure",
                      leg.title="Signature", cex.axis = 1,
                      dpnames=c("BRAF", "NRAS", "NF1"),
                      cex.names=0.8)
dev.off()

braf_group_ci <-comp_dp_distn(hdp_multi)$cred.int[[2]]
nras_group_ci <- comp_dp_distn(hdp_multi)$cred.int[[3]]
nf1_group_ci <-comp_dp_distn(hdp_multi)$cred.int[[4]]

write.table(braf_group_ci, file="braf_credibility_intervals.txt",
            row.names=FALSE, sep="\t")

write.table(nras_group_ci, file="nras_group_credibility_intervals.txt",
            row.names=FALSE, sep="\t")

write.table(nf1_group_ci, file="nf1_credibility_intervals.txt",
            row.names=FALSE, sep="\t")

# match to sigs in cosmic.
cosmic.sigs <- read.table('http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt', header=TRUE, sep='\t')
#  sort by Substitution Type and Trinucleotide
 cosmic.sigs <- cosmic.sigs[order(cosmic.sigs$Substitution.Type, cosmic.sigs$Trinucleotide),]
cosmic.sigs <- as.matrix(cosmic.sigs[,grep('Signature', colnames(cosmic.sigs))])

my.sigs <- t(comp_categ_distn(hdp_multi)$mean)

cosine_matches <- apply(my.sigs, 2, function(x){
                        apply(cosmic.sigs, 2, function(y){
                            lsa::cosine(x, y)
                        })
                    })

cosine_matches[which(cosine_matches<0.85)] <- NA
round(cosine_matches, 3)

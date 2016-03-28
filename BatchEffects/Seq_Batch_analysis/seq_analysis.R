# batch effect analysis in 1k genomes pilot 3 data
# chromosome 16

# utilities
source("./BatchEffects/Seq_Batch_analysis/funcs.R")

################
# preliminaries

# get information about samples
get.batch.data()

# preprocess the count data 
preprocess()

# get exon annotation
#do.exons()

######################
# do the analysis

# load data
load("./BatchEffects/Seq_Batch_analysis/chr16_exons.rda")
load("./BatchEffects/Seq_Batch_analysis/chr16.rda")

# 10K bins, data is in 1K
library(Matrix)
ids <- as.numeric(factor(start(bins)%/%10000 + 1))

mm <- sparseMatrix(i=ids, j=1:length(ids), x=1)
tt.10k <- mm %*% tt2

bins.10k <- start(bins)[!duplicated(ids)]
bins.10k <- IRanges(start=bins.10k, width=10000)

keep <- bins.10k %in% exons
mat.10k <- as.matrix(tt.10k)[keep,]

# do only CEPH
keep <- annot$population == "CEU"
png("./BatchEffects/Seq_Batch_analysis/10k_CEPH_pvals.png", width=432, height=432, type="cairo")
.junk <- runit(mat.10k[,keep], annot[keep,],noout=TRUE,plot=TRUE)
dev.off()

# do all groups
png("10k_all_pvals.png", width=864, height=864, type="cairo")
res <- runit(mat.10k,annot,plot=TRUE)
dev.off()


# plot fstats
png("fstats.png", width=432,height=432, type="cairo")
mypar(3,1)
plot(res$f1,type="h",main="Surrogate",ylab="F-stat",ylim=c(0,500),xlab="Regions",xaxt="n")
plot(res$f2,type="h",main="PC",ylab="F-stat",xlab="Regions",xaxt="n",ylim=c(0,500))
plot(res$f3,type="h",main="Outcome",ylab="F-stat",xlab="Regions",xaxt="n")
dev.off()

# set data up for plotting
dat <- normalize.quantiles(mat.10k)
day <- as.numeric(annot$date)
o <- order(day)
sorted.mat <- dat[,o]

# do the plot
do.plot(sorted.mat, start(bins.10k), day, "seq-batch-effect")



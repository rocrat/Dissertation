get.batch.data <- function() {
  cat("Getting bamfile info\n")
  # read alignment index file
  als.index <- url("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/pilot_data.alignment.index")
  als.index <- read.delim(als.index, header=FALSE, as.is=TRUE)

  dim(als.index)
  names(als.index) <- c("alignment.filename",
                        "alignment.file.md5",
                        "study",
                        "individual",
                        "index.filename",
                        "index.file.md5",
                        "statistics.filename",
                        "statistics.file.md5")
  
  als.info <- strsplit(als.index$alignment.filename, split="\\.")

  cat("Processing bamfile info\n")
  tmp <- sapply(als.info,length)
  als.info.unmapped <- as.data.frame(do.call(rbind, als.info[tmp==7]))
  als.info <- as.data.frame(do.call(rbind,als.info[tmp==6]))
  
                                        # these are all SOLID files for some reason
# we'll add a missing aligner column
  tmp <- (als.info[,5]=="unmapped")
  als.info.unmapped2 <- als.info[tmp,]
  als.info <- als.info[!tmp,]
  
# these alignments are split by chromosome
  tmp <- (als.info.unmapped[,6] != "unmapped")
  als.info2 <- als.info.unmapped[tmp,]
  als.info.unmapped <- als.info.unmapped[!tmp,]
  
# make alignment information data frame
  names(als.info) <- c("individual",
                       "platform",
                       "aligner",
                       "study",
                       "date",
                       "format")

# create empty chromosome column
  als.info$chromosome <- ""

# add the alignments split by chromosome
  als.info2 <- als.info2[,c(1,3:7,2)]
  names(als.info2) <- names(als.info)
  als.info <- rbind(als.info, als.info2)
  
# make the unaligned read info data frame
# drop the "unmapped" column
  als.info.unmapped <- als.info.unmapped[,-6]
  names(als.info.unmapped) <- names(als.info)[-ncol(als.info)]

# add an empty aligner column to those SOLID files
  als.info.unmapped2 <- als.info.unmapped2[,c(1:2,5,3:4,6)]
  names(als.info.unmapped2) <- names(als.info.unmapped)
  als.info.unmapped2$aligner <- ""
  als.info.unmapped <- rbind(als.info.unmapped, als.info.unmapped2)

# fix the individual columns
  als.info$individual <- gsub(".*/(.*)$", "\\1", als.info$individual)
  als.info.unmapped$individual <- gsub(".*/(.*)$","\\1", als.info.unmapped$individual)

# let's check the columns should assert to false
  any(!grepl("NA", als.info$individual))
  any(!grepl("NA", als.info.unmapped$individual))

# fix the levels
  for (j in 2:ncol(als.info))
    als.info[,j] <- factor(als.info[,j])

  for (j in 2:ncol(als.info.unmapped))
  als.info.unmapped[,j] <- factor(als.info.unmapped[,j])
  
  batch.als.info <- subset(als.info, study=="SRP000033" & platform=="SLX" & aligner=="maq")

  cat("Gettign sample info\n")
# load sequence.index
  seq.index <- url("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/SRP000033.sequence.index")
  seq.index <- read.delim(seq.index, header=FALSE, as.is=TRUE)

  names(seq.index) <- c("FASTQ_FILE",
                        "MD5",
                        "RUN_ID",
                        "STUDY_ID",
                        "STUDY_NAME",
                        "CENTER_NAME",
                        "SUBMISSION_ID",
                        "SUBMISSION_DATE",
                        "SAMPLE_ID",
                        "SAMPLE_NAME",
                        "POPULATION",
                        "EXPERIMENT_ID",
                        "INSTRUMENT_PLATFORM",
                        "INSTRUMENT_MODEL",
                        "LIBRARY_NAME",
                        "RUN_NAME",
                        "RUN_BLOCK_NAME",
                        "INSERT_SIZE",
                        "LIBRARY_LAYOUT",
                        "PAIRED_FASTQ",
                        "WITHDRAWN",
                        "WITHDRAWN_DATE",
                        "COMMENT",
                        "READ_COUNT",
                        "BASE_COUNT")

  cat("Processing sample info\n")
  batch.seq.index <- subset(seq.index, SAMPLE_NAME %in% batch.als.info$individual)
  batch.seq.index <- subset(batch.seq.index, CENTER_NAME=="BI" & INSTRUMENT_PLATFORM == "ILLUMINA" & WITHDRAWN==0)
  
  keep <- which(!duplicated(batch.seq.index$RUN_ID))
  batch.seq.index <- batch.seq.index[keep,]

  cat("Preparing batch information\n")
  dates <- gsub("^.*\\.(\\d+)_SL.*$", "\\1", batch.seq.index$RUN_NAME, perl=TRUE)

  library(chron)
  dates2 <- chron(dates, format="ymd")
  wks <- cut(dates2, "weeks")
  ds <- cut(dates2, "days")

  table(table(ds))

  tmp <- table(ds)
  keep2 <- which(ds %in% names(tmp)[which(tmp>=4)])
  
  batch.seq.index2 <- batch.seq.index[keep2,]
  batch.seq.index2 <- subset(batch.seq.index2, WITHDRAWN==0)
  
  batch.als.info2 <- subset(batch.als.info, individual %in% batch.seq.index2$SAMPLE_NAME)

  dates3 <- dates2[keep2]
  ds <- cut(dates3, "days")
  
  filenames <- do.call(paste,c(batch.als.info2[,-7],sep="."))
  write(filenames,file="./BatchEffects/Seq_Batch_analysis/broad_bamfiles.txt")

  pop <- batch.seq.index2$POPULATION
  pop <- sapply(strsplit(pop, split="-"), function(x) x[1])
  pop <- gsub("\\s+", "", pop, perl=TRUE)

#   pop[pop=="CEU_1"] <- "CEPH"
#   pop[pop=="YRI_1"] <- "Yoruba"
# 
  # pop <- factor(pop, levels=c("Japanese", "HanChinese", "Yoruba", "Tuscan", "CEPH", "DenverChinese"))
  pop <- factor(pop, levels = c("JPT", "CHB", "YRI", "TSI", "CEU", "CHD"))

  png("./BatchEffects/Seq_Batch_analysis/population.png", height=432, width=864, pointsize=10, type="cairo")
  par(mar=c(5,6,0.2,0.2)+0.1,mgp=c(4,0.5,0))
  plot(as.numeric(ds), jitter(as.numeric(pop)), xlab="days", ylab="population", yaxt="n",xaxt="n")
  axis(2, las=1, at=1:length(levels(pop)), labels=levels(pop))
  axis(1, las=1, at=1:length(levels(ds)), labels=levels(ds))
  dev.off()

  save(batch.als.info2, batch.seq.index2, dates3, filenames, pop, file="./BatchEffects/Seq_Batch_analysis/broad_batches.rda")
}

preprocess <- function() {
  # read in batch descriptions
  cat("Reading batch info\n")
  load("./BatchEffects/Seq_Batch_analysis/broad_batches.rda")

  # read in count data
  cat("Reading count info\n")
  tab <- read.table("./BatchEffects/Seq_Batch_analysis/chr16_1KG_pilot3.tab")

  # keep only the runs we care about and that we actually got data for
  rid <- factor(batch.seq.index2$RUN_ID)
  m <- match(colnames(tab), levels(rid))

  keep <- !is.na(m)
  tab2 <- tab[,keep]

  # keep runs we got data for
  keep <- levels(rid)[unique(sort(Filter(function(x) !is.na(x), m)))]
  keep <- rid %in% keep

  library(chron)
  ds <- cut(dates3,"days")
  ds <- ds[keep]
  ds <- factor(ds, levels=levels(ds)[table(ds)>0])

  pop2 <- pop[keep]
  batch.seq.index3 <- batch.seq.index2[keep,]

  annot <- data.frame(runid=batch.seq.index3$RUN_ID,
                      individual=batch.seq.index3$SAMPLE_NAME,
                      population=pop2,
                      date=ds)
  date.to.keep <- aggregate(as.numeric(annot$date), annot["individual"], max)
  keep <- as.numeric(annot$date) == date.to.keep$x[match(annot$individual, date.to.keep$individual)]

  annot2 <- annot[keep,]
  annot <- annot2

  for (i in ncol(annot)) {
    if (is.factor(annot[[i]])) {
      annot[[i]] <- factor(annot[[i]])
    }
  }

  m <- match(annot$runid, colnames(tab2))
  tab2 <- tab2[,m]

  mm <- model.matrix(~-1+annot$individual)
  tt <- as.matrix(tab2)
  tt2 <- tt %*% mm
  nms <- gsub("annot\\$individual","",colnames(tt2))
  colnames(tt2) <- nms

  keep2 <- !duplicated(annot$individual)
  annot <- annot[keep2,]
  m <- match(nms, annot$individual)
  annot <- annot[m,]

  # families
  cat("Take only one individual from pedigrees\n")
  peds <- read.table("./BatchEffects/Seq_Batch_analysis/relationships_w_pops_121708.txt", header=TRUE)
  keep <- !(annot$individual %in% peds$IID)
  peds <- peds[peds$IID %in% annot$individual,]
  keep2 <- !duplicated(peds$FID)

  nms.to.keep <- union(peds$IID[keep2], annot$individual[keep])
  m <- match(nms.to.keep, annot$individual)
  annot <- annot[m,]

  m <- match(nms.to.keep, colnames(tt2))
  tt2 <- tt2[,m]

  for (i in ncol(annot)) {
    if (is.factor(annot[[i]])) {
      annot[[i]] <- factor(annot[[i]])
    }
  }

  cat("Get starting points for bins\n")
  require(IRanges)
  bins <- as.integer(gsub("\\d+_(\\d+)","\\1",rownames(tt2)))
  bins <- IRanges(start=bins+1, width=1000)

  save(tt2,annot,bins,file="./BatchEffects/Seq_Batch_analysis/chr16.rda")
}

do.exons <- function()
  {
# let's grab exons
    library(biomaRt)
    cat("Query biomart for exon annotation\n")
    mart <- useMart("ensembl", "hsapiens_gene_ensembl")
    ens <- getBM(c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", "rank", "strand"),
             filters=c("chromosome_name"), values=list("16"), mart=mart)

    cat("Create exon IRanges object\n")
    exons <- IRanges(start=ens$exon_chrom_start, end=ens$exon_chrom_end)
    save(exons, file="./BatchEffects/Seq_Batch_analysis/chr16_exons.rda")
  }

########
# analysis

f.pvalue <- function(dat,mod,mod0){
  # This is a function for performing
  # parametric f-tests on the data matrix
  # dat comparing the null model mod0
  # to the alternative model mod. 
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0,m)
  Id <- diag(n)
  
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))
  
  rss1 <- resid^2 %*% rep(1,n)
  rss0 <- resid0^2 %*% rep(1,n)
  
  fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
  p <-  1-pf(fstats,df1=(df1-df0),df2=(n-df1))
  return(list(p=p,fstats=fstats))
}


library(corpcor)
library(FactoMineR)
library(preprocessCore)

runit <- function(mat, annot, noout=FALSE, plot=FALSE) {
  browser()
  dat <- normalize.quantiles(mat)

  # this is the outcome factor
  grp <- if (noout) NULL else factor(annot$population)

  # this is our "batch" factor
  dts <- factor(annot$date)
  

  # Calculate p-values for dts
  mod <- model.matrix( ~as.factor(dts))
  mod0 <- cbind(mod[,1])
  tmp1 <- f.pvalue(dat,mod,mod0)
  pp <- tmp1$p
  pp.adj <- p.adjust(pp,method="BH")

  if (plot) {
    if (noout) {
      mypar(2,1)
    } else {
      mypar(2,2)
    }
    
    hist(pp,nc=100, main="Surrogate")
    hist(pp.adj,nc=100,main="Surrogate adjusted")
  }
  
  cat("Surrogate Q:", mean(pp.adj < 0.05), "\n")

# RV coefficient between dts model matrix and
# grp model matrix

  if (!noout) {
    ymat <- model.matrix(~ as.factor(dts))
    gmat <- model.matrix(~ as.factor(grp))
    cat("RV: ", coeffRV(ymat,gmat)$rv, "\n")
  }

# Calculate SVD

  ss <- fast.svd(t(scale(t(dat),scale=F)),tol=0)

# Calculate p-values for svs 

  mod <- model.matrix( ~ ss$v[,1:5])
  mod0 <- cbind(mod[,1])
  tmp2 <- f.pvalue(dat,mod,mod0)
  pp <- tmp2$p
  pp.adj <- p.adjust(pp,method="BH")

  cat("PCA Q: ", mean(pp.adj < 0.05), "\n")


# Find adjusted multiple R^2 for each sv with dts (take the max)

  cc <- rep(0,5)
  
  for(i in 1:5){
    cc[i] <- summary(lm(ss$v[,i] ~ as.factor(dts)))$adj.r.squared
  }

  cat("PC Surrogate R^2 max: ", max(cc), " index: ", which.max(cc), "\n")

# Find adjusted multiple R^2 for each sv with outcome (take the max)

  if(!noout) {
    cc <- rep(0,5)
    
    for(i in 1:5){
      cc[i] <- summary(lm(ss$v[,i] ~ as.factor(grp)))$adj.r.squared
    }

    cat("PC Outcome R^2 max: ", max(cc), " index: ", which.max(cc), "\n")

}

# Calculate p-values for group

  if (!noout) {
    mod <- model.matrix( ~ as.factor(grp))
    mod0 <- cbind(mod[,1])
    tmp3 <- f.pvalue(dat,mod,mod0)
    pp <- tmp3$p
    pp.adj <- p.adjust(pp,method="BH")
    
    if (plot) {
      hist(pp,nc=100,main="Outcome")
      hist(pp.adj,nc=100,main="Outcome adjusted")
    }

    cat("Outcome R^2: ", mean(pp.adj < 0.05), "\n")
  } else {
    tmp3 <- list(fstats=NULL)
  }
  return(list(f1=tmp1$fstats, f2=tmp2$fstats,f3=tmp3$fstats))
}

mypar <- function(a=1,b=1,brewer.n=8,brewer.name="Dark2",...){
  require(RColorBrewer)
  par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
  par(mfrow=c(a,b),...)
  palette(brewer.pal(brewer.n,brewer.name))
}

do.plot <- function(mat, locations, day, title) {
  mat=sorted.mat
  mat[mat>65]=NA
  mm=apply(mat,1,median,na.rm=TRUE)
  ss=apply(mat,1,mad,na.rm=TRUE);ss[ss<3] <- 3
  mat=(mat-mm)/ss
  K=2.54
  mat[mat>K]<-K;mat[mat< -K]<- -K
  Index=1:300
  day=as.numeric(sort(annot[,4]))
  dayname=gsub("day ","",sort(annot[,4]))
  tmp=table(day);tmp=as.numeric(names(tmp))[tmp>3]
  cIndex=which(day%in%tmp)

  pdf(paste(title,".pdf", sep=""),height=6,width=9)
  mypar()
  par(mar=c(2.5,4.5,1.6,1.1),mgp=c(1.5,.5,0))
  image(locations[Index],1:ncol(mat[,cIndex]),mat[Index,cIndex],col=brewer.pal(11,"RdYlBu"),yaxt="n",xlab="Genome location",ylab="")
  where=which(diff(c(day[cIndex]))!=-0)
  abline(h=where)
  where=(c(0,where)+c(where,length(cIndex)))/2 
  axis(side=2,at=where,label=unique(dayname[cIndex]),tick=FALSE,las=1)
  axis(side=2,at=119/2,label="Sample ordered by date",line=2,tick=FALSE)
  dev.off()
  
  postscript(paste(title,".eps", sep=""),height=6, width=9)
  mypar()
  par(mar=c(2.5,4.5,1.6,1.1),mgp=c(1.5,.5,0))
  image(locations[Index],1:ncol(mat[,cIndex]),mat[Index,cIndex],col=brewer.pal(11,"RdYlBu"),yaxt="n",xlab="Genome location",ylab="")
  where=which(diff(c(day[cIndex]))!=-0)
  abline(h=where)
  where=(c(0,where)+c(where,length(cIndex)))/2 
  axis(side=2,at=where,label=unique(dayname[cIndex]),tick=FALSE,las=1)
  axis(side=2,at=119/2,label="Sample ordered by date",line=2,tick=FALSE)
  dev.off()
}

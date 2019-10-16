# heterozygous and keep flags of the SNPs
procSnps <- function(rcmat, ndepth=35, het.thresh=0.25, snp.nbhd=250, gbuild="hg19", unmatched=FALSE, ndepthmax=5000) {
    # keep only chromsomes 1-22 & X for humans and 1-19, X for mice
    if (gbuild %in% c("hg19", "hg38", "hg18")) {
        chromlevels <- c(1:22,"X")
    } else {
        chromlevels <- c(1:19,"X")
    }
    chr.keep <- rcmat$Chromosome %in% chromlevels
    # keep only snps with normal read depth between ndepth and 1000
    depthN.keep <- (rcmat$NOR.DP >= ndepth) & (rcmat$NOR.DP < ndepthmax)
    # reduce the data frame to these snps
    rcmat <- rcmat[chr.keep & depthN.keep,]
    # output data frame
    out <- list()
    out$chrom <- rcmat$Chromosome
    out$maploc <- rcmat$Position
    out$rCountT <- rcmat$TUM.DP
    out$rCountN <- rcmat$NOR.DP
    # if count matrix has unmatched normal as well include it
    if (hasName(rcmat, "UMN.DP")) out$umNrCount <- rcmat$UMN.DP
    out$vafT <- 1 - rcmat$TUM.RD/rcmat$TUM.DP
    out$vafN <- 1 - rcmat$NOR.RD/rcmat$NOR.DP
    # make chromosome ordered and numeric
    out$chrom <- as.numeric(ordered(out$chrom, levels=chromlevels))
    # call a snp heterozygous if min(vafN, 1-mafN) > het.thresh
    if (unmatched) {
        if (het.thresh == 0.25) het.thresh <- 0.1
        out$het <- 1*(pmin(out$vafT, 1-out$vafT) > het.thresh & out$rCountT >= 50)
    } else {
        out$het <- 1*(pmin(out$vafN, 1-out$vafN) > het.thresh)
    }
    # scan maploc for snps that are close to one another (within snp.nbhd bases)
    # heep all the hets (should change if too close) and only one from a nbhd
    out$keep <- scanSnp(out$maploc, out$het, snp.nbhd)
    as.data.frame(out)
}

#determine sex of sample and unmatched normals based on number of chrX het SNPs
#males should not have het X
procXSnps <- function(unorms, ndepth=35, het.thresh=0.25, snp.nbhd=250, gbuild="hg19", unmatched=FALSE, ndepthmax=5000, nhet=10, normCount=NULL) {
    
    chromlevels = "X"
    chr.keep <- unorms$Chromosome %in% chromlevels
    # keep only snps with normal read depth between ndepth and 1000
    depthN.keep <- (unorms$NOR.DP >= ndepth) & (unorms$NOR.DP < ndepthmax)
    # reduce the data frame to these snps
    rcmatX <- unorms[chr.keep & depthN.keep,]
    # output data frame
    out <- list()
    out$chrom <- rcmatX$Chromosome
    out$maploc <- rcmatX$Position
    
    
    out$rCountT <- rcmatX$TUM.DP
    out$rCountN <- rcmatX$NOR.DP
    out$vafT <- 1 - rcmatX$TUM.RD/rcmatX$TUM.DP
    out$vafN <- 1 - rcmatX$NOR.RD/rcmatX$NOR.DP
    out = as.data.frame(out)
    
    for(i in 3:(normCount+2)){
        tempVAF = paste('File', i, "VAF", sep="")
        tempR = paste("File", i, "R", sep="")
        tempDP = paste("File", i, "DP", sep="")
        tempHET = paste("File", i, "DPhet", sep="")
        out[,tempVAF] = 1 - (rcmatX[,tempR]/rcmatX[,tempDP])
        out[,tempHET] =  1*(pmin(out[,tempVAF], 1-out[,tempVAF]) > het.thresh )
    }
    out$NOR.DPhet <- 1*(pmin(out$vafN, 1-out$vafN) > het.thresh)
    
    out.hets = out[,grep("het", colnames(out))]
    
    out.hets = as.data.frame(colSums(out.hets, na.rm = T))
    colnames(out.hets) = "numHet"  
    out.hets$sampleSex = ifelse(out.hets$numHet>nhet,"Female", "Male") 
    rownames(out.hets) = gsub("het", "", rownames(out.hets))
    
    out.hets
}

scanSnp <- function(maploc, het, nbhd) {
    n <- length(maploc)
    zzz <- .Fortran("scansnp",
                    as.integer(n),
                    as.double(maploc),
                    as.double(het),
                    keep=double(n),
                    as.double(nbhd))
    zzz$keep
}

# obtain logR and logOR from read counts and GC-correct logR
counts2logROR <- function(mat, gbuild, unmatched=FALSE, MandUnormal=FALSE, f, spanT, spanA, spanX) {
    out <- mat[mat$keep==1,]
    #out$chrom = gsub('X', '23', out$chrom) #testing replace X with 23

    # gc percentage
    out$gcpct <- rep(NA_real_, nrow(out))
    # get GC percentages from pctGCdata package
    # loop thru chromosomes
    #nchr <- max(mat$chrom) # IMPACT doesn't have X so only 22

    for (i in c(1:23)) {
   # for (i in 1:nchr) {
        ii <- which(out$chrom==i)
        # allow for chromosomes with no SNPs i.e. not targeted
        if (length(ii) > 0) {
            out$gcpct[ii] <- getGCpct(i, out$maploc[ii], gbuild)
        }

    }
    out = out[which(!is.na(out$gcpct)),]
    x.idx <- grep("X|23",out$chrom)
    ##### log-ratio with gc correction and maf log-odds ratio steps
    chrom <- out$chrom
    maploc <- out$maploc
    if (hasName(mat, "umNrCount")) {
        rCountN <- out$umNrCount
    } else {
        rCountN <- out$rCountN
    }
    rCountT <- out$rCountT
    vafT <- out$vafT
    vafN <- out$vafN
    het <- out$het
    gcpct <- out$gcpct
    gcpct.auto = gcpct[-x.idx]
    gcpct.x = gcpct[x.idx]

    # compute gc bias
    ncount <- tapply(rCountN, gcpct, sum)
    tcount <- tapply(rCountT, gcpct, sum)
    pctgc <- as.numeric(names(ncount))
    tscl <- sum(ncount)/sum(tcount)
    gcb <- lowess(pctgc, log2(tcount*tscl)-log2(ncount), f=f)
    jj <- match(gcpct, gcb$x)
    gcbias <- gcb$y[jj]
    # compute cn log-ratio (gc corrected) and baf log odds-ratio
    #####################################
    #square root transform count vectors.
    tumor_sqrt = sqrt(rCountT)
    tumor_sqrt.auto = tumor_sqrt[-x.idx]
    tumor_sqrt.x    = tumor_sqrt[x.idx]

    normal_sqrt      = sqrt(rCountN)
    normal_sqrt.auto = normal_sqrt[-x.idx]
    normal_sqrt.x    = normal_sqrt[x.idx]

    #loess regression for lr tumor autosomes and X seperately.
    loess_tumor.auto = lowess(gcpct.auto,tumor_sqrt.auto, f=spanT) #need to change this to input value from span.fits
    jj=match(gcpct.auto, loess_tumor.auto$x)
    fit<-loess_tumor.auto$y[jj]
    #loess_tumor.auto <-loess(tumor_sqrt.auto~gcpct.auto,span=f);
   # temp<-predict(loess_tumor.auto);
    normalized_t.auto<-(tumor_sqrt.auto-fit+median(tumor_sqrt.auto))/(median(tumor_sqrt.auto[which(tumor_sqrt.auto != 0)]));

    loess_tumor.x <-lowess(gcpct.x,tumor_sqrt.x,f=spanT);
    jj=match(gcpct.x, loess_tumor.x$x)
    fit<-loess_tumor.x$y[jj]
    #temp<-predict(loess_tumor.x);
    normalized_t.x<-(tumor_sqrt.x-fit+median(tumor_sqrt.x))/(median(tumor_sqrt.x[which(tumor_sqrt.x != 0)]));

    tumor_rt = c(normalized_t.auto, normalized_t.x)

    #loess regression for lr normal autosomes and X seperately.
    loess_normal.auto <-lowess(gcpct.auto, normal_sqrt.auto,f=spanA);
   # temp<-predict(loess_normal.auto);
    jj=match(gcpct.auto, loess_normal.auto$x)
    fit<-loess_normal.auto$y[jj]
    normalized_n.auto<-(normal_sqrt.auto-fit+median(normal_sqrt.auto))/(median(normal_sqrt.auto[which(normal_sqrt.auto != 0)]));

    loess_normal.x <-lowess(gcpct.x,normal_sqrt.x,f=spanX);
    jj=match(gcpct.x, loess_normal.x$x)
    fit<-loess_normal.x$y[jj]
    #temp<-predict(loess_normal.x);
    normalized_n.x<-(normal_sqrt.x-fit+median(normal_sqrt.x))/(median(normal_sqrt.x[which(normal_sqrt.x != 0)]));

    normal_rt = c(normalized_n.auto, normalized_n.x)

    #calculate log2 ratios
    cnlr = log2(tumor_rt) - log2(normal_rt)

    #####################################
    #use old method of cnlr calc if matched normal
    if (!MandUnormal) cnlr <- log2(1+rCountT*tscl) - log2(1+rCountN) - gcbias
    # minor allele log-odds ratio and weights
    rCountN <- out$rCountN # reset normal depth in case umNrCount exists
    lorvar <- valor <- rep(NA_real_, length(maploc))
    if (unmatched) {
        # read count matrix for odds ratio etc
        rcmat <- round(cbind(vafT[het==1]*rCountT[het==1], (1-vafT[het==1])*rCountT[het==1]))
        # folded log of Tukey (with 1/6 correction)
        valor[het==1] <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
        # variance - approximation using delta method
        lorvar[het==1] <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
    } else {
        # read count matrix for odds ratio etc
        rcmat <- round(cbind(vafT[het==1]*rCountT[het==1], (1-vafT[het==1])*rCountT[het==1], vafN[het==1]*rCountN[het==1], (1-vafN[het==1])*rCountN[het==1]))
        # log-odds-ratio (Haldane correction)
        valor[het==1] <- log(rcmat[,1]+0.5) - log(rcmat[,2]+0.5) - log(rcmat[,3]+0.5) + log(rcmat[,4]+0.5)
        # variance of log-odds-ratio (Haldane; Gart & Zweifel Biometrika 1967)
        lorvar[het==1] <- (1/(rcmat[,1]+0.5) + 1/(rcmat[,2]+0.5) + 1/(rcmat[,3]+0.5) + 1/(rcmat[,4]+0.5))
    }
    # put them together
    out$lorvar <- out$valor <- out$cnlr <- out$gcbias <- rep(NA_real_, nrow(out))
    out$gcbias <- gcbias
    out$cnlr <- cnlr
    out$valor <- valor
    out$lorvar <- lorvar
    out
}

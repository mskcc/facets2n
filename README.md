# facets2n
Algorithm to implement Fraction and Allelic Copy number Estimate from Tumor/normal Sequencing. Package created to test FACETS with both matched and unmatched normal.

This implementation of FACETS requires a tumor sample, matched normal and user defined number of unmatched normal samples. 

- :sparkles: The normal sample with computed minimal noise is selected for copy number log2 ratio calculations, while the matched normal is always selected for VALOR calculation. 
      *At present, autosomes and chromosome X are normalized separately.*

- :construction: (beta) infer sex from matched normal sample, and select an unmatched normal with same sex for chrX normalization.

## Install R library

```
devtools::install_github("rptashkin/facets2n")
```

## snp-pileup command options used for IMPACT data
*hat tip kpjonsson: https://github.com/taylor-lab/facets-suite*
```
inst/extcode/snp-pileup-wrapper.R \
  --snp-pileup-path <optional, path to snp-pileup executable, defaults to snp-pileup in your PATH> \
  --vcf-file <path to SNP VCF> \
  --normal-bam normal.bam \
  --tumor-bam tumor.bam \
  --output-prefix <prefix for output file, preferrably tumorSample__normalSample>
  --unmatched-normal-BAMS <if using unmatched BAMs for logR normalization, e.g. "<some/path/>*-N.bam" >
```

### Example generation of the input counts file for Reference normals using snp-pileup. (Recommended)
*We suggest at least 5 male and 5 female reference normal BAMs, that were processed in the laboratory and with an analysis pipeline using the same parameters of the tumor sample you are analyzing. For hybridization capture sequencing, we have found that a pooled sample of non-neoplastic tissue (e.g. blood from healthy donors) from 10 individuals captured and sequenced together with the data being analyzed often outperforms the matched normal and individual reference normals for minimizing noise in log ratio plots.*
```
inst/extcode/snp-pileup-wrapper.R --output-prefix standard_normals_cv3heme  \
  --vcf-file dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT_NOCHR.vcf \
  --unmatched-normal-BAMS "<some/path_to_bam_directory>/*-N*.bam"
```

### Example generation of the input counts file for a tumor-normal pair + sequencing batch control
*Multiple BAM files can be suppplied as a quoted, space seperated, string to --unmatched-normal-BAMS*
```
inst/extcode/snp-pileup-wrapper.R --output-prefix P-0029502_NN_TH_PoolNormal  
  --vcf-file dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT_NOCHR.vcf
  --normal-bam normal.bam
  --tumor-bam tumor.bam
  --unmatched-normal-BAMS "<some/path>/PoolNormal.bam"
```

## Run FACETS with unmatched normal samples for CNLR

#### There are two ways to run the unmatched analyis. The reference normal loess data can be generated independently (as above, with non-pooled samples) and used for subsequent analysis as shown below. This will substantially decrease analysis time if you will be using the same set of unmatched normals for multiple analyses. 
*The runtime of readSnpMatrix() increases with number of unmatched normal samples*

#### Generating the reference loess files will take several minutes, but needs to be performed just once for a given sequencing panel.
```
library(facets2n)
MakeLoessObject(pileup = PreProcSnpPileup(filename = "tests/standard_normals_cv3heme.snp_pileup.gz", is.Reference = TRUE), write.loess = TRUE, outfilepath = "tests/standard_normals_cv3heme.loess.txt", is.Reference = TRUE)
```
#### including the argument refX=TRUE in call to readSnpMatrix() will force selection of a individual reference sample from chrX normalization (preffered). 
```
readu <- readSnpMatrix(filename = "tests/P-0029502_NN_TH_PoolNormal.snp_pileup.gz",
  MandUnormal = TRUE,
  ReferencePileupFile = "tests/standard_normals_cv3heme.snp_pileup.gz",
  ReferenceLoessFile = "tests/standard_normals_cv3heme.loess.txt",
  useMatchedX = FALSE,
  refX=TRUE)
```

#### Alternatively, the tumor, normal, and reference unmatched normals pileup data can all be provided in a single counts file and loess normalization for all samples will take place at run-time.
```
readu <- readSnpMatrix(
  filename = "tests/P-0029502_NN_TH_PoolNormal.snp_pileup.gz",
  MandUnormal = TRUE,
  useMatchedX = FALSE)

```

#### The following steps are applicable to either of approaches above.
```
xx <- preProcSample(readu$rcmat, unmatched = F,
  ndepth = 50,het.thresh = 0.25, ndepthmax = 5000,
  spanT = readu$spanT, spanA=readu$spanA, spanX = readu$spanX,
  MandUnormal = TRUE)
```
```
oo <- procSample(xx,min.nhet = 10, cval = 150)
```
```
dlr <- oo$dipLogR
```
```
oo <- procSample(xx,min.nhet = 10, cval = 50, dipLogR = dlr)
fit <- emcncf(oo, min.nhet = 10)
```
```
png(filename = "tests/P-0029502_unmatched_CNLR.png",width = 4, height = 6, units = "in",res = 500)
plotSample(x=oo,emfit=fit, plot.type = "both")
dev.off()
```

### User still has the ability to run FACETS with the matched normal, if preferred.
see readSnpMatrix() argument  ```MandUnormal = FALSE```
```
library(facets2n)
readm <-  readSnpMatrix(filename = "tests/P-0029502_NN_TH_PoolNormal.snp_pileup.gz",MandUnormal = FALSE)
```
```
xx <- preProcSample(readm, unmatched = F,ndepth = 50,het.thresh = 0.25,ndepthmax = 5000, MandUnormal = FALSE)
oo=procSample(xx,min.nhet = 10, cval = 150)
```
```
dlr <- oo$dipLogR
```
```
oo <- procSample(xx,min.nhet = 10, cval = 50, dipLogR = dlr)
fit <- emcncf(oo, min.nhet = 10)
```
```
png(filename = "tests/P-0029502_matched_CNLR.png",width = 4, height = 6, units = "in",res = 500)
plotSample(x=oo,emfit=fit, plot.type = "both")
dev.off()
```

Results using matched normal for CNLR                     |  Results using unmatched normal for CNLR
:--------------------------------------------------------:|:------------------------------------------------------------:
![matched normal cnlr](/tests/P-0029502_matched_CNLR.png) | ![unmatched normal cnlr](/tests/P-0029502_unmatched_CNLR.png)



### Starting with v0.3.0, allele specific copy number with transplant cases is possible by including a seperate counts matrix for donor sample(s). Requires a baseline host sample as matched normal (e.g. nails or other source of non-neoplastic cells). 
```
library(facets2n)
#read counts matrix for tumor and baseline host sample
readu = readSnpMatrix(filename = "tests/P-0023145-TH_NN.snp_pileup.gz", MandUnormal = TRUE, ReferencePileupFile ="tests/standard_normals_cv3heme.snp_pileup.gz", ReferenceLoessFile = "tests/standard_normals_cv3heme.loess.txt", useMatchedX = FALSE, refX = TRUE)
```
```
#read counts matrix for donor sample
readonor = readSnpMatrix(filename = "tests/P-0023145-ND.snp_pileup.gz", donorCounts = TRUE)
```

```
# preprocess sample and limit hets to those that are het in both host and donor
xx <- preProcSample(readu$rcmat, unmatched = F,
                    ndepth = 50,het.thresh = 0.3, ndepthmax = 5000,
                    spanT = readu$spanT, spanA=readu$spanA, spanX = readu$spanX,
                    MandUnormal = TRUE, donorCounts = readonor)
```

```
oo <- procSample(xx,min.nhet = 10, cval = 150)
dlr <- oo$dipLogR
oo <- procSample(xx,min.nhet = 10, cval = 50, dipLogR = dlr)
fit <- emcncf(oo, min.nhet = 10)

png(filename = "tests/P-0023145_host_donor_snps.png",width = 4, height = 6, units = "in",res = 500)
plotSample(x=oo,emfit=fit, plot.type = "both")
dev.off()
```

```
# running without baseline donor sample produces inaccurate logOR. allele specific copy number is not possible.
xx <- preProcSample(readu$rcmat, unmatched = F,
                    ndepth = 50,het.thresh = 0.25, ndepthmax = 5000,
                    spanT = readu$spanT, spanA=readu$spanA, spanX = readu$spanX,
                    MandUnormal = TRUE)

oo <- procSample(xx,min.nhet = 10, cval = 150)
dlr <- oo$dipLogR
oo <- procSample(xx,min.nhet = 10, cval = 50, dipLogR = dlr)
fit <- emcncf(oo, min.nhet = 10)

png(filename = "tests/P-0023145_only_host_snps.png",width = 4, height = 6, units = "in",res = 500)
plotSample(x=oo,emfit=fit, plot.type = "both")
dev.off()
```
Results using baseline host and donor samples             |  Results using baseline host sample only
:--------------------------------------------------------:|:------------------------------------------------------------:
![host and donor](/tests/P-0023145_host_donor_snps.png) | ![host only](/tests/P-0023145_only_host_snps.png)



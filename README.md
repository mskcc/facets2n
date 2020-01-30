# facets2n
Algorithm to implement Fraction and Allelic Copy number Estimate from Tumor/normal Sequencing. Package created to test FACETS with both matched and unmatched normal.

This implementation of FACETS requires a tumor sample, matched normal and user defined number of unmatched normal samples. 

- [x] The normal sample with computed minimal noise is selected for copy number log2 ratio calculations, while the matched normal is always selected for VALOR calculation. 

*At present, autosomes and chromosome X are normalized separately.*

- [x] Planned: infer sex from matched normal sample, and select an unmatched normal with same sex for chrX normalization.

## Install R library

```
devtools::install_github("rptashkin/facets2n")
```

## snp-pileup command options used for IMPACT data
*hat tip kpjonsson: https://github.com/taylor-lab/facets-suite*
```
R/snp-pileup-wrapper.R \
  --snp-pileup-path <optional, path to snp-pileup executable, defaults to snp-pileup in your PATH> \
  --vcf-file <path to SNP VCF> \
  --normal-bam normal.bam \
  --tumor-bam tumor.bam \
  --output-prefix <prefix for output file, preferrably tumorSample__normalSample>
  --unmatched_normal_BAMS <if using unmatched BAMs for logR normalization, e.g. "<some/path/>*-N.bam" >
```

### Example generation of the input counts file for reference normals using snp-pileup.
*We recommend at least 5 male and 5 female reference normal BAMs, that were processed in the laboratory and analysis pipeline using the same parameters of the tumor sample you are analyzing*
```
R/snp-pileup-wrapper.R --output-prefix standard_normal_heme_bams  \
  --vcf-file dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT_NOCHR.vcf \
  --unmatched-normal-BAMS "<some/path_to_bam_directory>/*-N*.bam"
```

### Example generation of the input counts file for a tumor-normal pair + sequencing batch control
*Multiple BAM files can be suppplied as a quoted, space seperated, string to --unmatched-normal-BAMS*
```
R/snp-pileup-wrapper.R --output-prefix standard_normal_bams  
  --vcf-file dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT_NOCHR.vcf
  --normal-bam normal.bam
  --tumor-bam tumor.bam
  --unmatched-normal-BAMS "<some/path>/PoolNormal.bam"
```

## Run FACETS with matched normal
```
library(facets2n)
readm <-  readSnpMatrix(filename = "tests/countsMerged_uNormals_P-0029502.dat.gz",MandUnormal = FALSE)
```
```
xx <- preProcSample(readm, unmatched = F,ndepth = 10,het.thresh = 0.25,ndepthmax = 5000, MandUnormal = FALSE)
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

## Run FACETS with unmatched normal samples for CNLR

#### There are two ways to run the unmatched analyis. The reference normal loess data can be generated independently (as above) and used for subsequent analysis as shown below. This will substantially decrease analysis time if you will be using the same set of unmatched normals for multiple analyses.
*The runtime of readSnpMatrix() increases with number of unmatched normal samples*

```
library(facets2n)
MakeLoessObject(pileup = PreProcSnpPileup(filename = "tests/standard_normal_heme_bams.snp_pileup.gz", is.TandMN = FALSE), write.loess = TRUE, outfilepath = "tests/standard_normal_heme_bams.loess.txt", is.TandMN = FALSE)
```
```
readu <- readSnpMatrix(filename = "tests/countsMerged_MatchedNormalandTumor_P-0029502.dat.gz",
  MandUnormal = TRUE,
  ReferencePileupFile = "tests/standard_normal_heme_bams.snp_pileup.gz",
  ReferenceLoessFile = "tests/standard_normal_heme_bams.loess.txt",
  useMatchedX = FALSE)
```

#### Alternatively, the tumor, normal, and reference unmatched normals pileup data can all be provided in a single file and loess normalization for all samples will take place at run-time.
```
readu <- readSnpMatrix(
  filename = "tests/countsMerged_uNormals_P-0029502.dat.gz",
  MandUnormal = TRUE,
  useMatchedX = FALSE)

```

#### The following steps are applicable to either of approaches above.
```
xx <- preProcSample(readu$rcmat, unmatched = F,
  ndepth = 10,het.thresh = 0.25, ndepthmax = 5000,
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

Results using matched normal for CNLR                     |  Results using unmatched normal for CNLR
:--------------------------------------------------------:|:------------------------------------------------------------:
![matched normal cnlr](/tests/P-0029502_matched_CNLR.png) | ![unmatched normal cnlr](/tests/P-0029502_unmatched_CNLR.png)

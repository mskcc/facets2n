# facets2n
Algorithm to implement Fraction and Allelic Copy number Estimate from Tumor/normal Sequencing. Package created to test FACETS with both matched and unmatched normal.

This implementation of FACETS requires a tumor sample, matched normal and user defined number of unmatched normal samples. 

- [x] The normal sample with computed minimal noise is selected for copy number log2 ratio calculations, while the matched normal is always selected for VALOR calculation. 

*At present, autosomes and chromosome X are normalized separately.*

- [ ] Planned: infer sex from matched normal sample, and select an unmatched normal with same sex for chrX normalization.

## Install R library

```
devtools::install_github("rptashkin/facets2n")
```

## snp-pileup command options used for IMPACT data:
```
snp-pileup
 -g
 -p
 -A
 -d 20000
 -r 10,0,10,...10 (10, for each unmatched normal)
 -q 0
 -Q 0
 -v <dbsnp VCF>

<counts file>
 <BAM files, in order: matched normal, tumor, unmatched normals>
```

## Example generation of the input counts file using snp-pileup (required) and 18 assay specific unmatched diploid normals and 1 batch Pooled Normal:
```
snp-pileup -g -p -A -d 20000 -r 10,0,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10
-q 0 -Q 0 -v dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT_NOCHR.vcf countsMerged_uNormals_P-0029502.dat.gz
P-0030397-NN.bam P-0029502-TH.bam <path_to_assay_specific_unmatched_diploid_normals>/*bam PoolNormal.bam
```

## Run FACETS with matched normal
```
library(facets2n)
readm =  readSnpMatrix(filename = "tests/countsMerged_uNormals_P-0029502.dat.gz",MandUnormal = FALSE)

xx = preProcSample(readm, unmatched = F,ndepth = 10,het.thresh = 0.25,ndepthmax = 5000, MandUnormal = FALSE)
oo=procSample(xx,min.nhet = 10, cval = 150)

dlr = oo$dipLogR

oo=procSample(xx,min.nhet = 10, cval = 50, dipLogR = dlr)
fit=emcncf(oo, min.nhet = 10)

png(filename = "tests/P-0029502_matched_CNLR.png",width = 4, height = 6, units = "in",res = 500)
plotSample(x=oo,emfit=fit, plot.type = "both")
dev.off()
```

## Run FACETS with unmatched normal samples for CNLR
*The runtime of readSnpMatrix() increases with number of unmatched normal samples*
```
library(facets2n)

# There are two ways to run the unmatched analyis. The reference normal loess data can be generated independently and used for subsequent analysis as shown below. This will substantially decrease analysis time if you will be using the same set of unmatched normals for multiple analyses.

MakeLoessObject(PreProcSnpPileup("./tests/countsMerged_uNormals.dat.gz", is.TandMN = FALSE), write.loess = TRUE, outfilepath = "./tests/countsMerged_uNormals.loess.txt", is.TandMN = FALSE)

readu <- readSnpMatrix("./tests/countsMerged_MatchedNormalandTumor_P-0029502.dat.gz", MandUnormal = TRUE, ReferencePileupFile = "./tests/countsMerged_uNormals.dat.gz", ReferenceLoessFile = "./countsMerged_uNormals.loess.txt", matchedNormalforX = FALSE)

# Alternatively, the tumor, normal, and reference unmatched normals pileuo data can all be provided in a single file for one-time analysis.

readu = readSnpMatrix(filename = "./tests/countsMerged_uNormals_P-0029502.dat.gz", MandUnormal = TRUE, matchedNormalforX = FALSE)

# The following steps are applicable to either of approaches above.

xx = preProcSample(readu$rcmat, unmatched = F,ndepth = 10,het.thresh = 0.25,ndepthmax = 5000,spanT = readu$spanT, spanA=readu$spanA, spanX = readu$spanX, MandUnormal = TRUE)
oo=procSample(xx,min.nhet = 10, cval = 150)

dlr = oo$dipLogR

oo=procSample(xx,min.nhet = 10, cval = 50, dipLogR = dlr)
fit=emcncf(oo, min.nhet = 10)

png(filename = "tests/P-0029502_unmatched_CNLR.png",width = 4, height = 6, units = "in",res = 500)
plotSample(x=oo,emfit=fit, plot.type = "both")
dev.off()
```

Results using matched normal for CNLR                     |  Results using unmatched normal for CNLR
:--------------------------------------------------------:|:------------------------------------------------------------:
![matched normal cnlr](/tests/P-0029502_matched_CNLR.png) | ![unmatched normal cnlr](/tests/P-0029502_unmatched_CNLR.png)

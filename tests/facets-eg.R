######################################################################
# Code to try out code coverage
######################################################################

library("facets2n")
# check if .Random.seed exists
seedexists <- exists(".Random.seed")
# save seed
if(seedexists) oldSeed <- .Random.seed
# Alway use the same random seed
set.seed(0xfade)

# running the stomach example in the vignette 
datafile = system.file("extdata", "stomach.csv.gz", package="facets2n")
# read the data
rcmat = readSnpMatrix(datafile)

# fit segmentation tree
xx = preProcSample(rcmat)
# estimate allele specific copy numbers
oo=procSample(xx,cval=150)
# EM fit version 1
fit=emcncf(oo)
# EM fit version 2
fit2=emcncf2(oo)
# finished

# Reset to previous random seed
if(seedexists) .Random.seed <- oldSeed
set.seed(0xfade)

#test transplant case with baseline donor sample
transplantFile = system.file("extdata", "transplant.csv.gz", package="facets2n")
normalsFile = system.file("extdata", "test_standard_normals.csv.gz", package="facets2n")
loessFile = system.file("extdata", "test_standard_normals.loess.txt", package="facets2n")

readu = readSnpMatrix(transplantFile, MandUnormal = TRUE, ReferencePileupFile = normalsFile, ReferenceLoessFile = loessFile, useMatchedX = FALSE, refX = TRUE)

#read counts matrix for donor sample
donorFile = system.file("extdata", "donor.csv.gz", package="facets2n")
readonor = readSnpMatrix(donorFile, donorCounts = TRUE)


# preprocess sample and limit hets to those that are het in both host and donor
xxx <- preProcSample(readu$rcmat, unmatched = F,
                    ndepth = 50,het.thresh = 0.3, ndepthmax = 5000,
                    spanT = readu$spanT, spanA=readu$spanA, spanX = readu$spanX,
                    MandUnormal = TRUE, donorCounts = readonor)

ooo <- procSample(xxx,min.nhet = 10, cval = 150)
dlr <- ooo$dipLogR
ooo <- procSample(xx,min.nhet = 10, cval = 50, dipLogR = dlr)
fit <- emcncf(ooo, min.nhet = 10)

# Reset to previous random seed
if(seedexists) .Random.seed <- oldSeed
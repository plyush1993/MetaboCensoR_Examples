# Scalability Testing Notes

For scalability testing [`mzrtsim`](https://github.com/yufree/mzrtsim/tree/master) R package was used to generate simulated datasets. <br>
Testing was performed in [`R 4.5.0`](https://cran.r-project.org/bin/windows/base/old/4.5.0/); [session info](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/scalability%20testing/session%20info%20for%20scalability.txt). 

### Laptop Testing
Code for generating simulated dataset:
```r
library(mzrtsim)
library(tidyverse)
data("monams1")

simdata <- mzrtsim(ncomp = 100, ncond = 2, ncpeaks = 0.05,
   nbatch = 10, nbpeaks = 0.1, npercond = 1000, nperbatch = rep(200, 10), seed = 42, batchtype = 'mb', db=monams1)
mtr <- cbind(Feature = 1:length(simdata[["mz"]]), mz = simdata[["mz"]], rt = simdata[["rt"]], simdata[["data"]]) %>% as.data.frame()
write_csv(mtr, "test 2k*2k.csv")
```
**Hardware:** 11th Gen Intel® Core™ i5-11300H @ 3.10 GHz, 4 cores / 8 logical processors, 16 GB RAM. <br>
**Input dataset:** 69 MB, 1921 peaks * 2000 samples.<br>
**Applied Filters:** MS filters (Isotopes, Adducts, NLs, ISFs), and Peak Filters (mz, rt, RMD, AMD). `Limit table previews` enabled.<br>
**Processing Time:** 1.30 min.<br>

### Desktop Testing
Code for generating simulated dataset:
```r
library(mzrtsim)
library(tidyverse)
data("monams1")

simdata <- mzrtsim(ncomp = 150, ncond = 2, ncpeaks = 0.05,
  nbatch = 10, nbpeaks = 0.1, npercond = 2500, nperbatch = rep(500, 10), seed = 42, batchtype = 'mb', db=monams1)
mtr <- cbind(Feature = 1:length(simdata[["mz"]]), mz = simdata[["mz"]], rt = simdata[["rt"]], simdata[["data"]]) %>% as.data.frame()
write_csv(mtr, "test 6k*5k.csv")
```
**Hardware:** Intel(R) Xeon(R) W-3225 CPU @ 3.70GHz, 8 cores / 16 logical processors, 64 GB RAM. <br>
**Input dataset:** 576 MB, 6426 peaks * 5000 samples.<br>
**Applied Filters:** MS filters (Isotopes, Adducts, NLs, ISFs), and Peak Filters (mz, rt, RMD, AMD). `Limit table previews` enabled.<br>
**Processing Time:** 30 min.<br>

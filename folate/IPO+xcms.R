#...........................................................
# IPO processing ----
#...........................................................

library(IPO)
library(xcms)
library(doParallel)
library(BiocParallel)
library(stringr)

setwd("C:/../") # main folder

# select files
wd_1 <- c("C:/../") # folder with mzXML files 
files <- list.files(wd_1, recursive = T, full.names = T, pattern = ".mzXML") 
files_QC <- files[which(str_detect(files, "QC_pool_n"))] # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.

# NEG 
val <- "NEG_"
files_QC <- files[which(str_detect(files, "QC_pool_n"))]
files_QC <- files_QC[which(str_detect(files_QC, val))]

nCore=detectCores()-1

# PeakPickingParameters (adjust to your data)
peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$noise=c(100,1000)
peakpickingParameters$value_of_prefilter=c(3,800)
peakpickingParameters$min_peakwidth<- c(3,15)
peakpickingParameters$max_peakwidth<- c(30,40)
peakpickingParameters$ppm<- c(3,20)
param=SnowParam(workers = nCore)
nSlaves=1 # or nCore or more

resultPeakpicking <-
  optimizeXcmsSet(files = files_QC,
                  params = peakpickingParameters,
                  BPPARAM = param,
                  nSlaves = nSlaves,
                  subdir = NULL,
                  plot = TRUE)

optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset
#
save(resultPeakpicking, file = "IPO_optimiz_xcms_7QC_NEG.RData")

# Retention Time Alignment Optimization (adjust to your data)
retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- c(0.33,1)
retcorGroupParameters$gapExtend <- c(2.0,3.0)
retcorGroupParameters$minfrac=c(0.2,0.9)
retcorGroupParameters$response=c(9.0,18.0)
retcorGroupParameters$gapInit=c(0.10, 0.50)
retcorGroupParameters$mzwid=c(0.0010, 0.050)
BiocParallel::register(BiocParallel::SerialParam())
nSlaves=1 # or nCore or more

resultRetcorGroup <-
  optimizeRetGroup(xset = optimizedXcmsSetObject,
                   params = retcorGroupParameters,
                   nSlaves = nSlaves,
                   subdir = NULL,
                   plot = TRUE)

resultRetcorGroup$best_settings
#
save(resultRetcorGroup, file = "IPO_optimiz_ret_align_7QC_NEG.RData")

# Get final results
# NOTE!!! It is better to optimize bw and minfrac manually to your final peak table after applying IPO optimized parameters
writeRScript(resultPeakpicking$best_settings$parameters, resultRetcorGroup$best_settings)
param <- c(resultPeakpicking$best_settings$parameters, resultRetcorGroup$best_settings)
#
save(param,file = 'all params IPO 7QC NEG.RData')

funs_params <- capture.output(writeRScript(resultPeakpicking$best_settings$parameters, resultRetcorGroup$best_settings), type = "message")
#
save(funs_params,file = 'funs params IPO 7QC NEG.RData')

#...........................................................
# xcms processing ----
#...........................................................

# Settings 
library(xcms)
library(data.table)
library(dplyr)
library(stringr)
library(BiocParallel)
library(doParallel)

setwd("C:/...")

# select files
wd_1 <- c("C:/.../") # folder with mzXML files 
files <- list.files(wd_1, recursive = T, full.names = T, pattern = ".mzXML") 
files_all <- files
rname <- files_all # obtain all info from rownames
rname <- str_remove(rname, wd_1) # remove folder name
rname <- str_remove(rname, ".mzXML") # remove some pattern from vendor-specific format
all_id <- base::lapply(1:length(rname), function(y) unlist(str_split(rname[y], "_"))) # split info from rownames
all_id <- as.data.frame(do.call(rbind, all_id))
ro_id <- all_id[,2] # obtain run order ID (every [2] element)
files_all_df <- as.data.frame(cbind(files_all, ro = as.numeric(ro_id)))
files_all_df <- files_all_df[order(as.numeric(files_all_df$ro), decreasing = F),] # sort by run order
files_all <- files_all_df[,-2]

# NEG 
val <- "NEG_"
files_all <- files_all[which(str_detect(files_all, val))]

wd_1 <- c("C:/.../") # folder with mzXML files 
files <- list.files(wd_1, recursive = T, full.names = T, pattern = ".mzXML") 
files_QC <- files[which(str_detect(files, "QC_pool_n"))] # adjust to your data format: ".CDF", ".mzXML", ".mzML", etc.
val <- "NEG_"
files_QC <- files[which(str_detect(files, "QC_pool_n"))]
files_QC <- files_QC[which(str_detect(files_QC, val))]

center_NEG <- files_QC[resultRetcorGroup$best_settings$center]
center_NEG <- which(files_all == center_NEG)

# Create a phenodata data.frame
pd <- data.frame(sample_name = sub(basename(files_all), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE), stringsAsFactors = FALSE) # download filenames

rname <- pd
all_id <- as.data.frame(sapply(1:nrow(rname), function(y) unlist(str_split(rname[y,], "_")))) # split info from rownames
all_id1 <- as.numeric(all_id[2,]) # as numeric run order in [[x]][1]

p1_id <- unlist(lapply(all_id, function(y) unlist(y[3]))) # obtain patient ID (every [3] element)
n_gr_t <- data.frame(p1_id)

vec_gr <- as.numeric(as.factor(n_gr_t$p1_id))
sample_gr <- unique(as.numeric(as.factor(n_gr_t$p1_id)))
n_gr <- sapply(1:length(sample_gr), function(y) length(vec_gr[vec_gr == y]))
min_frac_man <- min(round(n_gr/length(vec_gr), 2)) # calculate min_frac manually

# download files
raw_data <- readMSData(files = files_all, pdata = new("NAnnotatedDataFrame", n_gr_t), mode = "onDisk") # or use pd only as: pdata = new("NAnnotatedDataFrame", pd) or pdata = new("NAnnotatedDataFrame", n_gr_t)

# parallel processing
cores = detectCores()-1
register(bpstart(SnowParam(cores)))
BiocParallel::register(BiocParallel::SerialParam())

# load best parameters from IPO:
load("C:/.../IPO_optimiz_xcms_7QC_NEG.RData")
load("C:/.../IPO_optimiz_ret_align_7QC_NEG.RData")

# NOTE!!! You could also set all parameters manually
# feature detection
cwp <- xcms::CentWaveParam(ppm = resultPeakpicking$best_settings$parameters$ppm, # maximal tolerated m/z deviation in consecutive scans, in ppm
                           peakwidth = c(resultPeakpicking$best_settings$parameters$min_peakwidth, resultPeakpicking$best_settings$parameters$max_peakwidth), # min/max chromatographic peak width
                           snthresh = resultPeakpicking$best_settings$parameters$snthresh, # Signal/Noise threshold
                           prefilter = c(resultPeakpicking$best_settings$parameters$prefilter, resultPeakpicking$best_settings$parameters$value_of_prefilter), # Prefilter step for the first phase. Mass traces are only retained if they contain at least [prefilter peaks] peaks with intensity >= [prefilter intensity]
                           mzCenterFun = resultPeakpicking$best_settings$parameters$mzCenterFun,
                           integrate = resultPeakpicking$best_settings$parameters$integrate, # Integration method
                           mzdiff = resultPeakpicking$best_settings$parameters$mzdiff, # minimum difference in m/z for peaks with overlapping retention times, can be negative to allow overlap
                           fitgauss = resultPeakpicking$best_settings$parameters$fitgauss,
                           noise = resultPeakpicking$best_settings$parameters$noise) # optional argument which is useful for data that was centroided without any intensity threshold, centroids with intensity < noise are omitted from ROI detection

feat_det <- xcms::findChromPeaks(raw_data, param = cwp)

# retention time correction
BiocParallel::register(BiocParallel::SerialParam())
app <- xcms::ObiwarpParam(binSize = resultRetcorGroup$best_settings$profStep, # step size (in m/z) to use for profile generation from the raw data files
                          center = center_NEG, # or use 1st sample , Sample as center sample
                          response = resultRetcorGroup$best_settings$response,
                          distFun = resultRetcorGroup$best_settings$distFunc,
                          gapInit = resultRetcorGroup$best_settings$gapInit,
                          gapExtend = resultRetcorGroup$best_settings$gapExtend,
                          factorDiag = resultRetcorGroup$best_settings$factorDiag,
                          factorGap = resultRetcorGroup$best_settings$factorGap,
                          localAlignment = ifelse(resultRetcorGroup$best_settings$localAlignment==0, F,T))

ret_cor <- xcms::adjustRtime(feat_det, param = app)

# peak grouping
# NOTE!!! It is better to optimize bw and minfrac manually to your final peak table after applying IPO optimized parameters
pgp <- xcms::PeakDensityParam(sampleGroups = as.numeric(as.factor(n_gr_t$p1_id)), # rep(1, length(fileNames(feat_det))) or as.numeric(as.factor(n_gr_t$n_gr_t))
                              bw = resultRetcorGroup$best_settings$bw, # Allowable retention time deviations, in seconds. In more detail: bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram
                              minFraction =  0.9, # or resultRetcorGroup$best_settings$minfrac , minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group
                              minSamples = resultRetcorGroup$best_settings$minsamp, # minimum number of samples necessary in at least one of the sample groups for it to be a valid group
                              binSize = resultRetcorGroup$best_settings$mzwid, # width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples
                              maxFeatures = resultRetcorGroup$best_settings$max) # maximum number of groups to identify in a single m/z slice

pk_gr <- xcms::groupChromPeaks(ret_cor, param = pgp)

# peak filling
pk_fil <- xcms::fillChromPeaks(pk_gr) # see also params: "expandMz", "expandRt", "ppm"

# final feature table
ft_tbl <- featureValues(pk_fil, value = "into")  # "index" for height "into" for area

# final peak info table
ft_inf <- featureDefinitions(pk_fil)[, c("mzmed", "rtmed")]

# join peak table and save
fvals <- cbind(feature = rownames(ft_inf), ft_inf, ft_tbl)
fwrite(fvals, "xcms after IPO NEG.csv", row.names = T)

# save all xcms objects
save(feat_det, file = "xcms obj feat_det NEG.RData")
save(ret_cor, file = "xcms obj ret_cor NEG.RData")
save(pk_gr, file = "xcms obj pk_gr NEG.RData")
save(pk_fil, file = "xcms obj pk_fil NEG.RData")

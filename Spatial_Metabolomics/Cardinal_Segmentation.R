args <- commandArgs(trailingOnly = TRUE)


library(Cardinal)


DATA_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/Data/"
OUT_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/"

i <- as.integer(args[1])

print(args)
print(class(i))


data_list <- c("VLP94C/vlp94c_dhb","VLP94A/vlp94a_dhb","VLP94D/vlp94d_dhb","VLP94B/vlp94b_dhb")
sample_name <-c("T1","C1","T2","C2")

folder <- paste0("/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/",sample_name[i],"/Segmentation/")

data <- readImzML(data_list[i],folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)


data_mean <- summarizeFeatures(data, "mean")
saveRDS(data, paste0(folder,sample_name[i],"_mean_vignette.RDS"))

data_tic <- summarizePixels(data, c(tic="sum"))
saveRDS(data, paste0(folder,sample_name[i],"_tic_vignette.RDS"))

print("generating reference peaks")
data_ref_test <- data_mean %>%
  peakPick(SNR=3) %>%
  peakAlign(ref="mean",
            tolerance=0.5,
            units="mz") %>%
  peakFilter() %>%
  process()

saveRDS(data_ref_test, paste0(folder,sample_name[i],"_test_ref_peak_vignette.RDS"))

data_ref <- data_mean %>%
  peakPick(method = "mad") %>%
  peakAlign(ref="mean",
            tolerance=0.01,
            units="mz") %>%
  peakFilter() %>%
  process()



saveRDS(data_ref, paste0(folder,sample_name[i],"_ref_peak_vignette.RDS"))

print("generating binned peaks")
data_peaks_test <- data %>%
  normalize(method="rms") %>%
  peakBin(ref=mz(data_ref_test)) %>%
  process()

saveRDS(data_peaks_test, paste0(folder,sample_name[i],"_test_peak_vignette.RDS"))


data_peaks <- data %>%
  normalize(method="rms") %>%
  peakBin(ref=mz(data_ref_test)) %>%
  process()

saveRDS(data_peaks, paste0(folder,sample_name[i],"_peak_vignette.RDS"))


print("Running PCA ...")
data_pca <- PCA(data_peaks, ncomp=3)
test_pca <- PCA(data_peaks_test, ncomp=3)
saveRDS(data_pca, paste0(folder,sample_name[i],"_pca_vignette.RDS"))
saveRDS(test_pca, paste0(folder,sample_name[i],"_pca_test_vignette.RDS"))

print("Running SSC Test ...")
set.seed(1)
test_ssc <- spatialShrunkenCentroids(data_peaks_test, method="adaptive",
                                       r=2, s=c(0,5,10,15,20,25), k=10)

print("Running SSC ...")
set.seed(1)
data_scc <- spatialShrunkenCentroids(data_peaks, method="adaptive",
                                       r=2, s=c(0,5,10,15,20,25), k=10)


saveRDS(data_scc, paste0(folder,sample_name[i],"_ssc_vignette.RDS"))
saveRDS(test_ssc, paste0(folder,sample_name[i],"_ssc_test_vignette.RDS"))







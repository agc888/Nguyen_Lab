args <- commandArgs(trailingOnly = TRUE)

#### Import Library ####
library(Cardinal)
########################


#### Setup Directories ####
DATA_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/Data/"
OUT_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/"
###########################


#### Get Index Input from Bash/Job ####
i <- as.integer(args[1])
#######################################


#### Import Data ####
data_list <- c("VLP94C/vlp94c_dhb","VLP94A/vlp94a_dhb","VLP94D/vlp94d_dhb","VLP94B/vlp94b_dhb")
sample_name <-c("T1","C1","T2","C2")
folder <- paste0("/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/",sample_name[i],"/Refined_Binning/")
data <- readImzML(data_list[i],folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
#####################



#### Run Pre-Processing Analysis ####
print("Starting Pre-Processing Analysis ... ")

#data_mean <- summarizeFeatures(data, "mean")
#saveRDS(data, paste0(folder,sample_name[i],"_mean.RDS"))

#data_tic <- summarizePixels(data, c(tic="sum"))
#saveRDS(data, paste0(folder,sample_name[i],"_summarizedPixels.RDS"))

print("Normalising Data ...")
data_pre <- data %>%
  normalize(method="rms")

print("Generating Reference Peaks ...... ")

data_ref <- data_pre %>%
  peakPick(method = "mad") %>%
  peakFilter(freq.min=0.005) %>%
  process()

#saveRDS(data_ref, paste0(folder,runs[i],"_Reference_Peaks.RDS"))

print("Generating Binned Data .........")

data_peaks <- data_pre %>%
  peakBin(ref=mz(data_ref)) %>%
  process()

saveRDS(data_peaks, paste0(folder,sample_name[i],"_Binned_Data.RDS"))

#### Run PCA and SSC Analysis ####

print("Running PCA ............")
data_pca <- PCA(data_peaks, ncomp=3)
saveRDS(data_pca, paste0(folder,sample_name[i],"_pca.RDS"))

print("Running SSC ...............")
set.seed(1)
data_scc <- spatialShrunkenCentroids(data_peaks, method="adaptive",
                                       r=2, s=c(0,5,10,15,20,25), k=5)

saveRDS(data_scc, paste0(folder,sample_name[i],"_ssc.RDS"))


print("Done! - Files are saved in '/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/'")


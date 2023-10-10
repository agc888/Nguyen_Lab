args <- commandArgs(trailingOnly = TRUE)



library(Cardinal)


DATA_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/Data/"
OUT_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/"

i <- as.integer(args[1])

print(args)
print(class(i))


data_list <- c("VLP94C/vlp94c_dhb","VLP94A/vlp94a_dhb","VLP94D/vlp94d_dhb","VLP94B/vlp94b_dhb")
sample_name <-c("T1","C1","T2","C2")

folder <- paste0("/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/",sample_name[i],"/")

data <- readImzML(data_list[i],folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)


new_data <- normalize(data, method="rms")
new_data <- process(new_data)

saveRDS(new_data, paste0(folder,sample_name[i],"_pre.RDS"))



peaks <- peakPick(new_data, method="mad")
peaks <- peakFilter(new_data, freq.min=0.005)
peaked_data <- process(peaks)
  
###Peak Binninng
peak_bins <- peakBin(new_data, ref=mz(peaked_data))
binned_data <- peak_bins %>% process()
  
saveRDS(binned_data, paste0(folder,"binned_",sample_name[i],".RDS"))

### Run Pixel Summary
data_summary <- summarizePixels(binned_data, c(tic="sum"))
saveRDS(binned_data, paste0(folder,sample_name[i],"_summarizePixels.RDS"))


### Run PCA
data_pca <- PCA(binned_data, ncomp=3)
saveRDS(data_pca, paste0(folder,sample_name[i],"_PCA.RDS"))


### Run Segmenation
set.seed(1)
data_ssc <- spatialShrunkenCentroids(binned_data, method="adaptive",r=2, s=c(0,5,10,15,20,25), k=10)
saveRDS(data_ssc, paste0(folder,sample_name[i],"_SSC.RDS"))

### Run feature Segmenation
data_dgmm <- spatialDGMM(binned_data, r=1, k=5, method="adaptive")
saveRDS(data_dgmm, paste0(folder,sample_name[i],"_dgmm.RDS"))
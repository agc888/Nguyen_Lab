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
folder <- paste0("/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/All_Samples/Refined/")

C1 <- readImzML("VLP94A/vlp94a_dhb",folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
T1 <- readImzML("VLP94C/vlp94c_dhb",folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
C2 <- readImzML("VLP94B/vlp94b_dhb",folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
T2 <- readImzML("VLP94D/vlp94d_dhb",folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)

#####################


#### Combine Samples ####
print("Combining Samples ... ")
set_centroid_to_true <- function(data){
  centroided(data) <- TRUE
  return(data)
}

### Run 1
C1 <- set_centroid_to_true(C1)
T1 <- set_centroid_to_true(T1)
C2 <- set_centroid_to_true(C2)
T2 <- set_centroid_to_true(T2)


##### T1 #####
T1_add <- T1
T1_add_pixel_data <- pixelData(T1_add)
T1_add_pixel_data@coord$x <- T1_add_pixel_data@coord$x + 200 # Adds values to pixels so they dont overlap on the plots
pixelData(T1_add) <- T1_add_pixel_data 
##############


##### C2 #####
C2_add <- C2
C2_add_pixel_data <- pixelData(C2_add)
C2_add_pixel_data@coord$y <- C2_add_pixel_data@coord$y + 200 # Adds values to pixels so they dont overlap on the plots
pixelData(C2_add) <- C2_add_pixel_data 
##############


##### T1 #####
T2_add <- T2
T2_add_pixel_data <- pixelData(T2_add)
T2_add_pixel_data@coord$x <- T2_add_pixel_data@coord$x + 200 # Adds values to pixels so they dont overlap on the plots
T2_add_pixel_data@coord$y <- T2_add_pixel_data@coord$y + 200 # Adds values to pixels so they dont overlap on the plots
pixelData(T2_add) <- T2_add_pixel_data 
##############



##### Combine Samples #####
data <- Cardinal::combine(C1, T1_add, C2_add, T2_add)
##########################

print("Normalising Data ...")
data_pre <- data %>%
  normalize(method="rms") %>%
  process()

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

saveRDS(data_peaks, paste0(folder,"All_Samples_Binned_Data.RDS"))

#### Run PCA and SSC Analysis ####

print("Running PCA ............")
data_pca <- PCA(data_peaks, ncomp=3)
saveRDS(data_pca, paste0(folder,"All_Samples_pca.RDS"))


print("Running SSC ...............")
set.seed(1)
data_scc <- spatialShrunkenCentroids(data_peaks, method="adaptive",
                                       r=2, s=c(0,5,10,15,20,25), k=5)

saveRDS(data_scc, paste0(folder,"All_Samples_ssc.RDS"))

print("Done! - Files are saved in '/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/'")


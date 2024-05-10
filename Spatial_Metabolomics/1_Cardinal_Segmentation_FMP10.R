print("This script if for the generation of ssc segmentation of MALDI mass-spectrometry data")



#### Import Library ###############################################################
library(Cardinal)
#################################################################################







#### Setup Directories ###########################################################
DATA_DIR <- "/QRISdata/Q5291/VLP94_MALDI/2023_10_01_tissue_regions_imzml_file/" 
OUT_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/9aa/"

if (!(dir.exists(OUT_DIR))){
    print(paste0("Creating New Folder - ", OUT_DIR))
    print("All Results and Analyses will be saved here! ")
    dir.create(folder)
}
##################################################################################





generate_ssc <- function(file_path) {
        
        
    #### Import Data #################################################################
    print("Samples are input using a mass.range = (160,1500) and a resolution of 10ppm")
    C1 <- readImzML("vlp94a-9aa-tissue region",folder = DATA_DIR, mass.range = c(160,1500), resolution = 10) ## Change the filename for different matrix
    T1 <- readImzML("vlp94c-9aa-tissue region",folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
    C2 <- readImzML("vlp94b-9aa-tissue region",folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
    T2 <- readImzML("vlp94d-9aa-tissue region",folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
    ##################################################################################
    
    
    
    
    # we must sent the centroid of each sample = TRUE so that we can merge them 
    
    #### Combine Samples ############################################################
    print("Combining Samples ...... ")
    set_centroid_to_true <- function(data){
      centroided(data) <- TRUE
      return(data)
    }
    
    ### Run 1
    C1 <- set_centroid_to_true(C1)
    T1 <- set_centroid_to_true(T1)
    C2 <- set_centroid_to_true(C2)
    T2 <- set_centroid_to_true(T2)
    
    
    ##### Combine Samples #####
    data <- Cardinal::combine(C1, T1, C2, T2)
    #################################################################################
    
    
    
    
    
    
    #Now we can run ssc segmentation analysis. PCA and SSC cardinal objects will be saved in this process
    
    #### Run Pre-Processing Analysis ################################################
    print("Starting Pre-Processing Analysis ......... ")
    
    data_mean <- summarizeFeatures(data, "mean")
    
    
    print("Generating Reference Peaks ............ ")
    
    data_ref <- data_mean %>%
      peakPick(method = "mad") %>%
      peakFilter(freq.min=0.005) %>%
      process()
    
    print("Generating Binned Data ...............")
    
    data_peaks <- data %>%
      normalize(method="rms") %>%
      peakBin(ref=mz(data_ref)) %>%
      process()
    
    saveRDS(data_peaks, paste0(OUT_DIR,"All_Samples_Binned_Data.RDS")) #this saves the binned data
    
    #### Run PCA and SSC Analysis ####
    
    print("Running PCA ..................")
    data_pca <- PCA(data_peaks, ncomp=3)
    saveRDS(data_pca, paste0(OUT_DIR,"All_Samples_pca.RDS"))
    
    print("Running SSC .....................")
    set.seed(1)
    data_scc <- spatialShrunkenCentroids(data_peaks, method="adaptive",
                                           r=2, s=c(0,5,10,15,20,25,30), k=10) #can change these parameters if needed
    
    saveRDS(data_scc, paste0(OUT_DIR,"All_Samples_ssc.RDS"))
    #################################################################################
    
    
    print(paste0("Done! - Files are saved in '",OUT_DIR))


print("This script if for the generation of ssc segmentation of MALDI mass-spectrometry data")



#### Import Library ###############################################################
library(Cardinal)
#################################################################################

matricies <- c("DHB","FMP10", "9aa")

for (matrix in matricies) {

    if (matrix == "DHB") {
        mtx.name <- "dhb"
    } else if (matrix == "FMP10"){
        mtx.name <- "fmp10"
    } else {
        mtx.name <- "9aa"
    }

    
    #### Setup Directories ###########################################################
    DATA_DIR <- paste0("/QRISdata/Q5291/MS_Spatial_metabolomics/10. Data/raw_data/MALDI/",matrix)
    OUT_DIR <- paste0("/QRISdata/Q5291/MS_Spatial_metabolomics/10. Data/processed_data/", matrix, "/ssc/")
    
    if (!(dir.exists(OUT_DIR))){
        print(paste0("Creating New Folder - ", OUT_DIR))
        print("All Results and Analyses will be saved here! ")
        dir.create(folder)
    }
    ##################################################################################
    
    
    
    #### Import Data #################################################################
    print("Samples are input using a mass.range = (160,1500) and a resolution of 10ppm")
    C1 <- readImzML(paste0("vlp94a_",mtx.name),folder = DATA_DIR, mass.range = c(160,1500), resolution = 10) ## Change the filename for different matrix
    T1 <- readImzML(paste0("vlp94c_",mtx.name),folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
    C2 <- readImzML(paste0("vlp94b_",mtx.name),folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
    T2 <- readImzML(paste0("vlp94d_",mtx.name),folder = DATA_DIR, mass.range = c(160,1500), resolution = 10)
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
    data <- Cardinal::combine(C1, T1, C2, T2) #merge samples into one object
    #################################################################################

    
    
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
    
    saveRDS(data_peaks, paste0(OUT_DIR,mtx.name,"_Binned_Data.RDS")) #this saves the binned data
    
    #### Run PCA and SSC Analysis ####
    
    print("Running PCA ..................")
    data_pca <- PCA(data_peaks, ncomp=3)
    saveRDS(data_pca, paste0(OUT_DIR,mtx.name,"_pca.RDS"))
    
    print("Running SSC .....................")
    set.seed(1)
    data_scc <- spatialShrunkenCentroids(data_peaks, method="adaptive",
                                           r=2, s=c(0,5,10,15,20,25,30), k=10) #can change these parameters if needed
    
    saveRDS(data_scc, paste0(OUT_DIR,mtx.name,"_ssc.RDS"))
    #################################################################################
    rm(data)
    rm(data_mean)
    rm(data_ref)
    rm(data_peaks)
    rm(data_pca)
    rm(data_scc)
    
    print(paste0("Done! - Files are saved in '",OUT_DIR, "'"))
    
}
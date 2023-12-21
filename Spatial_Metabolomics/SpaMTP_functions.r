#### SpaMTP Cardinal Preprocessing Functions ############################################################################################################################################################################

### Calculates coordinate change for 
calculate_coordinates <- function(index, ncol, padding) {
        row <- ceiling(index / ncol)
        col <- (index - 1) %% ncol + 1
        x <- (col - 1) * padding
        y <- (row - 1) * padding
        return(list(x = x, y = y))
}

### This function sets the centroids to be = TRUE so that you can merge objects together
set_centroid_to_true <- function(data){
        centroided(data) <- TRUE
        return(data)
}

### merges datasets together so they are on the same plot
merge.cardinalData <- function(data_list, ncols = 2, padding = 200){
    
    ### NOTE: mass.range and resolution of each sample in the data list must be the same to merge
    
    message("Setting Centroids as TRUE ...... ")

    data_list <- lapply(data_list, set_centroid_to_true)

    ### Changes pixel coordinates so that they do not overlap on the same plot
    message("Shifting Pixel Coordinates to Combining Samples ...... ")
         
    n <- 1
    while ( n <= length(data_list)){
        
        data_copy <- data_list[[n]]
        data_copy_pixel_data <- pixelData(data_copy)
        coordnate_list <- calculate_coordinates(n,ncols, padding)
        data_copy_pixel_data@coord$x <- data_copy_pixel_data@coord$x + as.numeric(coordnate_list$x) # Adds values to pixels so they dont overlap on the plots
        data_copy_pixel_data@coord$x <- data_copy_pixel_data@coord$y + as.numeric(coordnate_list$y) 
        pixelData(data_list[[n]]) <- data_copy_pixel_data 
        n <- n+1
    }

    return(Cardinal::combine(data_list))
    
}


########################################################################################################################################################################################################################



#### SpaMTP Cardinal SSC clustering segmentation ###############################################################################################################################################################################

## THIS PROCESS IS VERY LONG/COMPUTATIONALLY EXPENSIVE
generate_SSCclusters <- function(data, output_files = NULL){ 
    message("Starting Pre-Processing Analysis ......... ")

    data_mean <- summarizeFeatures(data, "mean")

    message("Generating Reference Peaks ............ ")

    data_ref <- data_mean %>%
      peakPick(method = "mad") %>%
      peakFilter(freq.min=0.005) %>%
      process()
  

    message("Generating Binned Data ...............")

    data_peaks <- data %>%
      normalize(method="rms") %>%
      peakBin(ref=mz(data_ref)) %>%
      process()

    

    #### Run PCA and SSC Analysis ####

    message("Running PCA ..................")
    data_pca <- PCA(data_peaks, ncomp=3)
    

    message("Running SSC .....................")
    set.seed(1)
    data_ssc <- spatialShrunkenCentroids(data_peaks, method="adaptive",
                                       r=2, s=c(0,5,10,15,20,25), k=10)

    message("Done! - Files are saved in '/QRISdata/Q1851/Andrew_C/Metabolomics/Pipeline/Segmentation/'")

    if (!(is.null(output_files))){
    
        saveRDS(data_ref, paste0(output_files,"All_Samples_Reference_Peaks.RDS"))
        saveRDS(data_pca, paste0(output_files,"All_Samples_pca.RDS"))
        saveRDS(data_peaks, paste0(output_files,"All_Samples_Binned_Data.RDS"))
        saveRDS(data_ssc, paste0(output_files,"All_Samples_ssc.RDS"))
    }
    return(data_ssc)
}


### Adds ssc annotation to binned/normalised m/z data
add_ssc_annotation <- function(data_binned, data_ssc, resolution = 25){

    message(paste0("Getting cluster segments for resolution (s) = ", resolution))
    data_bin <- data_binned
    cluster_idx <- which(modelData(data_ssc)[["s"]] == resolution)

    classes <- resultData(data_ssc)[[cluster_idx]][[1]]
    pixel_data <- pixelData(data_bin)
  
  
    pixel_data[["ssc"]] <- classes
  
    pixelData(data_bin) <- pixel_data
    return(data_bin)

}

########################################################################################################################################################################################################################






#### SpaMTP Cardinal to Seurat Functions ###############################################################################################################################################################################
cardinal_to_seurat <- function(data,run_name, seurat.coord = NULL){

    message("Convering Cardinal object to Seurat object .... ")
    run_data <- subsetPixels(data, run(data) == paste0(run_name))
    sparse_matrix <- spectra(run_data)

    if (!(is.null(seurat.coord))){
        message("Convering Cardinal Coordinates to Seurat Visium Coordinates specified in the seurat.coord file .... ")
        pixel_data <- pixelData(run_data)
        pixel_data[["x_coord",]] <- seurat.coord$X_new # changes coordinates to matched Visium Object
        pixel_data[["y_coord",]] <- seurat.coord$Y_new
        pixelData(run_data) <- pixel_data 
    }


    message("Generating Seurat Barcode Labels from Pixel Coordinates .... ")
    spot_name <- c()
  
    for(idx in seq(1,length(pixelData(run_data)[[1]]))){
        x_coord <- pixelData(run_data)[["x_coord",]][idx]
        y_coord <- pixelData(run_data)[["y_coord",]][idx]
        name <- paste0(x_coord,"_",y_coord)
        spot_name <- c(spot_name, name)
    }
    
    
    colnames(sparse_matrix)<- spot_name
    rownames(sparse_matrix)<- paste("mz-", featureData(run_data)@mz, sep = "") 

    message("Constructing Seurat Object ....")
    mat <- as.matrix(sparse_matrix)
    
    seuratobj <- CreateSeuratObject(mat, assay = "Spatial")

    message("Adding Pixel Metadata ....")
    seuratobj <- AddMetaData(seuratobj,col.name = "sample", metadata = run(run_data))

    for (name in names(pixelData(run_data))){
        seuratobj <- AddMetaData(seuratobj,col.name = name, metadata = pixelData(run_data)[[name,]])
    }

    message("Creating Centroids for Spatial Seurat Object ....")
    ## Add spatial data
    cents <- CreateCentroids(data.frame(x = c(pixelData(run_data)[["x_coord",]]), y = c(pixelData(run_data)[["y_coord",]]), cell = c(spot_name)))


    segmentations.data <- list(
        "centroids" = cents,
        "segmentation" = NULL
      )

    coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = NULL,
        assay = "Spatial"
    )

    seuratobj[["fov"]] <- coords
  
    return(seuratobj)
}

########################################################################################################################################################################################################################




#### SpaMTP Seurat Plotting Functions #################################################################################################################################################################################
find_nearest <- function(data, target_mz){
  
  numbers <- as.numeric(gsub("mz-", "", Features(data)))
  closest_number <- numbers[which.min(abs(numbers - target_mz))]
  return(paste0("mz-",closest_number))
}


plusminus <- function(data, target_mz, plus_minus){
    feature_list <- c()
    center <- find_nearest(data, target_mz)
    
    center_value <- as.numeric(gsub("mz-", "", center))
    up_value <- center_value + as.numeric(plus_minus) 
    low_value <- center_value - as.numeric(plus_minus)

    upper <- find_nearest(data, up_value)
    lower <- find_nearest(data, low_value)

    feature_list <- c(feature_list, center)
    if (!(upper %in% feature_list)){
        feature_list <- c(feature_list, upper)
    }
    if (!(lower %in% feature_list)){
        feature_list <- c(feature_list, lower)
    }

    if (length(feature_list) == 1){
        print("No other m/z peaks in plusminus range -> increase to include more peaks")
    }
    return(feature_list)
}



bin.mz <- function(data, mz_list){
    data_copy <- data
    assay_counts <- data_copy@assays$Spatial$counts
    selected_genes <- assay_counts[mz_list, , drop = FALSE]
    #return(selected_genes)
    if (is.null(selected_genes)) {
        stop("One or more genes not found in the assay counts.")
    }
  
    #data_copy[["combined_mz"]]<- colSums(selected_genes)
    return(colSums(selected_genes))
}


ImageMZPlot <- function(object, 
                        mzs,
                        plusminus = NULL,
                        fov = NULL,
                        boundaries = NULL,
                        cols = if (isTRUE(x = blend)) {
                            c("lightgrey", "#ff0000", "#00ff00")
                        } else {
                            c("lightgrey", "firebrick1")
                        },
                        size = 0.5,
                        min.cutoff = NA,
                        max.cutoff = NA,
                        split.by = NULL,
                        molecules = NULL,
                        mols.size = 0.1,
                        mols.cols = NULL,
                        nmols = 1000,
                        alpha = 1,
                        border.color = "white",
                        border.size = NULL,
                        dark.background = TRUE,
                        blend = FALSE,
                        blend.threshold = 0.5,
                        crop = FALSE,
                        cells = NULL,
                        scale = c("feature", "all", "none"),
                        overlap = FALSE,
                        axes = FALSE,
                        combine = TRUE,
                        coord.fixed = TRUE
                       ){

    if (is.null(mzs)){
        stop("No mz values have been supplied")
    } else{
        
        mz_list <- c()
        for (target_mz in mzs){
            mz_string <- find_nearest(object, target_mz)
            mz_list <- c(mz_list, mz_string)
        }
    }

    if (!(is.null(plusminus))){
        
        data_copy <- object
        col_names_to_plot <- c()
        plot_titles <- c()
        
        for (target_mz in mz_list){
            mz_integer <- as.numeric(strsplit(target_mz, "-")[[1]][2])

            meta_col <- bin.mz(data_copy, plusminus(data_copy, mz_integer, plusminus))
            
            col_name <- paste0(target_mz,"_plusminus_", plusminus)
            plot_name <- paste0("mz: ", round(mz_integer, 3)," ± ", plusminus)
        
            col_names_to_plot <- c(col_names_to_plot,col_name)
            plot_titles <- c(plot_titles,plot_name)
            data_copy[[col_name]] = meta_col
        }

        plot <- ImageFeaturePlot(object = data_copy, 
                     features = col_names_to_plot, 
                     fov = fov,
                     boundaries = boundaries,
                     cols = cols,
                     size = size,
                     min.cutoff = min.cutoff,
                     max.cutoff = max.cutoff,
                     split.by = split.by,
                     molecules = molecules,
                     mols.size = mols.size,
                     mols.cols = mols.cols,
                     nmols = nmols,
                     alpha = alpha,
                     border.color = border.color,
                     border.size = border.size,
                     dark.background = dark.background,
                     crop = crop,
                     cells = cells,
                     scale = scale,
                     overlap = overlap,
                     axes = axes,
                     combine = combine,
                     coord.fixed = coord.fixed
                     )
        
        for (plot_idx in seq(1, length(plot_titles))){
            plot[[plot_idx]] <- plot[[plot_idx]] + ggtitle(plot_titles[[plot_idx]])+labs(fill = plot_titles[[plot_idx]])
        }
        
    } else {

        plot <- ImageFeaturePlot(object = object, 
                     features = mz_list, 
                     fov = fov,
                     boundaries = boundaries,
                     cols = cols,
                     size = size,
                     min.cutoff = min.cutoff,
                     max.cutoff = max.cutoff,
                     split.by = split.by,
                     molecules = molecules,
                     mols.size = mols.size,
                     mols.cols = mols.cols,
                     nmols = nmols,
                     alpha = alpha,
                     border.color = border.color,
                     border.size = border.size,
                     dark.background = dark.background,
                     crop = crop,
                     cells = cells,
                     scale = scale,
                     overlap = overlap,
                     axes = axes,
                     combine = combine,
                     coord.fixed = coord.fixed
                     ) 
        
        for (plot_idx in seq(1, length(mz_list))){
            mz_integer <- as.numeric(strsplit(mz_list[plot_idx], "-")[[1]][2])
            plot[[plot_idx]] <- plot[[plot_idx]] + ggtitle(paste0("mz: ",round(mz_integer,3)))+labs(fill = paste0("mz: ",round(mz_integer,3)))
        }
    }


    return(plot)

}
    
########################################################################################################################################################################################################################


#### SpaMTP Saving Data Objects ########################################################################################################################################################################################
saveSeuratData <- function(data, outdir, assay = "Spatial", slot = "counts"){
    
    message(paste0("Generating new directory to store output here: ", outdir))
    message(paste0("Writing ", slot," slot to matrix.mtx, barcode.tsv, genes.tsv"))
    write10xCounts(data[[assay]][slot], path = outdir, overwrite = TRUE)
    
    message("Writing @metadata slot to metadata.csv")
    fwrite(data@meta.data, paste0(outdir,"metadata.csv"))                   
}

########################################################################################################################################################################################################################














#!! ALL CODE BELOW WAS WRITEN BY Christopher Fitzgerald github https://github.com/ChemCharles !!#

#### SpaMTP Saving Data Objects ########################################################################################################################################################################################








########################################################################################################################################################################################################################





















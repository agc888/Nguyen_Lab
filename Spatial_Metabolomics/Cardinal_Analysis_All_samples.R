print("Starting Analysis Script ")

#### Import Library ####
library(Cardinal)
library(cowplot)
########################


#### Setup Directories #######################################################################################################
DATA_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/Data/"
OUT_DIR <- "/QRISdata/Q1851/Andrew_C/Metabolomics/"
##############################################################################################################################



#### Import Data #############################################################################################################
folder <- paste0("/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/All_Samples/Analysis/")

all_data_binned <- readRDS("/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/All_Samples/Refined/All_Samples_Binned_Data.RDS")
all_data_ssc <- readRDS("/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/All_Samples/All_Samples_ssc.RDS")

##############################################################################################################################





#### Plotting Segmentation Results ###########################################################################################

print("Plotting Segmentation Results ... ")
darkmode()

jpeg(file = paste0(folder,"All_Samples_Segmentation_Results.jpeg"), width = 2000, height = 1200)
image(all_data_ssc, model=list(s=c(0,5,10,15,20,25)))
dev.off()


print("Plotting Drug Intensities ... ")

## Drug 1
pdf(file = paste0(folder,"All_Samples_mz448.pdf"), width = 10, height = 8)
image(all_data_binned, mz=448.24, plusminus=0.005,  colorscale=c(rep("#161616",10),col.map("jet")[1:90]))
dev.off()

## Drug 2
pdf(file = paste0(folder,"All_Samples_mz450.pdf"), width = 10, height = 8)
image(all_data_binned, mz=450.25, plusminus=0.005,  colorscale=c(rep("#161616",10),col.map("jet")[1:90]))
dev.off()
##############################################################################################################################





#### Adding Segmentation #####################################################################################################

print("Splitting Sample by Tumour Section ...... ")

add_tumour_annotation <- function(data_binned, data_ssc, cancer_class){
  
  data_bin <- data_binned
  classes <- resultData(data_ssc)[[6]][[1]]
  pixel_data <- pixelData(data_bin)
  
  
  cancer_class_column <- c()
  for (row in classes){
    cancer <- FALSE
    if (row == cancer_class){
         cancer <- TRUE
    }
    cancer_class_column <- c(cancer_class_column,cancer)
  }
      
  pixel_data[["cancer_class"]] <- cancer_class_column
  
  treatment_list <- c()
  
  for (run in pixel_data@run){
    
    treated <- "Control"
    if (run == "VLP94C/vlp94c_dhb"|| run == "VLP94D/vlp94d_dhb"){
      treated <- "Treated"
    }
    treatment_list <- c(treatment_list, treated)
    
  }
  
  pixel_data[["treatment_class"]] <- treatment_list
  
  pixel_data[["region"]] <- classes
  
  pixelData(data_bin) <- pixel_data
  return(data_bin)

}

all_data <- add_tumour_annotation(all_data_binned, all_data_ssc, "3")
only_tumour <- subsetPixels(all_data, region==3)

saveRDS(all_data, paste0(folder, "All_Samples_annotated_data.RDS"))

print("Plotting Tumour Section ......... ")

pdf(file = paste0(folder,"All_Samples_tumour_section.pdf"), width = 10, height = 8)
image(only_tumour, mz=450.25, plusminus=0.005,  colorscale= c(rep("#161616",10),col.map("jet")[1:90]))
dev.off()

##############################################################################################################################






#### Add Region and Run Pixel Data ###########################################################################################
print("Adding Region and Run Pixel Data ............ ")

add_feature_data <- function(data_binned){
    
  data <- data_binned
  feature_data <- featureData(data)
    
  for (i in unique(data$region)){
    fdata <- data %>% subsetPixels(region == i) %>% summarizeFeatures("mean",  as = "DataFrame")
    feature_data[[paste0("class_",i)]] <- fdata$mean
  }
  
  for (i in unique(run(data))){
    fdata <- data %>% subsetPixels(run(data) == i) %>% summarizeFeatures("mean",  as = "DataFrame")
    feature_data[[paste0("run_",i)]] <- fdata$mean
  }
  
  for (i in unique(data$treatment_class)){
    fdata <- data %>% subsetPixels(treatment_class == i) %>% summarizeFeatures("mean",  as = "DataFrame")
    feature_data[[paste0(i, "_group")]] <- fdata$mean
  }
  
  featureData(data) <- feature_data
  return(data)
}

new_data <- add_feature_data(only_tumour)
e
##############################################################################################################################






#### Finding Significant Peaks  ###############################################################################################

print("Running Means Test to find top 500 Signficant Peaks between Treated and Control Samples ............... ")
mtest_tumour_tissue <- meansTest(new_data, ~ treatment_class, groups=run(only_tumour))
saveRDS(mtest_tumour_tissue, paste0(folder, "All_Samples_meansTest_tumour_region.RDS"))
#mtest_whole_tissue <- meansTest(all_data, ~ treatment_class, groups=run(all_data))

top_features_tumour_tissue <- topFeatures(mtest_tumour_tissue, p.adjust="fdr", AdjP < .05, n=500)
saveRDS(top_features_tumour_tissue, paste0(folder, "All_Samples_topFeatures_tumour_region.RDS"))
#top_features_whole_tissue <- topFeatures(mtest_whole_tissue, p.adjust="fdr", AdjP < .05, n=500)


## Add Fold Change info
print("Add Fold Change Values .................. ")
Treated_Fold_Change <- c()
for (i in rownames(top_features_tumour_tissue)){
  treated_mean <- top_features_tumour_tissue$Treated_group[as.integer(i)]
  control_mean <- top_features_tumour_tissue$Control_group[as.integer(i)]
  if (treated_mean > control_mean){
    change <- "UP"
  } else {
    change <- "DOWN"
  }
  Treated_Fold_Change <- c(Treated_Fold_Change, change)
}
top_features_tumour_tissue$Treated_Fold_Change <- Treated_Fold_Change


print("Writing Top Peaks To File ..................... ")
write.csv(top_features_tumour_tissue, paste0(folder,"All_Sample_tumour_features.csv"))

##############################################################################################################################






#### Generating top10 Significant Peak plots  #################################################################################

print("Generating Significnat Peak Plots ........................ ")


#### UP Regulated Peaks
up_features <- top_features_tumour_tissue[top_features_tumour_tissue$Treated_Fold_Change == "UP", "feature"][1:10]

plots <- list()
for (i in seq(1,length(up_features))){
    plot <- plot(mtest_tumour_tissue, ylab = "intensity", model = DataFrame(feature = c(up_features[i])), strip = FALSE)
    plots[[i]] <- plot
}

plots[[1]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[1],][["mz"]],digits = 3)))
p1 <- recordPlot()

plots[[2]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[2],][["mz"]],digits = 3)))
p2 <- recordPlot()

plots[[3]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[3],][["mz"]],digits = 3)))
p3 <- recordPlot()

plots[[4]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[4],][["mz"]],digits = 3)))
p4 <- recordPlot()

plots[[5]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[5],][["mz"]],digits = 3)))
p5 <- recordPlot()

plots[[6]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[6],][["mz"]],digits = 3)))
p6 <- recordPlot()

plots[[7]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[7],][["mz"]],digits = 3)))
p7 <- recordPlot()

plots[[8]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[8],][["mz"]],digits = 3)))
p8 <- recordPlot()

plots[[9]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[9],][["mz"]],digits = 3)))
p9 <- recordPlot()

plots[[10]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == up_features[10],][["mz"]],digits = 3)))
p10 <- recordPlot()


top10_plots <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow = 2, ncol = 5, scale = 0.8)

save_plot(filename = paste0(folder,"All_Samples_top10_UP_peaks.pdf"), top10_plots , base_height = 8, base_width = 10)


plots <- list()
for (i in seq(1,length(up_features))){
    plot <- image(only_tumour, mz = as.data.frame(fData(mtest_tumour_tissue)[as.integer(up_features[i]),0])[[1]], plusminus = 0.005, colorscale=c(rep("#161616",10),col.map("jet")[1:90]))
    plots[[i]] <- plot
}

plots[[1]]
p1 <- recordPlot()

plots[[2]]
p2 <- recordPlot()

plots[[3]]
p3 <- recordPlot()

plots[[4]]
p4 <- recordPlot()

plots[[5]]
p5 <- recordPlot()

plots[[6]]
p6 <- recordPlot()

plots[[7]]
p7 <- recordPlot()

plots[[8]]
p8 <- recordPlot()

plots[[9]]
p9 <- recordPlot()

plots[[10]]
p10 <- recordPlot()


top10_plots <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow = 4, ncol = 3, scale = 0.8)

save_plot(filename = paste0(folder,"All_Samples_top10_UP_spectra.pdf"),top10_plots, base_height = 30, base_width = 40)



#### Down Regulated Peaks
down_features <- top_features_tumour_tissue[top_features_tumour_tissue$Treated_Fold_Change == "DOWN", "feature"][1:10]


plots <- list()
for (i in seq(1,length(down_features))){
  plot <- plot(mtest_tumour_tissue, ylab = "intensity", model = DataFrame(feature = c(down_features[i])), strip = FALSE)
  plots[[i]] <- plot
}

plots[[1]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[1],][["mz"]],digits = 3)))
p1 <- recordPlot()

plots[[2]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[2],][["mz"]],digits = 3)))
p2 <- recordPlot()

plots[[3]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[3],][["mz"]],digits = 3)))
p3 <- recordPlot()

plots[[4]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[4],][["mz"]],digits = 3)))
p4 <- recordPlot()

plots[[5]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[5],][["mz"]],digits = 3)))
p5 <- recordPlot()

plots[[6]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[6],][["mz"]],digits = 3)))
p6 <- recordPlot()

plots[[7]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[7],][["mz"]],digits = 3)))
p7 <- recordPlot()

plots[[8]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[8],][["mz"]],digits = 3)))
p8 <- recordPlot()

plots[[9]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[9],][["mz"]],digits = 3)))
p9 <- recordPlot()

plots[[10]]
title(paste0("m/z = ",round(top_features_tumour_tissue[top_features_tumour_tissue$feature == down_features[10],][["mz"]],digits = 3)))
p10 <- recordPlot()



top10_plots <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow = 2, ncol = 5, scale = 0.8)

save_plot(filename = paste0(folder,"All_Samples_top10_DOWN_peaks.pdf"), top10_plots, base_height = 8, base_width = 10)




plots <- list()
for (i in seq(1,length(down_features))){
    plot <- image(only_tumour, mz = as.data.frame(fData(mtest_tumour_tissue)[as.integer(down_features[i]),0])[[1]], plusminus = 0.005, colorscale=c(rep("#161616",10),col.map("jet")[1:90]))
    plots[[i]] <- plot
}

plots[[1]]
p1 <- recordPlot()

plots[[2]]
p2 <- recordPlot()

plots[[3]]
p3 <- recordPlot()

plots[[4]]
p4 <- recordPlot()

plots[[5]]
p5 <- recordPlot()

plots[[6]]
p6 <- recordPlot()

plots[[7]]
p7 <- recordPlot()

plots[[8]]
p8 <- recordPlot()

plots[[9]]
p9 <- recordPlot()

plots[[10]]
p10 <- recordPlot()


top10_plots <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow = 4, ncol = 3, scale = 0.8)

save_plot(filename = paste0(folder,"All_Samples_top10_DOWN_spectra.pdf"), top10_plots, base_height = 30, base_width = 40)


##############################################################################################################################


print("Done! - Files are saved in '/QRISdata/Q1851/Andrew_C/Metabolomics/Analysis/All_Samples/Analysis/'")


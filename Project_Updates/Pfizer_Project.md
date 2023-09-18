# Pfizer Project

# Table of Contents
- [Pfizer Project](#pfizer-project)
- [Table of Contents](#table-of-contents)
- [Reference Dataset](#reference-dataset)
			- [UPDATE: 18/09/23](#update-180923)
	- [Single Cell Data](#single-cell-data)
	- [Visium Data](#visium-data)
	- [Xenium Data](#xenium-data)


# Reference Dataset

#### UPDATE: 18/09/23 #### 

1. Looked at single cell reference dataset (Wu et.al 2021):
	Performed label transfer on this dataset to our sample matched scRNA data
	-> initially the reference dataset object used only had 308 genes in it
	-> this was fixed using raw data and ~12,000 genes were used to perform label transfer again

2. Comparing results to before a lot more cells are labeled as cancer
	-> Sub clustering has improved with a lot more of the Luminal subtypes matching
	-> Data was normalised when running label transfer

NEXT STEP:

I will run label transfer again without SCTransform:
	-> also look at using clusters rather then the label transfer itself
	
REULTS:

![alt text](/~/Downloads/non_normalised_data_comparision.png "Title")



## Single Cell Data



## Visium Data



## Xenium Data


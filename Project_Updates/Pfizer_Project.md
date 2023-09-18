# Pfizer Project

## Table of Contents
1. [Reference Dataset](##Reference Dataset)
2. [Single Cell Data](##Single Cell Data)
3. [Visium Data](##Visium Data)
4. [Xenium Data](##Xenium Data)


## Reference Dataset

#### UPDATE: 18/09/23 #### 

1. Looked at single cell reference dataset (Wu et.al 2021):\
	Performed label transfer on this dataset to our sample matched scRNA data\
	-> initially the reference dataset object used only had 308 genes in it\
	-> this was fixed using raw data and ~12,000 genes were used to perform label transfer again\

2. Comparing results to before a lot more cells are labeled as cancer\
	-> Sub clustering has improved with a lot more of the Luminal subtypes matching\
	-> Data was normalised when running label transfer\

NEXT STEP:

I will run label transfer again without SCTransform:\
	-> also look at using clusters rather then the label transfer itself\

## Single Cell Data




## Visium Data






## Xenium Data

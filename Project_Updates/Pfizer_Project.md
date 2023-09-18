############################### Pfizer Project ############################### 

# Table of Contents
1. [Example](#example)
2. [Example2](#example2)
3. [Third Example](#third-example)
4. [Fourth Example](#fourth-examplehttpwwwfourthexamplecom)




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


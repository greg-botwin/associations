## Purpose
The purpose of this Shiny application is to interactively and quickly review previously found associations identified by the Cedars-Sinai F. Widjaja Foundation Inflammatory Bowel and Immunobiology Research Institute.

## Functionality
### Search
Users can search associations by **Gene** (e.g. NOD2), **Position** (e.g. Chr. 16, between 50731050-50766987), or by **RS ID** (e.g. rs5743293). Searching by position queries by the position of the SNP not by the position of the Gene. 

Users can search by multiple Genes or SNPs at a time by selecting multiple Genes or SNPs respectively. If a user has a large list of Genes or SNPs they can copy and paste them into the search box but the Genes/SNPs must be in comma separated format (e.g. IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10).

### Filtering
Users can filter their query to only include SNP associations that reach user defined thresholds. Users can filter by:

Users can filter their query to only include SNP associations that reach user defined thresholds. Users can filter by: 

1. P Value
2. Minor Allele Frequency
3. SNP Location (e.g. intron, coding, complex...)
4. Overall Missingness
5. Overall Hardy-Weinberg Assumptions

Once a table of results is returned, a user can sub-filter the table by the output to further narrow down the results. 

## Data 

PUT INFO HERE ON THE NUMBER AND TYPE OF PHENOTYPES TESTED WITH THE SCINETICS N of PATIENTS and  

### Guidance for Analysts Submitting Association Data
The backbone of the application is genetic association tests performed by scientists at the Cedars-Sinai Inflammatory Bowel and Immunobiology Research Institute. To have your association stored in this application the following information is needed in a CSV file per association performed.

The backbone of the application is genetic association tests performed by scientists at the Cedars-Sinai Inflammatory Bowel and Immunobiology Research Institute. To have your association stored in this application the following information is needed in a CSV file per association test. 

1. Illumina Chip ID
2. RS ID
3. P Value
4. Odds Ratio or B Value 
5. Confidence Intervals (if performed)

In addition to the above the analyst will need to submit information on:

1. The name of the analyst who performed the association
2. The year the analysis was run
3. The sample of patients the association was performed on (e.g. Disease Sub-Type (CD or UC), Race, Ethnicity, Jewish, other unique criteria of your analysis)
4. The number of patients in the association (total number in association, both cases and controls)
5. The phenotype tested


### Annotations
Specific information provided by Talin will go here. Plan is to use hg19 genome locations and UCSC gene information.
Hardy-Weinberg, minor allele frequency, and missingness will be performed globally on all runs. Additional details to follow. 

## Created By:
The Cedars-Sinai F. Widjaja Foundation Inflammatory Bowel and Immunobiology Research Institute Translational Genomics Group.





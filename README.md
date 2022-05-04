# Materials for The HLA-Associated Type-1-Diabetes in World Populations Data Analysis Session at the 18th International HLA & Immunogenetics Workshop

This repository contains materials for the _HLA-Associated Type-1-Diabetes in World Populations_ session at the 2022 18th IHIW Meeting.

These materials consist of a dataset for data analysis (DRB1.T1D.Azerbaijan.txt) and an R script (18IHIW.R) for use in the session.

To use these files, first make a new folder called "18IHIW" in your Documents folder. Then do the following:

Click the green Code button above, and click 'Download ZIP'.
Unzip the "18IHIW.T1D-main.zip" file in your Downloads folder (either by double-clicking it, or by right-clicking, depending on your system). This will create a new folder called "18IHIW.T1D-main" in your Downloads folder.
Move the DRB1.T1D.Azerbaijan.txt and 18IHIW.R files in the Downloads/18IHIW.T1D-main folder to your Documents/18IHIW folder.
To use these files, Workshop participants should have R, Rstudio, and BIGDAWG installed on their computers.

To install R on your computer, follow the instructions for linux, macOS, or Windows here.

Once R is installed, install Rstudio on your computer using the instructions for your variety of linux, macOS or Windows here.

Once Rstudio has been installed, open Rstudio (by double clicking the icon).

Install BIGDAWG in the Rstudio environment by doing the following:

1. Go to the Tools pulldown menu, and select "Install Packages ...".
2. In the popup window, select "Repository (CRAN)" in the "Install from:" menu, and enter "BIGDAWG" in the "Packages (separate multiple with space or comma):" menu.
3. Click the "Install" button.
4. Installing BIGDAWG will automatically install all of the other R packages that BIGDAWG requires. This may take a few minutes, and may generate a lot of text in your Rstudio console. This is normal.

---
# Azerbaijan T1D Dataset
Analyses of the associated dataset (AZ_18IHIW.txt) were published in [2018 by Ahdamov et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6384092/), as part of a study of Type-1-Diabetes (T1D) in Azerbaijan.

The dataset is a tab-delimtied text file that consists of four columns, *Sample_ID*, *Status*, *DRB1_1* and *DRB1_2*, containing data for 315 individuals, with each line of the file containing all of the data for one individual. 

The 106 Azerbaijani individuals with T1D are indicated by the value of 1 in the *Status* column, and the 209 Azerbaijani individuals without T1D (controls) are indicated by the value of 0 in the *Status* column.

The DRB1 genotypes are represented in two columns, *DRB1_1* and *DRB1_2*, with one two-field DRB1 allele name appearing in each column. 

If genotype data for other loci were available for this dataset, they would appear in other pairs of columns. For example, if DQB1 data were available, the DQB1 allele names for each individual would appear in *DQB1_1* and *DQB1_2* columns.

If HLA genotype data were missing for any individual, the columns for that individual would contain either '****' or 'NA' values.

When analyzing HLA data like these with BIGDAWG, it is important to ensure that there are no extraneous spaces associated with any of the data. Allele names with leading or trailing spaces (' 02:01 ') may be treated as distinct from allele names without any associated spaces ('02:01'). 

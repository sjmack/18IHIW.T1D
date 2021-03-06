# Materials for The HLA-Associated Type-1-Diabetes in World Populations Data Analysis Session at the 18th International HLA & Immunogenetics Workshop

This repository contains materials for the _HLA-Associated Type-1-Diabetes in World Populations_ session at the 2022 18th IHIW Meeting.

These materials consist of a data set for data analysis (AZ_18IHIW.txt), an R script (18IHIW.R) and four R Data Files (AZ.Locus.rda, AZ.HWE.rda, AZ.DRB1_0301_negative.rda and BDStrat.rda) for use in the session.

To use these files, first make a new folder called "18IHIW" in your Documents folder. Then do the following:

Click the green Code button above, and click 'Download ZIP'.
Unzip the "18IHIW.T1D-main.zip" file in your Downloads folder (either by double-clicking it, or by right-clicking, depending on your operating system). This will create a new folder called "18IHIW.T1D-main" in your Downloads folder.
Move the AZ_18IHIW.txt, 18IHIW.R, AZ.Locus.rda, AZ.HWE.rda, AZ.DRB1_0301_negative.rda and BDStrat.rda files in the Downloads/18IHIW.T1D-main folder to your Documents/18IHIW folder.

To use these files, Workshop participants should have R, Rstudio Desktop, and BIGDAWG installed on their computers.

**However, there is an issue with BIGDAWG in the current version if R (v4.2.0). If you have just (or recently) installed R and R studio, you will not need to run _BIGDAWG_.**
**Instead, you can use the AZ.Locus.rda, AZ.HWE.rda and AZ.DRB1_0301_negative.rda BIGDAWG output objects bundled with this repository.**

To install R on your computer, follow the instructions for linux, macOS, or Windows [here](https://cran.r-project.org).

Once R is installed, install Rstudio Desktop on your computer using the instructions for your variety of linux, macOS or Windows [here](https://www.rstudio.com/products/rstudio/download/).

Once Rstudio Desktop has been installed, open Rstudio desktop (by double clicking the icon).

Install BIGDAWG in the Rstudio Desktop environment by doing the following:

1. Go to the Rstudio Tools pulldown menu (at the top of your screen), and select "Install Packages ...".
2. In the popup window, select "Repository (CRAN)" in the "Install from:" menu, and enter "BIGDAWG" in the "Packages (separate multiple with space or comma):" menu.
3. Click the "Install" button.
4. Installing BIGDAWG will automatically install all of the other R packages that BIGDAWG requires. This may take a few minutes, and may generate a lot of text in your Rstudio console. This is normal. If you are asked 'Do you want to install from sources the package which needs compilation?', enter 'no'.

---
# Azerbaijan T1D Dataset
Analyses of the associated dataset (AZ_18IHIW.txt) were published in [2018 by Ahdamov et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6384092/), as part of a study of Type-1-Diabetes (T1D) in Azerbaijan.

The dataset is a tab-delimtied text file that consists of four columns, *Sample_ID*, *Status*, *DRB1_1* and *DRB1_2*, containing data for 315 individuals, with each line of the file containing all of the data for one individual. 

The 106 Azerbaijani individuals with T1D are indicated by the value of 1 in the *Status* column, and the 209 Azerbaijani individuals without T1D (controls) are indicated by the value of 0 in the *Status* column.

The DRB1 genotypes are represented in two columns, *DRB1_1* and *DRB1_2*, with one two-field DRB1 allele name appearing in each column. 

If genotype data for other loci were available for this dataset, they would appear in other pairs of columns. For example, if DQB1 data were available, the DQB1 allele names for each individual would appear in *DQB1_1* and *DQB1_2* columns.

If HLA genotype data were missing for any individual, the columns for that individual would contain either '****' or 'NA' values.

When analyzing HLA data like these with BIGDAWG, it is important to ensure that there are no extraneous spaces associated with any of the data. Allele names with leading or trailing spaces (' 02:01 ') may be treated as distinct from allele names without any associated spaces ('02:01'). 

# RMMLR
A rapid gene-based genome-wide association test in nuclear families.
The R code RMMLR_nuclearfamilies.r implements a rapid gene-based genomewide association test in nuclear families. The current version can handle families with two offspring. The 
offspring could MZ twins, DZ twins or biological offspring, adoptees. The program can also handle unrelated individuals. The input files for RMMLR_nuclearfamilies.r: 
 A file listing phenotypes and covariates (Phenotypefile.csv) 
 A genetic data file with individual level genotype data (Genotypefile.csv)
A file listing the gene names and the corresponding genetic variants within each gene (genlist.txt)
 The RMMLR_nuclearfamilies.r requires three functions that are listed in manDep.r,manova_func.r and useRFGLS.r
The code also requires R package RFGLS (https://cran.r-project.org/web/packages/RFGLS/index.html)

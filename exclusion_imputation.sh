#!/bin/bash -l
#$ -pe smp 2 -R y -binding linear:2
#$ -l h_rt=1:00:00
#$ -l s_rt=1:00:00
#$ -l h_vmem=8G
#$ -j y


source ~/.my.bashrc

cdpr
Rscript /medpop/esp/mkaminen/CAD_prs_proteomics_DM/2_primary_association_analyses.R
#Rscript /medpop/esp/mkaminen/CAD_prs_proteomics_DM/1_exclusion_imputation.R
#Rscript /medpop/esp/mkaminen/CAD_prs_proteomics_DM/check_baseline_assoc.R


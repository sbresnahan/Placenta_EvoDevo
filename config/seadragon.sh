#!/bin/bash
# Trusted infrastructure settings, following Placenta_xQTL_atlas / Placenta_POE.
# Data/reference/output paths are in the explicit cohort YAML, not here.
export POE_MODULE_INIT="${POE_MODULE_INIT:-/etc/profile.d/modules.sh}"
export POE_CONDA="${POE_CONDA:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
export POE_SAMTOOLS_ENV="${POE_SAMTOOLS_ENV:-samtools-1.16.1}"
export POE_BCFTOOLS_ENV="${POE_BCFTOOLS_ENV:-bcftools-1.16}"
export POE_ORNAMENTS_ACTIVATE="${POE_ORNAMENTS_ACTIVATE:-/rsrch5/home/epi/bhattacharya_lab/software/ornaments/bin/activate}"
export POE_R_SIF="${POE_R_SIF:-/risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif}"
export POE_R_LIBRARY="${POE_R_LIBRARY:-/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}"
export POE_SINGULARITY="${POE_SINGULARITY:-singularity}"

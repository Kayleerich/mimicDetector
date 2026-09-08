# ===============================================================================
# filt_lcrs.py : call functions to filter LCRs for mimicDetector
# Copyright (C) 2024 Kaylee D. Rich
# Read the COPYING file for license information.
# ================================================================================

from psychoscope import GreaterMimicDetection

config_dict = snakemake.config
config_dict['patho_name'] = snakemake.wildcards.pathogen
gmd = GreaterMimicDetection(**config_dict)
filt_lcr_file = f"{config_dict['outdir']}/{config_dict['patho_name']}/{config_dict['patho_name']}.{config_dict['runidIII']}_paired_mimics.tsv"
gmd.lcr_filter(snakemake.input.qfilt, 
               snakemake.input.segbed, 
               snakemake.input.hfa, 
               snakemake.input.pfa,
               filt_lcr_file) 

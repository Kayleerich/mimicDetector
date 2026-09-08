# ===============================================================================
# cnvrt_lcrs.py : convert LCR files for mimicDetector
# Copyright (C) 2024 Kaylee D. Rich
# Read the COPYING file for license information.
# ================================================================================

from psychoscope import Chipset

config_dict = snakemake.config
chip = Chipset(**config_dict)
seg_bed = chip.make_bed(str(snakemake.input.seg), "segmasker") 
try:
    assert(str(seg_bed) == str(snakemake.output.segbed))
except AssertionError:
    err=f"File names {seg_bed} and {snakemake.output.segbed} are not the same"
    raise AssertionError(err)

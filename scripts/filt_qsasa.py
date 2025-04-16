from psychoscope import MimicDetectionI
from psychoscope import MimicDetectionII

config_dict = snakemake.config
config_dict['patho_name'] = snakemake.wildcards.pathogen
mdI = MimicDetectionI(**config_dict)

filtfile = mdI.filter_blast_bitscore(str(snakemake.input.hblast), str(snakemake.input.cblast))
merge_dict = mdI.merge_ranges(filtfile)

mdII = MimicDetectionII(**config_dict)
all_qsasa_df = mdII.region_avg_qsasa(merge_dict, snakemake.params.ppops, snakemake.params.hpops) 
all_qsasa_file, filt_qsasa_file = mdII.qsasa_filt(all_qsasa_df) 

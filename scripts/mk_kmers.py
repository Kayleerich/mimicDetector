from psychoscope import Chipset

config_dict = snakemake.config
config_dict['patho_name'] = snakemake.wildcards.pathogen
patho_obj = Chipset(**config_dict)
kfa=f"{config_dict['outdir']}/{config_dict['patho_name']}/{config_dict['patho_name']}-{config_dict['k']}mers.fasta" 
patho_obj.fragment_fasta(snakemake.input.fa, kfa)
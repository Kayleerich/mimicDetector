from pathlib import Path
import yaml
import os, inspect
import glob
from time import strftime
import pandas as pd

def dir_empty(dir_path):
        return not next(os.scandir(dir_path), None)

def write_log(log_file, log_msg):
    with open(log_file, "a") as lg:
        lg.write(log_msg)
    return

def load_config(yml):
    if Path(yml).exists():
        with open(Path(yml), 'r') as f:
            config_dict = yaml.safe_load(f)
    else:
        err_msg = ''.join(['Cannot find the configuration file ',
                           str(Path(yml)), 
                           '\n\tCheck that you specified the correct path or, ', 
                           'to make the configuration file, see mimicDetector/mimic_configuration.py --help  '])
        raise FileNotFoundError(err_msg)
    log_file = f"{Path(yml).parent}/{config_dict['fileid']}{config_dict['date']}.log" 
    return config_dict, log_file

def make_config_file(args):
    def _check_output(outdir, fileid, force):
        if not outdir:
            outdir = str(Path(f'./mimicDetector_results/{str(strftime("%b%d_%Y")).strip("_")}').resolve())
        else:
            outdir = str(Path(outdir).resolve())

        if not Path(f'{outdir}/config.yaml').exists(): # ! json
            try:
                Path(outdir).mkdir(parents=True)
            except FileExistsError:
                pass
            except PermissionError:
                err_msg = f'Unable to create "{outdir}", permission denied'
                raise PermissionError()
        elif not force:
            err_msg = ''.join(['You are trying to overwrite currently existing output files (config file: ', 
                            str(Path(f'{outdir}/config.yaml')), 
                            ')\n\tUse --force to confirm you want to overwrite these files'])
            raise FileExistsError(err_msg)
        config_file = f'{outdir}/config.yaml'
        log_file = f'{outdir}/{fileid}{str(strftime("%b%d_%Y")).strip("_")}.log'
        return outdir, config_file, log_file
    
    config_dict = {}
    ## to get values from or update a given config file: 
    ## (to update the same file do not give an output directory and use --force)
    if args.yaml:
        with open(Path(args.yaml), 'r') as f:
            config_file_dict = yaml.safe_load(f) #json.load(f)
        config_dict = {k:'' for k in config_file_dict.keys()}
        match_args_dict = { ## 'argparse_key': 'yml_key'
            'indir':'indir',
            'outdir':'outdir',
            'pathogen_species':'patho_species', 
            'control_species':'control_name',
            'host_species':'host_name',
            'fileid':'fileid',
            'log_file':'log_file',
            'k_size':'k',
            'min_bitscore':'bit_min',
            'bitscore_diff':'bit_diff',
            'min_evalue':'min_e',
            'min_qsasa':'qsasa',
            'max_lcr':'lcr',
            'max_threads': 'max_threads',
            'min_unmasked':'n', 
            'mask':'mask',
            }
        for arg_key, j_key in match_args_dict.items(): 
            if arg_key != 'outdir' and arg_key != 'fileid':
                if not vars(args)[arg_key]:
                    config_dict[j_key] = config_file_dict[j_key]
                else:
                    config_dict[j_key] = vars(args)[arg_key]
            if arg_key == 'outdir':
                outdir = args.outdir if args.outdir else config_file_dict[j_key]
                if args.force:
                    config_dict['fileid'] = config_file_dict['fileid']
                else:
                    try:
                        config_dict['fileid'] = f"{args.fileid}.{config_dict['k']}mers_" if args.fileid else f"{config_dict['host_name']}.{config_dict['k']}mers_"
                    except KeyError:
                        config_dict['fileid'] = f"{args.fileid}.{config_file_dict['k']}mers_" if args.fileid else f"{config_file_dict['host_name']}.{config_file_dict['k']}mers_" # !! check if this can still throw exception
                config_dict['date'] = str(strftime("%b%d_%Y"))
                config_dict['log_file'] = f"{outdir}/{config_dict['fileid']}{config_dict['date']}.log"
                config_dict['outdir'], config_file, config_dict['log_file'] = _check_output(outdir, config_dict['fileid'], args.force)
    else:
        config_dict['host_name'] = args.host_species if args.host_species else [d.split('/')[-1].split('.')[0] for d in glob.glob(f"{config_dict['indir']}/host/*.fa*")][0]
        config_dict['k'] = int(args.k_size) if args.k_size else int(12)
        config_dict['fileid'] = f"{args.fileid}.{config_dict['k']}mers_" if args.fileid else f"{config_dict['host_name']}.{config_dict['k']}mers_"
        config_dict['date'] = str(strftime("%b%d_%Y"))
        config_dict['outdir'], config_file, config_dict['log_file'] = _check_output(args.outdir, config_dict['fileid'], args.force)
        config_dict['indir'] = str(Path(args.indir).resolve()) if args.indir else Path.cwd().resolve()
        config_dict['indir_rel'] = str(os.path.relpath(Path(args.indir).resolve(), start=Path(inspect.getfile(inspect.currentframe())).resolve().parent))
        config_dict['mdpath'] = str(Path(inspect.getfile(inspect.currentframe())).parent)
        config_dict['n'] = int(round(args.min_unmasked)) if args.min_unmasked else int(round(config_dict['k']/2))
        config_dict['mask'] = args.mask if args.mask else 'lowercase'
        config_dict['bit_min'] = int(args.min_bitscore) if args.min_bitscore else int(30)
        config_dict['bit_diff'] = int(args.bitscore_diff) if args.bitscore_diff else int(2)
        config_dict ['min_e'] = args.min_evalue if args.min_evalue else 0.01
        config_dict['qsasa'] = args.min_qsasa if args.min_qsasa else 0.75
        config_dict['lcr'] = args.max_lcr if args.max_lcr else 0.50
        config_dict['control_species'] = args.control_species if args.control_species else [d.split('/')[-1].split('.')[0] for d in glob.glob(f"{config_dict['indir']}/controls/*.fa*")]
        config_dict['control_name'] = args.control_name if args.control_name else f"{len(args.control_species)}CtrlSpp" if len(args.control_species) > 1 else args.control_species[0]
        config_dict['patho_species'] = args.pathogen_species if args.pathogen_species else [d for d in filter(os.path.isdir, 
                                                                                                              os.listdir(config_dict['indir'])) if d not in ['host', 'controls']]
    config_dict['e_str'] = str('{:,g}'.format(config_dict ['min_e'])).split('.')[1]
    config_dict['b_str'] = str(config_dict['bit_diff'])
    config_dict['q_str'] = str('{:.2f}'.format(round(config_dict['qsasa'], ndigits=2))).split('.')[1]
    config_dict['l_str'] = str('{:.2f}'.format(round(config_dict['lcr'], ndigits=2))).split('.')[1]
    config_dict['runidI'] = f"{config_dict['fileid']}b{str(config_dict['bit_diff'])}_e{config_dict['e_str']}" 
    config_dict['runidII'] = f"{config_dict['runidI']}_q{config_dict['q_str']}"
    config_dict['runidIII'] = f"{config_dict['runidII']}_l{config_dict['l_str']}"
    config_dict['max_threads'] = args.max_threads if args.max_threads else 8
    characters = 0

    def _fasta_reader(file):
        fasta_df = pd.read_csv(file, sep='>', lineterminator='>', header=None)
        fasta_df[['Accession', 'Sequence']] = fasta_df[0].str.split('\n', n=1, expand=True)
        fasta_df.drop(0, axis=1, inplace=True)
        fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True)
        return fasta_df
    
    # combine all control species fasta files into single file
    controls_fasta = [Path(f"{config_dict['indir']}/controls", f'{f}.fasta') for f in config_dict['control_species']]
    df = pd.concat(_fasta_reader(i) for i in controls_fasta)
    df['Accession'] = '>' + df['Accession']
    df.to_csv(f"{config_dict['indir']}/controls/{config_dict['control_name']}.fasta", sep='\n', index=None, header=None)

    # check that host and control fasta files exist and get combined database size
    for db in [('host', 'host_name'), ('controls', 'control_name')]:
        db_file=f"{config_dict['indir']}/{db[0]}/{config_dict[db[1]]}.fasta"
        if Path(db_file).exists():
            with open(db_file, "r") as f:
                for line in f:
                    if not line.startswith(">"):
                        line = line.strip(os.linesep)
                        characters=characters+len(line)
        else:
            err_msg = ''.join([f'Cannot find the {db[0]} fasta file "{db_file}"',
                               f"{config_dict['indir']}/{db[0]}/{config_dict[db[1]]}.fasta",
                                '\n\tCheck that you specified the correct path'])
            raise FileNotFoundError(err_msg)
    config_dict['dbsize'] = characters
    
    with open(Path(config_file), "w") as f:
        yaml.dump(config_dict, f) 
    return config_file

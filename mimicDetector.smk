import os
import glob
import pandas as pd

if not workflow.configfiles:
    err_msg = ['Please specify the configuration file using --configfile path/to/outdir/config.yaml', 
                '\nTo make the configuration file, see mimicDetector/mimic_configuration.py --help']
    raise FileNotFoundError(''.join(err_msg))

PATHOGENS = config['patho_species'] 

rule all:
    input: 
        expand("{outdir}/{pathogen}/{pathogen}.{f}b{b}_e{e}_q{q}_l{l}_paired_mimics.tsv", 
                pathogen=PATHOGENS, outdir=config['outdir'], f=config['fileid'], b=config["b_str"], 
                e=config["e_str"], q=config["q_str"], l=config["l_str"])


# make pathogen k-mers: Chipset
def make_kmers_input(wildcards):
    fa_pattern = f"{config['indir']}/{wildcards.pathogen}/{wildcards.pathogen}*.fasta"
    if glob.glob(fa_pattern): 
        try:
            assert(len(glob.glob(fa_pattern)) == 1)
            return glob.glob(fa_pattern)[0]
        except AssertionError:
            err_msg = f'More than one fasta file for {wildcards.pathogen} found in {config["indir"]}/{wildcards.pathogen}/'
            raise LookupError(err_msg)
    else:
        err_msg = f'No fasta files for {wildcards.pathogen} found in {config["indir"]}/{wildcards.pathogen}/'
        raise FileNotFoundError(err_msg)

rule make_kmers:
    threads: 1
    input:
        fa=make_kmers_input
    output:
        kfa="{outdir}/{pathogen}/{pathogen}-{k}mers.fasta" 
    script:
        "scripts/mk_kmers.py"

# run blastp with pathogen k-mers as query and either host or control proteins as databse
rule control_blast_db:
    input:
        fa=f"{config['indir']}/controls/{config['control_name']}.fasta" 
    output:
        multiext(f"{config['indir']}/controls/{config['control_name']}.fasta", ".phr", ".pdb", ".pin", ".pjs", ".pog", ".pos", ".pot", ".psq", ".ptf", ".pto") 
    shell:
        "makeblastdb -in {input.fa} -dbtype prot -parse_seqids"

rule host_blast_db:
    input:
        fa=f"{config['indir']}/host/{config['host_name']}.fasta" 
    output:
        multiext(f"{config['indir']}/host/{config['host_name']}.fasta", ".phr", ".pdb", ".pin", ".pjs", ".pog", ".pos", ".pot", ".psq", ".ptf", ".pto") 
    shell:
        "makeblastdb -in {input.fa} -dbtype prot -parse_seqids"

rule run_host_kmer_blastp:
    threads: config['max_threads']
    resources:
        runtime="18000",
        mem_mb=150000,
        cpus_per_task=config['max_threads'],
    input:
        rules.host_blast_db.output,
        kmers="{outdir}/{pathogen}/{pathogen}-{k}mers.fasta"
    output:
        bph="{outdir}/{pathogen}/{pathogen}.{id}.{k}mers_blastp.out"
    params:
        db=f"{config['indir']}/host/{config['host_name']}.fasta",
        dbsize=config['dbsize'], 
    shell:
        "blastp -db {params.db} -query {input.kmers} -num_threads {threads} -outfmt 6 -out {output.bph} -evalue 1 -comp_based_stats 0 -matrix PAM30 -word_size 2 -ungapped -dbsize {params.dbsize}"

rule run_control_kmer_blastp:
    threads: config['max_threads']
    resources:
        runtime="18000",
        mem_mb=150000,
        cpus_per_task=config['max_threads'],
    input:
        rules.control_blast_db.output,
        kmers="{outdir}/{pathogen}/{pathogen}-{k}mers.fasta"
    output:
        bpc="{outdir}/{pathogen}/{pathogen}." + config['control_name'] + ".{k}mers_blastp.out"
    params:
        db=f"{config['indir']}/controls/{config['control_name']}.fasta",
        dbsize=config['dbsize'], 
    shell:
        "blastp -db {params.db} -query {input.kmers} -num_threads {threads} -outfmt 6 -out {output.bpc} -evalue 1 -comp_based_stats 0 -matrix PAM30 -word_size 2 -ungapped -dbsize {params.dbsize}"


# get qsasa for all pathogen and host proteins  
def get_pstruct_files(wildcards):
    indir=f'{config["indir"]}/{wildcards.pathogen}/structures/'
    pstruct = [Path(''.join([indir, f])) for f in os.listdir(indir) if Path(''.join([indir, f])).exists()] 
    return pstruct

rule run_pathogen_popscomp: 
    input:
        files=get_pstruct_files
    output:
        # flag file made only after popscomp has run on all pathogen structures
        flag='{outdir}/{pathogen}/ppops_flag.out'
    params: 
        indir=config["indir"]+'/{pathogen}/structures',
        outdir=config["outdir"]+'/{pathogen}/pops',
    run:
        Path(params.outdir).mkdir(parents=True, exist_ok=True)
        outfile_str = '/'.join([str(params.outdir), 'pops_'])
        for file in input.files:
            z = ('--zipped', '.gz') if Path(file).suffix == '.gz' else ('', '')
            outfile = ''.join([outfile_str, str(Path(file).name).replace(f'.pdb{z[1]}', '.out')])
            if os.stat(file).st_size:
                shell(f'~/POPScomp/POPSC/src/pops --pdb {file} --residueOut --popsOut {outfile} {z[0]}') 
            else:
                with open(outfile, "w") as o:
                    o.write('')
        struct_prots=[str(Path(p).stem).split('.')[0] for p in os.listdir(params.indir)] 
        pops_prots=[str(Path(p).stem).split('.')[0].replace('pops_', '') for p in os.listdir(params.outdir)]
        if sorted(pops_prots) == sorted(struct_prots):
            with open(output.flag, "w") as f:
                f.write(f"Number of protein structures processed: {len(pops_prots)}")

def get_hstruct_files(wildcards):
    indir=f'{config["indir"]}/host/structures/'
    hstruct = [Path(''.join([indir, f])) for f in os.listdir(indir) if Path(''.join([indir, f])).exists()] 
    return hstruct

rule run_host_popscomp:
    input:
        files=get_hstruct_files
    output:
        # flag file made only after popscomp has run on all host structures
        flag=f'{config["outdir"]}/{config["host_name"]}/hpops_flag.out'
    params: 
        indir=f'{config["indir"]}/host/structures',
        outdir=f'{config["outdir"]}/{config["host_name"]}/pops',
    run:
        Path(params.outdir).mkdir(parents=True, exist_ok=True)
        outfile_str = '/'.join([str(params.outdir), 'pops_'])
        for file in input.files:
            z = ('--zipped', '.gz') if Path(file).suffix == '.gz' else ('', '')
            outfile = ''.join([outfile_str, str(Path(file).name).replace(f'.pdb{z[1]}', '.out')])
            if os.stat(file).st_size:
                shell(f'~/POPScomp/POPSC/src/pops --pdb {file} --residueOut --popsOut {outfile} {z[0]}') 
            else:
                with open(outfile, "w") as o:
                    o.write('')
        struct_prots=[str(Path(p).stem).split('.')[0] for p in os.listdir(params.indir)] 
        pops_prots=[str(Path(p).stem).split('.')[0].replace('pops_', '') for p in os.listdir(params.outdir)]
        if sorted(pops_prots) == sorted(struct_prots):
            with open(output.flag, "w") as f:
                f.write(f"Number of protein structures processed: {len(pops_prots)}")
        

# blast minimum bitscore difference and qsasa filtering: MimicDetectionI and MimicDetectionII 
def get_ppops_files(wildcards):
    indir=f'{config["indir"]}/{wildcards.pathogen}/structures/'
    outstr=f'{wildcards.outdir}/{wildcards.pathogen}/pops/pops_'
    pstruct=[''.join([indir, f]) for f in os.listdir(indir) if Path(''.join([indir, f])).exists()]
    ppops=[p.replace(indir, outstr).replace('.pdb', '.out') for p in pstruct]
    return ppops

def get_hpops_files(wildcards):
    indir=f'{config["indir"]}/host/structures/'
    outstr=f'{wildcards.outdir}/host/pops/pops_'
    hstruct=[''.join([indir, f]) for f in os.listdir(indir) if Path(''.join([indir, f])).exists()]
    hpops=[p.replace(indir, outstr).replace('.pdb', '.out') for p in hstruct]
    return hpops

rule filter_qsasa:
    threads: 1
    input:
        pflag=rules.run_pathogen_popscomp.output, 
        hflag=rules.run_host_popscomp.output,
        hblast=rules.run_host_kmer_blastp.output,
        cblast=rules.run_control_kmer_blastp.output,
    output:
        bfilt="{outdir}/{pathogen}/{pathogen}.{id}.{k}mers_b{b}_e{e}_hcfiltered_blastp.out",
        qsasa="{outdir}/{pathogen}/{pathogen}.{id}.{k}mers_b{b}_e{e}_qsasa_averages.tsv",
        qfilt="{outdir}/{pathogen}/{pathogen}.{id}.{k}mers_b{b}_e{e}_q"+f"{config['q_str']}_filtered.tsv"
    params:
        ppops=get_ppops_files, 
        hpops=get_hpops_files, 
    script:
        "scripts/filt_qsasa.py"


# low complexity sequence filtering: GreaterMimicDetection
rule get_host_lcrs:
    threads: 1
    input:
        f"{config['indir']}/host/{config['host_name']}.fasta"
    output:
        f"{config['indir']}/host/{config['host_name']}-segmask-intervals.txt"
    shell:
        """
        segmasker -in {input} -infmt fasta -parse_seqids -out {output}
        """

rule convert_host_lcrs:
    threads: 1
    input:
        seg=rules.get_host_lcrs.output
    output:
        segbed=f"{config['indir']}/host/{config['host_name']}-segmask-intervals.bed"
    script:
        "scripts/cnvrt_lcrs.py"

rule filter_lcrs:
    threads: 1
    input:
        qfilt=rules.filter_qsasa.output.qfilt, 
        hfa=f"{config['indir']}/host/{config['host_name']}.fasta", 
        pfa=config['indir']+"/{pathogen}/{pathogen}.fasta",
        segbed= rules.convert_host_lcrs.output.segbed 
    output:
        hseqbed=temp("{outdir}/{pathogen}/{pathogen}.{id}.{k}mers_b{b}_e{e}_q{q}_l{l}_seqs_host"),
        pseqbed=temp("{outdir}/{pathogen}/{pathogen}.{id}.{k}mers_b{b}_e{e}_q{q}_l{l}_seqs_pathogen"),
        pairs="{outdir}/{pathogen}/{pathogen}.{id}.{k}mers_b{b}_e{e}_q{q}_l{l}_paired_mimics.tsv" 
    script:
        "scripts/filt_lcr.py"

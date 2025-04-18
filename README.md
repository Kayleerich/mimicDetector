# mimicDetector
mimicDetector is a Snakemake pipeline for identifying short regions of pathogen proteins that share high sequence similarity to proteins of their host.  

## Running the pipeline
The pipeline requires >=Python3.6 as well as pre-installation of [POPScomp](https://github.com/Fraternalilab/POPScomp). The other programs/packages can be installed into a new conda environment using the [snake_mimics.yml](env/snake_mimics.yml) file: `conda env create -f snake_mimics.yml`

### Basic usage
The pipeline is contained in a Snakemake workflow, and requires a custom configuration file to be passed via `--configfile`:  
`snakemake -s mimicDetector.smk --configfile path/to/outdir/config.yaml`

### Create the configuration file for mimicDetector
All of the parameters required to run the pipeline are found in the config.yaml file in the output directory, which needs to be created using [mimic_configuration.py](mimic_configuration.py). The only parameters that are absolutely required are a single species for `--host_species` and at least one species for each of `--pathogen_species` and `--control_species`: `python mimic_configuration.py -p PATHOGEN(S) -s HOST -c CONTROL(S)`  

Example: specify multiple pathogen species. Each of these pathogen species will be run through the pipeline separately  
`python mimic_configuration.py -p pathogen_speciesA pathogen_speciesB pathogen_speciesC -s host_species -c control_species`  

Example: specify multiple control species with a specific control group name. This is the name that will be used in certain output files and the name of the concatenated FASTA controls file used for all pathogens  
`python mimic_configuration.py -p pathogen_species -s host_species -c control_speciesA control_speciesB --control_name controlGroup1`  

Example: specify input and output directories (see below for input directory organization)  
`python mimic_configuration.py -p pathogen_species -s host_species -c control_species -i ../path/to/input/directory -o path/to/output/directory`

All parameter arguments for mimic_configuration.py:
```
  -h, --help                show this help message and exit
  -p, --pathogen_species    List of pathogen species names (space delim)
  -s , --host_species       Host species name
  -c, --control_species     List of control species names (space delim)
  --control_name            Abbreviated name for controls set
  -i, --indir               Directory containing data (default is current directory)
  -o, --outdir              Output directory (default is <cwd>/mimicDetector_results/<date>/)
  -f, --fileid              Name/identifier for output files (default is host species name)
  -k, --k_size              INT, Length of k-mer, i.e. sequence fragment, to use (default is k=12)
  -b, --bitscore_diff       NUM >=0, Minimum difference in bitscore between pathogen-host and
                             pathogen-control blastp hits (default is b=2)
  --min_bitscore            NUM >=0, Minimum bitscore allowed for pathogen-host blastp hits 
                             (default is 30)
  -e, --min_evalue          FLOAT between 0 and 1, Maximum E-value allowed for pathogen-host blastp 
                             hits (default is e=0.01)
  -q, --min_qsasa           FLOAT between 0 and 1, Minimum average solvent accessibility required for 
                             mimicry candidates (default is q=0.75)
  -l, --max_lcr             FLOAT between 0 and 1, Maximum low-complexity region overlap allowed for 
                             mimicry candidates (default is 50%, l=0.50)
  --mask                    Single character string (i.e. X), default is lowercase
  --min_unmasked            INT < k_size. Minimum number of unmasked residues in k-mer 
                             (default: k_size/2)
  -t, --max_threads         Maximum number of threads to use during workflow (default: 8)
  --force                   Overwrite previous configuration file
```

### Input directory and file organization
Ensure that sequences in FASTA files have unique names. Sequence names should not contain spaces, `\` or `:`. Uniprot identifiers work well (for example: >sp|P04004|VTNC_HUMAN Vitronectin to >P04004), and are easily renamed using bioawk:  
`bioawk -c fastx '{split($name, a, "|"); print ">"a[2]"\n"$seq}' sequences.fasta`  
Structure files must also be named the same as the proteins in the corresponding FASTA file (e.g. the protein with the FASTA header ">P04004" corresponds to the structure file "structures/P04004.pdb(.gz)"). 

Since there are a large number of files required for multiple species, file names and organization are important.  
In the specified input directory, there should be at minimum three directories: two named "host" and "control", and one for each of the pathogen species specified by `-p`.  
Each of the pathogen directories and the host directory should contain a FASTA file named to match the species and a directory named "structures" containing `.pdb` or `.pdb.gz` files for each protein structure. 
The "controls" directory should contain a FASTA file for each of the control species specified by `-c`. 
```
├── indir  
│   ├── host  
|   │   ├── host_name.fasta  
|   │   ├── structures  
|   │   │   ├── prot1.pdb  
|   │   │   ├── prot2.pdb  
|   │   │   ├── ...  
│   ├── pathogenA_name
|   │   ├── pathogenA_name.fasta
|   │   ├── structures
|   │   │   ├── protA1.pdb
|   │   │   ├── protA2.pdb
|   │   │   ├── ...
│   ├── pathogenB_name
|   │   ├── pathogenB_name.fasta
|   │   ├── structures
|   │   │   ├── protB1.pdb
|   │   │   ├── protB2.pdb
|   │   │   ├── ...
:
│   ├── controls
|   │   ├── control1_name.fasta
|   │   ├── control2_name.fasta
|   │   ├── ...
```
  
After the pipeline has finished, there will be a directory for each of the pathogen species as well as for the host species. 
Files are named according to the parameters specified in the config.yaml file.  
The final output (`*_paired_mimics.tsv`) contains 10 columns: `pathogen_prot, pathogen_start, pathogen_end, pathogen_sequence, pathogen_qsasa, host_prot, host_start, host_end, host_sequence, host_qsasa` where pathogen_start/pathogen_end and host_start/host_end are the start and end protein coordinates (inclusive) for each mimicry candidate, and pathogen_qsasa/host_qsasa are the average QSASA for all residues in those sequence ranges.
```
├── outdir
│   ├── config.yaml
│   ├── host_name
|   │   ├── hpops_flag.out
|   │   ├── pops
|   │   │   ├── pops_prot1.out
|   │   │   ├── pops_prot2.out
|   │   │   ├── ...
│   ├── pathogenA_name
|   │   ├── ppops_flag.out
|   │   ├── pathogenA_name-12mers.fasta
|   │   ├── pathogenA_name.control_name.12mers_blastp.out 
|   │   ├── pathogenA_name.fileID.12mers_blastp.out
|   │   ├── pathogenA_name.fileID.12mers_b2_e01_hcfiltered_blastp.out
|   │   ├── pathogenA_name.fileID.12mers_b2_e01_qsasa_averages.tsv
|   │   ├── pathogenA_name.fileID.12mers_b2_e01_q75_filtered.tsv
|   │   ├── pathogenA_name.fileID.12mers_b2_e01_q75_l50_paired_mimics.tsv
|   │   ├── pops
|   │   │   ├── pops_protA1.out
|   │   │   ├── pops_protA2.out
|   │   │   ├── ...
│   ├── pathogenB_name
|   │   ├── ...
```

## Data selection recommendations and ramblings
First, data selection is important. As the adage goes, “garbage in, garbage out.” Appropriate controls include species that are more closely related to the pathogen than to the host, since if the control set includes species that are more closely related to the host than to the pathogen, it is likely that this will remove host taxon-specific homologous sequences. 
After data selection, the pipeline filtering thresholds can be modified to further tailor the analysis. Any adjustment to the filtering thresholds will either increase or decrease the final number of mimicry candidates, the desirability of which will likely depend on the species of interest. For example, if the pathogen species is only distantly related to its host, decreasing stringency in certain filtering steps to increase the number of final mimicry candidates may be desired. Conversely, if the pathogen and host are closely related (i.e. helminths and humans), or if the number of mimicry candidates is unmanageable, an increase in stringency during certain filtering steps may be desired.
 
To increase the number of mimicry candidates obtained, the first parameter change recommended is the solvent accessibility threshold (q). Decreasing this value will increase the number of candidates in the final results without compromising quality (i.e. increased low complexity sequences). A q threshold greater than 0.5 indicates that the majority of the mimicry region is accessible for protein-protein interactions. However, low mean solvent accessibility for a given mimicry candidate may indicate that fewer residues are involved in mimicry interactions, so these candidates should be further screened for similarity between the accessible residues. Changing either E-value (e) or LCR content (l) thresholds will affect the quality of the final mimicry candidates. Decreasing l will allow for more mimicry candidates that overlap low complexity sequences in the final results. Increasing e will have a similar effect by permitting sequences with reduced statistical relevance. 
Decreasing either bits difference threshold (b) or k-mer size (k) is not recommended. Lowering b will result in mimicry candidates closely related to a control, and decreasing k will eliminate any statistical relevance of the BLASTP E-values. However, k may be increased at the expense of increased time and memory usage, though this should not be done without also adjusting b and e. To decrease in the number of final mimicry candidates, increasing b and/or decreasing e is recommended. If the pathogen is closely related to host, we recommend increasing b to ensure greater distinction between the potential mimic and a control sequence. 

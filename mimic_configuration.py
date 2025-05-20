import argparse
from scripts.utils import make_config_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='%(prog)s [-h] -p PATHOGENS -s HOST -c CONTROLS [-i DIR -o DIR -f FILEID -k INT -b INT -e FLOAT -q FLOAT -l FLOAT]', 
                                     description='Create configuration file for mimicDetector')
    parser.add_argument('-p', '--pathogen_species', metavar='', nargs='+', required=True, 
                        help='List of pathogen species names (space delim)')
    parser.add_argument('-s', '--host_species', metavar='', type=str, required=True, 
                        help='Host species name')
    parser.add_argument('-c', '--control_species', metavar='', nargs='+', required=True, 
                        help='List of control species names (space delim)')
    ## group=input/output control
    parser.add_argument('-i', '--indir', metavar='INDIR', required=False, 
                        help='Directory containing all FASTA files (default is current directory)') 
    parser.add_argument('-o', '--outdir', metavar='', type=str, required=False, 
                        help='Directory to save output files (default is <cwd>/mimicDetector/<date>/)') 
    parser.add_argument('-f', '--fileid', metavar='', type=str, required=False, 
                        help='Name/identifier for output files')
    parser.add_argument('--control_name', metavar='', type=str, required=False, 
                        help='Abbreviated name for controls set')
    ## group=parameters
    parser.add_argument('-k', '--k_size', metavar='', type=int, required=False, 
                        help='INT, Length of k-mer, i.e. sequence fragment, to use (default is k=12)')
    parser.add_argument('-b', '--bitscore_diff', metavar='', type=int, required=False, 
                        help='NUM >=0, Minimum difference in bitscore between pathogen-host and pathogen-control blastp hits (default is b=2)')
    # parser.add_argument('--min_bitscore', metavar='', type=int, required=False, 
    #                     help='NUM >=0, Minimum bitscore allowed for pathogen-host blastp hits (default is 30)') 
    parser.add_argument('-e', '--min_evalue', metavar='', type=float, required=False, 
                        help='FLOAT between 0 and 1, Maximum E-value allowed for pathogen-host blastp hits (default is e=0.01)')
    parser.add_argument('-q', '--min_qsasa', metavar='', type=float, required=False, 
                        help='FLOAT between 0 and 1, Minimum average solvent accessibility required for mimicry candidates (default is q=0.75)')
    parser.add_argument('-l', '--max_lcr', metavar='', type=float, required=False, 
                        help='FLOAT between 0 and 1, Maximum low-complexity region overlap allowed for mimicry candidates (default is l=0.50)')
    parser.add_argument('--mask', metavar='', type=str, required=False, 
                        help='Single character string (i.e. X, default is lowercase)')
    parser.add_argument('--min_unmasked', metavar='', type=int, required=False, 
                        help='INT < k_size, Minimum number of unmasked aa in k-mer (default: k_size/2)')
    parser.add_argument('-t', '--max_threads', metavar='', type=int, required=False, 
                        help='Maximum number of threads to use during workflow (default: 8)')
    parser.add_argument('--force', action='store_true', required=False, 
                        help='Overwrite previous output files')
    parser.add_argument('-y', '--yaml', metavar='', type=str, required=False, 
                        help='input configuration file') 

    args = parser.parse_args()
    config_file = make_config_file(args)
    print(f'Configuration file saved to {config_file}')

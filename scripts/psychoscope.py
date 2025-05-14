import numpy as np
from utils import write_log
import re
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
import glob
from pybedtools import BedTool
import sys
import os

class Chipset():
    __chipset_vars = ('outdir', 'patho_name', 'host_name', 'control_name', 'k', 
                      'n', 'mask', 'bit_min', 'bit_diff', 'b_str', 'min_e', 
                      'e_str', 'qsasa', 'q_str', 'lcr', 'l_str', 'runid')
    """
    mimicDetector basic functions to prepare data for filtering steps
    
    fragment_fasta:     creates k-mer FASTA file or returns k-mer dictionary, input: FASTA file
    make_bed:           creates BED file from BLAST tabular output, SEGmasker interval output or Phobius short output
                            options are: 'blast-query', 'blast-subject', 'phobius', 'segmasker'
    merge_ranges_dict:  returns dictionary with pathogen-host predicted mimicry sites, input: BLASTP tabular file
    """

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            try:
                assert(k in self.__class__.__chipset_vars)
                setattr(self, k, v)
            except AssertionError:
                pass

    ### make sure these fns are the most recent/up to date 
    def fragment_fasta(self, fastafile, k_file):
        """
        Read fasta file and create k-mers for each sequence.
        Can take soft- or hard-masked sequences and remove k-mers that exceed the threshold given for number of masked redisues.
        By default k-mers containing more than half masked residues will be removed, assumes soft-masking to be lowercase letters.
        """
        if self.n == False:
            self.n = int(self.k/2)
        REGEX = re.compile(r'[a-z]')
        def _count_lower(s):
            return len(REGEX.findall(s))
        k_dict = {}
        with open(k_file, 'w') as fragfile:
            records = SeqIO.parse(Path(fastafile), 'fasta')
            for record in records:
                sequence = record.seq
                for i in range(0, len(sequence)):
                    seq_fragment = sequence[i:i+self.k]
                    header = f"{record.description}/{i}:{i+self.k}"
                    record.id = header
                    #remove fragments that have too many masked characters
                    if self.mask != 'lowercase':
                        if seq_fragment.count(self.mask) > self.n:
                            continue
                    elif _count_lower(str(seq_fragment)) > self.n:
                        continue
                    #remove fragments that are less than k amino acids in length
                    if len(seq_fragment) != self.k:
                        continue
                    else:
                        k_dict[header] = seq_fragment.upper()
                        fragfile.write(f">{header}\n{seq_fragment}\n")
        return
    
    def make_bed(self, infile, intype):
        """
        Make .bed file from BLAST tabular output, SEGmasker interval output or Phobius short output
        """
        bedfile=infile.replace('.txt', '.bed') 
        with open(Path(bedfile), 'w') as f:
            if "blast" in intype:
                with open(infile, "r") as blastfile:
                    for ln in blastfile: 
                        ln = ln.split()
                        if intype == "blast-subject":
                            bed = f"{ln[1]}\t{int(ln[8]) - 1}\t{int(ln[9])}"
                        else:
                            bed = f"{ln[0]}\t{int(ln[6]) - 1}\t{int(ln[7])}"
                        f.write(f"{bed}\n")
                    blastfile.close()

            elif intype == "phobius":
                with open(infile, 'r') as phobiusout:
                    for ln in phobiusout:
                        ln = ln.split()
                        if ln[2] == "Y":
                            start = re.search(r"n\d+-", ln[3]) 
                            end = re.search(r"-\d+c", ln[3])
                            start = start.group()
                            end = end.group()
                            start = int(start[1:-1]) - 1
                            end = int(end[1:-1]) - 1
                            signal = f"{ln[0]}\t{start}\t{end}" 
                            f.write(f"{signal}\n")
                    phobiusout.close()
            
            elif intype == "segmasker":
                with open(infile, 'r') as segintervals:
                    for ln in segintervals:
                        if ln.startswith(">"): 
                            name = ln.split("|")[1]
                            name = name.split(" ")[0]
                        else:
                            start = ln.split(" ")[0]
                            end = int(ln.split(" ")[2]) + 1
                            f.write(f"{name}\t{start}\t{end}\n")
                    segintervals.close()
            f.close()
        return bedfile
    

class MimicDetectionI():
    __md1_vars = ('fileid', 'log_file', 'patho_name', 'host_name', 'control_name', 
                  'k', 'bit_min', 'bit_diff', 'b_str', 'min_e', 'e_str', 'qsasa', 
                  'q_str', 'lcr', 'l_str', 'runidI', 'outdir', 'indir') 
    """
    mimicDetector k-mer comparison functions: bitscore and e-value filtering 
    filter_blast_bitscore(): returns filtered host blastp results from input pathogen-host and pathogen-control BLASTP tabular files
                             if control BLASTP file is too large: chunk=True
                             to write max control bitscores to file: write_max=True
                             to return full dataframe instead of filtered file name: outdf=True
    """
    def __init__(self, **kwargs): 
        for k, v in kwargs.items():
            try:
                assert(k in self.__class__.__md1_vars)
                setattr(self, k, v)
            except AssertionError:
                pass
    
    def filter_blast_bitscore(self, host_blast, control_blast, write_all=False, chunk=False, write_max=False, outdf=False):
        def _max_bitscore(blastfile, chunk=False, write=write_max):
            max_bits_file = Path(self.outdir, f'{self.patho_name}.controls_max_bitscores_blastp.out')
            def _max(df):
                df.columns = ['query', 'subject', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
                df = df.loc[df.groupby('query')['bitscore'].idxmax()]
                return df
            if chunk:
                df_list = []
                chunksize = 10 ** 6
                with pd.read_csv(blastfile, chunksize=chunksize, sep='\s+', header=None) as reader:
                    for chunk in reader:
                        chunkdf = _max(chunk) 
                        df_list.append(chunkdf)
                max_blast = pd.concat([df_list])
                max_blast = _max(max_blast)
            else:
                blastdf = pd.read_csv(blastfile, sep='\s+', header=None)
                max_blast = _max(blastdf)
            if write:
                max_blast.to_csv(max_bits_file, sep='\t', header=None, index=False)
                log_msg = f'Top control BLASTP results saved to {max_bits_file}\n'
                # write_log(self.log_file, log_msg)
            return max_blast
        
        cont_df = _max_bitscore(Path(control_blast), chunk=chunk, write=False)
        cont_df = cont_df.loc[cont_df.groupby('query')['bitscore'].idxmax()]
        host_df = pd.read_csv(Path(host_blast), sep='\s+', header=None)
        host_df.columns = ['query', 'subject', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        host_df = host_df.loc[host_df.groupby('query')['bitscore'].idxmax()]

        topbits_df = host_df.set_index('query').join(cont_df.set_index('query'), lsuffix='_host', rsuffix='_control')
        topbits_df = topbits_df.fillna(0)
        topbits_df['bits_diff'] = (topbits_df['bitscore_host'] - topbits_df['bitscore_control'])
        topbits_df['numIDs_diff'] = ((topbits_df['length_host'] - topbits_df['mismatch_host']) - (topbits_df['length_control'] - topbits_df['mismatch_control']))
        if self.bit_diff == 'stdev':
            self.bit_diff = round(topbits_df[(topbits_df['bitscore_control'] > 0) & (topbits_df['bitscore_host'] > 0)]['bits_diff'].std())
            self.b_str = str(self.bit_diff).split('.')[0]
            self.runidI = f'{self.fileid}b{self.b_str}_e{self.e_str}'
        filtered_df = topbits_df[(topbits_df['bits_diff'] >= float(self.bit_diff)) & (topbits_df['evalue_host'] <= float(self.min_e))]
        if write_all:
            unfilt_bits_file=Path(self.outdir, self.patho_name, f'{self.patho_name}.{self.fileid}_top_hits_unfiltered.tsv')
            topbits_df[(topbits_df['bitscore_host'] > 0)].to_csv(unfilt_bits_file, sep='\t', header=True, index=True)
        hcols = filtered_df.columns.str.contains(r'_host')
        bfilt_file = Path(self.outdir, self.patho_name, f'{self.patho_name}.{self.runidI}_hcfiltered_blastp.out')
        filtered_df.iloc[:, hcols].to_csv(bfilt_file, sep='\t', header=False, index=True) 
        if outdf:
            return filtered_df, self.bit_diff
        else:
            log_msg = f'Filtered BLASTP results saved to {bfilt_file}\n'
            # write_log(self.log_file, log_msg)
            return bfilt_file
        
    def merge_ranges(self, tsv_file):
        """
        Merge overlapping k-mers (BLASTP tab-separated results) into predicted mimicry sites. 
        Returns a dictionary.
        """
        def _is_overlaping(a, b):
            if (b[0] >= a[0] and b[0] <= a[1]):
                return True
            elif (a[0] >= b[0] and a[0] <= b[1]):
                return True
            elif (a[1] >= b[0] and a[1] <= b[1]):
                return True
            elif (b[1] >= a[0] and b[1] <= a[1]):
                return True
            else:
                return False
            
        def _merge_range(range1, range2):
            range_start = min(min(range1), min(range2))
            range_end = max(max(range1), max(range2))
            new_range = (range_start, range_end)
            return new_range
        
        coords_dict = {}
        with open(tsv_file, 'r') as file:
            for ln in file:
                ln = ln.split()
                pair = ln[0].split('/')[0] + "." + ln[1] 
                kcoords = ln[0].split('/')[1]
                kstart = int(kcoords.split(':')[0]) - 1
                kend = int(kcoords.split(':')[1])
                qstart = (kstart + int(ln[6]))
                qend = kend - (self.k - int(ln[7])) + 1
                query_range = [qstart, qend]
                try:
                    coords_dict[pair]['query'].append(query_range)
                except KeyError:
                    coords_dict[pair] = {'query':[], 'host':[]}
                    coords_dict[pair]['query'] = [query_range]

                hstart = (int(ln[8]) - 1)
                hend = int(ln[9]) + 1
                host_range = [hstart, hend]
                try:
                    coords_dict[pair]['host'].append(host_range)
                except KeyError:
                    coords_dict[pair]['host'] = [host_range]
        for pair in coords_dict.keys():
            coords_dict[pair]['num_aln'] = len(coords_dict[pair]['query'])
        
        merged_coords_dict = {}
        for pair in coords_dict.keys():

            query_lists = coords_dict[pair]['query']
            host_lists = coords_dict[pair]['host']
            range_num = coords_dict[pair]['num_aln']

            ## assign variables and add first range to range_lists
            q_range = tuple(query_lists[0])
            h_range = tuple(host_lists[0])

            query_range_list = [q_range]
            host_range_list = [h_range]
            len_query_range_list = len(query_range_list)

            i = 0 ## position of range to compare from coords_dict
            while i < range_num:
                j = 0 ## position of range to compare from range_lists

                while j < len_query_range_list: ## check ranges added to range_lists
                    if (q_range == query_range_list[j]) and (h_range == host_range_list[j]):
                        ## exact same range, already in range_lists -- go to next in query_lists
                        if i < range_num:
                            q_range = tuple(query_lists[i])
                            h_range = tuple(host_lists[i])
                            i += 1
                        break
                    else:
                        ## if current range is not equal to previous, check range_lists for overlaps
                        any_overlap = []
                        for prev in range(len(query_range_list)):
                            tf = _is_overlaping(q_range, query_range_list[prev]) and _is_overlaping(h_range, host_range_list[prev])
                            any_overlap.append(tf)
                        if any(any_overlap):
                            ## if any overlapping ranges are found, find index from query_range_list for overlaps
                            idx = sorted([n for n, o in enumerate(any_overlap) if o == True], reverse=True)

                            ## get overlapping ranges and remove them from range_lists
                            q_ovrlp = [query_range_list.pop(n) for n in idx]
                            h_ovrlp = [host_range_list.pop(n) for n in idx]

                            ## get start and stop values from overlapping ranges
                            qmin = min(q_ovrlp)[0]
                            qmax = max(q_ovrlp, key = lambda t: t[1])[1]
                            q_ovrlp = (qmin, qmax)

                            hmin = min(h_ovrlp)[0]
                            hmax = max(h_ovrlp, key = lambda t: t[1])[1]
                            h_ovrlp = (hmin, hmax)

                            ## overwrite range variables with new start/stop coordinates and add to range_lists
                            q_range = _merge_range(q_range, q_ovrlp)
                            h_range = _merge_range(h_range, h_ovrlp)
                            query_range_list.append(q_range)
                            host_range_list.append(h_range)
                            ## update length variable of range_lists (j)
                            len_query_range_list = len(query_range_list)
                        else:
                            ## no overlaps in current range_lists, append current ranges
                            query_range_list.append(q_range)
                            host_range_list.append(h_range)
                        if i < range_num:
                            q_range = tuple(query_lists[i])
                            h_range = tuple(host_lists[i])
                            i += 1
                        j += 1
                len_query_range_list = len(query_range_list)
            paired_coords_set = set()

            for item in zip(query_range_list, host_range_list):
                paired_coords_set.add(item)
            for item in paired_coords_set:
                query_range = [item[0][0], item[0][1]]
                host_range = [item[1][0], item[1][1]]
                try:
                    merged_coords_dict[pair]['query'].append(query_range)
                    merged_coords_dict[pair]['host'].append(host_range)
                except KeyError:
                    merged_coords_dict[pair] = {'query':[], 'host':[]}
                    merged_coords_dict[pair]['query'] = [query_range]
                    merged_coords_dict[pair]['host'] = [host_range]
        for pair in merged_coords_dict.keys():
            merged_coords_dict[pair]['num_aln'] = len(merged_coords_dict[pair]['query'])
        return merged_coords_dict


class MimicDetectionII():
    __md2_vars = ('fileid', 'log_file', 'patho_name', 'host_name', 'control_name', 'k', 'bit_min', 
                  'bit_diff', 'b_str', 'min_e', 'e_str', 'qsasa', 'q_str', 'lcr', 'l_str', 
                  'runidI', 'runidII', 'outdir', 'indir')
    """
    mimicDetector amino acid solvent accessibility assesment functions: QSASA filtering
    kmer_coords_dict is a nested dictionary made with merge_ranges(filt_blast) from MimicDetectionI()
    region_avg_qsasa:   returns list of lists with merged region coordinates and mean qsasa
                            [prot_name, q_region_avg, q_region_len, q_range[0], (q_range[1] - 1),
                             h_prot_name, h_region_avg, h_region_len, h_range[0], (h_range[1] - 1)]
    """
    def __init__(self, **kwargs): 
        for k, v in kwargs.items():
            try:
                assert(k in self.__class__.__md2_vars)
                setattr(self, k, v)
            except AssertionError:
                pass

    def pops_files(self, ppops_path, hpops_path):
        # add checks for proper paths
        ppops_files = glob.glob(f"{ppops_path}/pops_*.out")
        if not ppops_files:
            err_msg = f"No pops files found for pathogen in {ppops_path}"
            # write_log(self.log_file, f'{err_msg}\n')
            raise FileNotFoundError(err_msg)
        hpops_files = glob.glob(f"{hpops_path}/pops_*.out")
        if not hpops_files:
            err_msg = f"No pops files found for host in {hpops_path}"
            # write_log(self.log_file, f'{err_msg}\n')
            raise FileNotFoundError(err_msg)
        return ppops_files, hpops_files
    
    def region_avg_qsasa(self, kmer_coords_dict, ppops_files, hpops_files):
        colnames = ['ResidNe', 'Chain', 'ResidNr', 'iCode', 'Phob/A^2', 'Phil/A^2', 'SASA/A^2', 'Q(SASA)', 'N(overl)', 'Surf/A^2']
        region_vals_lists = []
        for pfile in ppops_files:
            ## check each pathogen pops file, get protein name
            prot_name = pfile.split('pops_')[1]
            prot_name = prot_name.split('.')[0]
            prot_name = prot_name.split('-')[0]
            ## get all pairs that include pathogen protein name
            all_pair_keys = [pair_key for pair_key in kmer_coords_dict.keys() if prot_name in pair_key] 

            if len(all_pair_keys) > 0:
                ## if there is at least one pathogen-host pair for that protein, check number of pops files
                pfiles = [file for file in ppops_files if prot_name in file]
                if len(pfiles) > 1:
                    log_msg = f"Multiple pops files for pathogen protein {prot_name}\n"
                    write_log(self.log_file, log_msg)

                ## open pops file, make dataframe
                try: 
                    prot_sasa = pd.read_csv(pfile, sep="\s+", header=None, engine='python', skiprows=3, skipfooter=3)
                except pd.errors.EmptyDataError:
                    continue
                prot_sasa.columns = colnames
                for pair in all_pair_keys:
                    q_region_vals_lists = []
                    h_region_vals_lists = []
                    ## get host protein file, need to add check if exists
                    h_prot_name = pair.split('.')[1]
                    ## find all pops files
                    hfiles = [file for file in hpops_files if h_prot_name in file]
                    if len(hfiles) > 1:
                        log_msg = f"Multiple pops files for host protein {h_prot_name}\n"
                        write_log(self.log_file, log_msg)

                    ## assign variables
                    query_ranges = kmer_coords_dict[pair]['query'] 
                    host_ranges = kmer_coords_dict[pair]['host']
                    for i in range(len(query_ranges)):
                        q_range = query_ranges[i]
                        h_range = host_ranges[i]
                        q_region_vals = prot_sasa.iloc[q_range[0]:q_range[1]]['Q(SASA)'].to_list()
                        q_region_len = (q_range[1] - q_range[0])
                        if not q_region_vals:
                            log_msg = f"Pathogen protein {prot_name} region {q_range[0]}:{q_range[1]} out of pops file range\n"
                            write_log(self.log_file, log_msg)
                            p += 1
                            q_region_vals = 0
                            q_region_avg = "NaN"
                        else:
                            q_region_len = len(q_region_vals) 
                            q_region_avg = (sum(q_region_vals) / q_region_len)
                        q_region_vals_lists.append(q_region_vals)
                        h_region_len = (h_range[1] - h_range[0])
                        if hfiles:
                            hfile = f"{self.outdir}/{self.host_name}/pops/pops_{h_prot_name}.out"
                            try: 
                                h_prot_sasa = pd.read_csv(f'{hfile}', sep="\s+", header=None, engine='python', skiprows=3, skipfooter=3)
                            except pd.errors.EmptyDataError:
                                h_region_vals = 0
                                h_region_vals_lists.append(h_region_vals)
                                h_region_avg = "NaN"
                                continue
                            h_prot_sasa.columns = colnames
                            h_region_vals = h_prot_sasa.iloc[h_range[0]:h_range[1]]['Q(SASA)'].to_list()
                            if not h_region_vals:
                                log_msg = f"Host protein {h_prot_name} region {h_range[0]}:{h_range[1]} out of pops file range\n"
                                write_log(self.log_file, log_msg)
                                h += 1
                                h_region_vals = 0
                                h_region_avg = "NaN"
                            else:
                                h_region_len = len(h_region_vals)
                                h_region_avg = (sum(h_region_vals) / h_region_len)                    
                            h_region_vals_lists.append(h_region_vals)
                        else: 
                            h_region_vals = 0
                            h_region_vals_lists.append(h_region_vals)
                            h_region_avg = "NaN"
                        region_vals_lists.append([prot_name, q_range[0], (q_range[1] - 1), q_region_avg, q_region_len, 
                            h_prot_name, h_range[0], (h_range[1] - 1), h_region_avg, h_region_len])
        all_qsasa_df = pd.DataFrame(region_vals_lists, columns = ['q_prot', 'q_start', 'q_end', 
                                                                  'q_avg', 'q_length', 'h_prot', 
                                                                  'h_start', 'h_end', 'h_avg', 
                                                                  'h_length'])
        return all_qsasa_df 
    
    def qsasa_filt(self, all_qsasa_df, write=True):
        qsasa_df_filtered = all_qsasa_df[(all_qsasa_df['h_avg'] != 'NaN') & (all_qsasa_df['q_avg'] != 'NaN')]
        qsasa_df_filtered = qsasa_df_filtered[(qsasa_df_filtered['q_avg'] > self.qsasa) & (qsasa_df_filtered['h_avg'] > self.qsasa)] 
        if write:
            all_qsasa_file=Path(self.outdir, self.patho_name, f"{self.patho_name}.{self.runidI}_qsasa_averages.tsv")
            all_qsasa_df.to_csv(all_qsasa_file, sep="\t", header=False, index=False)
            filt_qsasa_file=Path(self.outdir, self.patho_name, f"{self.patho_name}.{self.runidII}_filtered.tsv")
            qsasa_df_filtered.to_csv(filt_qsasa_file, sep="\t", header=False, index=False) 

            log_msg = ''.join([all_qsasa_df.shape[0], 
                               " paired regions before QSASA filtering",
                               " saved to: \n\t", 
                               str(all_qsasa_file),  
                               qsasa_df_filtered.shape[0],
                               " paired regions with minimum mean QSASA of 0.", 
                               self.q_str, 
                               " saved to: \n\t", 
                               str(filt_qsasa_file), 
                               "\n"])
            # write_log(self.log_file, log_msg)
            return all_qsasa_file, filt_qsasa_file
        else:
            log_msg = ''.join([all_qsasa_df.shape[0], 
                               " paired regions before QSASA filtering\n", 
                               qsasa_df_filtered.shape[0],
                               " paired regions with minimum mean QSASA of 0.", 
                               self.q_str, 
                               "\n"])
            # write_log(self.log_file, log_msg)
            return qsasa_df_filtered
        
    
class GreaterMimicDetection():
    __gmd_vars = ('fileid', 'log_file', 'patho_name', 'host_name', 'control_name', 'k', 'bit_min', 
                  'bit_diff', 'b_str', 'min_e', 'e_str', 'qsasa', 'q_str', 'lcr', 'l_str', 'dbsize',
                  'runidI', 'runidII', 'runidIII', 'outdir', 'indir')
    """
    mimicDetector final mimicry candidate filtering: LCR filtering 

    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            try:
                assert(k in self.__class__.__gmd_vars)
                setattr(self, k, v)
            except AssertionError:
                pass

    def lcr_filter(self, qsasa_filtered, lcr_coords, host_fasta, pathogen_fasta, outfile):
        qfilt_df = pd.read_csv(qsasa_filtered, sep='\t', skip_blank_lines=True,
                               header=None, names=['pathogen_prot', 'pathogen_start',
                                                   'pathogen_end', 'pathogen_qsasa', 
                                                   'pathogen_length', 'host_prot', 
                                                   'host_start', 'host_end', 
                                                   'host_qsasa', 'host_length'])

        seg_bed = BedTool(lcr_coords)
        seg_df = seg_bed.to_dataframe()
        seg_df.columns = ['host_prot', 'host_start', 'host_end']

        qhfilt_bed = BedTool.from_dataframe(qfilt_df[['host_prot', 'host_start', 'host_end']])
        hostlcr_bed = qhfilt_bed.intersect(seg_bed, wo=True, f=self.lcr, v=True)
        hostlcr_df = hostlcr_bed.to_dataframe()
        hostlcr_df.columns = ['host_prot', 'host_start', 'host_end']

        paired_df = hostlcr_df.merge(qfilt_df, on=['host_prot', 'host_start', 'host_end'], how='right') 

        pbed = BedTool.from_dataframe(paired_df[['pathogen_prot', 'pathogen_start', 'pathogen_end']])
        hbed = BedTool.from_dataframe(paired_df[['host_prot', 'host_start', 'host_end']])
                
        hseqbed=f"{self.outdir}/{self.patho_name}/{self.patho_name}.{self.runidIII}_seqs_host"
        pseqbed=f"{self.outdir}/{self.patho_name}/{self.patho_name}.{self.runidIII}_seqs_pathogen"

        pseq_bed = pbed.sequence(fi=BedTool(pathogen_fasta), tab=True, fo=pseqbed) 
        pseq_df= pd.read_csv(pseqbed, sep='\t', header=None, names=['chr', 'pathogen_sequence'])
        pseq_df[['pathogen_prot','coords']] = pseq_df['chr'].str.split(':', expand=True)
        pseq_df[['pathogen_start', 'pathogen_end']] = pseq_df['coords'].str.split('-', expand=True)
        pseq_df = pseq_df[['pathogen_prot', 'pathogen_start', 'pathogen_end', 'pathogen_sequence']]
        pseq_df[['pathogen_start', 'pathogen_end']] = pseq_df[['pathogen_start', 'pathogen_end']].astype(int)

        hseq_bed = hbed.sequence(fi=BedTool(host_fasta), tab=True, fo=hseqbed)
        hseq_df= pd.read_csv(hseqbed, sep='\t', header=None, names=['chr', 'host_sequence'])
        hseq_df[['host_prot','coords']] = hseq_df['chr'].str.split(':', expand=True)
        hseq_df[['host_start', 'host_end']] = hseq_df['coords'].str.split('-', expand=True)
        hseq_df = hseq_df[['host_prot', 'host_start', 'host_end', 'host_sequence']]
        hseq_df[['host_start', 'host_end']] = hseq_df[['host_start', 'host_end']].astype(int)
        paired_df = paired_df.infer_objects()

        hpaired_df = hseq_df.merge(paired_df, how='outer', on=['host_prot', 'host_start', 'host_end'])
        filt_paired_df = pseq_df.merge(hpaired_df, how='inner', on=['pathogen_prot', 'pathogen_start', 
                                                                    'pathogen_end']).round({'pathogen_qsasa': 3, 
                                                                                            'host_qsasa': 3})
        filt_paired_df = filt_paired_df[['pathogen_prot', 'pathogen_start', 'pathogen_end', 
                                         'pathogen_sequence', 'pathogen_qsasa', 
                                         'host_prot', 'host_start', 'host_end', 
                                         'host_sequence', 'host_qsasa']]
        
        def _aln(seq1, seq2, matrix="PAM30", lamb=0.339, kconst=0.28, gopen=-15, gext=-3):  
            """LAMBDA and K values retrieved from BLASTP v2.16.0 stats table for PAM30"""
            aligner = Align.PairwiseAligner()
            aligner.substitution_matrix = substitution_matrices.load(matrix)
            aligner.open_gap_score = gopen
            aligner.extend_gap_score = gext
            bits = aligner.score(seq1, seq2)
            bitsc =  ((lamb*bits) - np.log(kconst))/ np.log(2)
            evalue = kconst*len(seq1)*self.dbsize*(2.71828**(-lamb*bits))
            return bitsc, evalue
        
        filt_paired_df[['bitscore', 'e_value']] = filt_paired_df.apply(lambda x: _aln(x.pathogen_sequence, 
                                                                                      x.host_sequence), 
                                                                                      axis = 1, 
                                                                                      result_type='expand'
                                                                                      ).round({'bitscore':1, 'e_value':5})

        # log_msg = ''.join([filt_paired_df.shape[0], 
        #                    " paired regions with maximum LCR overlap of 0.", 
        #                    self.l_str, 
        #                    " saved to: \n\t", str(outfile), 
        #                    "\n"])
        # write_log(self.log_file, log_msg)
        filt_paired_df.to_csv(Path(outfile), sep='\t', header=True, index=False)
        return 


class Mimesis():
    __mim_vars = ()
    """
    mimicDetector functions for metagenomic commensal analysis 
    Takes dictionary containing all species to be compared: symb_dict['symbiotic_species_name'] = 'control_name'
    Host is assumed to be the same for all 
    """
    def __init__(self, symb_dict={}, **kwargs):
        err_str = "Mimesis: mimicry analysis in metagenomic commensal data has not been implemented yet, sorry!"
        raise NotImplementedError(err_str)
        for k, v in kwargs.items():
            try:
                assert(k in self.__class__.__mim_vars)
                setattr(self, k, v)
            except AssertionError:
                pass
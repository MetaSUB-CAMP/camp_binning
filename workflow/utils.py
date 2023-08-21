'''Utilities.'''


# --- Workflow setup --- #


import glob
import gzip
import os
from os import makedirs, symlink
from os.path import abspath, basename, exists, join
import pandas as pd
import shutil


def extract_from_gzip(ap, out):
    if open(ap, 'rb').read(2) == b'\x1f\x8b': # If the input is gzipped
        with gzip.open(ap, 'rb') as f_in, open(out, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else: # Otherwise, symlink
        symlink(ap, out)


def scrub_fasta_tags(fi, fo):
    with open(fi, 'r') as f_in, open(fo, 'w') as f_out:
        for l in f_in:
            if l.startswith('>'):
                f_out.write(l.replace(' ', '_')) # VAMB doesn't allow spaces
            else:
                f_out.write(l)


def ingest_samples(samples, tmp):
    df = pd.read_csv(samples, header = 0, index_col = 0) # name, ctgs, fwd, rev
    s = list(df.index)
    lst = [str(l) for l in df.values.tolist()]
    for i,l in enumerate(lst):
        if not exists(join(tmp, s[i] + '.fasta')):
            scrub_fasta_tags(l[0], join(tmp, s[i] + '.fasta'))
            # symlink(abspath(l[0]), join(tmp, s[i] + '.fasta'))
        if not exists(join(tmp, s[i] + '_1.fastq.gz')):
            extract_from_gzip(abspath(l[1]), join(tmp, s[i] + '_1.fastq'))
            extract_from_gzip(abspath(l[2]), join(tmp, s[i] + '_2.fastq'))
    return s


def check_make(d):
    if not exists(d):
        makedirs(d)


class Workflow_Dirs:
    '''Management of the working directory tree.'''
    OUT = ''
    TMP = ''
    LOG = ''

    def __init__(self, work_dir, module):
        self.OUT = join(work_dir, module)
        self.TMP = join(work_dir, 'tmp')
        self.LOG = join(work_dir, 'logs')
        check_make(self.OUT)
        out_dirs = ['0_contig_coverage', '1_metabat2', '2_concoct', '3_semibin', '4_maxbin2', '5_dastool', 'final_reports']
        for d in out_dirs:
            check_make(join(self.OUT, d))
        # Add a subdirectory for symlinked-in input files
        check_make(self.TMP)
        # Add custom subdirectories to organize rule logs
        check_make(self.LOG)
        log_dirs = ['map_sort', 'calculate_depth', 'metabat2_binning', 'concoct_binning', 'semibin_binning', 'maxbin2_binning', 'dastool_refinement']
        for d in log_dirs:
            check_make(join(self.LOG, d))


def cleanup_files(work_dir, df):
    smps = list(df.index)
    for s in smps:
        for f in glob.glob(join(work_dir, 'binning', '0_contig_coverage', s, 'coverage.*')):
            os.remove(abspath(f))
        for f in glob.glob(join(work_dir, 'binning', '0_contig_coverage', s, '*.bt*')):
            os.remove(abspath(f))
        for f in glob.glob(join(work_dir, 'binning', '2_concoct', s, '*.fasta')):
            os.remove(abspath(f))
        for f in glob.glob(join(work_dir, 'binning', '2_concoct', s, 'original_data_gt*')):
            os.remove(abspath(f))
        shutil.rmtree(join(work_dir, 'binning', '3_semibin', s, 'output_prerecluster_bins'))


def print_cmds(f):
    # fo = basename(log).split('.')[0] + '.cmds'
    # lines = open(log, 'r').read().split('\n')
    fi = [l for l in f.split('\n') if l != '']
    write = False
    with open('commands.sh', 'w') as f_out:
        for l in fi:
            if 'rule' in l:
                f_out.write('# ' + l.strip().replace('rule ', '').replace(':', '') + '\n')
                write = False
            if 'wildcards' in l:
                f_out.write('# ' + l.strip().replace('wildcards: ', '') + '\n')
            if 'resources' in l:
                write = True
                l = ''
            if write:
                f_out.write(l.strip() + '\n')
            if 'rule make_config' in l:
                break


# --- Workflow functions --- #


from Bio import SeqIO
from collections import Counter
from io import StringIO
from os.path import basename, splitext
import re
import subprocess
import sys


def chunks(l, n, o):
    # Yield successive n-sized chunks from l with given overlap o between the chunks.
    assert n > o
    for i in range(0, len(l) - n + 1, n - o):
        yield l[i:i + n] if i + n + n - o <= len(l) else l[i:]


def cut_up_fasta(ctg, frag_size, olap_size, out_dir, out_fa, out_bed):
    '''Split up assembly contigs into windows of 'frag_size' bp 
    with 'olap_size' overlap.'''
    if not exists(out_dir):
        makedirs(out_dir)
    with open(out_fa, 'w') as of, open(out_bed, 'w') as ob:
        for r in SeqIO.parse(ctg, "fasta"):
            if len(r.seq) >= 2 * frag_size:
                i = 0
                for split_seq in chunks(r.seq, frag_size, olap_size):
                    print(">{}.concoct_part_{}\n{}".\
                        format(r.id, i, split_seq), file = of)
                    print("{0}\t{2}\t{3}\t{0}.concoct_part_{1}".\
                        format(r.id, i, frag_size*i, frag_size*i+len(split_seq)),\
                        file = ob)
                    i += 1
            else:
                print(">{}.concoct_part_0\n{}".\
                        format(r.id, r.seq), file = of)
                print("{0}\t0\t{1}\t{0}.concoct_part_0".format(r.id, len(r.seq)),\
                    file = ob)


def make_concoct_table(in_bed, in_bam, output):
    '''Reads input files into dictionaries then prints everything in the table
    format required for running CONCOCT.'''
    p = subprocess.Popen(["samtools", "bedcov", in_bed, in_bam], \
        stdout = subprocess.PIPE)
    out, err = p.communicate()
    # Header
    col_names = [splitext(basename(in_bam))[0]] # Use index if no sample names given in header
    header = ["cov_mean_sample_{}".format(n) for n in col_names]
    # Content
    fh = StringIO(out.decode('utf-8'))
    df = pd.read_table(fh, header=None)
    avg_coverage_depth = df[df.columns[4:]].divide((df[2]-df[1]), axis=0)
    avg_coverage_depth.index = df[3]
    avg_coverage_depth.columns = header
    avg_coverage_depth.to_csv(output, index_label='contig', sep='\t', \
        float_format='%.3f')


def extract_bin_id(contig_id):
    CONTIG_PART_EXPR = re.compile("(.*)\.concoct_part_([0-9]*)")
    n = contig_id.split('.')[-1]
    try:
        original_contig_id, part_id = CONTIG_PART_EXPR.match(contig_id).group(1,2)
        return [original_contig_id, part_id]
    except AttributeError: # No matches for concoct_part regex
        return contig_id, 0


def split_concoct_output(concoct, ctg, out_dir):
    # Match CONCOCT bin labels to the original contig IDs
    df = pd.read_csv(concoct, header=0)  # contig_id,cluster_id
    new_cols = df.apply(lambda row: extract_bin_id(row['contig_id']), axis=1)
    df['original_contig_id'] = [i[0] for i in new_cols]
    df['part_id'] = [i[1] for i in new_cols]
    # Find the best bin label for each original contig
    cluster_mapping = {}
    cluster_fastas = {}
    for curr_ctg in df.original_contig_id.unique():  # For each original contig...
        sub_df = df[df['original_contig_id'] == curr_ctg]
        if sub_df.shape[1] > 1:  # If there are multiple assignments
            c = Counter(list(sub_df['cluster_id']))
            majority_vote = c.most_common(1)[0][0]
            possible_bins = [(a, b) for a, b in c.items()]
            if len(c.values()) > 1:
                sys.stderr.write('No consensus cluster for \
                    contig {}: {}\t CONCOCT cluster: {}\n' \
                                 .format(curr_ctg, possible_bins, majority_vote))
        else:
            majority_vote = list(sub_df['cluster_id'])[0]
        cluster_mapping[curr_ctg] = majority_vote
        cluster_fastas[majority_vote] = []
    # Split the assembly FastA into separate bin FastAs
    cluster_fastas['unbinned'] = []
    curr_bin = ''
    curr_contig = ''
    tmp_lines = []
    for line in open(ctg, 'r'):
        if line.startswith('>'):
            if curr_contig != '':
                cluster_fastas[curr_bin].extend(tmp_lines)
                tmp_lines = []
            curr_contig = line[1:-1].split('.')[0].split()[0]
            curr_bin = cluster_mapping[curr_contig] if curr_contig in cluster_mapping else 'unbinned'
            # print('{} {}'.format(curr_contig, curr_bin))
        tmp_lines.append(line.strip())
    # Write each separate bin FastA
    for k, v in cluster_fastas.items():
        with open('{}/bin.{}.fa'.format(out_dir, k), 'w') as f_out:
            for line in v:
                f_out.write(line + '\n')


def get_dastool_unbinned(ctg, tsv, out_dir): # No header, just ctg\tbin
    df = pd.read_csv(tsv, sep = '\t', header = None)
    unbinned_df = df[df[1] == 'unbinned']
    unbinned_ctgs = list(unbinned_df[0])
    with open(ctg, 'r') as fi, open(join(out_dir, 'bin.unbinned.fa'), 'w') as fo:
        print_ctg = False
        for l in fi:
            if '>' in l: 
                print_ctg = False
                for c in unbinned_ctgs:
                    if c in l:
                        print_ctg = True
                        break
            if print_ctg:
                fo.write(l)

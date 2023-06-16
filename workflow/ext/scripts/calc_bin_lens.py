import argparse
import glob
from os.path import abspath, basename, join
import pandas as pd
import numpy as np


def main(args): # sample_name,binners,bin_num,num_ctgs,total_size,mean_bin_size,stdev_bin_size
    parts = args.in_dir.split('/')
    sample = parts[-2]
    binner = parts[-3]
    bin_lst = []
    for fi in glob.glob(join(args.in_dir, '*.fa')):
        bin_num = basename(fi).split('.')[1]
        with open(abspath(fi), 'r') as f:
            ctg_lens = []
            curr_len = 0
            for l in f:
                if l.startswith('>'):
                    ctg_lens.append(curr_len)
                    curr_len = 0
                else:
                    curr_len += len(l.strip())
            ctg_lens.append(curr_len)
            ctg_lens.pop(0)  # Remove the first 0
            bin_lst.append([sample, binner, bin_num, len(ctg_lens), sum(ctg_lens), sum(ctg_lens) / len(ctg_lens), np.std(ctg_lens)])
    pd.DataFrame(bin_lst).to_csv(args.output, header = False, index = False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_dir", help="MAG directory")
    parser.add_argument("output", help="Summary statistics of MAGs inferred from a sample by a binner")
    args = parser.parse_args()
    main(args)

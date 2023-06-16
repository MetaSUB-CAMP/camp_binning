from vamb.vambtools import Reader, read_clusters, loadfasta, write_bins
from os import listdir
from os.path import join
import shutil
import sys

# sys.argv[1]: VAMB-generated cluster TSV
# sys.argv[2]: De novo assembled contig FastA
# sys.argv[3]: Input directory
# sys.argv[4]: Output directory

bins = read_clusters(open(sys.argv[1], 'r'))
fasta = loadfasta(Reader(sys.argv[2], 'rb'))
write_bins(sys.argv[3], bins, fasta)
for f in listdir(sys.argv[3]):
    if f.endswith('.fna'): # 1.fna, 2.fna, etc.
        bin_num = f.split('.')[0]
        shutil.copy(join(sys.argv[3], f), join(sys.argv[4], 'bin.' + bin_num + '.fa'))
open(sys.argv[3] + '_done.txt', 'w').write('')

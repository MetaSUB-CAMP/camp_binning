from vamb.vambtools import Reader, read_clusters, loadfasta, write_bins
from os import listdir
from os.path import join
import shutil
import sys

# sys.argv[1]: VAMB-generated cluster TSV
# sys.argv[2]: De novo assembled contig FastA
# sys.argv[3]: Input directory
# sys.argv[4]: Output directory

bins = read_clusters(open(sys.argv[1], 'r'), min_size = 2)
fasta = loadfasta(Reader(sys.argv[2], 'rb'))
write_bins(sys.argv[4], bins, fasta, maxbins=10000, minsize=int(sys.argv[3]))
for f in listdir(sys.argv[4]):
    if f.endswith('.fna'): # 1.fna, 2.fna, etc.
        bin_num = f.split('.')[0]
        shutil.copy(join(sys.argv[4], f), join(sys.argv[5], 'bin.' + bin_num + '.fa'))
open(sys.argv[4] + '_done.txt', 'w').write('')

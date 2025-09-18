from .dbg import DBG
from utils import read_data, n50
import sys
import os

sys.setrecursionlimit(1000000)


if __name__ == "__main__":
    argv = sys.argv
    short1, short2, long1 = read_data(os.path.join('./', argv[1]))

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    contigs = []
    # dbg.show_count_distribution()
    with open(os.path.join('./', argv[1], 'contig.fasta'), 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            # print(i, len(c))
            contigs.append(c)
            f.write('>contig_%d\n' % i)
            f.write(c + '\n')

    contig_lengths = [len(c) for c in contigs]
    print(f"N50: {n50(contig_lengths)}")

from dbg_codon import DBG
from utils_codon import read_data, n50
import sys

if __name__ == "__main__":
    argv = sys.argv
    short1, short2, long1 = read_data( argv[1])

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    contigs = []
    # dbg.show_count_distribution()
    with open(argv[1] + '/contig.fasta', 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            contigs.append(c)
            f.write(f">contig_{i}\n")
            f.write(c + '\n')

    contig_lengths = [len(c) for c in contigs]
    print(f"N50: {n50(contig_lengths)}")
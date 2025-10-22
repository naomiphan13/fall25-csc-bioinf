import Alignment_codon as Alignment
from util_codon import read_fasta, data_dir, join_path
import time
import re

def load_map(path):
    return {h.split()[0]: s for h, s in read_fasta(path)}

def pairs():
    q = load_map(join_path(data_dir(), "q1.fa"))
    t = load_map(join_path(data_dir(), "t1.fa"))
    human = load_map(join_path(data_dir(), "MT-human.fa"))
    orang = load_map(join_path(data_dir(), "MT-orang.fa"))
    for n in range(1, 6):
        qid = f"q{n}"
        tid = f"t{n}"
        yield qid, tid, q[qid], t[tid]

    hid, hseq = next(iter(human.items()))
    oid, oseq = next(iter(orang.items()))
    yield hid, oid, hseq, oseq

def norm_id(s: str) -> str:
    # lower-case and turn non-alnum into underscores for tidy tags
    s = s.lower()
    return re.sub(r'[^0-9a-z]+', '_', s).strip('_')

def time_call(tag: str, fn):
    t0 = time.time()
    fn()
    ms = (time.time() - t0) * 1000
    print(f"{tag}\t\tcodon\t{ms}ms")

for seq1_id, seq2_id, seq1, seq2 in pairs():
    s1 = seq1.casefold()
    s2 = seq2.casefold()
    tag_base = f"{norm_id(seq1_id)}_{norm_id(seq2_id)}"

    aln = Alignment.Alignment(s1, s2)
    time_call(f"global-{tag_base}",  lambda: aln.global_alignment())
    time_call(f"local-{tag_base}",   lambda: aln.local_alignment())
    time_call(f"fitting-{tag_base}", lambda: aln.fitting_alignment())
    time_call(f"affine-{tag_base}",  lambda: aln.affine_alignment())




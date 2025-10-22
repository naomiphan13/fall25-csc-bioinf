import os

def data_dir():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    week4_root = os.path.abspath(os.path.join(current_dir, "../"))
    data_dir = os.path.join(week4_root, "data")
    return data_dir

def read_fasta(path):
    with open(path, "r") as f:
        header = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]                  # full header (no '>')
                seq_chunks = []
            else:
                seq_chunks.append(line)           # lines may be wrapped
        if header is not None:
            yield header, "".join(seq_chunks)

    return header, seq_chunks



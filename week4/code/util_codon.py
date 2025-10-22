from typing import List, Optional

def data_dir() -> str:
    return "../data"

def join_path(dir_path: str, filename: str) -> str:
    return dir_path + "/" + filename

def read_fasta(path: str):
    header: Optional[str] = None
    seq_chunks: List[str] = []

    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
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



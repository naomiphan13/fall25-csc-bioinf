def read_fasta(path, name):
    data = []
    with open('/'.join([path, name]), 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] != '>':
                data.append(line)
    # print(name, len(data), len(data[0]))
    # print('Sample:', data[0])
    return data


def read_data(path):
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1, short2, long1

def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    lengths = sorted([L for L in lengths if L > 0], reverse=True)
    total = sum(lengths)
    half = total / 2
    cumm = 0
    for L in lengths:
        cumm += L
        if cumm >= half:
            return L
    return 0

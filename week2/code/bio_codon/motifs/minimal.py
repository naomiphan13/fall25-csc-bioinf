# Copyright 2018 by Ariel Aptekmann.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Module for the support of MEME minimal motif format."""

from .__init__ import create_motif_from_counts
from .__init__ import Motif
from typing import Optional, Dict, List

def read(handle):
    motif_number = 0
    record = Record()
    _read_version(record, handle)
    _read_alphabet(record, handle)
    _read_background(record, handle)

    while True:
        for line in handle:
            if line.startswith("MOTIF"):
                break
        else:
            return record
        name = line.split()[1]
        motif_number += 1

        length, num_occurrences, evalue = _read_motif_statistics(handle)

        counts = _read_lpm(record, handle, length, num_occurrences)

        # {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25}
        motif = create_motif_from_counts(alphabet=record.alphabet, counts=counts)
        motif.background = record.background
        motif.length = motif.counts.length
        motif.num_occurrences = num_occurrences
        motif.evalue = evalue
        motif.name = name
        record.append(motif)
        assert len(record) == motif_number
    return record


class Record():
    """Class for holding the results of a minimal MEME run."""

    version: str
    datafile: str
    command: str
    alphabet: Optional[str]
    background: Dict[str, float]
    sequences: List[str]
    data: List[Motif]

    def __init__(self):
        """Initialize record class values."""
        self.version = ""
        self.datafile = ""
        self.command = ""
        self.alphabet = None
        self.background = {}
        self.sequences = []

    def append(self, item: Motif):
        self.data.append(item)

    def __len__(self):
        return len(self.data)
    
    def __iter__(self):
        return self.data.__iter__()

    def __getitem__(self, key: str):
        """Return the motif of index key."""
        if isinstance(key, str):
            for motif in self:
                if motif.name == key:
                    return motif
        else:
            return self.data.__getitem__(self, key)


# Everything below is private


def _read_background(record, handle):
    """Read background letter frequencies (PRIVATE)."""
    for line in handle:
        if line.startswith("Background letter frequencies"):
            background_freqs = []
            for line in handle:
                line = line.rstrip()
                if line:
                    background_freqs.extend(
                        [
                            float(freq)
                            for i, freq in enumerate(line.split(" "))
                            if i % 2 == 1
                        ]
                    )
                else:
                    break
            if not background_freqs:
                raise Exception(
                    "Unexpected end of stream: Expected to find line starting background frequencies."
                )
            break
    else:
        raise Exception(
            "Improper input file. File should contain a line starting background frequencies."
        )
    record.background = dict(zip(record.alphabet, background_freqs))


def _read_version(record, handle):
    """Read MEME version (PRIVATE)."""
    for line in handle:
        if line.startswith("MEME version"):
            break
    else:
        raise Exception(
            "Improper input file. File should contain a line starting MEME version."
        )
    line = line.strip()
    ls = line.split()
    record.version = ls[2]


def _read_alphabet(record, handle):
    """Read alphabet (PRIVATE)."""
    for line in handle:
        if line.startswith("ALPHABET"):
            break
    else:
        raise Exception(
            "Unexpected end of stream: Expected to find line starting with 'ALPHABET'"
        )
    if not line.startswith("ALPHABET= "):
        raise Exception(f"Line does not start with 'ALPHABET': {line}")
    line = line.strip().replace("ALPHABET= ", "")
    if line == "ACGT":
        al = "ACGT"
    elif line == "ACGU":
        al = "ACGU"
    else:
        # al = "ACDEFGHIKLMNPQRSTVWY"
        raise Exception("Only parsing of DNA and RNA motifs is implemented")
    record.alphabet = al


def _read_lpm(record, handle, length, num_occurrences):
    """Read letter probability matrix (PRIVATE)."""
    counts = [[], [], [], []]
    for line in handle:
        freqs = line.split()
        if len(freqs) != 4:
            break
        counts[0].append(round(float(freqs[0]) * num_occurrences))
        counts[1].append(round(float(freqs[1]) * num_occurrences))
        counts[2].append(round(float(freqs[2]) * num_occurrences))
        counts[3].append(round(float(freqs[3]) * num_occurrences))
        if length and len(counts[0]) == length:
            break
    c = dict(zip(record.alphabet, counts))
    return c


def _read_motif_statistics(handle):
    """Read motif statistics (PRIVATE)."""
    for line in handle:
        if line.startswith("letter-probability matrix:"):
            break

    # The "nsites= source sites" will default to 20 if it is not provided.
    num_occurrences = (
        int(line.split("nsites=")[1].split()[0]) if line.find("nsites=") != -1 else 20
    )
    # Length can be infered later if it is not provided.
    length = int(line.split("w=")[1].split()[0]) if line.find("w=") != -1 else None
    # E-value will default to zero if it is not provided.
    evalue = float(line.split("E=")[1].split()[0]) if line.find("E=") != -1 else 0.0
    return length, num_occurrences, evalue


def _read_motif_name(handle):
    """Read motif name (PRIVATE)."""
    for line in handle:
        if "sorted by position p-value" in line:
            break
    else:
        raise Exception("Unexpected end of stream: Failed to find motif name")
    line = line.strip()
    words = line.split()
    name = " ".join(words[0:2])
    return name

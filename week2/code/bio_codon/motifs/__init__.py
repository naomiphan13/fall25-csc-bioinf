from python import importlib
from python import Bio
from typing import Dict, List, Optional

Alignment = importlib.import_module("Bio.Align").Alignment

try:
    from python import numpy as np
except Exception:
    raise Exception(
        "Install NumPy if you want to use Bio.motifs."
    )

import matrix

class Motif:
    name: str
    counts: Optional[matrix.FrequencyPositionMatrix]
    length: int
    alignment: Optional[pyobj]
    alphabet: str
    _pseudocounts: Dict[str, float]
    _background: Optional[Dict[str, float]]
    _mask: Optional[List[int]]
    num_occurrences: int
    evalue: float

    """A class representing sequence motifs."""
    def __init__(self, alphabet="ACGT", alignment: Optional[pyobj]=None, counts: Optional[matrix.FrequencyPositionMatrix]=None):
        """Initialize the class."""
        self.name = ""
        self.num_occurrences = 0
        self.evalue = 0.0

        if counts is not None and alignment is not None:
            raise Exception("Specify either counts or an alignment, don't specify both")

        elif counts is not None:
            self.alignment = None
            self.counts = counts
            self.length = self.counts.length

        elif alignment is not None:
            length = alignment.length
            frequencies: Dict[str, List[float]] = {}
            for letter in alphabet:
                if letter not in list(alignment.frequencies.keys()):
                    frequencies[letter] = [0.0] * int(length)
                else:
                    frequencies[letter] = [float(alignment.frequencies[letter][i]) for i in range(length)]
            self.counts = matrix.FrequencyPositionMatrix(alphabet, frequencies)
            self.alignment = alignment
            self.length = length
            
        else:
            self.counts = None
            self.alignment = None
            self.length = 0
        self.alphabet = alphabet
        self._pseudocounts = dict.fromkeys(self.alphabet, 0.0)
        self._background: Optional[Dict[str, float]] = None
        self._mask = None
        
    
    @property
    def mask(self):
        return self.__mask
    
    @mask.setter
    def mask(self, mask):     
        if self.length is None or self.length == 0:
            self.__mask = []
        elif mask is None:
            self.__mask = (1,) * self.length
        elif len(mask) != self.length:
            raise Exception(
                f"The length {len(mask)} of the mask is inconsistent with the length {self.length} of the motif"
            )
        elif isinstance(mask, str):
            _mask = []
            for char in mask:
                if char == "*":
                    _mask.append(1)
                elif char == " ":
                    _mask.append(0)
                else:
                    raise Exception(
                        f"Mask should contain only '*' or ' ' and not a '{char}'"
                    )
            self.__mask = _mask
        else:
            self.__mask = [int(bool(c)) for c in mask]

    @property
    def pseudocounts(self):
        return self._pseudocounts
    
    @pseudocounts.setter
    def pseudocounts(self, value):
        self._pseudocounts = {}
        if isinstance(value, Dict[str, float]):
            self._pseudocounts = {letter: float(value.get(letter, 0.0)) for letter in self.alphabet}
        else:
            if value is None:
                value = 0.0
            self._pseudocounts = dict.fromkeys(self.alphabet, value)

    @property
    def background(self):
        return self._background
    
    @background.setter
    def background(self, value = None):
        if isinstance(value, Dict[str, float]):
            self._background: dict[str, float] = {letter: value[letter] for letter in self.alphabet}
        elif value is None:
            self._background = dict.fromkeys(self.alphabet, 1.0)
        else:
            if not self._has_dna_alphabet() and not self._has_rna_alphabet():
                raise Exception(
                    "Setting the background to a single value only works for DNA and RNA"
                    "motifs (in which case the value is interpreted as the GC content)"
                )
            value = float(value)
            T_or_U = "T" if self._has_dna_alphabet() else "U"
            self._background["A"] = (1.0 - value) / 2.0
            self._background["C"] = value / 2.0
            self._background["G"] = value / 2.0
            self._background[T_or_U] = (1.0 - value) / 2.0
        total = sum(self._background.values())
        for letter in self.alphabet:
            self._background[letter] /= total

    def __getitem__(self, key: slice):
        alphabet = self.alphabet
        if self.alignment is None:
            alignment = None
            if self.counts is None:
                counts = None
            else:
                counts: Optional[Dict[str, float]] = {letter: float(self.counts[letter][key]) for letter in alphabet}
        else:
            alignment = self.alignment[:, key]
            counts = None
        motif = Motif(alphabet=alphabet, alignment=alignment, counts=counts)
        motif.mask = self.mask[key]
        if alignment is None and counts is None:
            try:
                length = self.length
            except Exception:
                pass
            else:
                motif.length = len(range(*key.indices(length)))
        motif.pseudocounts = self.pseudocounts.copy()
        motif.background = self.background.copy()
        return motif

    @property
    def pwm(self):
        """Calculate and return the position weight matrix for this motif."""
        return self.counts.normalize(self._pseudocounts)

    @property
    def pssm(self):
        """Calculate and return the position specific scoring matrix for this motif."""
        return self.pwm.log_odds(self._background)

    def __str__(self, masked=False):
        """Return string representation of a motif."""
        text = ""
        if self.alignment is not None:
            text += "\n".join(self.alignment)

        if masked:
            for i in range(self.length):
                if self.__mask[i]:
                    text += "*"
                else:
                    text += " "
            text += "\n"
        return text

    def __len__(self):
        """Return the length of a motif.

        Please use this method (i.e. invoke len(m)) instead of referring to m.length directly.
        """
        return 0 if self.length is None else self.length

    def _has_dna_alphabet(self):
        return sorted(self.alphabet) == ["A", "C", "G", "T"]

    def _has_rna_alphabet(self):
        return sorted(self.alphabet) == ["A", "C", "G", "U"]

    def reverse_complement(self):
        """Return the reverse complement of the motif as a new motif."""
        alphabet = self.alphabet
        if not self._has_dna_alphabet() and not self._has_rna_alphabet():
            raise Exception(
                "Calculating reverse complement only works for DNA and RNA motifs"
            )
        T_or_U = "T" if self._has_dna_alphabet() else "U"
        if self.alignment is not None:
            alignment = self.alignment.reverse_complement()
            sequences = alignment.sequences
            if T_or_U == "U":
                seqs = []
                for s in sequences: seqs.append(s.replace("T", "U"))
                res = Motif(alphabet=alphabet, alignment=Alignment(sequences=seqs))
            else:
                res = Motif(alphabet=alphabet, alignment=alignment)
        else:  # has counts
            counts = {
                "A": self.counts[T_or_U][::-1],
                "C": self.counts["G"][::-1],
                "G": self.counts["C"][::-1],
                T_or_U: self.counts["A"][::-1],
            }
            res = Motif(alphabet=alphabet, counts=counts)
        res.__mask = self.__mask[::-1]
        res.background = {
            "A": self.background[T_or_U],
            "C": self.background["G"],
            "G": self.background["C"],
            T_or_U: self.background["A"],
        }
        res.pseudocounts = {
            "A": self.pseudocounts[T_or_U],
            "C": self.pseudocounts["G"],
            "G": self.pseudocounts["C"],
            T_or_U: self.pseudocounts["A"],
        }
        return res

    @property
    def consensus(self):
        """Return the consensus sequence."""
        return self.counts.consensus

    @property
    def anticonsensus(self):
        """Return the least probable pattern to be generated from this motif."""
        return self.counts.anticonsensus

    @property
    def degenerate_consensus(self):
        """Return the degenerate consensus sequence.

        Following the rules adapted from
        D. R. Cavener: "Comparison of the consensus sequence flanking
        translational start sites in Drosophila and vertebrates."
        Nucleic Acids Research 15(4): 1353-1361. (1987).

        The same rules are used by TRANSFAC.
        """
        return self.counts.degenerate_consensus

    @property
    def relative_entropy(self):
        """Return an array with the relative entropy for each column of the motif."""
        background = self.background
        pseudocounts = self.pseudocounts
        alphabet = self.alphabet
        counts = self.counts
        length = self.length
        values = np.zeros(length)
        if self.alignment is None:
            total = np.array([
                sum(counts[c][i] + pseudocounts[c] for c in alphabet)
                for i in range(length)
            ])
        else:
            total = np.zeros(length)
            for letter in alphabet:
                total += np.array([counts[letter][i] for i in range(length)]) + pseudocounts[letter]

        for letter in alphabet:
            frequencies = np.array([counts[letter][i] for i in range(length)]) + pseudocounts[letter]
            mask = frequencies > 0
            probs = frequencies[mask] / total[mask]
            values[mask] += probs * np.log2(probs / background[letter])

        return values
    #############################################################

    # def __format__(self, format_spec: Optional[str], **kwargs):
    #     """Return a string representation of the Motif in the given format.

    #     Currently supported formats:
    #      - clusterbuster: Cluster Buster position frequency matrix format
    #      - pfm : JASPAR single Position Frequency Matrix
    #      - jaspar : JASPAR multiple Position Frequency Matrix
    #      - transfac : TRANSFAC like files

    #     """
    #     if format_spec in ("pfm", "jaspar"):
    #         jaspar = importlib.import_module("Bio.motifs").jaspar

    #         motifs = [self]
    #         return str(jaspar.write(motifs, format_spec))
    #     elif format_spec == "transfac":
    #         transfac = importlib.import_module("Bio.motifs").transfac

    #         motifs = [self]
    #         return str(transfac.write(motifs))
    #     elif format_spec == "clusterbuster":
    #         clusterbuster = importlib.import_module("Bio.motifs").clusterbuster

    #         motifs = [self]
    #         return str(clusterbuster.write(motifs, **kwargs))
    #     elif not format_spec:
    #         # Follow python convention and default to using __str__
    #         return str(self)
    #     else:
    #         raise ValueError("Unknown format type %s" % format_spec)

    # def format(self, format_spec) -> str:
    #     """Return a string representation of the Motif in the given format.

    #     Currently supported formats:
    #      - clusterbuster: Cluster Buster position frequency matrix format
    #      - pfm : JASPAR single Position Frequency Matrix
    #      - jaspar : JASPAR multiple Position Frequency Matrix
    #      - transfac : TRANSFAC like files

    #     """
    #     return self.__format__(format_spec)




def create(instances: list[pyobj], alphabet: str = "ACGT") -> Motif:
    """Create a Motif object."""
    alignment = Alignment(instances)
    return Motif(alignment=alignment, alphabet=alphabet)

def create_motif_from_counts(alphabet: str, counts: dict[str, list[float]]) -> Motif:
    fpm = matrix.FrequencyPositionMatrix(alphabet, counts)
    return Motif(alphabet=alphabet, counts=fpm)


def parse(handle, fmt, strict=True):
    fmt = fmt.lower()
    if fmt == "minimal":
        import minimal
        return minimal.read(handle)
    
    else:
        raise Exception(f"Unknown format {fmt}")


# def read(handle, fmt, strict=True):
#     """Read a motif from a handle using the specified file-format.

#     This supports the same formats as Bio.motifs.parse(), but
#     only for files containing exactly one motif.  For example,
#     reading a JASPAR-style pfm file:

#     >>> from Bio import motifs
#     >>> with open("motifs/SRF.pfm") as handle:
#     ...     m = motifs.read(handle, "pfm")
#     >>> m.consensus
#     Seq('GCCCATATATGG')

#     Or a single-motif MEME file,

#     >>> from Bio import motifs
#     >>> with open("motifs/meme.psp_test.classic.zoops.xml") as handle:
#     ...     m = motifs.read(handle, "meme")
#     >>> m.consensus
#     Seq('GCTTATGTAA')

#     If the handle contains no records, or more than one record,
#     an exception is raised:

#     >>> from Bio import motifs
#     >>> with open("motifs/alignace.out") as handle:
#     ...     motif = motifs.read(handle, "AlignAce")
#     Traceback (most recent call last):
#         ...
#     ValueError: More than one motif found in handle

#     If however you want the first motif from a file containing
#     multiple motifs this function would raise an exception (as
#     shown in the example above).  Instead use:

#     >>> from Bio import motifs
#     >>> with open("motifs/alignace.out") as handle:
#     ...     record = motifs.parse(handle, "alignace")
#     >>> motif = record[0]
#     >>> motif.consensus
#     Seq('TCTACGATTGAG')

#     Use the Bio.motifs.parse(handle, fmt) function if you want
#     to read multiple records from the handle.

#     If strict is True (default), the parser will raise a ValueError if the
#     file contents does not strictly comply with the specified file format.
#     """
#     fmt = fmt.lower()
#     motifs = parse(handle, fmt, strict)
#     if len(motifs) == 0:
#         raise ValueError("No motifs found in handle")
#     if len(motifs) > 1:
#         raise ValueError("More than one motif found in handle")
#     motif = motifs[0]
#     return motif

# def write(motifs, fmt, **kwargs):
#     """Return a string representation of motifs in the given format.

#     Currently supported formats (case is ignored):
#      - clusterbuster: Cluster Buster position frequency matrix format
#      - pfm : JASPAR simple single Position Frequency Matrix
#      - jaspar : JASPAR multiple PFM format
#      - transfac : TRANSFAC like files

#     """
#     fmt = fmt.lower()
#     if fmt in ("pfm", "jaspar"):
#         jaspar = importlib.import_module("Bio.motifs").jaspar
#         jaspar_string = jaspar.write(motifs, fmt)
#         return str(jaspar_string)
    
#     elif fmt == "transfac":
#         transfac = importlib.import_module("Bio.motifs").transfac
#         transfac_string = transfac.write(motifs)
#         return str(transfac_string)
#     elif fmt == "clusterbuster":
#         clusterbuster = importlib.import_module("Bio.motifs").clusterbuster
#         clusterbuster_string = clusterbuster.write(motifs, **kwargs)
#         return str(clusterbuster_string)
#     else:
#         raise ValueError("Unknown format type %s" % fmt)

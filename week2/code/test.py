import math
import unittest
from unittest import test
import time
from typing import Tuple, NoneType

from python import importlib, pathlib, sys, os
from python import Bio

testCase = unittest.TestCase()

sys.path.append(os.getcwd())

try:
    from python import numpy as np
except Exception:
    raise Exception(
        "Install numpy if you want to use Bio.motifs."
    )

Seq = importlib.import_module("Bio.Seq").Seq
    
import bio_codon.motifs as motifs

# @test
# def test_format():
#     m = motifs.create([Seq("ATATA")])
#     m.name = "Foo"
    
#     expected_pfm = """  1.00   0.00   1.00   0.00  1.00
#     0.00   0.00   0.00   0.00  0.00
#     0.00   0.00   0.00   0.00  0.00
#     0.00   1.00   0.00   1.00  0.00
#     """

#     expected_jaspar = """>None Foo
#     A [  1.00   0.00   1.00   0.00   1.00]
#     C [  0.00   0.00   0.00   0.00   0.00]
#     G [  0.00   0.00   0.00   0.00   0.00]
#     T [  0.00   1.00   0.00   1.00   0.00]
#     """
   
#     expected_transfac = """P0      A      C      G      T
#     01      1      0      0      0      A
#     02      0      0      0      1      T
#     03      1      0      0      0      A
#     04      0      0      0      1      T
#     05      1      0      0      0      A
#     XX
#     //
#     """

#     s1 = motifs.Motif.format(m, "pfm")
#     s2 = motifs.Motif.format(m, "jaspar")
#     s3 = motifs.Motif.format(m, "transfac")

#     assert s1 == expected_pfm

#     assert s2 == expected_jaspar

#     assert s3 == expected_transfac

#     try:
#         motifs.Motif.format(m, "foo_bar")

#         assert False, "Expected ValueError but none was raised"
#     except ValueError:
#         pass
# test_format()

@test
def test_relative_entropy():
    m = motifs.create([Seq("ATATA"), Seq("ATCTA"), Seq("TTGTA")])

    length = len(m.alignment)
    assert length == 3

    exp_background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    testCase.assertEqual(m.background, exp_background)

    exp_pseudocounts = {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}
    assert m.pseudocounts == exp_pseudocounts

    exp_entropy = np.array([
        1.0817041659455104,
        2.0,
        0.4150374992788437,
        2.0,
        2.0
        ])
    testCase.assertTrue(
        np.allclose(m.relative_entropy, exp_entropy)
    )

    m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
    exp_entropy = np.array([
                    0.8186697601117167,
                    1.7369655941662063,
                    0.5419780939258206,
                    1.7369655941662063,
                    1.7369655941662063,
                ])
    testCase.assertTrue(np.allclose(m.relative_entropy, exp_entropy))
    

    m.background = None
    exp_background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    assert m.background == exp_background

    pseudocounts = math.sqrt(len(m.alignment))
    m.pseudocounts = {
        letter: m.background[letter] * pseudocounts for letter in "ACGT"
    }

    exp_entropy = np.array([
                    0.3532586861097656,
                    0.7170228827697498,
                    0.11859369972847714,
                    0.7170228827697498,
                    0.7170228827697499,
                ])
    
    testCase.assertTrue(
        np.allclose(m.relative_entropy, exp_entropy)
    )
        
    m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
    exp_entropy = np.array([
                    0.19727984803857979,
                    0.561044044698564,
                    0.20984910512125132,
                    0.561044044698564,
                    0.5610440446985638,
                ])
    
    testCase.assertTrue(np.allclose(m.relative_entropy, exp_entropy))
test_relative_entropy()

# @test
# def test_reverse_complement():
#     """Test if motifs can be reverse-complemented."""
#     background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
#     pseudocounts = 0.5
#     m = motifs.create([Seq("ATATA")])
#     m.background = background
#     m.pseudocounts = pseudocounts
#     received_forward = format(m, "transfac")
#     expected_forward = """\
#     P0      A      C      G      T
#     01      1      0      0      0      A
#     02      0      0      0      1      T
#     03      1      0      0      0      A
#     04      0      0      0      1      T
#     05      1      0      0      0      A
#     XX
#     //
#     """
#     assert received_forward == expected_forward

#     expected_forward_pwm = """\
#     0      1      2      3      4
#     A:   0.50   0.17   0.50   0.17   0.50
#     C:   0.17   0.17   0.17   0.17   0.17
#     G:   0.17   0.17   0.17   0.17   0.17
#     T:   0.17   0.50   0.17   0.50   0.17
#     """
#     assert str(m.pwm) == expected_forward_pwm

#     m = m.reverse_complement()
#     received_reverse = format(m, "transfac")
#     expected_reverse = """\
#     P0      A      C      G      T
#     01      0      0      0      1      T
#     02      1      0      0      0      A
#     03      0      0      0      1      T
#     04      1      0      0      0      A
#     05      0      0      0      1      T
#     XX
#     //
#     """
#     assert received_reverse == expected_reverse

#     expected_reverse_pwm = """\
#     0      1      2      3      4
#     A:   0.17   0.50   0.17   0.50   0.17
#     C:   0.17   0.17   0.17   0.17   0.17
#     G:   0.17   0.17   0.17   0.17   0.17
#     T:   0.50   0.17   0.50   0.17   0.50
#     """
#     assert str(m.pwm) == expected_reverse_pwm

#     # Same but for RNA motif.
#     background_rna = {"A": 0.3, "C": 0.2, "G": 0.2, "U": 0.3}
#     pseudocounts = 0.5
#     m_rna = motifs.create([Seq("AUAUA")], alphabet="ACGU")
#     m_rna.background = background_rna
#     m_rna.pseudocounts = pseudocounts
#     expected_forward_rna_counts = """\
#     0      1      2      3      4
#     A:   1.00   0.00   1.00   0.00   1.00
#     C:   0.00   0.00   0.00   0.00   0.00
#     G:   0.00   0.00   0.00   0.00   0.00
#     U:   0.00   1.00   0.00   1.00   0.00
#     """
#     assert str(m_rna.counts) == expected_forward_rna_counts

#     expected_forward_rna_pwm = """\
#     0      1      2      3      4
#     A:   0.50   0.17   0.50   0.17   0.50
#     C:   0.17   0.17   0.17   0.17   0.17
#     G:   0.17   0.17   0.17   0.17   0.17
#     U:   0.17   0.50   0.17   0.50   0.17
#     """
#     assert str(m_rna.pwm) == expected_forward_rna_pwm

#     expected_reverse_rna_counts = """\
#     0      1      2      3      4
#     A:   0.00   1.00   0.00   1.00   0.00
#     C:   0.00   0.00   0.00   0.00   0.00
#     G:   0.00   0.00   0.00   0.00   0.00
#     U:   1.00   0.00   1.00   0.00   1.00
#     """
#     assert str(m_rna.reverse_complement().counts) == expected_reverse_rna_counts

#     expected_reverse_rna_pwm = """\
#     0      1      2      3      4
#     A:   0.17   0.50   0.17   0.50   0.17
#     C:   0.17   0.17   0.17   0.17   0.17
#     G:   0.17   0.17   0.17   0.17   0.17
#     U:   0.50   0.17   0.50   0.17   0.50
#     """
#     assert str(m_rna.reverse_complement().pwm) == expected_reverse_rna_pwm

#     # Same thing, but now start with a motif calculated from a count matrix
#     m = motifs.create([Seq("ATATA")])
#     counts = m.counts
#     m = motifs.Motif(counts=counts)
#     m.background = background
#     m.pseudocounts = pseudocounts
#     received_forward = format(m, "transfac")
#     assert received_forward == expected_forward

#     assert str(m.pwm) == expected_forward_pwm

#     m = m.reverse_complement()
#     received_reverse = format(m, "transfac")
#     assert received_reverse == expected_reverse

#     assert str(m.pwm) == expected_reverse_pwm

#     # Same, but for RNA count matrix
#     m_rna = motifs.create([Seq("AUAUA")], alphabet="ACGU")
#     counts = m_rna.counts
#     m_rna = motifs.Motif(counts=counts, alphabet="ACGU")
#     m_rna.background = background_rna
#     m_rna.pseudocounts = pseudocounts
#     assert str(m_rna.counts) == expected_forward_rna_counts
#     assert str(m_rna.pwm) == expected_forward_rna_pwm
#     assert str(m_rna.reverse_complement().counts) == expected_reverse_rna_counts
#     assert str(m_rna.reverse_complement().pwm) == expected_reverse_rna_pwm
# test_reverse_complement()

@test
def test_minimal_meme_parser():
    """Parse data/minimal_test.meme file."""
    path = pathlib.Path(sys.argv[0]).resolve()
    current_dir = path.parent
    minimal_test_file = current_dir.parent / "data" / "minimal_test.meme"
    with open(minimal_test_file) as stream:
        record = motifs.parse(stream, "minimal")

    assert record.version == "4"
    assert record.alphabet == "ACGT"
    assert len(record.sequences) == 0
    assert record.command == ""
    assert len(record) == 3
    motif = next(iter(record))
    assert motif.name == "KRP"
    # testCase.assertEqual(record["KRP"], motif)
    assert motif.num_occurrences == 17
    assert motif.length == 19

    background = motif.background["A"]
    exp_value = 0.30269730269730266
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    background = motif.background["C"]
    exp_value = 0.1828171828171828
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    background = motif.background["G"]
    exp_value = 0.20879120879120877
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    background = motif.background["T"]
    exp_value = 0.30569430569430567
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    evalue = motif.evalue
    exp_value = 4.1e-09
    assert math.isclose(
        evalue, 
        exp_value, 
        rel_tol=1e-10,
        abs_tol=0.0
    )

    assert motif.alphabet == "ACGT"

    assert motif.alignment is None

    assert motif.consensus == "TGTGATCGAGGTCACACTT"

    assert motif.degenerate_consensus == "TGTGANNNWGNTCACAYWW"
    
    exp_entropy = np.array([
                    1.1684297174927525,
                    0.9432809925744818,
                    1.4307101633876265,
                    1.1549413780465179,
                    0.9308256303218774,
                    0.009164393966550805,
                    0.20124190687894253,
                    0.17618542656995528,
                    0.36777933103380855,
                    0.6635834532368525,
                    0.07729943368061855,
                    0.9838293592717438,
                    1.72489868427398,
                    0.8397561713453014,
                    1.72489868427398,
                    0.8455332015343343,
                    0.3106481207768122,
                    0.7382733641762232,
                    0.537435993300495,
                ])
    testCase.assertTrue(np.allclose(motif.relative_entropy, exp_entropy))

    # assert motif[2:9].consensus == "TGATCGA"

    motif = next(item for i, item in enumerate(record) if i == 1)
    assert motif.name == "IFXA"
    # assert record["IFXA"] == motif
    assert motif.num_occurrences, 14
    assert motif.length, 18

    background = motif.background["A"]
    exp_value = 0.30269730269730266
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )
    
    background = motif.background["C"]
    exp_value = 0.1828171828171828
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    background = motif.background["G"]
    exp_value = 0.20879120879120877
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    background = motif.background["T"]
    exp_value = 0.30569430569430567
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    evalue = motif.evalue
    exp_value = 3.2e-35
    assert math.isclose(
        evalue, 
        exp_value, 
        rel_tol=1e-36,
        abs_tol=0.0
    )

    assert motif.alphabet == "ACGT"

    assert motif.alignment is None

    assert motif.consensus == "TACTGTATATATATCCAG"

    assert motif.degenerate_consensus == "TACTGTATATAHAWMCAG"

    exp_entropy = np.array([
                    0.9632889858595118,
                    1.02677956765017,
                    2.451526420551951,
                    1.7098384161433415,
                    2.2598671267551107,
                    1.7098384161433415,
                    1.02677956765017,
                    1.391583804103081,
                    1.02677956765017,
                    1.1201961888781142,
                    0.27822438781180836,
                    0.36915366971717867,
                    1.7240522753630425,
                    0.3802185945622609,
                    0.790937683007783,
                    2.451526420551951,
                    1.7240522753630425,
                    1.3924085743645374,
                ])
    testCase.assertTrue(np.allclose(motif.relative_entropy, exp_entropy))

    # assert motif[2:9].consensus == "CTGTATA"
    

    motif = motif = next(item for i, item in enumerate(record) if i == 2)
    assert motif.name == "IFXA_no_nsites_no_evalue"
    # assert record["IFXA_no_nsites_no_evalue"] == motif
    assert motif.num_occurrences == 20
    assert motif.length == 18

    background = motif.background["A"]
    exp_value = 0.30269730269730266
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    background = motif.background["C"]
    exp_value = 0.1828171828171828
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    background = motif.background["G"]
    exp_value = 0.20879120879120877
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    background = motif.background["T"]
    exp_value = 0.30569430569430567
    assert math.isclose(
        background, 
        exp_value, 
        rel_tol=1e-7,
        abs_tol=0.0
    )

    evalue = motif.evalue
    exp_value = 0.0
    assert math.isclose(
        evalue, 
        exp_value, 
        rel_tol=1e-36,
        abs_tol=0.0
    )

    assert motif.alphabet == "ACGT"
    assert motif.alignment is None
    assert motif.consensus == "TACTGTATATATATCCAG"
    assert motif.degenerate_consensus == "TACTGTATATAHAWMCAG"

    exp_entropy = np.array([
                    0.99075309,
                    1.16078104,
                    2.45152642,
                    1.70983842,
                    2.25986713,
                    1.70983842,
                    1.16078104,
                    1.46052586,
                    1.16078104,
                    1.10213019,
                    0.29911041,
                    0.36915367,
                    1.72405228,
                    0.37696488,
                    0.85258086,
                    2.45152642,
                    1.72405228,
                    1.42793329,
                ])
    testCase.assertTrue(np.allclose(motif.relative_entropy, exp_entropy))

    # assert motif[2:9].consensus == "CTGTATA"
test_minimal_meme_parser()

    

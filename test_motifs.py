from Bio.Seq import Seq

from enhancers import Enhancer
from motifs import get_motif
from motifs import find_bs_motifs


def test_get_motif():
    m = get_motif("test_motif_elf1.pcm")

    assert m is not None
    assert len(m) == 14
    assert m.name == 'ELF1_HUMAN.H11MO.0.A'
    assert m.consensus == 'GAACCCGGAAGTGG'


def test_get_motif_no_header():
    m = get_motif("test_motif_no_header.pcm")
    assert m is not None
    assert len(m) == 1
    assert m.name is None
    assert m.consensus == 'G'


def test_find_bs_motifs():
    enhancer = Enhancer("1", 10, 55, "ensembleId", True)
    enhancer.sequence = "AGTCCCTCGCAACAGGGACTCTGCGAACCCGGAAGTGGACGTGGC"
    motif_elf1 = get_motif("test_motif_elf1.pcm")  # should match on + strand
    motif_hsf1 = get_motif("test_motif_hsf1.pcm")  # should not match
    motif_revs = get_motif("test_motif_reverse.pcm")  # should match on both strands
    hits = find_bs_motifs(enhancer, [motif_elf1, motif_hsf1, motif_revs])

    assert hits is not None
    assert len(hits) == 3
    assert hits[0] == (24, Seq('GAACCCGGAAGTGG'), '+')
    assert hits[1] == (13, Seq('AGGGACT'), '+')
    assert hits[2] == (38, Seq('AGGGACT'), '-')

    #     0 1 2 3 4 5 6 7 8 9 10 11 12   coordinates
    # 5'  G G G A C G C T T A A A A +   3'  -->
    # 3'  C C C T G C G A A T T T T     5'  <--
    #        <--  G C G A A  <---


def test_find_bs_motifs_multiple_times():
    enhancer = Enhancer("1", 10, 55, "ensembleId", True)
    enhancer.sequence = "AGTCGAACCCGGAAGTGGACTCACGAACCCGGAAGTGGACGTGGC"
    motif_elf1 = get_motif("test_motif_elf1.pcm")  # should match twice on + strand
    hits = find_bs_motifs(enhancer, [motif_elf1])

    assert hits is not None
    assert len(hits) == 2
    assert hits[0] == (4, Seq('GAACCCGGAAGTGG'), '+')
    assert hits[1] == (24, Seq('GAACCCGGAAGTGG'), '+')


def test_find_bs_motifs_no_match():
    enhancer = Enhancer("1", 10, 55, "ensembleId", True)
    enhancer.sequence = "AGTCGAACACCGGAAGTGGACTCACGAACCCGGGAAGTGGACGGC"
    motif_elf1 = get_motif("test_motif_elf1.pcm")  # should match twice on + strand
    hits = find_bs_motifs(enhancer, [motif_elf1])

    assert hits is not None
    assert len(hits) == 0


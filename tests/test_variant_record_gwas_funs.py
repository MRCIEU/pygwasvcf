import pygwasvcf
import os
import pytest
from math import log10

CHROM = "1"
START = 49297
STOP = 49298
TRAIT = "UKB-b:13008"
FILE = os.path.join(os.path.dirname(__file__), "data", "case.control.example.vcf.gz")


def test_get_pval():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_pval(rec, TRAIT) == pytest.approx(0.9599999)


def test_get_se():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_se(rec, TRAIT) == pytest.approx(1.62715e-05)


def test_get_beta():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_beta(rec, TRAIT) == pytest.approx(-8.975e-07)


def test_get_af():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_af(rec, TRAIT) == pytest.approx(0.623765)


def test_get_id():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_id(rec, TRAIT) == "rs10399793"


def test_get_ss():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            with pytest.raises(KeyError):
                assert pygwasvcf.VariantRecordGwasFuns.get_ss(rec, TRAIT) == 0


def test_transform_pval():
    p = 0.01
    logp = -log10(p)
    assert logp == 2
    assert pygwasvcf.VariantRecordGwasFuns.transform_logpval(2) == 0.01

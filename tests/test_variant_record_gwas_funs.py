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


def test_get_id_rsid():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_id(rec, TRAIT) == "rs10399793"


def test_get_id_chrpos():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_id(rec, TRAIT, create_if_missing=False) is not None
            del rec.samples[TRAIT]['ID']
            with pytest.raises(KeyError):
                assert pygwasvcf.VariantRecordGwasFuns.get_id(rec, TRAIT, create_if_missing=False)
            assert pygwasvcf.VariantRecordGwasFuns.get_id(rec, TRAIT, create_if_missing=True) == "1-49298-T-C"


def test_get_ss():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            with pytest.raises(KeyError):
                assert pygwasvcf.VariantRecordGwasFuns.get_ss(rec, TRAIT) == 0


def test_get_ss_from_metadata():
    with pygwasvcf.GwasVcf(FILE) as g:
        metadata = g.get_metadata()
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_ss(rec, TRAIT, metadata) == (463001 + 9)
        del metadata[TRAIT]['TotalCases']
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_ss(rec, TRAIT, metadata) == 463001


def test_get_nc():
    with pygwasvcf.GwasVcf(FILE) as g:
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            with pytest.raises(KeyError):
                assert pygwasvcf.VariantRecordGwasFuns.get_nc(rec, TRAIT)


def test_get_nc_from_metadata():
    with pygwasvcf.GwasVcf(FILE) as g:
        metadata = g.get_metadata()
        for rec in g.query(contig=CHROM, start=START, stop=STOP):
            assert pygwasvcf.VariantRecordGwasFuns.get_nc(rec, TRAIT, metadata) == 9
        del metadata[TRAIT]['TotalCases']
        with pytest.raises(KeyError):
            for rec in g.query(contig=CHROM, start=START, stop=STOP):
                pygwasvcf.VariantRecordGwasFuns.get_nc(rec, TRAIT, metadata)


def test_transform_pval():
    p = 0.01
    logp = -log10(p)
    assert logp == 2
    assert pygwasvcf.VariantRecordGwasFuns.transform_logpval(2) == 0.01

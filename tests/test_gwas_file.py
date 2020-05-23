import pygwasvcf
import pysam
import os
import sqlite3
import pytest

FILE = os.path.join(os.path.dirname(__file__), "data", "case.control.example.vcf.gz")


def test_get_sample_metadata():
    g = pygwasvcf.GwasVcf(FILE)
    recs = g.get_sample_metadata()
    assert "UKB-b:13008" in recs
    assert "TotalVariants" in recs["UKB-b:13008"]
    assert "VariantsNotRead" in recs["UKB-b:13008"]
    assert "HarmonisedVariants" in recs["UKB-b:13008"]
    assert "VariantsNotHarmonised" in recs["UKB-b:13008"]
    assert "SwitchedAlleles" in recs["UKB-b:13008"]
    assert "TotalControls" in recs["UKB-b:13008"]
    assert "TotalCases" in recs["UKB-b:13008"]
    assert "StudyType" in recs["UKB-b:13008"]
    g.close()


def test_format_variant_record_for_rsidx():
    vcf = pysam.VariantFile(FILE)
    for rec in vcf.fetch():
        for rsid, chrom, pos in pygwasvcf.GwasVcf.format_variant_record_for_rsidx(rec):
            assert chrom is not None
            assert isinstance(chrom, str)
            assert pos is not None
            assert isinstance(pos, int)
            assert rsid is not None
            assert isinstance(rsid, int)
    vcf.close()


def test_index_rsid():
    # index GWAS-VCF using rsidx
    g = pygwasvcf.GwasVcf(FILE)
    g.index_rsid()
    g.close()
    # check index exists
    assert os.path.exists(FILE + ".rsidx")
    # check contents of index
    with sqlite3.connect(FILE + ".rsidx") as dbconn:
        cur = dbconn.cursor()
        cur.execute("SELECT * FROM rsid_to_coord")

        for rec in cur.fetchall():
            assert rec[0] is not None
            assert isinstance(rec[0], int)
            assert rec[1] is not None
            assert isinstance(rec[1], str)
            assert rec[2] is not None
            assert isinstance(rec[2], int)


def test_get_location_from_rsid():
    g = pygwasvcf.GwasVcf(FILE)
    g.index_rsid()
    chrom, start, end = g.get_location_from_rsid("rs10399793")
    g.close()
    assert chrom == "1"
    assert start == 49298
    assert end == 49298
    g = pygwasvcf.GwasVcf(FILE, rsidx_path=FILE + ".rsidx")
    chrom, start, end = g.get_location_from_rsid("rs10399793")
    g.close()
    assert chrom == "1"
    assert start == 49298
    assert end == 49298


def check_first_row(row):
    assert row.chrom == "1"
    assert row.pos == 49298
    assert row.get_pval("UKB-b:13008") == pytest.approx(0.9599999)


def test_query():
    g = pygwasvcf.GwasVcf(FILE)
    g.index_rsid()
    for num, row in enumerate(g.query(chrom="1", start=49298, end=49298)):
        assert num == 1
        check_first_row(row)
    for num, row in enumerate(g.query(variant_id="rs10399793")):
        assert num == 1
        check_first_row(row)
    g.close()

import pygwasvcf
import pysam
import os
import sqlite3
import pytest

FILE = os.path.join(os.path.dirname(__file__), "data", "case.control.example.vcf.gz")


def test_close():
    with pygwasvcf.GwasVcf(FILE) as g:
        assert not g.is_closed()
    assert g.is_closed()
    g = pygwasvcf.GwasVcf(FILE)
    assert g.is_closed()


def test_get_metadata():
    with pygwasvcf.GwasVcf(FILE) as g:
        recs = g.get_metadata()
        assert "UKB-b:13008" in recs
        assert "TotalVariants" in recs["UKB-b:13008"]
        assert "VariantsNotRead" in recs["UKB-b:13008"]
        assert "HarmonisedVariants" in recs["UKB-b:13008"]
        assert "VariantsNotHarmonised" in recs["UKB-b:13008"]
        assert "SwitchedAlleles" in recs["UKB-b:13008"]
        assert "TotalControls" in recs["UKB-b:13008"]
        assert "TotalCases" in recs["UKB-b:13008"]
        assert "StudyType" in recs["UKB-b:13008"]


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
    # delete old index if present
    if os.path.exists(FILE + ".rsidx"):
        os.remove(FILE + ".rsidx")

    # index GWAS-VCF
    with pygwasvcf.GwasVcf(FILE) as g:
        g.index_rsid()

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
    with pygwasvcf.GwasVcf(FILE) as g:
        g.index_rsid()
        chrom, pos = g.get_location_from_rsid("rs10399793")
        assert chrom == "1"
        assert pos == 49298

    with pygwasvcf.GwasVcf(FILE, rsidx_path=FILE + ".rsidx") as g:
        chrom, pos = g.get_location_from_rsid("rs10399793")
        assert chrom == "1"
        assert pos == 49298


def check_first_row(row):
    assert row.chrom == "1"
    assert row.pos == 49298
    assert row.get_pval("UKB-b:13008") == pytest.approx(0.9599999)


def test_query():
    with pygwasvcf.GwasVcf(FILE) as g:
        g.index_rsid()
        for num, row in enumerate(g.query(contig="1", start=49297, stop=49298)):
            assert num == 0
            check_first_row(row)
        for num, row in enumerate(g.query(variant_id="rs10399793")):
            assert num == 0
            check_first_row(row)
        j = 0
        for num, row in enumerate(g.query()):
            j = num
            if num == 0:
                check_first_row(row)
        assert j > 0


def test_vcf_query_can_be_null():
    vcf = pysam.VariantFile(FILE)
    j = 0
    for num, row in enumerate(vcf.fetch(contig=None, start=None, stop=None)):
        j = num
        if num == 0:
            check_first_row(row)
    assert j > 0
    vcf.close()

def test_class_assignment():
    vcf = pysam.VariantFile(FILE)
    for rec in vcf.fetch():
        print(type(rec), rec)
        rec.__class__ = pygwasvcf.VariantRecordGwas
    vcf.close()

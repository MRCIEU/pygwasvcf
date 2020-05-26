import pygwasvcf
import pysam
import os
import sqlite3

TRAIT = "UKB-b:13008"
FILE = os.path.join(os.path.dirname(__file__), "data", "case.control.example.vcf.gz")


def check_first_row(chrom, pos):
    assert chrom == "1"
    assert pos == 49298


def test_close():
    with pygwasvcf.GwasVcf(FILE) as g:
        assert not g.is_closed()
    assert g.is_closed()
    g = pygwasvcf.GwasVcf(FILE)
    assert g.is_closed()


def test_get_metadata():
    with pygwasvcf.GwasVcf(FILE) as g:
        recs = g.get_metadata()
        assert TRAIT in recs
        assert "TotalVariants" in recs[TRAIT]
        assert "VariantsNotRead" in recs[TRAIT]
        assert "HarmonisedVariants" in recs[TRAIT]
        assert "VariantsNotHarmonised" in recs[TRAIT]
        assert "SwitchedAlleles" in recs[TRAIT]
        assert "TotalControls" in recs[TRAIT]
        assert "TotalCases" in recs[TRAIT]
        assert "StudyType" in recs[TRAIT]


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
        check_first_row("1", 49298)

    with pygwasvcf.GwasVcf(FILE, rsidx_path=FILE + ".rsidx") as g:
        chrom, pos = g.get_location_from_rsid("rs10399793")
        assert chrom == "1"
        assert pos == 49298


def test_query_by_chr_pos():
    with pygwasvcf.GwasVcf(FILE) as g:
        for num, row in enumerate(g.query(contig="1", start=49297, stop=49298)):
            check_first_row(row.chrom, row.pos)
        assert num == 0


def test_query_by_rsid():
    with pygwasvcf.GwasVcf(FILE) as g:
        g.index_rsid()
        for num, row in enumerate(g.query(variant_id="rs10399793")):
            check_first_row(row.chrom, row.pos)
        assert num == 0


def test_query_all():
    with pygwasvcf.GwasVcf(FILE) as g:
        for num, row in enumerate(g.query()):
            if num == 0:
                check_first_row(row.chrom, row.pos)
        assert num > 0

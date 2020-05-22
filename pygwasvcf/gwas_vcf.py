from pysam import VariantFile
from variant_record_gwas import VariantRecordGwas
import sqlite3
import os

"""Class to parse GWAS-VCF file
"""


class GwasVcf:
    def __init__(self, file_path, rsidx_path=None):
        """
        Constructor for GwasVcf class
        :param file_path: Path to GWAS-VCF
        :param rsidx_path: Path to RSIDX (optional)
        """
        self.file_path = file_path
        self.vcf = VariantFile(file_path)
        self.rsidx_path = rsidx_path

    def close(self):
        """
        Close GWAS-VCF file handle
        """
        self.vcf.close()

    @staticmethod
    def format_variant_record_for_rsidx(rec):
        """
        Function to extract data for RSIDX query
        :param rec:
        :return: rsid: dbSNP identifier as integer (rs removed)
        :return: chrom: chromosome for association
        :return: pos: base-position for association
        """
        for assoc in rec.samples:
            var_id = rec.samples[assoc]['ID']
            if var_id is not None:
                for rsid in var_id.split(";"):
                    if rsid[0:2] == "rs":
                        yield int(rsid[2:]), rec.chrom, rec.pos

    def index_rsid(self):
        """
        Index GWAS-VCF using rsID and pysam adapted from [rsidx](https://github.com/bioforensics/rsidx)
        """
        idx_path = self.file_path + ".rsidx"
        if os.path.exists(idx_path):
            os.remove(idx_path)

        with sqlite3.connect(idx_path) as dbconn:
            # prepare database
            c = dbconn.cursor()
            c.execute(
                'CREATE TABLE rsid_to_coord ('
                'rsid INTEGER PRIMARY KEY, '
                'chrom TEXT NULL DEFAULT NULL, '
                'coord INTEGER NOT NULL DEFAULT 0)'
            )
            dbconn.commit()

            # add variant records to SQLite DB
            for rec in self.vcf.fetch():
                c.executemany('INSERT OR IGNORE INTO rsid_to_coord VALUES (?,?,?)',
                              GwasVcf.format_variant_record_for_rsidx(rec))
                dbconn.commit()

        self.rsidx_path = idx_path

    def get_sample_metadata(self):
        """
        Extract metadata about the GWAS trait
        :return: res: Dict of Dict containing a key=value pairs for each trait in the GWAS-VCF
        """
        res = dict()
        for rec in self.vcf.header.records:
            if rec.key == "SAMPLE":
                res[rec['ID']] = dict()
                for k in rec:
                    if k != "ID":
                        res[rec['ID']][k] = rec[k]
        return res

    def get_location_from_rsid(self, rsid):
        """
        Helper function to convert rsID to chromosome and position using [rsidx](https://github.com/bioforensics/rsidx)
        :param rsid: dbsnp indentifier
        :return: res: chromosome
        :return: res: start
        :return: res: end
        """
        if not rsid.startswith("rs"):
            raise ValueError("Variant ID query must be an rsID")
        if self.rsidx_path is None:
            raise ValueError("Cannot query by variant identifier without providing an rsidx")
        q = (int(rsid[2:]),)
        with sqlite3.connect(self.rsidx_path) as dbconn:
            cur = dbconn.cursor()
            cur.execute('SELECT DISTINCT chrom,coord FROM rsid_to_coord WHERE rsid =?', q)
            res = cur.fetchone()
        return res[0], res[1], res[1]

    def query(self, chrom=None, start=None, end=None, variant_id=None, exclude_filtered=True):
        """
        Variant-trait association query function
        :param chrom: Chromosome to query
        :param start: Start position of interval (1-based)
        :param end: End position of interval (1-based)
        :param variant_id: rsID to query using [rsidx](https://github.com/bioforensics/rsidx)
        :param exclude_filtered: Boolean flag to remove record that do not meet QC
        :return: rec: VariantRecordGwas object containing chromosome, position, alleles, association statistics
        """
        if variant_id is None:
            if chrom is None or start is None or end is None:
                raise ValueError(
                    "You must provide a genomic location range or variant identifier to perform the query.")
        else:
            if chrom is not None or start is not None or end is not None:
                raise ValueError("Cannot provide chromosome, start or end with variant ID. Choose one query.")
            chrom, start, end = self.get_location_from_rsid(variant_id)

        # extract variant(s) from GWAS-VCF
        for rec in self.vcf.fetch(chrom, start, end):
            # Extend to VariantRecord to provide useful funcs for GWAS assoc
            rec.__class__ = VariantRecordGwas
            rec.check_biallelic()

            # skip variants not meeting filter requirements
            if exclude_filtered and rec.filter != "PASS":
                continue

            # lazy return record
            yield rec

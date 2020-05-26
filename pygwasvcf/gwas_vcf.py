from pysam import VariantFile
import pygwasvcf
import sqlite3
import os

"""
Class to parse GWAS-VCF file using pysam
"""


class GwasVcf:
    def __init__(self, file_path, rsidx_path=None):
        """
        Constructor for GwasVcf class
        :param file_path: Path to GWAS-VCF
        :param rsidx_path: Path to RSIDX (optional)
        """
        self.__file_path = file_path
        self.__rsidx_path = rsidx_path
        self.__vcf = None

    def __enter__(self):
        """
        Open VariantFile using with resources
        """
        self.__vcf = VariantFile(self.__file_path)
        return self

    def __exit__(self, type, value, traceback):
        """
        Close GWAS-VCF file handle
        """
        self.__vcf.close()
        self.__vcf = None

    def is_closed(self):
        return self.__vcf is None

    @staticmethod
    def format_variant_record_for_rsidx(rec):
        """
        Function to extract data for RSIDX query
        :param rec: pysam.VariantRecord
        :return: rsid: dbSNP identifier as integer (rs removed)
        :return: chrom: chromosome for association
        :return: pos: position for association (1-based)
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
        if self.__vcf is None:
            raise ValueError("Cannot use methods on the VCF object before the file is open.")

        idx_path = self.__file_path + ".rsidx"
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
            for rec in self.__vcf.fetch():
                c.executemany('INSERT OR IGNORE INTO rsid_to_coord VALUES (?,?,?)',
                              GwasVcf.format_variant_record_for_rsidx(rec))
                dbconn.commit()

        self.__rsidx_path = idx_path

    def get_metadata(self):
        """
        Extract metadata about the GWAS trait(s)
        :return: res: Dict of Dict containing a key=value pairs for each trait in the GWAS-VCF
        """
        if self.__vcf is None:
            raise ValueError("Cannot use methods on the VCF object before the file is open.")
        res = dict()
        for rec in self.__vcf.header.records:
            if rec.key == "SAMPLE":
                res[rec['ID']] = dict()
                for k in rec:
                    if k != "ID":
                        res[rec['ID']][k] = rec[k]
        return res

    def get_traits(self):
        """
        Extract list of traits in the GWAS-VCF
        :return: traits: List of traits
        """
        return list(self.__vcf.header.samples)

    def get_location_from_rsid(self, rsid):
        """
        Helper function to convert rsID to chromosome and position using [rsidx](https://github.com/bioforensics/rsidx)
        :param rsid: dbsnp indentifier
        :return: res: chromosome
        :return: res: position (1-based)
        """
        if not rsid.startswith("rs"):
            raise ValueError("Variant ID query must be an rsID")
        if self.__rsidx_path is None:
            raise ValueError("Cannot query by variant identifier without providing an rsidx")
        q = (int(rsid[2:]),)
        with sqlite3.connect(self.__rsidx_path) as dbconn:
            cur = dbconn.cursor()
            cur.execute('SELECT DISTINCT chrom,coord FROM rsid_to_coord WHERE rsid =?', q)
            res = cur.fetchone()
        return res[0], res[1]

    def query(self, contig=None, start=None, stop=None, variant_id=None, exclude_filtered=True):
        """
        Variant-trait association query function
        :param contig: Chromosome to query
        :param start: Start position of interval (0-based)
        :param stop: End position of interval (0-based)
        :param variant_id: rsID to query using [rsidx](https://github.com/bioforensics/rsidx)
        :param exclude_filtered: Boolean flag to remove record that do not meet QC
        :return: rec: pysam.VariantRecord object containing chromosome, position, alleles, association statistics
        """
        if self.__vcf is None:
            raise ValueError("Cannot use methods on the VCF object before the file is open.")

        if variant_id is not None:
            if contig is not None or start is not None or stop is not None:
                raise ValueError("Cannot provide chromosome, start or end with variant ID. Choose one query.")
            contig, pos = self.get_location_from_rsid(variant_id)
            start = pos - 1
            stop = pos

        # extract variant(s) from GWAS-VCF
        for rec in self.__vcf.fetch(contig=contig, start=start, stop=stop):
            # check multiallelics are on separate rows which is required for functions
            pygwasvcf.VariantRecordGwasFuns.check_biallelic(rec)

            # skip variants not meeting filter requirements
            if exclude_filtered and rec.filter.keys()[0] != "PASS":
                continue

            # lazy return record
            yield rec

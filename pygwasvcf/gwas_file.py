from pysam import VariantFile


class GwasFile:
    def __init__(self, file_path, rsidx_path=None):
        self.vcf = VariantFile(file_path)
        self.rsidx_path = rsidx_path

    def close(self):
        """Close GWAS-VCF file handle"""
        self.vcf.close()

    def get_sample_metadata(self):
        """Extract metadata about the GWAS trait"""
        res = []
        for rec in self.vcf.header.records:
            if rec.type == "SAMPLE":
                res.append(rec)
        return res

    @staticmethod
    def parse(rec):
        """Helper function to extract GWAS related variables from the VCF record"""
        return rec

    def get_location_for_rsid(self, rsid):
        """Helper function to convert rsID to chromosome and position using [rsidx](https://github.com/bioforensics/rsidx)"""
        if self.rsidx_path is None:
            raise ValueError("Cannot query by variant identifier without providing an rsidx")
        return None, None, None

    def query(self, chrom=None, start=None, end=None, variant_id=None, studies_to_include=None, exclude_filtered=True,
              pval_threshold=None):
        """Variant-trait association query function"""

        if variant_id is not None:
            if chrom is not None or start is not None or end is not None:
                raise ValueError("Cannot provide chromosome, start or end with variant ID. Choose one query.")
            chrom, start, end = self.get_location_for_rsid(variant_id)
        else:
            if chrom is None or start is None or end is None:
                raise ValueError(
                    "You must provide a genomic location range or variant identifier to perform the query.")

        # extract variant(s) from GWAS-VCF
        for rec in self.vcf.fetch(chrom, start, end):

            # skip variants not meeting filter requirements
            if exclude_filtered and rec.filter != "PASS":
                continue

            # lazy return parsed record
            yield GwasFile.parse(rec)

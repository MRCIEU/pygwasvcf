"""
Class to provide helper functions on the pysam.VariantRecord object
"""


# subclassing pysam.VariantRecord does not seem to work (I think this must be done in Cython)
class VariantRecordGwasFuns:
    @staticmethod
    def transform_logpval(p):
        return 10 ** -p

    @staticmethod
    def get_pval(variant_record, trait):
        p = variant_record.samples[trait]['LP'][0]
        if p == 0:
            return 1
        elif p == 999:
            return 0
        else:
            return VariantRecordGwasFuns.transform_logpval(p)

    @staticmethod
    def get_se(variant_record, trait):
        return variant_record.samples[trait]['SE'][0]

    @staticmethod
    def get_beta(variant_record, trait):
        return variant_record.samples[trait]['ES'][0]

    @staticmethod
    def get_af(variant_record, trait):
        return variant_record.samples[trait]['AF'][0]

    @staticmethod
    def get_id(variant_record, trait):
        return variant_record.samples[trait]['ID'][0]

    @staticmethod
    def get_ss(variant_record, trait):
        return variant_record.samples[trait]['SS'][0]

    @staticmethod
    def check_biallelic(variant_record):
        assert len(variant_record.alts) == 1

"""
Class to provide helper functions on the pysam.VariantRecord object
"""


# subclassing pysam.VariantRecord does not seem to work (I think this must be done in Cython)
class VariantRecordGwasFuns:
    """
    Transforms -log10 P back to P value
    :param p: -log10 P value
    :return P value in the range 0,1
    """

    @staticmethod
    def transform_logpval(p):
        return 10 ** -p

    """
    Getter for the variant-trait association P value
    :param variant_record: pysam.VariantRecord object for the VCF row
    :param trait: Name of the trait
    :return P value in the range 0,1
    """

    @staticmethod
    def get_pval(variant_record, trait):
        p = variant_record.samples[trait]['LP'][0]
        if p == 0:
            return 1
        elif p == 999:
            return 0
        else:
            return VariantRecordGwasFuns.transform_logpval(p)

    """
    Getter for the variant-trait association standard error value
    :param variant_record: pysam.VariantRecord object for the VCF row
    :param trait: Name of the trait
    :return Association standard error
    """

    @staticmethod
    def get_se(variant_record, trait):
        return variant_record.samples[trait]['SE'][0]

    """
    Getter for the variant-trait beta value
    :param variant_record: pysam.VariantRecord object for the VCF row
    :param trait: Name of the trait
    :return effect size coefficient
    """

    @staticmethod
    def get_beta(variant_record, trait):
        return variant_record.samples[trait]['ES'][0]

    """
    Getter for the variant-trait allele-frequency in the study
    :param variant_record: pysam.VariantRecord object for the VCF row
    :param trait: Name of the trait
    :return alternative (effect) allele frequency
    """

    @staticmethod
    def get_af(variant_record, trait):
        return variant_record.samples[trait]['AF'][0]

    """
    Getter for the variant-trait variant ID
    :param variant_record: pysam.VariantRecord object for the VCF row
    :param trait: Name of the trait
    :param create_if_missing: Create ID in the format chrom-pos-ref-alt if no ID is available
    :return Variant/marker identifier
    """

    @staticmethod
    def get_id(variant_record, trait, create_if_missing=False):
        if "ID" in variant_record.samples[trait] and variant_record.samples[trait]['ID'] is not None:
            return variant_record.samples[trait]['ID']
        elif create_if_missing:
            return variant_record.chrom + "-" + str(variant_record.pos) + "-" + variant_record.ref + "-" + \
                   variant_record.alts[0]
        else:
            raise KeyError("No ID available for this record")

    """
    Getter for the variant-trait sample size used to estimate the effect
    :param variant_record: pysam.VariantRecord object for the VCF row
    :param trait: Name of the trait
    :param metadata: If the per-variant sample size if missing then taken from global metadata (optional)
    :return Total sample size used to estimate the association effect size
    """

    @staticmethod
    def get_ss(variant_record, trait, metadata=None):
        if 'SS' in variant_record.samples[trait]:
            return variant_record.samples[trait]['SS'][0]
        elif metadata is not None and 'TotalControls' in metadata[trait]:
            if 'TotalCases' in metadata[trait]:
                return int(metadata[trait]['TotalControls']) + int(metadata[trait]['TotalCases'])
            else:
                return int(metadata[trait]['TotalControls'])
        else:
            raise KeyError("No sample size available")

    """
    Getter for the variant-trait number of cases size used to estimate the effect
    :param variant_record: pysam.VariantRecord object for the VCF row
    :param trait: Name of the trait
    :param metadata: If the per-variant sample size if missing then taken from global metadata (optional)
    :return Number of cases used to estimate the association effect size
    """

    @staticmethod
    def get_nc(variant_record, trait, metadata=None):
        if 'NC' in variant_record.samples[trait]:
            return variant_record.samples[trait]['NC'][0]
        elif metadata is not None and 'TotalCases' in metadata[trait]:
            return int(metadata[trait]['TotalCases'])
        else:
            raise KeyError("No sample size available")

    """
    Check the VCF record only contains a bi-allelic change; multi-allelic changes should be on separate rows 
    :param variant_record: pysam.VariantRecord object for the VCF row
    """

    @staticmethod
    def check_biallelic(variant_record):
        assert len(variant_record.alts) == 1

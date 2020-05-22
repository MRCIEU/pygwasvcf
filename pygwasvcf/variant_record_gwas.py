from pysam.libcbcf import VariantRecord


class VariantRecordGwas(VariantRecord):
    @staticmethod
    def transform_pval(p):
        return 10 ** -p

    def get_pval(self, trait):
        p = self.samples[trait]['LP'][0]
        if p == 0:
            return 1
        elif p == 999:
            return 0
        else:
            return VariantRecordGwas.transform_pval(p)

    def check_biallelic(self):
        assert len(self.alts) == 1

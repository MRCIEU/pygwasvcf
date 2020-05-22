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

    def get_se(self, trait):
        return (self.samples[trait]['SE'][0])

    def get_beta(self, trait):
        return (self.samples[trait]['ES'][0])

    def get_af(self, trait):
        return (self.samples[trait]['AF'][0])

    def check_biallelic(self):
        assert len(self.alts) == 1

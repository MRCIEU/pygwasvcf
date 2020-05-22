from variant_record_gwas import VariantRecordGwas
from math import log10


def test_transform_pval():
    p = 0.01
    logp = -log10(p)
    assert logp == 2
    assert VariantRecordGwas.transform_pval(2) == 0.01

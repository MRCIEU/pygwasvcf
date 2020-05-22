from gwas_file import GwasFile


def test_gwas_file():
    g = GwasFile("/Users/ml/GitLab/pygwasvcftools/pygwasvcf/test/data/case.control.example.vcf.gz")
    g.get_sample_metadata()
    g.close()

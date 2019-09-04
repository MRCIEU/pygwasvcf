from gwas_file import GwasFile


def test_gwas_file():
    g = GwasFile("/Users/ml/GitLab/pygwasvcftools/pygwasvcftools/test/data/case.control.example.vcf.gz",
                 '/Users/ml/Desktop/human_g1k_v37.fasta')
    f = g.read(exclude_filtered=False)
    f.close()

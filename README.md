# GWAS-VCF Python parser

<!-- badges: start -->
[![Build Status](https://travis-ci.org/MRCIEU/pygwasvcf.svg?branch=master)](https://travis-ci.org/MRCIEU/pygwasvcf)
<!-- badges: end -->

The package provides a thin wrapper around [pysam](https://pysam.readthedocs.io/en/latest/index.html) and [rsidx](https://github.com/bioforensics/rsidx) to support extraction of GWAS metadata.

## Install

```shell script
pip install git+https://github.com/mrcieu/pygwasvcf
```

## GWAS-VCF files

Download over 10,000 GWAS-VCF files from the [IEU GWAS database](https://gwas.mrcieu.ac.uk/)

## Examples

Read GWAS trait/study metadata

```python
from pygwasvcf.gwas_vcf import GwasVcf
g = GwasVcf("/path/to/gwasvcf.vcf.gz")

# print dictionary of GWAS metadata
print(g.get_sample_metadata())

# always close when done to release resources
g.close()
```

Query variant-trait association(s) by genomic location

```python
from pygwasvcf.gwas_vcf import GwasVcf
g = GwasVcf("/path/to/gwasvcf.vcf.gz")

# query by chromosome and position interval
for variant in g.query(chrom="1", start=1, end=1):
    print(variant)

# index on dbSNP identifier
# based on [rsidx](https://github.com/bioforensics/rsidx)
# only need to do this once and then provide the index path to the constructor
# i.e. GwasVcf("/path/to/gwasvcf.vcf.gz", rsidx_path="/path/to/gwasvcf.vcf.gz.rsidx")
g.index_rsid()

# rapid query by rsID  
for variant in g.query(variant_id="rs1245"):
    print(variant)

# always close when done to release resources
g.close()
```

Extract summary statistics for variant

```python
from pygwasvcf.gwas_vcf import GwasVcf
g = GwasVcf("/path/to/gwasvcf.vcf.gz")

# query by chromosome and position interval
for variant in g.query(chrom="1", start=1, end=1):
    # print variant-trait P value
    print(variant.get_pval("trait_name"))
    # print variant-trait SE
    print(variant.get_se("trait_name"))
    # print variant-trait beta
    print(variant.get_beta("trait_name"))
    # print variant-trait allele frequency in study
    print(variant.get_af("trait_name"))
# always close when done to release resources
g.close()
```
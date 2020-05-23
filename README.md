# GWAS-VCF Python parser

<!-- badges: start -->
[![Build Status](https://travis-ci.org/MRCIEU/pygwasvcf.svg?branch=master)](https://travis-ci.org/MRCIEU/pygwasvcf)
<!-- badges: end -->

The package provides a thin wrapper around [pysam](https://pysam.readthedocs.io/en/latest/index.html) and [rsidx](https://github.com/bioforensics/rsidx) to parse VCF files containing GWAS summary statistics and trait metadata.

Parses GWAS-VCF with version 1.0 of the [specification](https://github.com/MRCIEU/gwas-vcf-specification/releases/tag/1.0.0)

## Install

```shell script
pip install git+https://github.com/mrcieu/pygwasvcf
```

## Summary statistics in GWAS-VCF

Download over 10,000 GWAS-VCF files contain full summary statistics from the [IEU GWAS database](https://gwas.mrcieu.ac.uk/)

## Parser examples

Read GWAS trait/study metadata

```python
import pygwasvcf
g = pygwasvcf.GwasVcf("/path/to/gwas.vcf.gz")

# print dictionary of GWAS metadata
print(g.get_sample_metadata())

# always close when done to release resources
g.close()
```

Query variant-trait association(s) by chromosome and position location

```python
import pygwasvcf
g = pygwasvcf.GwasVcf("/path/to/gwas.vcf.gz")

# query by chromosome and position interval
for variant in g.query(chrom="1", start=1, end=1):
    print(variant)

# always close when done to release resources
g.close()
```

Query variant-trait association(s) by dbSNP rsID

```python
import pygwasvcf
g = pygwasvcf.GwasVcf("/path/to/gwas.vcf.gz")

# index on dbSNP identifier
# based on [rsidx](https://github.com/bioforensics/rsidx)
# only need to do this once and then provide the index path to the constructor
# i.e. GwasVcf("/path/to/gwas.vcf.gz", rsidx_path="/path/to/gwas.vcf.gz.rsidx")
g.index_rsid()

# rapid query by rsID  
for variant in g.query(variant_id="rs1245"):
    print(variant)

# always close when done to release resources
g.close()
```

Extract summary statistics from a variant object

```python
import pygwasvcf
g = pygwasvcf.GwasVcf("/path/to/gwas.vcf.gz")

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

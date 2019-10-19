import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pygwasvcf",
    version="0.0.1",
    author="Matt Lyon",
    author_email="ml18692@bristol.ac.uk",
    description="A package for reading GWAS summary statistics stored in VCF/BCF format",
    long_description="Genome wide association studies produce variant-trait association metrics across the genome. Often data are stored in plain text files but this leads to data non-standardisation and poor query performance. VCF/BCF format is a specific format for storing and querying / manipulating genetic data with high performance. This package can read such files.",
    long_description_content_type="text/markdown",
    url="https://github.com/MRCIEU/pygwasvcf",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

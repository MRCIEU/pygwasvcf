import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pygwasvcf",
    version="0.0.4",
    author="Matt Lyon",
    author_email="ml18692@bristol.ac.uk",
    description="A package for reading GWAS summary statistics stored in VCF/BCF format",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MRCIEU/pygwasvcf",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["pysam", "pytest"]
)

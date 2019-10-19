import subprocess
from pysam import VariantFile
import numpy as np
import logging


class GwasFile:
    # TODO
    BCFTOOLS_BIN = "/Users/ml/GitLab/pygwasvcf/pygwasvcf/bcftools-1.9/bcftools"

    def __init__(self, file_path, genome_path):
        self.file_path = file_path
        self.genome_path = genome_path

    def read(self, chrom=None, start=None, end=None, variant_id=None, studies_to_include=None, exclude_filtered=True,
             pval_threshold=None):

        view_args = [
            GwasFile.BCFTOOLS_BIN,
            "view",
            "-O",
            "u"
        ]

        if exclude_filtered:
            view_args.append("-i")
            view_args.append("'FILTER=PASS'")

        if chrom is not None:
            if start is None:
                view_args.append("-r")
                view_args.append("{}".format(chrom))
            elif end is None:
                view_args.append("-r")
                view_args.append("{}:{}".format(chrom, start))
            else:
                view_args.append("-r")
                view_args.append("{}:{}-{}".format(chrom, start, end))

        if studies_to_include is not None:
            view_args.append("-s")
            view_args.append(",".join(studies_to_include))

        if pval_threshold is not None:
            view_args.append("-i")
            view_args.append("'FORMAT/LP > {}'".format(-np.log10(pval_threshold)))

        if variant_id is not None:
            view_args.append("-i")
            view_args.append("'ID={}'".format(variant_id))

        view_args.append(self.file_path)

        norm_args = [
            GwasFile.BCFTOOLS_BIN,
            "norm",
            "-c", "e",
            "-f", self.genome_path,
            "-m", "-any",
            "-O", "u"
        ]

        logging.debug("bcftool cmd: {}".format(view_args))
        logging.debug("bcftool cmd: {}".format(norm_args))

        # execute bcftools command and parse with pysam
        # TODO use tmp file?
        with subprocess.Popen(view_args, stdout=subprocess.PIPE) as view_proc:
            with subprocess.Popen(norm_args, stdout=subprocess.PIPE, stdin=view_proc.stdout) as norm_proc:
                return VariantFile(norm_proc.stdout)

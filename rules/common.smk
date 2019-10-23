import shutil
import sys
import os
import os.path
import re
import pandas as pd
import numpy as np
import subprocess
from collections import defaultdict
import yaml


def message(mes):
    sys.stderr.write("*** " + mes + " ***\n")

#def get_fastq(wildcards):
#    return reads.loc[wildcards.sample, ["read1", "read2"]].dropna()

def get_fastq1(wildcards):
    return reads.loc[wildcards.sample, ["read1"]].dropna()

def get_fastq2(wildcards):
    return reads.loc[wildcards.sample, ["read2"]].dropna()

def get_gtf(wildcards):
    if wildcards.annot == "ref":
        return GTF
    elif wildcards.annot == "new":
        return os.path.join(WORKDIR, "new/assembling/merged.gtf")

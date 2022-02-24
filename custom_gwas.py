#!/usr/bin/env python

import sys
import argparse
from cyvcf2 import VCF
from os import path
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description = "Code to run custom GWAS with VCF files")
parser.add_argument('--vcf', type = str, help = "Path to tabix'ed VCF", required = True)
parser.add_argument('--pheno', type = str, help = "Path to phenotype file", required = True)
parser.add_argument('--covar', type = str, help = "Path to covariate file", required = True)
parser.add_argument('--model', type = str, help = "Path to association test function", required = True)
parser.add_argument('--out', type = str, help = "Output prefix", required = True)

np.set_printoptions(precision = 3)
args = parser.parse_args()

if path.exists(args.vcf):
    vcf_file = VCF(args.vcf)
    samples = vcf_file.samples
    for variant in vcf_file:
        dosages = variant.format('DS')
        print(samples)
        print(dosages.ravel())
        print(dosages.shape)
        print(dosages.ndim)
        print(dosages.size)
        print(dosages.dtype)
        sys.exit()
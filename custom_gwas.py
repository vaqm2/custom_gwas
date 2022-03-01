#!/usr/bin/env python

import sys
import argparse
from tarfile import PAX_NAME_FIELDS
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

if path.exists(args.pheno):
    pheno_df = pd.read_csv(args.pheno, sep = ' ')
else:
    print('ERROR: Unable to read PHENOTYPE file: ')
    print(args.pheno)
    sys.exit()

if path.exists(args.covar):
    covar_df = pd.read_csv(args.covar, sep = ' ')
else:
    print('ERROR: Unable to read COVARIATE file: ')
    print(args.covar)
    sys.exit()

pheno_cov_df = pd.merge(pheno_df, covar_df, on = 'IID')

if path.exists(args.vcf):
    vcf_file = VCF(args.vcf)
    samples_df = pd.DataFrame(vcf_file.samples, columns = ['IID'])
    
    for variant in vcf_file:
        dosages = variant.format('DS')
        dosages_df = pd.DataFrame(dosages, columns = ['DOSAGE'])
        samples_ds_df = pd.concat([samples_df.reset_index(drop = True), dosages_df.reset_index(drop = True)], axis = 1)
        to_regress_df = pd.merge(pheno_cov_df, samples_ds_df, on = 'IID')
        print(to_regress_df.head())
        sys.exit()
#!/usr/bin/env python

import sys
import argparse
from cyvcf2 import VCF
from os import path
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def file_not_found_error(file_path):
    print("ERROR: Unable to locate file: " + file_path + "\n")
    sys.exit()

parser = argparse.ArgumentParser(description = "Code to run custom GWAS with VCF files")
parser.add_argument('--vcf', type = str, help = "Path to TABIX'ed VCF", required = True)
parser.add_argument('--pheno', type = str, help = "Path to PHENOTYPE file", required = True)
parser.add_argument('--covar', type = str, help = "Path to COVARIATE file", required = True)
parser.add_argument('--rscript', type = str, help = "Path to the Rscript that implements the desired model", required = True)
parser.add_argument('--model', type = str, help = "Name of the model within the Rscript that you wish to fit", required = True)
parser.add_argument('--out', type = str, help = "Output prefix", required = True)

np.set_printoptions(precision = 3)
args = parser.parse_args()

if path.exists(args.pheno):
    pheno_df = pd.read_csv(args.pheno, sep = ' ')
else:
    file_not_found_error(args.pheno)

if path.exists(args.covar):
    covar_df = pd.read_csv(args.covar, sep = ' ')
else:
    file_not_found_error(args.covar)

if path.exists(args.rscript):
    robjects.r['source'](args.rscript)
    model_to_fit = robjects.globalenv[args.model]
else:
    file_not_found_error(args.rscript)

pheno_cov_df = pd.merge(pheno_df, covar_df, on = 'IID')

if path.exists(args.vcf):
    vcf_file   = VCF(args.vcf)
    samples_df = pd.DataFrame(vcf_file.samples, columns = ['IID'])
    out_fh     = open(args.out, "w")
    out_fh.writelines("CHR POS SNP REF ALT ESTIMATE SE P\n")
    
    for variant in vcf_file:
        refAllele = variant.REF
        altAllele = ''.join(variant.ALT)

        if(len(refAllele) > 1 or len(altAllele) > 1): # Ignore multi-allelic loci
            continue
        else:
            dosages       = variant.format('DS')
            dosages_df    = pd.DataFrame(dosages, columns = ['DOSAGE'])
            samples_ds_df = pd.concat([samples_df.reset_index(drop = True), dosages_df.reset_index(drop = True)], axis = 1)
            to_regress_df = pd.merge(pheno_cov_df, samples_ds_df, on = 'IID')
            with localconverter(robjects.default_converter + pandas2ri.converter): # Convert pandas dataframe to R dataframe
                to_regress_df_r = robjects.conversion.py2rpy(to_regress_df)
            result_df_r = model_to_fit(to_regress_df_r)
            with localconverter(robjects.default_converter + pandas2ri.converter): # Convert R dataframe back to pandas dataframe
                result_df_pd = robjects.conversion.rpy2py(result_df_r)
            out_fh.write(variant.CHROM + " ")
            out_fh.write(str(variant.start + 1) + " ") # CyVCF2 returns zero-based start positions
            out_fh.write(variant.ID + " ")
            out_fh.write(refAllele + " ")
            out_fh.write(altAllele + " ") # Accounts for multi-allelic sites by merging CyVCF2 returned character array
            out_fh.write(result_df_pd.estimate.to_string(index = False) + " ")
            out_fh.write(result_df_pd.se.to_string(index = False) + " ")
            out_fh.write(result_df_pd.p.to_string(index = False) + "\n")
    out_fh.close()
else:
    print("ERROR: Unable to read VCF file: ")
    print(args.vcf)
    print("\n")
    sys.exit()
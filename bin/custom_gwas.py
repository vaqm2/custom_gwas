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
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description = "Code to run custom GWAS with VCF files")
    parser.add_argument('--vcf', 
        type = str, 
        help = "Path to TABIX'ed VCF", 
        required = True)
    parser.add_argument('--pheno', 
        type = str, 
        help = "Path to PHENOTYPE file", 
        required = True)
    parser.add_argument('--covar',
        type = str,
        help = "Path to COVARIATE file",
        required = True)
    parser.add_argument('--rscript', 
        type = str, 
        help = "Path to the Rscript that implements the desired model", 
        required = True)
    parser.add_argument('--model', 
        type = str,
        help = "Name of the model within the Rscript that you wish to fit", 
        required = True)
    parser.add_argument('--out', 
        type = str, 
        help = "Output prefix",
        required = True)

    np.set_printoptions(precision = 4)
    args = parser.parse_args()

    try:
        pheno_df = pd.read_csv(args.pheno, sep = "\s+")
        pheno_df.IID = pheno_df.IID.astype(str)
    except:
        print("ERROR: When opening PHENOTYPE file: ", sys.exc_info()[0], "occurred!")
        sys.exit()

    try:
        covar_df = pd.read_csv(args.covar, sep = "\s+")
        covar_df.IID = covar_df.IID.astype(str)     
    except:
        print("ERROR: When opening COVARIATE file: ", sys.exc_info()[0], "occurred!")
        sys.exit()

    try:
        robjects.r['source'](args.rscript)
        model_to_fit = robjects.globalenv[args.model]
    except:
        print("ERROR: When opening rscript: ", sys.exc_info()[0], "occurred!")
        sys.exit()

    try:
        pheno_cov_df = pd.merge(pheno_df, covar_df, on = 'IID')
    except:
        print("ERROR: When intersecting PHENOTYPE & COVARIATE files: ", 
            sys.exc_info()[0], 
            "occurred!")
        sys.exit()

    try:
        vcf_file = VCF(args.vcf)
    except:
        print("ERROR: ", sys.exc_info()[0], "occurred!")
    else:
        samples_df = pd.DataFrame(vcf_file.samples, columns = ['IID'])
        samples_df.IID = samples_df.IID.astype(str)
        try:
            out_fh = open(args.out, "w")
        except:
            print("ERROR: ", sys.exc_info()[0], "occurred!\n")
        else:
            outputBuffer = "CHR POS SNP REF ALT ESTIMATE SE P N\n"
            linesProcessed = 0

            for variant in vcf_file:
                chromosome = variant.CHROM
                position = str(variant.start + 1)
                rsId = variant.ID
                refAllele = variant.REF
                altAllele = ''.join(variant.ALT)
                
                # Ignore multi-allelic loci
                if(len(refAllele) > 1 or len(altAllele) > 1):
                    print("INFO: Skipping INDEL/MULTI-ALLELIC variant: ", 
                        rsId, " ",
                        chromosome, " ", 
                        position, " ", 
                        refAllele, " ", 
                        altAllele)
                    continue
                else:
                    try:
                        dosages = variant.format('DS')
                    except:
                        print("INFO: Missing DS format field at: ",
                            rsId, " ",
                            chromosome, " ",
                            position, " ",
                            refAllele, " ",
                            altAllele, " ",
                            "Skipping..")
                    else:
                        dosages_df = pd.DataFrame(dosages, columns = ['DOSAGE'])
                        samples_ds_df = pd.concat([samples_df.reset_index(drop = True), 
                            dosages_df.reset_index(drop = True)], 
                            axis = 1)
                        to_regress_df = pd.merge(pheno_cov_df, samples_ds_df, how = 'inner', on = 'IID')
                        N = str(to_regress_df.shape[0])

                    # Convert pandas dataframe to R dataframe
                    with localconverter(robjects.default_converter + pandas2ri.converter): 
                        to_regress_df_r = robjects.conversion.py2rpy(to_regress_df)

                    # Call the user-specified R function to perform the association test
                    try:
                        result_df_r = model_to_fit(to_regress_df_r)
                    except:
                        print("ERROR: When fitting user-specified model: ", 
                            sys.exc_info()[0], 
                            "occurred!")

                    # Convert R dataframe back to pandas dataframe
                    with localconverter(robjects.default_converter + pandas2ri.converter):
                        result_df_pd = robjects.conversion.rpy2py(result_df_r)

                    outputBuffer += chromosome + ""
                    outputBuffer += position + " "
                    outputBuffer += rsId + " "
                    outputBuffer += refAllele + " "
                    outputBuffer += altAllele + " "
                    outputBuffer += result_df_pd.estimate.to_string(index = False) + " "
                    outputBuffer += result_df_pd.se.to_string(index = False) + " "
                    outputBuffer += result_df_pd.p.to_string(index = False) + " "
                    outputBuffer += N + "\n"
                    linesProcessed += 1

                    if(linesProcessed > 0 and linesProcessed % 10000 == 0):
                        print ("PROGRESS: Processed ", linesProcessed, " ", "variants at ", datetime.now())

            # Empty buffer to output all at once
            out_fh.write(outputBuffer)
            out_fh.close()
      
if __name__ == "__main__":
    main()
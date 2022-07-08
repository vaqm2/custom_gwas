#!/usr/bin/env python

from multiprocessing.dummy import freeze_support
import sys
import argparse
import multiprocessing as mp
import itertools as its
from cyvcf2 import VCF
from sgkit.io.vcf import partition_into_regions
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

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
    num_cpus = mp.cpu_count()

    try:
        pheno_df = pd.read_csv(args.pheno, sep = "\s+")
    except:
        print("ERROR: When opening PHENOTYPE file: ", sys.exc_info()[0], "occurred!")
        sys.exit()
    else:
        pheno_df.IID = pheno_df.IID.astype(str)
        pheno_df.set_index('IID', inplace = True)

    try:
        covar_df = pd.read_csv(args.covar, sep = "\s+")
    except:
        print("ERROR: When opening COVARIATE file: ", sys.exc_info()[0], "occurred!")
        sys.exit()
    else:
        covar_df.IID = covar_df.IID.astype(str)
        covar_df.set_index('IID', inplace = True)

    try:
        pheno_cov_df = pheno_df.join(covar_df, how = 'inner')
    except:
        print("ERROR: When intersecting PHENOTYPE & COVARIATE files: ", sys.exc_info()[0], "occurred!")
        sys.exit()

    try:
        regions = partition_into_regions(args.vcf, num_parts = num_cpus)
        print(regions)
        num_chunks = len(regions)
    except:
        print("ERROR: When partitioning VCF file ", sys.exc_info()[0], " occurred!")
        sys.exit()
    else:
        try:
            out_fh = open(args.out, "w")
        except:
            print("ERROR: When creating output file ", sys.exc_info()[0], " occurred!")
            sys.exit()
        else:
            try:
                vcf_file = VCF(args.vcf)
            except:
                print("ERROR: When reading VCF file ", sys.exc_info()[0], " occurred!")
                sys.exit()
            else:
                vcf_samples_df = pd.DataFrame()
                vcf_samples_df['IID'] = vcf_file.samples
                vcf_samples_df.IID = vcf_samples_df.IID.astype(str)
                vcf_samples_df.set_index('IID', inplace = True)
                analysis_df = vcf_samples_df.join(pheno_cov_df, how = 'left')
                out_fh.write("CHR POS SNP REF ALT ESTIMATE SE P N\n")
                vcf_file.close()

                # Start multiple processes for different VCF chunk regions
                p = mp.Pool(num_chunks)
                model_args = []
                for i in range(0, num_chunks):
                    model_args.append((args.vcf, regions[i], analysis_df, pheno_cov_df.columns, args.rscript, args.model))
                results = p.starmap(run_model, model_args)

                # Write results to output file
                for chunk_result in results:
                    out_fh.write(chunk_result)
                out_fh.close()
        
def run_model(vcf_path, region, analysis_df, covariate_names, r_script, r_model):
    try:
        robjects.r['source'](r_script)
    except:
        print("ERROR: When opening rscript: ", sys.exc_info()[0], " occurred!")
        sys.exit()
    else:
        model_to_fit = robjects.globalenv[r_model]

    try:
        vcf_file = VCF(vcf_path)
    except:
        print("ERROR: When reading VCF file ", sys.exc_info()[0], " occurred!")
    else:
        chunk_result = ""
        for variant in vcf_file(region):
            chromosome = variant.CHROM
            position = str(variant.start + 1)
            rsId = variant.ID
            refAllele = variant.REF
            altAllele = ''.join(variant.ALT)

            # Ignore multi-allelic loci
            if(len(refAllele) > 1 or len(altAllele) > 1):
                print("INFO: Skipping INDEL/MULTI-ALLELIC variant: ", chromosome, " ", position, " ", refAllele, " ", altAllele)
                continue
            else:
                try:
                    dosages = variant.format('DS')
                except:
                    print("INFO: Missing DS format field at: ", chromosome, " ", position, " ", refAllele, " ", altAllele, " ", "Skipping..")
                else:
                    analysis_df['DOSAGE'] = dosages
                    analysis_df.dropna(subset = covariate_names)
                    N = str(analysis_df.shape[0])

                    # Convert pandas dataframe to R dataframe
                    with localconverter(robjects.default_converter + pandas2ri.converter): 
                        analysis_df_r = robjects.conversion.py2rpy(analysis_df)

                    # Call the user-specified R function to perform the association test
                    try:
                        result_df_r = model_to_fit(analysis_df_r)
                    except:
                        print("ERROR: When fitting user-specified model: ", sys.exc_info()[0], "occurred!")
                        sys.exit()

                    # Convert R dataframe back to pandas dataframe
                    with localconverter(robjects.default_converter + pandas2ri.converter):
                        result_df_pd = robjects.conversion.rpy2py(result_df_r)

                    chunk_result += chromosome + " "
                    chunk_result += position + " "
                    chunk_result += rsId + " "
                    chunk_result += refAllele + " "
                    chunk_result += altAllele + " "
                    chunk_result += result_df_pd.estimate.to_string(index = False) + " "
                    chunk_result += result_df_pd.se.to_string(index = False) + " "
                    chunk_result += result_df_pd.p.to_string(index = False) + " "
                    chunk_result += N + "\n"
        vcf_file.close()
        return chunk_result
      
if __name__ == "__main__":
    freeze_support()
    main()
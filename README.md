# custom_gwas
Python Code using CyVCF2 to run custom GWAS models on imputed VCF files

## AUTHOR

[Vivek Appadurai](https://github.com/vaqm2) 

Post Doc at [Institute of Biological Psychiatry](https://github.com/biopsyk)
## PURPOSE

This is code written in python, that uses the [CyVCF2](https://brentp.github.io/cyvcf2/) library to parse imputed, bgzip'ed and tabix indexed VCF files and accepts custom models written in R as a command line option to the script to run a genome wide association study, using [rpy2](https://rpy2.github.io/doc/v3.4.x/html/index.html), an interface, that facilitates r-functions to be seamlessly run in python code.

## USAGE

Usage can be obtained by running the following command:

    python custom_gwas.py -h

which prints:

    usage: custom_gwas.py [-h] --vcf VCF --pheno PHENO --covar COVAR --rscript RSCRIPT --model MODEL --out OUT

where

VCF => The bgzip, tab indexed VCF file containing dosage values in the format column 'DS'

PHENO => File of phenotypes with samples listed under the column name IID and a phenotype name of your choice

COVAR => File of covariates with samples listed under the column name IID and your choice of covariates, weights etc.

RSCRIPT => some file.R that contains one or more functions pertaining to the association models you wish to fit.

MODEL => Name of the particular function in the file.R that you wish to fit

OUT => Output file name

An example for the VCF, PHENO, COVAR and OUT files are provided in the test/ folder.

An example of a custom R-script and function is provided at /bin/custom_function.R

Ultimately the freedom is given to the user to specify columns of their choice in the phenotype and covariate files and thereby using them in the Rscript they provide in specifying the model. The only major requirements are that:

1. VCF contains DS column in the FORMAT fields
2. The pheno and covariate files can be interesected using an IID column
3. The custom gwas Rscript returns a dataframe with the column names:
    estimate
    se
    p

## RUNNING A TEST

The script can be tested on the test data using the following command:

    python /Users/vapp0002/Documents/postdoc_ibp_work/custom_gwas/bin/custom_gwas.py \
    --vcf test/test.vcf.gz \
    --pheno test/test_pheno.txt \
    --covar test/test_covar.txt \
    --rscript bin/custom_function.R \ 
    --model custom_fit \
    --out test_out.txt

## DEPENDENCIES

The code has the following dependencies:

    python
    cyvcf2
    numpy
    pandas
    rpy2
    R
    dplyr
    broom

All dependencies can be installed using conda or the requisite environment can be created using the environment.yml file with the following command:

    conda create --file environment.yml

For any questions or bug reports, please create an issue on github.
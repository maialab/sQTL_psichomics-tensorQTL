import pandas as pd
import numpy as np
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import matplotlib.pyplot as plt

# define paths to data
plink_prefix_path = 'phenotype.bed'
expression_bed = 'genotype.bed.gz'
covariates_file = 'covariates.txt'
prefix = 'results'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# map all cis-associations (results for each chromosome are written to file)

# nominal p-values for all genes
cis = cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, prefix)
cis.to_csv('cis.csv')

#empirical p-values for  all genes
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
cis_df.to_csv('cis_results.csv')

# to limit output size, only associations with p-value <= 1e-5 are returned
trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=20000,
                           return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
trans_df.to_csv('trans_results.csv')
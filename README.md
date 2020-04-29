# Alternative Splicing in Breast Cancer Risk
Breast cancer has a strong genetic risk component, and recently,
cis-regulatory single nucleotide polymorphisms (SNP) have been
associated with it. Current efforts are focused on the full
understanding of the cis-regulatory mechanisms involved.
Most studies functionally characterising GWAS loci for breast
cancer, have focused solely on the effect of regulatory SNPs
(rSNPs) on transcription factor binding at promoters and enhancers.
However, sequence changes can also have potential effects on, for
example, splicing, microRNA (miRNA) activity and epigenetic
regulation. Here we have initiated the study of genetic variants affecting
alternative splicing in breast cancer, by identifying splicing QTLs which are associated with risk.
We have used psichomics to quantify alternative splicing isoforms in normal breast RNA-seq data (phs000424.v8.p2 NHGRI GTEx), and mapped the sQTL using tensorQTL. Then we compared the list of significant sQTLs with breast cancer known risk variants, by using gwarapidd to retrieve data from the GWAS Catalog.
Results presented are preliminary.

This repository has all code developed for this end divided into 3 parts:

1) Psi calculation using psichomics;

2) Mapping using tensorQTL;

3) co-localization analysis with GWAS risk loci using Gwasrapidd and Ensemblr.

## Original Articles:
###[Gwasrapidd](https://github.com/ramiromagno/gwasrapidd)
Magno, R. and Maia, A.-T. (2019) ‘gwasrapidd: an R package to query, download and wrangle GWAS catalog data’, Bioinformatics. Oxford University Press (OUP). doi: 10.1093/bioinformatics/btz605.

###[Ensemblr](https://github.com/ramiromagno/ensemblr)
Magno, R. Unpublished Work

### [Psichomics](https://github.com/nuno-agostinho/psichomics)
Saraiva-Agostinho, N. and Barbosa-Morais, N. L. (2019) ‘psichomics: graphical application for alternative splicing quantification and analysis’, Nucleic Acids Research. Oxford University Press, 47(2), pp. e7–e7. doi: 10.1093/nar/gky888.


### [TensorQTL](https://github.com/broadinstitute/tensorqtl)
Taylor-Weiner, A. et al. (2019) ‘Scaling computational genomics to millions of individuals with GPUs’, Genome Biology. BioMed Central Ltd., 20(1), p. 228. doi: 10.1186/s13059-019-1836-7.

## Data
All data is available from [GTEx](https://gtexportal.org/home/).

## Contact info
For any questions mail me through
aarbduarte[at] gmail[dot]com

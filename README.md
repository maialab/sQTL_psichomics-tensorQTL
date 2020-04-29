# Alternative Splicing in Breast Cancer Risk
Breast cancer has a strong genetic risk component, and recently, cis-regulatory single nucleotide polymorphisms (SNP) have been associated with it. Current efforts are focused on the full understanding of the cis-regulatory mechanisms involved.
Most studies functionally characterising GWAS loci for breast cancer, have focused solely on the effect of regulatory SNPs (rSNPs) on transcription factor binding at promoters and enhancers.
However, sequence changes can also have potential effects on, for example, splicing, microRNA (miRNA) activity and epigenetic regulation. We have initiated studying genetic variants affecting alternative splicing in breast cancer, by identifying splicing QTLs. However, the best approach to studying cis-regulation is by measuring allelic effects, as such the most powerful method to apply in our studies is Allele Specific Alternative Splicing analysis. Here we aim to extend the in silico analysis that we have been conducting to determine the possible contribution of cis-regulatory SNPs influencing splicing to breast cancer susceptibility. We will (1) identify breast cancer risk variants that are associated with alternative splice isoforms (using bioinformatics tools such as PAIRADISE and ASARP, and publicly available data from GTEx project) and (2) perform in-silico and in-vitro functional validation of the best candidates identified.

This repository has all code developed for this end divided into 3 parts:
1) Psi calculation using psichomics;
2) Mapping using tensorQTL;
3) co-localization analysis with GWAS risk loci using GwasRapidd and Ensemblr.

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

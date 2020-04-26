# TensorQTL
TensorQTL implements phenotype mapping with GPU accelaration for a faster run

It needs 3 different inputs:
* Phenotype
* Genotype
* Covariates


The phenotype file needs to be BED format and later zipped and indexed
```
bgzip phenotypes.bed && tabix -p bed phenotypes.bed.gz
```

The genotype needs to be in PLINK format
```
plink2 --make-bed --output-chr chrM --vcf ${plink_prefix_path}.vcf.gz --out ${plink_prefix_path}
```

The covariates is a simple covariates x samples .txt file
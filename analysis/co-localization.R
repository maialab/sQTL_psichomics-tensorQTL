library(dplyr)
library(ensemblr)
library(gwasrapidd)
library(stringr)
library(biomaRt)

#genomic window size(in kb):
window_size <- 100L
#LD query, R^2 threshold:
LD <- 0.8

#alternative presented by joana: 1) see what phenotypes were included into breast cancer; 2) select relevant traits
studies <- get_studies(efo_trait = "breast carcinoma")
table (studies@studies$reported_trait) #chooses which traits to keep.
#choose which traits to include
phenotype <- c("Breast cancer",
               "Breast cancer (early onset)", 
               "Breast cancer (estrogen-receptor negative)", 
               "Breast cancer (male)", 
               "Breast Cancer in BRCA1 mutation carriers", 
               "Breast cancer in BRCA2 mutation carriers", 
               "breast cancer male", 
               "Breast cancer and/or colorectal cancer")
rows_to_keep <-which(studies@studies$reported_trait %in% phenotype) #dá-me as linhas com os estudos de interesse
studies_to_keep <- studies@studies$study_id[rows_to_keep]
studies@ancestral_groups

for(i in 1:length(studies_to_keep)){
  assign(studies_to_keep[i], get_variants(study_id = studies_to_keep[i]))
}
###
populations <- get_populations()
pop_loc <- str_locate_all(populations$population,":")
populations$pop_acr <- rep(NA,nrow(populations))

for(i in 1: nrow(populations)){
  populations$pop_acr[i] <- str_sub(populations$population[i], unlist(pop_loc[[i]][2,2])+1, nchar(populations$population[i]))
}
###

#the following list was made taking into account information available from Table 1, Morales et al 2018 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5815218/table/Tab1/?report=objectonly)
dictionary <- list(European = c("CEU","FIN","GBR","IBS","TSI"),
                   Other = c(),
                   `East Asian` =c("CDX", "CHB", "CHS", "JPT"),
                   `African American or Afro-Caribbean` = c("ACB", "ASW"),
                   `Hispanic or Latin American`=c("CLM", "MXL", "PEL", "PUR" ),
                   `NA` = c(),
                   `Asian unspecified` = c("CDX", "CHB", "CHS", "JPT", "KHV","BEB", "GIH", "ITU", "PJL", "STU" ),
                   `Sub-Saharan African` = c("ESN", "LWK", "GWD", "MSL", "MKK", "YRI"),
                   `African unspecified` = c("ACB", "ASW", "ESN", "LWK", "GWD", "MSL", "MKK", "YRI"))


in_LD <- data.frame(matrix(ncol = 6))
colnames(in_LD)<- c("study","population","variant_id1", "variant_id2", "r_squared","d_prime")

#function to retrieve population names to input on ensemblr
dic_search <- function(study){
  z <- sum(studies@ancestral_groups$study_id==study)
  all_pops <-c()
  for(j in 1:z){
      all_pops <- append(all_pops, dictionary[[studies@ancestral_groups$ancestral_group[1]]])
  }
  all_pops <- populations$population[which(populations$pop_acr %in% unique(all_pops))]
  return(all_pops)
}


#Retrieve all snps in LD for the specific populations
for(i in 1:length(studies_to_keep)){
  pops <-dic_search(studies_to_keep[i])
  num_pops <- length(pops)
  for(k in 1:num_pops){
    z <- get_ld_variants_by_window(get(studies_to_keep[i])@variants$variant_id, 
                                                           genomic_window_size = window_size,
                                                           r_squared = LD,
                                                           population =pops[k])
    if(nrow(z)>0){
      in_LD<- bind_rows(in_LD, z)
      in_LD[(nrow(in_LD)-nrow(z)+1):nrow(in_LD),1] <- rep(studies_to_keep[i],nrow(z))
    }
  }
}

#reorganize data
in_LD <- in_LD[2:nrow(in_LD),c(1,2,8,9,5,6)]
unique(in_LD$variant_id2)#14788 unique snps

#retrieve genomic positions for snps
gwas_snps <- getBM(attributes = c('refsnp_id','chr_name','chrom_start','allele','minor_allele_freq','synonym_name'),
                   filters = c('snp_filter'),
                   values = list(unique(in_LD$variant_id2)),
                   mart = snpmart)
#load sQTL results and genotype positions and ids
sQTL <- read.csv("~/psi_tens/permutations.csv", as.is = T)
positions <- read.delim("~/files_map/geno_pos.txt",as.is = T, sep = " ")
sQTL$chr <- sQTL$pos <- rep(NA,nrow(sQTL))

#cross sQTL snp id with position on genotype
for(i in 1:nrow(sQTL)){
  posi <- which(sQTL$variant_id[i]==positions$ID)
  sQTL$chr[i] <- str_sub(positions$X.CHROM[posi],4L,nchar(posi))
  sQTL$pos[i] <- positions$POS[posi]
  print(paste0(i,"/",28603))
}

connection <- data.frame(sQTL=sQTL$phenotype_id, gwas=rep(NA,nrow(sQTL)))

for(j in 1:nrow(sQTL)){
  if(sQTL$pos[j] %in% gwas_snps$chrom_start){
    q <- which(gwas_snps$chrom_start == sQTL$pos[j])
    if(sQTL$pos[j] %in% gwas_snps$chrom_start[q]){
      connection$gwas[j] <- gwas_snps$refsnp_id[q]
    }else{
      connection$gwas[j] <- F
    }
  }else{
    connection$gwas[j] <- F
  }
}  
coloco <- sQTL[which(connection$gwas != F),]
coloco <- coloco[coloco$pval_perm<0.01,]

gwas_snps[which(gwas_snps$chrom_start %in% coloco$pos),]

write.csv(coloco[coloco$pval_perm<0.01,],file = "coloco")
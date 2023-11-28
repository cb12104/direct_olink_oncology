##Code to conduct Mendelian randomization and colocalization analyses in "Impact of weight loss on cancer-related proteins in serm:results for a randomized-controlled trial"

#proteins: GPNMB, FURIN, WIF-1, TLR3, PPY, ERBB2, ESM1, HGF, RET
#cancers: colorectal, kidney, ovarian, and liver cancer

rm(list=ls(all=TRUE)) 

#load packages 
pacman::p_load(devtools, ggplot2, ggrepel, igraph, gridExtra, plyr, dplyr, TwoSampleMR, writexl, tidyverse, plyr, dplyr, ieugwasr, ggplot2, ggforestplot)
               
#legacy version 2smr
devtools::install_github("MRCIEU/TwoSampleMR@0.4.26")

#Sun UKBB pQTLs, downloaded from https://www.biorxiv.org/content/10.1101/2022.06.17.496443v1.supplementary-material
ukbb_pqtls <- readxl::read_xlsx(("Sun2022_media-2.xlsx"), sheet = 11, skip = 2)

#filter intervention associated proteins
prt_list <- c("GPNMB", "FURIN", "WIF1", "TLR3", "PPY", "ERBB2", "HGF", "ESM1", "RET")
ukbb_instruments <- filter(ukbb_pqtls, str_detect(ukbb_pqtls$`Assay Target`, paste(prt_list, collapse="|")))

#remove RETN
ukbb_instruments <- ukbb_instruments[!grepl("RETN", ukbb_instruments$`Assay Target`),]

#drop trans SNPs
ukbb_instruments <- ukbb_instruments[!grepl("trans", ukbb_instruments$`cis/trans`),]

#extract effect and other alleles
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
ukbb_instruments$A0 <- substrRight((ukbb_instruments$`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`),10)
ukbb_instruments$other_allele.exposure <- substr(ukbb_instruments$A0, 1,1) 
ukbb_instruments$A1 <- substrRight((ukbb_instruments$`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`),8)
ukbb_instruments$effect_allele.exposure <- substr(ukbb_instruments$A1, 1,1) 
ukbb_instruments$beta.exposure <- ukbb_instruments$`BETA (discovery, wrt. A1)`
ukbb_instruments$se.exposure <- ukbb_instruments$`SE (discovery)`
ukbb_instruments$eaf.exposure <- ukbb_instruments$`A1FREQ (discovery)`
ukbb_instruments$SNP <- ukbb_instruments$rsID
ukbb_instruments$exposure <- ukbb_instruments$`Assay Target`
ukbb_instruments$id.exposure <- ukbb_instruments$`Assay Target`
ukbb_instruments$samplesize.exposure <- 10708
ukbb_instruments$pval.exposure <- 10^-(ukbb_instruments$`log10(p) (discovery)`)

#compute f stats for pQTLs
ukbb_instruments$t_stat <- (ukbb_instruments$beta.exposure /ukbb_instruments$se.exposure )
ukbb_instruments$f_stat <- (ukbb_instruments$t_stat)^2

#compute r2
ukbb_instruments$num <- 2*(ukbb_instruments$beta.exposure^2)*ukbb_instruments$eaf.exposure*(1-ukbb_instruments$eaf.exposure)
ukbb_instruments$den <- 2*(ukbb_instruments$beta.exposure^2)*ukbb_instruments$eaf.exposure*(1-ukbb_instruments$eaf.exposure) + (ukbb_instruments$se.exposure^2)*2*ukbb_instruments$samplesize.exposure*ukbb_instruments$eaf.exposure*(1-ukbb_instruments$eaf.exposure)
ukbb_instruments$pve <- ukbb_instruments$num/ukbb_instruments$den

#extract outcome data - ieu open gwas summary stats
#Endometrial
#cancer_ids <- c("ebi-a-GCST006464")
#ukbb_instruments<-c("rs3764354","rs4242051","rs6227","rs75801644","rs5745687","rs2795507","rs3775291","rs462010")
#endometrial <- extract_outcome_data(
#  snps = ukbb_instruments,
#  outcomes = cancer_ids,
#  proxies = TRUE,
#  rsq = 0.8,
#  align_alleles = 1,
#  palindromes = 1,
#  maf_threshold = 0.01,
#  access_token = ieugwasr::check_access_token(),
#  splitsize = 10000,
#  proxy_splitsize = 500
#)
outfilepath<-"endometrial_snps.txt"
endometrial<-outcome_dat <- read_outcome_data(
  snps = ukbb_instruments$rsID,
  filename = outfilepath,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "eaf.outcome"
)
endometrial$outcome<-"Endometrial_cancer"
endometrial<-endometrial[,-c(2,3)]


#summary stats
#crc
outfilepath<-"overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
crc<-outcome_dat <- read_outcome_data(
  snps = ukbb_instruments$rsID,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
crc$outcome<-"Colorectal_cancer"

#breast
outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
breast <- read_outcome_data(
  snps = ukbb_instruments$rsID,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
breast$outcome<-"Breast_cancer"

#Looking for proxy for WIF1 in breast cancer dataset
proxy <- LDlinkR::LDproxy("rs462010", pop = "EUR", r2d = "r2", token = "a563a6f57fcc")
proxy <- proxy %>% dplyr::filter(R2 > 0.8)

outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
wif1_breast <- read_outcome_data(
  snps = proxy$RS_Number,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
wif1_breast$outcome<-"wif1_Breast_cancer"
wif1_breast$SNP<-"rs462010"
wif1_breast$effect_allele.outcome<-"A"
wif1_breast$other_allele.outcome<-"G"

breast<-rbind(breast,wif1_breast)

#liver
outfilepath<-"finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC"
liver <- read_outcome_data(
  snps = ukbb_instruments$rsID,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
liver$outcome<-"Liver_cancer"
liver<-liver[,-1]

#gallbladder
outfilepath<-"finngen_R8_C3_GALLBLADDER_EXALLC"
gallbladder <- read_outcome_data(
  snps = ukbb_instruments$rsID,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
gallbladder$outcome<-"Gallbladder_cancer"
gallbladder<-gallbladder[,-1]

#pancreas
outfilepath<-"finngen_R8_C3_PANCREAS_EXALLC"
pancreas <- read_outcome_data(
  snps = ukbb_instruments$rsID,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
pancreas$outcome<-"Pancreas_cancer"
pancreas<-pancreas[,-1]

cancer_all<-rbind(endometrial,crc,breast,liver,gallbladder,pancreas)

##run MR
#ieu open gwas summary stats 
ieu_harm_dat <- harmonise_data(ukbb_instruments, cancer_all, action = 2)
ieu_mr_results <- mr(ieu_harm_dat, method_list=c("mr_wald_ratio"))

ieu_mr_results$OR<-exp(ieu_mr_results$b)
ieu_mr_results$LCI<-exp(ieu_mr_results$b-1.96*ieu_mr_results$se)
ieu_mr_results$UCI<-exp(ieu_mr_results$b+1.96*ieu_mr_results$se)
ieu_mr_results$col<-paste(round(ieu_mr_results$OR,digits=2)," (",round(ieu_mr_results$LCI,digits=2),"-",round(ieu_mr_results$UCI,digits=2),")",sep="")

write.table(ieu_mr_results, row.names = FALSE, col.names = TRUE, quote = F, sep = '\t', file = "pqtl_mr_results.txt")


# Figure 2 ----------------------------------------------------------------
df<-data.table::fread("pqtl_mr_results.txt")
df$outcome<-gsub("_"," ",df$outcome)
df$outcome<-gsub("wif1 ","",df$outcome)

p1<- ggforestplot::forestplot(
  df = df,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 1,
  xlab = "Odds ratio (95% CI) cancer per 1-SD protein",
  logodds=TRUE,
  se=se,
  colour = outcome
)+
  scale_color_manual(values=rev(c("#D88E00","#59A4E7","#2D8F5E","#235DA5","#C24900","#BA6197")))+ guides(colour=guide_legend(title="",reverse=TRUE))

p1

ggsave(filename="Figure 2.png", plot=last_plot(),width = 220, height = 200, units = "mm",bg = "white")



# Colocalization analysis -------------------------------------------------

rm(list=ls(all=TRUE)) 

#use coloc susie package - basic coloc
# if(!require("remotes"))
#   install.packages("remotes") # if necessary
# library(remotes)
# install_github("chr1swallace/coloc",build_vignettes=TRUE)

# devtools::install_github("mrcieu/gwasglue")

#load coloc package
library(coloc)
library(dplyr)
library(gwasglue)
library(data.table)
library(TwoSampleMR)

##format protein summary stats - Pietzner et al https://www.science.org/doi/10.1126/science.abj1541
#read in protein and cancer summary stats +/- 10Mb around pqtls

#ESM1
headers <- c("snp", "MarkerName", "Allele1",  "Allele2",  "Freq1",  "FreqSE", "MinFreq",  "MaxFreq",  "beta", "se", "Pvalue", "Direction",  "HetISq", "HetChiSq", "HetDf",  "HetPVal",  "TotalSampleSize",  "chr", "position")
esm1 <- read.delim("esm1_coloc_10mb.txt", header = FALSE, sep = "\t")
names(esm1) <- headers

#FURIN
furin <- read.delim("furin_coloc_10mb.txt", header = FALSE, sep = "\t")
names(furin) <- headers

#GPNMB
gpnmb <- read.delim("gpnmb_coloc_10mb.txt", header = FALSE, sep = "\t")
names(gpnmb) <- headers

#HGF
hgf <- read.delim("hgf_coloc_10mb.txt", header = FALSE, sep = "\t")
names(hgf) <- headers

#RET
ret <- read.delim("ret_coloc_10mb.txt", header = FALSE, sep = "\t")
names(ret) <- headers

#TLR3
tlr3 <- read.delim("tlr3_coloc_10mb.txt", header = FALSE, sep = "\t")
names(tlr3) <- headers

#WIF1
wif1 <- read.delim("wif1_coloc_10mb.txt", header = FALSE, sep = "\t")
names(wif1) <- headers

#Outcome datasets
#endometrial
outfilepath<-"GCST006464.txt"
endometrial<-outcome_dat <- read_outcome_data(
  snps = esm1$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  eaf_col = "hm_effect_allele_frequency"
)

endometrial$samplesize.outcome<-121885

#crc
outfilepath<-"overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
crc<-outcome_dat <- read_outcome_data(
  snps = esm1$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
crc$samplesize.outcome<-98715
  
#breast
outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
breast <- read_outcome_data(
  snps = esm1$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
breast$samplesize.outcome<-247173

#liver
outfilepath<-"finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC"
liver <- read_outcome_data(
  snps = esm1$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
liver$samplesize.outcome<-260231
  
#gallbladder
outfilepath<-"finngen_R8_C3_GALLBLADDER_EXALLC"
gallbladder <- read_outcome_data(
  snps = esm1$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
gallbladder$samplesize.outcome<-259668
  
#pancreas
outfilepath<-"finngen_R8_C3_PANCREAS_EXALLC"
pancreas <- read_outcome_data(
  snps = esm1$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
pancreas$samplesize.outcome<-260832
  
#ESM1
esm1<-format_data(esm1,type="exposure",
                  snp_col="snp",
                  effect_allele_col = "Allele1",
                  other_allele_col = "Allele2",
                  eaf_col = "Freq1",
                  beta_col = "beta",
                  se_col = "se",
                  samplesize_col = "TotalSampleSize")

esm1_ec<-harmonise_data(esm1,endometrial)
df1 <- list(beta=esm1_ec$beta.exposure, varbeta = esm1_ec$se.exposure^2, MAF=esm1_ec$eaf.exposure, type = "quant", N = esm1_ec$samplesize.exposure,sdY=1)
df2 <- list(beta=esm1_ec$beta.outcome, varbeta = esm1_ec$se.outcome^2, MAF=esm1_ec$eaf.outcome, type = "cc", N = esm1_ec$samplesize.outcome,sdY=1)
esm1_ec_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

esm1_crc<-harmonise_data(esm1,crc)
df1 <- list(beta=esm1_crc$beta.exposure, varbeta = esm1_crc$se.exposure^2, MAF=esm1_crc$eaf.exposure, type = "quant", N = esm1_crc$samplesize.exposure,sdY=1)
df2 <- list(beta=esm1_crc$beta.outcome, varbeta = esm1_crc$se.outcome^2, MAF=esm1_crc$eaf.outcome, type = "cc", N = esm1_crc$samplesize.outcome,sdY=1)
esm1_crc_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

esm1_breast<-harmonise_data(esm1,breast)
df1 <- list(beta=esm1_breast$beta.exposure, varbeta = esm1_breast$se.exposure^2, MAF=esm1_breast$eaf.exposure, type = "quant", N = esm1_breast$samplesize.exposure,sdY=1)
df2 <- list(beta=esm1_breast$beta.outcome, varbeta = esm1_breast$se.outcome^2, MAF=esm1_breast$eaf.outcome, type = "cc", N = esm1_breast$samplesize.outcome,sdY=1)
esm1_breast_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

esm1_liver<-harmonise_data(esm1,liver)
df1 <- list(beta=esm1_liver$beta.exposure, varbeta = esm1_liver$se.exposure^2, MAF=esm1_liver$eaf.exposure, type = "quant", N = esm1_liver$samplesize.exposure,sdY=1)
df2 <- list(beta=esm1_liver$beta.outcome, varbeta = esm1_liver$se.outcome^2, MAF=esm1_liver$eaf.outcome, type = "cc", N = esm1_liver$samplesize.outcome,sdY=1)
esm1_liver_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

esm1_gallbladder<-harmonise_data(esm1,gallbladder)
df1 <- list(beta=esm1_gallbladder$beta.exposure, varbeta = esm1_gallbladder$se.exposure^2, MAF=esm1_gallbladder$eaf.exposure, type = "quant", N = esm1_gallbladder$samplesize.exposure,sdY=1)
df2 <- list(beta=esm1_gallbladder$beta.outcome, varbeta = esm1_gallbladder$se.outcome^2, MAF=esm1_gallbladder$eaf.outcome, type = "cc", N = esm1_gallbladder$samplesize.outcome,sdY=1)
esm1_gallbladder_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

esm1_pancreas<-harmonise_data(esm1,pancreas)
df1 <- list(beta=esm1_pancreas$beta.exposure, varbeta = esm1_pancreas$se.exposure^2, MAF=esm1_pancreas$eaf.exposure, type = "quant", N = esm1_pancreas$samplesize.exposure,sdY=1)
df2 <- list(beta=esm1_pancreas$beta.outcome, varbeta = esm1_pancreas$se.outcome^2, MAF=esm1_pancreas$eaf.outcome, type = "cc", N = esm1_pancreas$samplesize.outcome,sdY=1)
esm1_pancreas_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

#furin
furin<-format_data(furin,type="exposure",
                   snp_col="snp",
                   effect_allele_col = "Allele1",
                   other_allele_col = "Allele2",
                   eaf_col = "Freq1",
                   beta_col = "beta",
                   se_col = "se",
                   samplesize_col = "TotalSampleSize")

#Outcome datasets
#endometrial
outfilepath<-"GCST006464.txt"
endometrial<-outcome_dat <- read_outcome_data(
  snps = furin$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  eaf_col = "hm_effect_allele_frequency"
)

endometrial$samplesize.outcome<-121885

#crc
outfilepath<-"overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
crc<-outcome_dat <- read_outcome_data(
  snps = furin$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
crc$samplesize.outcome<-98715

#breast
outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
breast <- read_outcome_data(
  snps = furin$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
breast$samplesize.outcome<-247173

#liver
outfilepath<-"finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC"
liver <- read_outcome_data(
  snps = furin$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
liver$samplesize.outcome<-260231

#gallbladder
outfilepath<-"finngen_R8_C3_GALLBLADDER_EXALLC"
gallbladder <- read_outcome_data(
  snps = furin$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
gallbladder$samplesize.outcome<-259668

#pancreas
outfilepath<-"finngen_R8_C3_PANCREAS_EXALLC"
pancreas <- read_outcome_data(
  snps = furin$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
pancreas$samplesize.outcome<-260832


furin_ec<-harmonise_data(furin,endometrial)
df1 <- list(beta=furin_ec$beta.exposure, varbeta = furin_ec$se.exposure^2, MAF=furin_ec$eaf.exposure, type = "quant", N = furin_ec$samplesize.exposure,sdY=1)
df2 <- list(beta=furin_ec$beta.outcome, varbeta = furin_ec$se.outcome^2, MAF=furin_ec$eaf.outcome, type = "cc", N = furin_ec$samplesize.outcome,sdY=1)
furin_ec_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

furin_crc<-harmonise_data(furin,crc)
df1 <- list(beta=furin_crc$beta.exposure, varbeta = furin_crc$se.exposure^2, MAF=furin_crc$eaf.exposure, type = "quant", N = furin_crc$samplesize.exposure,sdY=1)
df2 <- list(beta=furin_crc$beta.outcome, varbeta = furin_crc$se.outcome^2, MAF=furin_crc$eaf.outcome, type = "cc", N = furin_crc$samplesize.outcome,sdY=1)
furin_crc_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

furin_breast<-harmonise_data(furin,breast)
df1 <- list(beta=furin_breast$beta.exposure, varbeta = furin_breast$se.exposure^2, MAF=furin_breast$eaf.exposure, type = "quant", N = furin_breast$samplesize.exposure,sdY=1)
df2 <- list(beta=furin_breast$beta.outcome, varbeta = furin_breast$se.outcome^2, MAF=furin_breast$eaf.outcome, type = "cc", N = furin_breast$samplesize.outcome,sdY=1)
furin_breast_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

furin_liver<-harmonise_data(furin,liver)
df1 <- list(beta=furin_liver$beta.exposure, varbeta = furin_liver$se.exposure^2, MAF=furin_liver$eaf.exposure, type = "quant", N = furin_liver$samplesize.exposure,sdY=1)
df2 <- list(beta=furin_liver$beta.outcome, varbeta = furin_liver$se.outcome^2, MAF=furin_liver$eaf.outcome, type = "cc", N = furin_liver$samplesize.outcome,sdY=1)
furin_liver_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

furin_gallbladder<-harmonise_data(furin,gallbladder)
df1 <- list(beta=furin_gallbladder$beta.exposure, varbeta = furin_gallbladder$se.exposure^2, MAF=furin_gallbladder$eaf.exposure, type = "quant", N = furin_gallbladder$samplesize.exposure,sdY=1)
df2 <- list(beta=furin_gallbladder$beta.outcome, varbeta = furin_gallbladder$se.outcome^2, MAF=furin_gallbladder$eaf.outcome, type = "cc", N = furin_gallbladder$samplesize.outcome,sdY=1)
furin_gallbladder_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

furin_pancreas<-harmonise_data(furin,pancreas)
df1 <- list(beta=furin_pancreas$beta.exposure, varbeta = furin_pancreas$se.exposure^2, MAF=furin_pancreas$eaf.exposure, type = "quant", N = furin_pancreas$samplesize.exposure,sdY=1)
df2 <- list(beta=furin_pancreas$beta.outcome, varbeta = furin_pancreas$se.outcome^2, MAF=furin_pancreas$eaf.outcome, type = "cc", N = furin_pancreas$samplesize.outcome,sdY=1)
furin_pancreas_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

#gpnmb
gpnmb<-format_data(gpnmb,type="exposure",
                   snp_col="snp",
                   effect_allele_col = "Allele1",
                   other_allele_col = "Allele2",
                   eaf_col = "Freq1",
                   beta_col = "beta",
                   se_col = "se",
                   samplesize_col = "TotalSampleSize")

#Outcome datasets
#endometrial
outfilepath<-"GCST006464.txt"
endometrial<-outcome_dat <- read_outcome_data(
  snps = gpnmb$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  eaf_col = "hm_effect_allele_frequency"
)

endometrial$samplesize.outcome<-121885

#crc
outfilepath<-"overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
crc<-outcome_dat <- read_outcome_data(
  snps = gpnmb$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
crc$samplesize.outcome<-98715

#breast
outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
breast <- read_outcome_data(
  snps = gpnmb$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
breast$samplesize.outcome<-247173

#liver
outfilepath<-"finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC"
liver <- read_outcome_data(
  snps = gpnmb$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
liver$samplesize.outcome<-260231

#gallbladder
outfilepath<-"finngen_R8_C3_GALLBLADDER_EXALLC"
gallbladder <- read_outcome_data(
  snps = gpnmb$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
gallbladder$samplesize.outcome<-259668

#pancreas
outfilepath<-"finngen_R8_C3_PANCREAS_EXALLC"
pancreas <- read_outcome_data(
  snps = gpnmb$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
pancreas$samplesize.outcome<-260832


gpnmb_ec<-harmonise_data(gpnmb,endometrial)
df1 <- list(beta=gpnmb_ec$beta.exposure, varbeta = gpnmb_ec$se.exposure^2, MAF=gpnmb_ec$eaf.exposure, type = "quant", N = gpnmb_ec$samplesize.exposure,sdY=1)
df2 <- list(beta=gpnmb_ec$beta.outcome, varbeta = gpnmb_ec$se.outcome^2, MAF=gpnmb_ec$eaf.outcome, type = "cc", N = gpnmb_ec$samplesize.outcome,sdY=1)
gpnmb_ec_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

gpnmb_crc<-harmonise_data(gpnmb,crc)
df1 <- list(beta=gpnmb_crc$beta.exposure, varbeta = gpnmb_crc$se.exposure^2, MAF=gpnmb_crc$eaf.exposure, type = "quant", N = gpnmb_crc$samplesize.exposure,sdY=1)
df2 <- list(beta=gpnmb_crc$beta.outcome, varbeta = gpnmb_crc$se.outcome^2, MAF=gpnmb_crc$eaf.outcome, type = "cc", N = gpnmb_crc$samplesize.outcome,sdY=1)
gpnmb_crc_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

gpnmb_breast<-harmonise_data(gpnmb,breast)
df1 <- list(beta=gpnmb_breast$beta.exposure, varbeta = gpnmb_breast$se.exposure^2, MAF=gpnmb_breast$eaf.exposure, type = "quant", N = gpnmb_breast$samplesize.exposure,sdY=1)
df2 <- list(beta=gpnmb_breast$beta.outcome, varbeta = gpnmb_breast$se.outcome^2, MAF=gpnmb_breast$eaf.outcome, type = "cc", N = gpnmb_breast$samplesize.outcome,sdY=1)
gpnmb_breast_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

gpnmb_liver<-harmonise_data(gpnmb,liver)
df1 <- list(beta=gpnmb_liver$beta.exposure, varbeta = gpnmb_liver$se.exposure^2, MAF=gpnmb_liver$eaf.exposure, type = "quant", N = gpnmb_liver$samplesize.exposure,sdY=1)
df2 <- list(beta=gpnmb_liver$beta.outcome, varbeta = gpnmb_liver$se.outcome^2, MAF=gpnmb_liver$eaf.outcome, type = "cc", N = gpnmb_liver$samplesize.outcome,sdY=1)
gpnmb_liver_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

gpnmb_gallbladder<-harmonise_data(gpnmb,gallbladder)
df1 <- list(beta=gpnmb_gallbladder$beta.exposure, varbeta = gpnmb_gallbladder$se.exposure^2, MAF=gpnmb_gallbladder$eaf.exposure, type = "quant", N = gpnmb_gallbladder$samplesize.exposure,sdY=1)
df2 <- list(beta=gpnmb_gallbladder$beta.outcome, varbeta = gpnmb_gallbladder$se.outcome^2, MAF=gpnmb_gallbladder$eaf.outcome, type = "cc", N = gpnmb_gallbladder$samplesize.outcome,sdY=1)
gpnmb_gallbladder_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

gpnmb_pancreas<-harmonise_data(gpnmb,pancreas)
df1 <- list(beta=gpnmb_pancreas$beta.exposure, varbeta = gpnmb_pancreas$se.exposure^2, MAF=gpnmb_pancreas$eaf.exposure, type = "quant", N = gpnmb_pancreas$samplesize.exposure,sdY=1)
df2 <- list(beta=gpnmb_pancreas$beta.outcome, varbeta = gpnmb_pancreas$se.outcome^2, MAF=gpnmb_pancreas$eaf.outcome, type = "cc", N = gpnmb_pancreas$samplesize.outcome,sdY=1)
gpnmb_pancreas_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

#hgf
hgf<-format_data(hgf,type="exposure",
                 snp_col="snp",
                 effect_allele_col = "Allele1",
                 other_allele_col = "Allele2",
                 eaf_col = "Freq1",
                 beta_col = "beta",
                 se_col = "se",
                 samplesize_col = "TotalSampleSize")

#Outcome datasets
#endometrial
outfilepath<-"GCST006464.txt"
endometrial<-outcome_dat <- read_outcome_data(
  snps = hgf$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  eaf_col = "hm_effect_allele_frequency"
)

endometrial$samplesize.outcome<-121885

#crc
outfilepath<-"overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
crc<-outcome_dat <- read_outcome_data(
  snps = hgf$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
crc$samplesize.outcome<-98715

#breast
outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
breast <- read_outcome_data(
  snps = hgf$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
breast$samplesize.outcome<-247173

#liver
outfilepath<-"finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC"
liver <- read_outcome_data(
  snps = hgf$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
liver$samplesize.outcome<-260231

#gallbladder
outfilepath<-"finngen_R8_C3_GALLBLADDER_EXALLC"
gallbladder <- read_outcome_data(
  snps = hgf$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
gallbladder$samplesize.outcome<-259668

#pancreas
outfilepath<-"finngen_R8_C3_PANCREAS_EXALLC"
pancreas <- read_outcome_data(
  snps = hgf$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
pancreas$samplesize.outcome<-260832


hgf_ec<-harmonise_data(hgf,endometrial)
df1 <- list(beta=hgf_ec$beta.exposure, varbeta = hgf_ec$se.exposure^2, MAF=hgf_ec$eaf.exposure, type = "quant", N = hgf_ec$samplesize.exposure,sdY=1)
df2 <- list(beta=hgf_ec$beta.outcome, varbeta = hgf_ec$se.outcome^2, MAF=hgf_ec$eaf.outcome, type = "cc", N = hgf_ec$samplesize.outcome,sdY=1)
hgf_ec_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

hgf_crc<-harmonise_data(hgf,crc)
df1 <- list(beta=hgf_crc$beta.exposure, varbeta = hgf_crc$se.exposure^2, MAF=hgf_crc$eaf.exposure, type = "quant", N = hgf_crc$samplesize.exposure,sdY=1)
df2 <- list(beta=hgf_crc$beta.outcome, varbeta = hgf_crc$se.outcome^2, MAF=hgf_crc$eaf.outcome, type = "cc", N = hgf_crc$samplesize.outcome,sdY=1)
hgf_crc_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

hgf_breast<-harmonise_data(hgf,breast)
df1 <- list(beta=hgf_breast$beta.exposure, varbeta = hgf_breast$se.exposure^2, MAF=hgf_breast$eaf.exposure, type = "quant", N = hgf_breast$samplesize.exposure,sdY=1)
df2 <- list(beta=hgf_breast$beta.outcome, varbeta = hgf_breast$se.outcome^2, MAF=hgf_breast$eaf.outcome, type = "cc", N = hgf_breast$samplesize.outcome,sdY=1)
hgf_breast_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

hgf_liver<-harmonise_data(hgf,liver)
df1 <- list(beta=hgf_liver$beta.exposure, varbeta = hgf_liver$se.exposure^2, MAF=hgf_liver$eaf.exposure, type = "quant", N = hgf_liver$samplesize.exposure,sdY=1)
df2 <- list(beta=hgf_liver$beta.outcome, varbeta = hgf_liver$se.outcome^2, MAF=hgf_liver$eaf.outcome, type = "cc", N = hgf_liver$samplesize.outcome,sdY=1)
hgf_liver_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

hgf_gallbladder<-harmonise_data(hgf,gallbladder)
df1 <- list(beta=hgf_gallbladder$beta.exposure, varbeta = hgf_gallbladder$se.exposure^2, MAF=hgf_gallbladder$eaf.exposure, type = "quant", N = hgf_gallbladder$samplesize.exposure,sdY=1)
df2 <- list(beta=hgf_gallbladder$beta.outcome, varbeta = hgf_gallbladder$se.outcome^2, MAF=hgf_gallbladder$eaf.outcome, type = "cc", N = hgf_gallbladder$samplesize.outcome,sdY=1)
hgf_gallbladder_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

hgf_pancreas<-harmonise_data(hgf,pancreas)
df1 <- list(beta=hgf_pancreas$beta.exposure, varbeta = hgf_pancreas$se.exposure^2, MAF=hgf_pancreas$eaf.exposure, type = "quant", N = hgf_pancreas$samplesize.exposure,sdY=1)
df2 <- list(beta=hgf_pancreas$beta.outcome, varbeta = hgf_pancreas$se.outcome^2, MAF=hgf_pancreas$eaf.outcome, type = "cc", N = hgf_pancreas$samplesize.outcome,sdY=1)
hgf_pancreas_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

#ret
ret<-format_data(ret,type="exposure",
                 snp_col="snp",
                 effect_allele_col = "Allele1",
                 other_allele_col = "Allele2",
                 eaf_col = "Freq1",
                 beta_col = "beta",
                 se_col = "se",
                 samplesize_col = "TotalSampleSize")

#Outcome datasets
#endometrial
outfilepath<-"GCST006464.txt"
endometrial<-outcome_dat <- read_outcome_data(
  snps = ret$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  eaf_col = "hm_effect_allele_frequency"
)

endometrial$samplesize.outcome<-121885

#crc
outfilepath<-"overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
crc<-outcome_dat <- read_outcome_data(
  snps = ret$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
crc$samplesize.outcome<-98715

#breast
outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
breast <- read_outcome_data(
  snps = ret$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
breast$samplesize.outcome<-247173

#liver
outfilepath<-"finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC"
liver <- read_outcome_data(
  snps = ret$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
liver$samplesize.outcome<-260231

#gallbladder
outfilepath<-"finngen_R8_C3_GALLBLADDER_EXALLC"
gallbladder <- read_outcome_data(
  snps = ret$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
gallbladder$samplesize.outcome<-259668

#pancreas
outfilepath<-"finngen_R8_C3_PANCREAS_EXALLC"
pancreas <- read_outcome_data(
  snps = ret$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
pancreas$samplesize.outcome<-260832


ret_ec<-harmonise_data(ret,endometrial)
df1 <- list(beta=ret_ec$beta.exposure, varbeta = ret_ec$se.exposure^2, MAF=ret_ec$eaf.exposure, type = "quant", N = ret_ec$samplesize.exposure,sdY=1)
df2 <- list(beta=ret_ec$beta.outcome, varbeta = ret_ec$se.outcome^2, MAF=ret_ec$eaf.outcome, type = "cc", N = ret_ec$samplesize.outcome,sdY=1)
ret_ec_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

ret_crc<-harmonise_data(ret,crc)
df1 <- list(beta=ret_crc$beta.exposure, varbeta = ret_crc$se.exposure^2, MAF=ret_crc$eaf.exposure, type = "quant", N = ret_crc$samplesize.exposure,sdY=1)
df2 <- list(beta=ret_crc$beta.outcome, varbeta = ret_crc$se.outcome^2, MAF=ret_crc$eaf.outcome, type = "cc", N = ret_crc$samplesize.outcome,sdY=1)
ret_crc_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

ret_breast<-harmonise_data(ret,breast)
df1 <- list(beta=ret_breast$beta.exposure, varbeta = ret_breast$se.exposure^2, MAF=ret_breast$eaf.exposure, type = "quant", N = ret_breast$samplesize.exposure,sdY=1)
df2 <- list(beta=ret_breast$beta.outcome, varbeta = ret_breast$se.outcome^2, MAF=ret_breast$eaf.outcome, type = "cc", N = ret_breast$samplesize.outcome,sdY=1)
ret_breast_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

ret_liver<-harmonise_data(ret,liver)
df1 <- list(beta=ret_liver$beta.exposure, varbeta = ret_liver$se.exposure^2, MAF=ret_liver$eaf.exposure, type = "quant", N = ret_liver$samplesize.exposure,sdY=1)
df2 <- list(beta=ret_liver$beta.outcome, varbeta = ret_liver$se.outcome^2, MAF=ret_liver$eaf.outcome, type = "cc", N = ret_liver$samplesize.outcome,sdY=1)
ret_liver_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

ret_gallbladder<-harmonise_data(ret,gallbladder)
df1 <- list(beta=ret_gallbladder$beta.exposure, varbeta = ret_gallbladder$se.exposure^2, MAF=ret_gallbladder$eaf.exposure, type = "quant", N = ret_gallbladder$samplesize.exposure,sdY=1)
df2 <- list(beta=ret_gallbladder$beta.outcome, varbeta = ret_gallbladder$se.outcome^2, MAF=ret_gallbladder$eaf.outcome, type = "cc", N = ret_gallbladder$samplesize.outcome,sdY=1)
ret_gallbladder_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

ret_pancreas<-harmonise_data(ret,pancreas)
df1 <- list(beta=ret_pancreas$beta.exposure, varbeta = ret_pancreas$se.exposure^2, MAF=ret_pancreas$eaf.exposure, type = "quant", N = ret_pancreas$samplesize.exposure,sdY=1)
df2 <- list(beta=ret_pancreas$beta.outcome, varbeta = ret_pancreas$se.outcome^2, MAF=ret_pancreas$eaf.outcome, type = "cc", N = ret_pancreas$samplesize.outcome,sdY=1)
ret_pancreas_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

#tlr3
tlr3<-format_data(tlr3,type="exposure",
                  snp_col="snp",
                  effect_allele_col = "Allele1",
                  other_allele_col = "Allele2",
                  eaf_col = "Freq1",
                  beta_col = "beta",
                  se_col = "se",
                  samplesize_col = "TotalSampleSize")

#Outcome datasets
#endometrial
outfilepath<-"GCST006464.txt"
endometrial<-outcome_dat <- read_outcome_data(
  snps = tlr3$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  eaf_col = "hm_effect_allele_frequency"
)

endometrial$samplesize.outcome<-121885

#crc
outfilepath<-"overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
crc<-outcome_dat <- read_outcome_data(
  snps = tlr3$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
crc$samplesize.outcome<-98715

#breast
outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
breast <- read_outcome_data(
  snps = tlr3$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
breast$samplesize.outcome<-247173

#liver
outfilepath<-"finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC"
liver <- read_outcome_data(
  snps = tlr3$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
liver$samplesize.outcome<-260231

#gallbladder
outfilepath<-"finngen_R8_C3_GALLBLADDER_EXALLC"
gallbladder <- read_outcome_data(
  snps = tlr3$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
gallbladder$samplesize.outcome<-259668

#pancreas
outfilepath<-"finngen_R8_C3_PANCREAS_EXALLC"
pancreas <- read_outcome_data(
  snps = tlr3$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
pancreas$samplesize.outcome<-260832


tlr3_ec<-harmonise_data(tlr3,endometrial)
df1 <- list(beta=tlr3_ec$beta.exposure, varbeta = tlr3_ec$se.exposure^2, MAF=tlr3_ec$eaf.exposure, type = "quant", N = tlr3_ec$samplesize.exposure,sdY=1)
df2 <- list(beta=tlr3_ec$beta.outcome, varbeta = tlr3_ec$se.outcome^2, MAF=tlr3_ec$eaf.outcome, type = "cc", N = tlr3_ec$samplesize.outcome,sdY=1)
tlr3_ec_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

tlr3_crc<-harmonise_data(tlr3,crc)
df1 <- list(beta=tlr3_crc$beta.exposure, varbeta = tlr3_crc$se.exposure^2, MAF=tlr3_crc$eaf.exposure, type = "quant", N = tlr3_crc$samplesize.exposure,sdY=1)
df2 <- list(beta=tlr3_crc$beta.outcome, varbeta = tlr3_crc$se.outcome^2, MAF=tlr3_crc$eaf.outcome, type = "cc", N = tlr3_crc$samplesize.outcome,sdY=1)
tlr3_crc_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

tlr3_breast<-harmonise_data(tlr3,breast)
df1 <- list(beta=tlr3_breast$beta.exposure, varbeta = tlr3_breast$se.exposure^2, MAF=tlr3_breast$eaf.exposure, type = "quant", N = tlr3_breast$samplesize.exposure,sdY=1)
df2 <- list(beta=tlr3_breast$beta.outcome, varbeta = tlr3_breast$se.outcome^2, MAF=tlr3_breast$eaf.outcome, type = "cc", N = tlr3_breast$samplesize.outcome,sdY=1)
tlr3_breast_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

tlr3_liver<-harmonise_data(tlr3,liver)
df1 <- list(beta=tlr3_liver$beta.exposure, varbeta = tlr3_liver$se.exposure^2, MAF=tlr3_liver$eaf.exposure, type = "quant", N = tlr3_liver$samplesize.exposure,sdY=1)
df2 <- list(beta=tlr3_liver$beta.outcome, varbeta = tlr3_liver$se.outcome^2, MAF=tlr3_liver$eaf.outcome, type = "cc", N = tlr3_liver$samplesize.outcome,sdY=1)
tlr3_liver_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

tlr3_gallbladder<-harmonise_data(tlr3,gallbladder)
df1 <- list(beta=tlr3_gallbladder$beta.exposure, varbeta = tlr3_gallbladder$se.exposure^2, MAF=tlr3_gallbladder$eaf.exposure, type = "quant", N = tlr3_gallbladder$samplesize.exposure,sdY=1)
df2 <- list(beta=tlr3_gallbladder$beta.outcome, varbeta = tlr3_gallbladder$se.outcome^2, MAF=tlr3_gallbladder$eaf.outcome, type = "cc", N = tlr3_gallbladder$samplesize.outcome,sdY=1)
tlr3_gallbladder_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

tlr3_pancreas<-harmonise_data(tlr3,pancreas)
df1 <- list(beta=tlr3_pancreas$beta.exposure, varbeta = tlr3_pancreas$se.exposure^2, MAF=tlr3_pancreas$eaf.exposure, type = "quant", N = tlr3_pancreas$samplesize.exposure,sdY=1)
df2 <- list(beta=tlr3_pancreas$beta.outcome, varbeta = tlr3_pancreas$se.outcome^2, MAF=tlr3_pancreas$eaf.outcome, type = "cc", N = tlr3_pancreas$samplesize.outcome,sdY=1)
tlr3_pancreas_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

#wif1
wif1<-format_data(wif1,type="exposure",
                  snp_col="snp",
                  effect_allele_col = "Allele1",
                  other_allele_col = "Allele2",
                  eaf_col = "Freq1",
                  beta_col = "beta",
                  se_col = "se",
                  samplesize_col = "TotalSampleSize")

#Outcome datasets
#endometrial
outfilepath<-"GCST006464.txt"
endometrial<-outcome_dat <- read_outcome_data(
  snps = wif1$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  eaf_col = "hm_effect_allele_frequency"
)

endometrial$samplesize.outcome<-121885

#crc
outfilepath<-"overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
crc<-outcome_dat <- read_outcome_data(
  snps = wif1$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
crc$samplesize.outcome<-98715

#breast
outfilepath<-"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
breast <- read_outcome_data(
  snps = wif1$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
breast$samplesize.outcome<-247173

#liver
outfilepath<-"finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC"
liver <- read_outcome_data(
  snps = wif1$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
liver$samplesize.outcome<-260231

#gallbladder
outfilepath<-"finngen_R8_C3_GALLBLADDER_EXALLC"
gallbladder <- read_outcome_data(
  snps = wif1$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
gallbladder$samplesize.outcome<-259668

#pancreas
outfilepath<-"finngen_R8_C3_PANCREAS_EXALLC"
pancreas <- read_outcome_data(
  snps = wif1$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
pancreas$samplesize.outcome<-260832


wif1_ec<-harmonise_data(wif1,endometrial)
df1 <- list(beta=wif1_ec$beta.exposure, varbeta = wif1_ec$se.exposure^2, MAF=wif1_ec$eaf.exposure, type = "quant", N = wif1_ec$samplesize.exposure,sdY=1)
df2 <- list(beta=wif1_ec$beta.outcome, varbeta = wif1_ec$se.outcome^2, MAF=wif1_ec$eaf.outcome, type = "cc", N = wif1_ec$samplesize.outcome,sdY=1)
wif1_ec_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

wif1_crc<-harmonise_data(wif1,crc)
df1 <- list(beta=wif1_crc$beta.exposure, varbeta = wif1_crc$se.exposure^2, MAF=wif1_crc$eaf.exposure, type = "quant", N = wif1_crc$samplesize.exposure,sdY=1)
df2 <- list(beta=wif1_crc$beta.outcome, varbeta = wif1_crc$se.outcome^2, MAF=wif1_crc$eaf.outcome, type = "cc", N = wif1_crc$samplesize.outcome,sdY=1)
wif1_crc_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

wif1_breast<-harmonise_data(wif1,breast)
df1 <- list(beta=wif1_breast$beta.exposure, varbeta = wif1_breast$se.exposure^2, MAF=wif1_breast$eaf.exposure, type = "quant", N = wif1_breast$samplesize.exposure,sdY=1)
df2 <- list(beta=wif1_breast$beta.outcome, varbeta = wif1_breast$se.outcome^2, MAF=wif1_breast$eaf.outcome, type = "cc", N = wif1_breast$samplesize.outcome,sdY=1)
wif1_breast_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

wif1_liver<-harmonise_data(wif1,liver)
df1 <- list(beta=wif1_liver$beta.exposure, varbeta = wif1_liver$se.exposure^2, MAF=wif1_liver$eaf.exposure, type = "quant", N = wif1_liver$samplesize.exposure,sdY=1)
df2 <- list(beta=wif1_liver$beta.outcome, varbeta = wif1_liver$se.outcome^2, MAF=wif1_liver$eaf.outcome, type = "cc", N = wif1_liver$samplesize.outcome,sdY=1)
wif1_liver_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

wif1_gallbladder<-harmonise_data(wif1,gallbladder)
df1 <- list(beta=wif1_gallbladder$beta.exposure, varbeta = wif1_gallbladder$se.exposure^2, MAF=wif1_gallbladder$eaf.exposure, type = "quant", N = wif1_gallbladder$samplesize.exposure,sdY=1)
df2 <- list(beta=wif1_gallbladder$beta.outcome, varbeta = wif1_gallbladder$se.outcome^2, MAF=wif1_gallbladder$eaf.outcome, type = "cc", N = wif1_gallbladder$samplesize.outcome,sdY=1)
wif1_gallbladder_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)

wif1_pancreas<-harmonise_data(wif1,pancreas)
df1 <- list(beta=wif1_pancreas$beta.exposure, varbeta = wif1_pancreas$se.exposure^2, MAF=wif1_pancreas$eaf.exposure, type = "quant", N = wif1_pancreas$samplesize.exposure,sdY=1)
df2 <- list(beta=wif1_pancreas$beta.outcome, varbeta = wif1_pancreas$se.outcome^2, MAF=wif1_pancreas$eaf.outcome, type = "cc", N = wif1_pancreas$samplesize.outcome,sdY=1)
wif1_pancreas_res=coloc.abf(p1=0.000001,p2=0.000001,p12=0.000001,df1,df2)



#BIG LIST OF RESULTS TO EXPORT
dfs <- lapply(ls(pattern="_res"), function(x) get(x))
summaries <- sapply(dfs, "[[", 1)
colheaders <- ls(pattern = "_res")
colnames(summaries) <- colheaders
write.table(summaries, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t", file = "coloc_results161023.txt")



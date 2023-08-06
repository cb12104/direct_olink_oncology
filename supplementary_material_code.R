##Code to conduct Mendelian randomization and colocalization analyses in "Impact of weight loss on cancer-related proteins in serm:results for a randomized-controlled trial"

#proteins: GPNMB, FURIN, WIF-1, TLR3, PPY, ERBB2, ESM1, HGF, RET
#cancers: colorectal, kidney, ovarian, and liver cancer

rm(list=ls(all=TRUE)) 

#load packages 
pacman::p_load(devtools, ggplot2, ggrepel, igraph, gridExtra, plyr, dplyr, TwoSampleMR, writexl, tidyverse, plyr, dplyr, ieugwasr, ggplot2, ggforestplot)
               
#legacy version 2smr
devtools::install_github("MRCIEU/TwoSampleMR@0.4.26")

ao <- available_outcomes()
head(ao)

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
cancer_ids <- c("ukb-b-1316", "ieu-b-4953", "ieu-a-1120")
cancer_all <- extract_outcome_data(
  snps = ukbb_instruments$rsID,
  outcomes = cancer_ids,
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.01,
  access_token = ieugwasr::check_access_token(),
  splitsize = 10000,
  proxy_splitsize = 500
)

wif1_proxies <- c("rs466776", "rs466273", "rs465530", "rs458449", "rs466104", "rs581475", "rs478380", "rs477646", "rs457919", "rs456576", "rs455429", "rs465308", "rs456500", "rs462059", "rs455668", "rs661179", "rs661169", "rs626852", "rs371493", "rs461829" )
proxies <- extract_outcome_data(
  snps = proxies,
  outcomes = "ieu-a-1120",
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.01,
  access_token = ieugwasr::check_access_token(),
  splitsize = 10000,
  proxy_splitsize = 500
)

#no wif1 proxies in oca summary stats

#crc summary stats
crc_headers <- c("chr_pos", "Chromosome", "Position", "MarkerName", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "FreqSE", "MinFreq", "MaxFreq", "beta.outcome", "se.outcome", "pval.outcome", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "RSq", "Neff", "chrom_name", "chrom_start", "SNP", "Ref", "ALT", "V6", "V7", "V8")
colon <- read.delim("ukbb_sun_pqtls_clumped_colon.txt", header = FALSE, sep = " ")
distal <- read.delim("ukbb_sun_pqtls_clumped_distal.txt", header = FALSE, sep = " ")
female_crc <- read.delim("ukbb_sun_pqtls_clumped_female.txt", header = FALSE, sep = " ")
left_crc <- read.delim("ukbb_sun_pqtls_clumped_left.txt", header = FALSE, sep = " ")
male_crc <- read.delim("ukbb_sun_pqtls_clumped_male.txt", header = FALSE, sep = " ")
crc <- read.delim("ukbb_sun_pqtls_clumped_overall.txt", header = FALSE, sep = " ")
proximal <- read.delim("ukbb_sun_pqtls_clumped_proximal.txt", header = FALSE, sep = " ")
rectal <- read.delim("ukbb_sun_pqtls_clumped_rectal.txt", header = FALSE, sep = " ")

names(colon) <- crc_headers
names(distal) <- crc_headers
names(female_crc) <- crc_headers
names(left_crc) <- crc_headers
names(male_crc) <- crc_headers
names(crc) <- crc_headers
names(proximal) <- crc_headers
names(rectal) <- crc_headers
colon$outcome = "4 colon"
distal$outcome = "2 distal"
female_crc$outcome = "female crc"
left_crc$outcome = "left crc"
male_crc$outcome = "male crc"
crc$outcome = "5 crc"
proximal$outcome = "3 proximal"
rectal$outcome = "1 rectal"

crc_stack <- bind_rows(list(colon, distal, female_crc, left_crc, male_crc, crc, proximal, rectal))

#reduce to snplist snps
crc_stack <- subset(crc_stack, SNP %in% ukbb_instruments$SNP)
crc_stack$id.outcome = crc_stack$outcome

#rcc summary stats
rcc <- read.delim("RCC.csv", header = TRUE, sep = ",")
rcc <- rcc %>% 
  rename(
    SNP = rsid_original,
    effect_allele.outcome = a1,
    other_allele.outcome = a0,
    beta.outcome = beta,
    se.outcome = beta_se,
    pval.outcome = p,
    eaf.outcome = AF_1000g_Allpop
  )

#reduce to snplist snps
rcc <- subset(rcc, SNP %in% ukbb_instruments$SNP)
# gpnmb rs75801644 snp not avail, no proxies >r2 0.8 identified using ldproxy tool
ukbb_instruments_rcc <- subset(ukbb_instruments, SNP %in% rcc$SNP)
rcc$id.outcome = "Renal cell carcinoma"
rcc$outcome = "Renal cell carcinoma"

##run MR
#ieu open gwas summary stats 
ieu_harm_dat <- harmonise_data(ukbb_instruments, cancer_all, action = 2)
ieu_mr_results <- mr(ieu_harm_dat, method_list=c("mr_wald_ratio"))

#gecco summary stats
crc_harm_dat <- harmonise_data(ukbb_instruments, crc_stack, action = 2)
crc_harm_dat$exposure = crc_harm_dat$id.exposure #add row to rcc outcome dat to permit harmonisation
crc_mr_results <- mr(crc_harm_dat)

#rcc summary stats
rcc_harm_dat <- harmonise_data(ukbb_instruments_rcc, rcc, action = 2)
rcc_harm_dat$exposure = rcc_harm_dat$id.exposure #add row to rcc outcome dat to permit harmonisation
rcc_mr_results <- mr(rcc_harm_dat)

#harmomised dfs of summary stats:
common_cols <- intersect(colnames(ieu_harm_dat), colnames(crc_harm_dat))
analysis_summary_stats <- rbind(
                            ieu_harm_dat[, common_cols], 
                            crc_harm_dat[, common_cols],
                            rcc_harm_dat[, common_cols]
                            )
write.table(analysis_summary_stats, row.names = FALSE, col.names = TRUE, quote = F, sep = '\t', file = "mr_analysis_summary_stats.txt")

#bind mr results
mr_results <- rbind(ieu_mr_results, crc_mr_results, rcc_mr_results)
mr_results[mr_results == "Type of cancer: ICD10: C64 Malignant neoplasm of kidney, except renal pelvis || id:ukb-b-1316"] <- "3 Kidney cancer (FIN)"
mr_results[mr_results == "Liver cell carcinoma || id:ieu-b-4953"] <- "3 Liver cancer (Burrows et al.)"
mr_results[mr_results == "Ovarian cancer || id:ieu-a-1120"] <- "2 Ovarian cancer (Phelan et al.)"
mr_results[mr_results == "5 crc"] <- "4 Colorectal cancer (Huyghe et al.)"
mr_results[mr_results == "Renal cell carcinoma"] <- "1 Renal cell carcinoma (Scelo et al.)"
write.table(mr_results, row.names = FALSE, col.names = TRUE, quote = F, sep = '\t', file = "pqtl_mr_results.txt")


##Colocalization analysis
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

##format protein summary stats - Pietzner et al https://www.science.org/doi/10.1126/science.abj1541
#read in protein and cancer summary stats +/- 10Mb around pqtls

#ESM1
headers <- c("snp", "MarkerName", "Allele1",  "Allele2",  "Freq1",  "FreqSE", "MinFreq",  "MaxFreq",  "beta", "se", "Pvalue", "Direction",  "HetISq", "HetChiSq", "HetDf",  "HetPVal",  "TotalSampleSize",  "chr", "position")
esm1 <- read.delim("esm1_coloc_10mb.txt", header = FALSE, sep = "\t")
names(esm1) <- headers
esm1$sdY = 1
esm1$type = "quant"
esm1$varbeta = esm1$se^2

#FURIN
furin <- read.delim("furin_coloc_10mb.txt", header = FALSE, sep = "\t")
names(furin) <- headers
furin$sdY = 1
furin$type = "quant"
furin$varbeta = furin$se^2

#GPNMB
gpnmb <- read.delim("gpnmb_coloc_10mb.txt", header = FALSE, sep = "\t")
names(gpnmb) <- headers
gpnmb$sdY = 1
gpnmb$type = "quant"
gpnmb$varbeta = gpnmb$se^2

#HGF
hgf <- read.delim("hgf_coloc_10mb.txt", header = FALSE, sep = "\t")
names(hgf) <- headers
hgf$sdY = 1
hgf$type = "quant"
hgf$varbeta = hgf$se^2

#RET
ret <- read.delim("ret_coloc_10mb.txt", header = FALSE, sep = "\t")
names(ret) <- headers
ret$sdY = 1
ret$type = "quant"
ret$varbeta = ret$se^2

#TLR3
tlr3 <- read.delim("tlr3_coloc_10mb.txt", header = FALSE, sep = "\t")
names(tlr3) <- headers
tlr3$sdY = 1
tlr3$type = "quant"
tlr3$varbeta = tlr3$se^2

#WIF1
wif1 <- read.delim("wif1_coloc_10mb.txt", header = FALSE, sep = "\t")
names(wif1) <- headers
wif1$sdY = 1
wif1$type = "quant"
wif1$varbeta = wif1$se^2

##format crc summary stats
crc_headers <- c("chr_pos", "Chromosome", "position", "MarkerName", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "FreqSE", "MinFreq", "MaxFreq", "beta", "se.outcome", "pval.outcome", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "RSq", "Neff", "chrom_name", "chrom_start", "snp", "Ref", "ALT", "V6", "V7", "V8")

#ESM1
crc_esm1 <- read.delim("esm1_coloc_gecco_10mb.txt", header = FALSE , sep = " ")
names(crc_esm1) <- crc_headers
crc_esm1$sdY = 1
crc_esm1$type = "cc"
crc_esm1$varbeta = crc_esm1$se.outcome^2
#remove duplicated snps
crc_esm1 <- crc_esm1 %>% distinct(snp, .keep_all = TRUE)
which(is.na(crc_esm1), arr.ind=TRUE)

colon_esm1 <- read.delim("esm1_coloc_gecco_10mb_colon.txt", header = FALSE , sep = " ")
names(colon_esm1) <- crc_headers
colon_esm1$sdY = 1
colon_esm1$type = "cc"
colon_esm1$varbeta = colon_esm1$se.outcome^2
#remove duplicated snps
colon_esm1 <- colon_esm1 %>% distinct(snp, .keep_all = TRUE)

rectal_esm1 <- read.delim("esm1_coloc_gecco_10mb_rectal.txt", header = FALSE , sep = " ")
names(rectal_esm1) <- crc_headers
rectal_esm1$sdY = 1
rectal_esm1$type = "cc"
rectal_esm1$varbeta = rectal_esm1$se.outcome^2
#remove duplicated snps
rectal_esm1 <- rectal_esm1 %>% distinct(snp, .keep_all = TRUE)

proximal_esm1 <- read.delim("esm1_coloc_gecco_10mb_proximal.txt", header = FALSE , sep = " ")
names(proximal_esm1) <- crc_headers
proximal_esm1$sdY = 1
proximal_esm1$type = "cc"
proximal_esm1$varbeta = proximal_esm1$se.outcome^2
#remove duplicated snps
proximal_esm1 <- proximal_esm1 %>% distinct(snp, .keep_all = TRUE)

distal_esm1 <- read.delim("esm1_coloc_gecco_10mb_distal.txt", header = FALSE , sep = " ")
names(distal_esm1) <- crc_headers
distal_esm1$sdY = 1
distal_esm1$type = "cc"
distal_esm1$varbeta = distal_esm1$se.outcome^2
#remove duplicated snps
distal_esm1 <- distal_esm1 %>% distinct(snp, .keep_all = TRUE)

#FURIN
crc_furin <- read.delim("furin_coloc_gecco_10mb.txt", header = FALSE , sep = " ")
names(crc_furin) <- crc_headers
crc_furin$sdY = 1
crc_furin$type = "cc"
crc_furin$varbeta = crc_furin$se.outcome^2
#remove duplicated snps
crc_furin <- crc_furin %>% distinct(snp, .keep_all = TRUE)

colon_furin <- read.delim("furin_coloc_gecco_10mb_colon.txt", header = FALSE , sep = " ")
names(colon_furin) <- crc_headers
colon_furin$sdY = 1
colon_furin$type = "cc"
colon_furin$varbeta = colon_furin$se.outcome^2
#remove duplicated snps
colon_furin <- colon_furin %>% distinct(snp, .keep_all = TRUE)

rectal_furin <- read.delim("furin_coloc_gecco_10mb_rectal.txt", header = FALSE , sep = " ")
names(rectal_furin) <- crc_headers
rectal_furin$sdY = 1
rectal_furin$type = "cc"
rectal_furin$varbeta = rectal_furin$se.outcome^2
#remove duplicated snps
rectal_furin <- rectal_furin %>% distinct(snp, .keep_all = TRUE)

proximal_furin <- read.delim("furin_coloc_gecco_10mb_proximal.txt", header = FALSE , sep = " ")
names(proximal_furin) <- crc_headers
proximal_furin$sdY = 1
proximal_furin$type = "cc"
proximal_furin$varbeta = proximal_furin$se.outcome^2
#remove duplicated snps
proximal_furin <- proximal_furin %>% distinct(snp, .keep_all = TRUE)

distal_furin <- read.delim("furin_coloc_gecco_10mb_distal.txt", header = FALSE , sep = " ")
names(distal_furin) <- crc_headers
distal_furin$sdY = 1
distal_furin$type = "cc"
distal_furin$varbeta = distal_furin$se.outcome^2
#remove duplicated snps
distal_furin <- distal_furin %>% distinct(snp, .keep_all = TRUE)

#GPNMB
crc_gpnmb <- read.delim("gpnmb_coloc_gecco_10mb.txt", header = FALSE , sep = " ")
names(crc_gpnmb) <- crc_headers
crc_gpnmb$sdY = 1
crc_gpnmb$type = "cc"
crc_gpnmb$varbeta = crc_gpnmb$se.outcome^2
#remove duplicated snps
crc_gpnmb <- crc_gpnmb %>% distinct(snp, .keep_all = TRUE)

colon_gpnmb <- read.delim("gpnmb_coloc_gecco_10mb_colon.txt", header = FALSE , sep = " ")
names(colon_gpnmb) <- crc_headers
colon_gpnmb$sdY = 1
colon_gpnmb$type = "cc"
colon_gpnmb$varbeta = colon_gpnmb$se.outcome^2
#remove duplicated snps
colon_gpnmb <- colon_gpnmb %>% distinct(snp, .keep_all = TRUE)

rectal_gpnmb <- read.delim("gpnmb_coloc_gecco_10mb_rectal.txt", header = FALSE , sep = " ")
names(rectal_gpnmb) <- crc_headers
rectal_gpnmb$sdY = 1
rectal_gpnmb$type = "cc"
rectal_gpnmb$varbeta = rectal_gpnmb$se.outcome^2
#remove duplicated snps
rectal_gpnmb <- rectal_gpnmb %>% distinct(snp, .keep_all = TRUE)

proximal_gpnmb <- read.delim("gpnmb_coloc_gecco_10mb_proximal.txt", header = FALSE , sep = " ")
names(proximal_gpnmb) <- crc_headers
proximal_gpnmb$sdY = 1
proximal_gpnmb$type = "cc"
proximal_gpnmb$varbeta = proximal_gpnmb$se.outcome^2
#remove duplicated snps
proximal_gpnmb <- proximal_gpnmb %>% distinct(snp, .keep_all = TRUE)

distal_gpnmb <- read.delim("gpnmb_coloc_gecco_10mb_distal.txt", header = FALSE , sep = " ")
names(distal_gpnmb) <- crc_headers
distal_gpnmb$sdY = 1
distal_gpnmb$type = "cc"
distal_gpnmb$varbeta = distal_gpnmb$se.outcome^2
#remove duplicated snps
distal_gpnmb <- distal_gpnmb %>% distinct(snp, .keep_all = TRUE)

#HGF
crc_hgf <- read.delim("hgf_coloc_gecco_10mb.txt", header = FALSE , sep = " ")
names(crc_hgf) <- crc_headers
crc_hgf$sdY = 1
crc_hgf$type = "cc"
crc_hgf$varbeta = crc_hgf$se.outcome^2
#remove duplicated snps
crc_hgf <- crc_hgf %>% distinct(snp, .keep_all = TRUE)

colon_hgf <- read.delim("hgf_coloc_gecco_10mb_colon.txt", header = FALSE , sep = " ")
names(colon_hgf) <- crc_headers
colon_hgf$sdY = 1
colon_hgf$type = "cc"
colon_hgf$varbeta = colon_hgf$se.outcome^2
#remove duplicated snps
colon_hgf <- colon_hgf %>% distinct(snp, .keep_all = TRUE)

rectal_hgf <- read.delim("hgf_coloc_gecco_10mb_rectal.txt", header = FALSE , sep = " ")
names(rectal_hgf) <- crc_headers
rectal_hgf$sdY = 1
rectal_hgf$type = "cc"
rectal_hgf$varbeta = rectal_hgf$se.outcome^2
#remove duplicated snps
rectal_hgf <- rectal_hgf %>% distinct(snp, .keep_all = TRUE)

proximal_hgf <- read.delim("hgf_coloc_gecco_10mb_proximal.txt", header = FALSE , sep = " ")
names(proximal_hgf) <- crc_headers
proximal_hgf$sdY = 1
proximal_hgf$type = "cc"
proximal_hgf$varbeta = proximal_hgf$se.outcome^2
#remove duplicated snps
proximal_hgf <- proximal_hgf %>% distinct(snp, .keep_all = TRUE)

distal_hgf <- read.delim("hgf_coloc_gecco_10mb_distal.txt", header = FALSE , sep = " ")
names(distal_hgf) <- crc_headers
distal_hgf$sdY = 1
distal_hgf$type = "cc"
distal_hgf$varbeta = distal_hgf$se.outcome^2
#remove duplicated snps
distal_hgf <- distal_hgf %>% distinct(snp, .keep_all = TRUE)

#RET
crc_ret <- read.delim("ret_coloc_gecco_10mb.txt", header = FALSE , sep = " ")
names(crc_ret) <- crc_headers
crc_ret$sdY = 1
crc_ret$type = "cc"
crc_ret$varbeta = crc_ret$se.outcome^2
#remove duplicated snps
crc_ret <- crc_ret %>% distinct(snp, .keep_all = TRUE)

colon_ret <- read.delim("ret_coloc_gecco_10mb_colon.txt", header = FALSE , sep = " ")
names(colon_ret) <- crc_headers
colon_ret$sdY = 1
colon_ret$type = "cc"
colon_ret$varbeta = colon_ret$se.outcome^2
#remove duplicated snps
colon_ret <- colon_ret %>% distinct(snp, .keep_all = TRUE)

rectal_ret <- read.delim("ret_coloc_gecco_10mb_rectal.txt", header = FALSE , sep = " ")
names(rectal_ret) <- crc_headers
rectal_ret$sdY = 1
rectal_ret$type = "cc"
rectal_ret$varbeta = rectal_ret$se.outcome^2
#remove duplicated snps
rectal_ret <- rectal_ret %>% distinct(snp, .keep_all = TRUE)

proximal_ret <- read.delim("ret_coloc_gecco_10mb_proximal.txt", header = FALSE , sep = " ")
names(proximal_ret) <- crc_headers
proximal_ret$sdY = 1
proximal_ret$type = "cc"
proximal_ret$varbeta = proximal_ret$se.outcome^2
#remove duplicated snps
proximal_ret <- proximal_ret %>% distinct(snp, .keep_all = TRUE)

distal_ret <- read.delim("ret_coloc_gecco_10mb_distal.txt", header = FALSE , sep = " ")
names(distal_ret) <- crc_headers
distal_ret$sdY = 1
distal_ret$type = "cc"
distal_ret$varbeta = distal_ret$se.outcome^2
#remove duplicated snps
distal_ret <- distal_ret %>% distinct(snp, .keep_all = TRUE)

#TLR3
crc_tlr3 <- read.delim("tlr3_coloc_gecco_10mb.txt", header = FALSE , sep = " ")
names(crc_tlr3) <- crc_headers
crc_tlr3$sdY = 1
crc_tlr3$type = "cc"
crc_tlr3$varbeta = crc_tlr3$se.outcome^2
#remove duplicated snps
crc_tlr3 <- crc_tlr3 %>% distinct(snp, .keep_all = TRUE)

colon_tlr3 <- read.delim("tlr3_coloc_gecco_10mb_colon.txt", header = FALSE , sep = " ")
names(colon_tlr3) <- crc_headers
colon_tlr3$sdY = 1
colon_tlr3$type = "cc"
colon_tlr3$varbeta = colon_tlr3$se.outcome^2
#remove duplicated snps
colon_tlr3 <- colon_tlr3 %>% distinct(snp, .keep_all = TRUE)

rectal_tlr3 <- read.delim("tlr3_coloc_gecco_10mb_rectal.txt", header = FALSE , sep = " ")
names(rectal_tlr3) <- crc_headers
rectal_tlr3$sdY = 1
rectal_tlr3$type = "cc"
rectal_tlr3$varbeta = rectal_tlr3$se.outcome^2
#remove duplicated snps
rectal_tlr3 <- rectal_tlr3 %>% distinct(snp, .keep_all = TRUE)

proximal_tlr3 <- read.delim("tlr3_coloc_gecco_10mb_proximal.txt", header = FALSE , sep = " ")
names(proximal_tlr3) <- crc_headers
proximal_tlr3$sdY = 1
proximal_tlr3$type = "cc"
proximal_tlr3$varbeta = proximal_tlr3$se.outcome^2
#remove duplicated snps
proximal_tlr3 <- proximal_tlr3 %>% distinct(snp, .keep_all = TRUE)

distal_tlr3 <- read.delim("tlr3_coloc_gecco_10mb_distal.txt", header = FALSE , sep = " ")
names(distal_tlr3) <- crc_headers
distal_tlr3$sdY = 1
distal_tlr3$type = "cc"
distal_tlr3$varbeta = distal_tlr3$se.outcome^2
#remove duplicated snps
distal_tlr3 <- distal_tlr3 %>% distinct(snp, .keep_all = TRUE)

#WIF1
crc_wif1 <- read.delim("wif1_coloc_gecco_10mb.txt", header = FALSE , sep = " ")
names(crc_wif1) <- crc_headers
crc_wif1$sdY = 1
crc_wif1$type = "cc"
crc_wif1$varbeta = crc_wif1$se.outcome^2
#remove duplicated snps
crc_wif1 <- crc_wif1 %>% distinct(snp, .keep_all = TRUE)

colon_wif1 <- read.delim("wif1_coloc_gecco_10mb_colon.txt", header = FALSE , sep = " ")
names(colon_wif1) <- crc_headers
colon_wif1$sdY = 1
colon_wif1$type = "cc"
colon_wif1$varbeta = colon_wif1$se.outcome^2
#remove duplicated snps
colon_wif1 <- colon_wif1 %>% distinct(snp, .keep_all = TRUE)

rectal_wif1 <- read.delim("wif1_coloc_gecco_10mb_rectal.txt", header = FALSE , sep = " ")
names(rectal_wif1) <- crc_headers
rectal_wif1$sdY = 1
rectal_wif1$type = "cc"
rectal_wif1$varbeta = rectal_wif1$se.outcome^2
#remove duplicated snps
rectal_wif1 <- rectal_wif1 %>% distinct(snp, .keep_all = TRUE)

proximal_wif1 <- read.delim("wif1_coloc_gecco_10mb_proximal.txt", header = FALSE , sep = " ")
names(proximal_wif1) <- crc_headers
proximal_wif1$sdY = 1
proximal_wif1$type = "cc"
proximal_wif1$varbeta = proximal_wif1$se.outcome^2
#remove duplicated snps
proximal_wif1 <- proximal_wif1 %>% distinct(snp, .keep_all = TRUE)

distal_wif1 <- read.delim("wif1_coloc_gecco_10mb_distal.txt", header = FALSE , sep = " ")
names(distal_wif1) <- crc_headers
distal_wif1$sdY = 1
distal_wif1$type = "cc"
distal_wif1$varbeta = distal_wif1$se.outcome^2
#remove duplicated snps
distal_wif1 <- distal_wif1 %>% distinct(snp, .keep_all = TRUE)

#read in all rcc snps
##format rcc summary stats
rcc_headers <- c("snp", "effect_allele.outcome", "other_allele.outcome", "beta", "se.outcome", "pval.outcome", "Direction", "chr_pos_hg19", "rsid_new_dbsnp150", "INFO_avg", "eaf_1000g_allpop", "chromosome", "position") 
rcc <- read.delim("for_Caroline_Bull_RCC.csv", header = TRUE , sep = ",")
names(rcc) <- rcc_headers
rcc$sdY = 1
rcc$type = "cc"
rcc$varbeta = rcc$se.outcome^2
#remove duplicated snps
rcc <- rcc %>% distinct(snp, .keep_all = TRUE)

##RUN COLOC##
#minimum info for coloc = "beta","varbeta","snp","position","type","sdY"
#ESM1/crc
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(crc_esm1, main="Dataset CRC")
crc_esm1 <- crc_esm1[complete.cases(crc_esm1[ , 21]),] #this line removes rows with missing SNP info required for coloc 
vres_esm1_crc <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=esm1, dataset2=crc_esm1)

par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(colon_esm1, main="Dataset Colon")
colon_esm1 <- colon_esm1[complete.cases(colon_esm1[ , 21]),]
vres_esm1_colon <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=esm1, dataset2=colon_esm1)

par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESDM1")
plot_dataset(rectal_esm1, main="Dataset Rectal")
rectal_esm1 <- rectal_esm1[complete.cases(rectal_esm1[ , 21]),]
vres_esm1_rectal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=esm1, dataset2=rectal_esm1)

par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(proximal_esm1, main="Dataset Proximal")
proximal_esm1 <- proximal_esm1[complete.cases(proximal_esm1[ , 21]),]
vres_esm1_proximal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=esm1, dataset2=proximal_esm1)

par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(distal_esm1, main="Dataset Distal")
distal_esm1 <- distal_esm1[complete.cases(distal_esm1[ , 21]),]
vres_esm1_distal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=esm1, dataset2=distal_esm1)

##rcc analysis
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
#subset chr5 for esm1 analysis
rcc_esm1 <- rcc %>% filter (chromosome == 5) 
plot_dataset(rcc_esm1, main="Dataset RCC")
vres_esm1_rcc <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=esm1, dataset2=rcc_esm1)

#FURIN/crc
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset FURIN")
plot_dataset(crc_furin, main="Dataset CRC")
#line 2068 doesn't have snp info data - remove to run coloc
crc_furin <- crc_furin[complete.cases(crc_furin[ , 21]),]
vres_furin <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=furin, dataset2=crc_furin)

par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(colon_furin, main="Dataset Colon")
colon_furin <- colon_furin[complete.cases(colon_furin[ , 21]),]
vres_furin_colon <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=furin, dataset2=colon_furin)

par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(rectal_furin, main="Dataset Rectal")
rectal_furin <- rectal_furin[complete.cases(rectal_furin[ , 21]),]
vres_furin_rectal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=furin, dataset2=rectal_furin)

par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(proximal_furin, main="Dataset Proximal")
proximal_furin <- proximal_furin[complete.cases(proximal_furin[ , 21]),]
vres_furin_proximal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=furin, dataset2=proximal_furin)

par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(distal_furin, main="Dataset Distal")
distal_furin <- distal_furin[complete.cases(distal_furin[ , 21]),]
vres_furin_distal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=furin, dataset2=distal_furin)

##rcc analysis
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
#subset chr15 for furin analysis
rcc_furin <- rcc %>% filter (chromosome == 15) 
plot_dataset(rcc_furin, main="Dataset RCC")
vres_furin_rcc <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=furin, dataset2=rcc_furin)


#GPNMB/crc
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset GPNMB")
plot_dataset(crc_gpnmb, main="Dataset CRC")
crc_gpnmb <- crc_gpnmb[complete.cases(crc_gpnmb[ , 21]),]
vres_gpnmb <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=gpnmb, dataset2=crc_gpnmb)

par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(colon_gpnmb, main="Dataset Colon")
colon_gpnmb <- colon_gpnmb[complete.cases(colon_gpnmb[ , 21]),]
vres_gpnmb_colon <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=gpnmb, dataset2=colon_gpnmb)

par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(rectal_gpnmb, main="Dataset Rectal")
rectal_gpnmb <- rectal_gpnmb[complete.cases(rectal_gpnmb[ , 21]),]
vres_gpnmb_rectal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=gpnmb, dataset2=rectal_gpnmb)

par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(proximal_gpnmb, main="Dataset Proximal")
proximal_gpnmb <- proximal_gpnmb[complete.cases(proximal_gpnmb[ , 21]),]
vres_gpnmb_proximal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=gpnmb, dataset2=proximal_gpnmb)

par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(distal_gpnmb, main="Dataset Distal")
distal_gpnmb <- distal_gpnmb[complete.cases(distal_gpnmb[ , 21]),]
vres_gpnmb_distal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=gpnmb, dataset2=distal_gpnmb)

##rcc analysis
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
#subset chr7 for gpnmb analysis
rcc_gpnmb <- rcc %>% filter (chromosome == 7) 
rcc_gpnmb <- rcc_gpnmb %>% filter (position > 21766522, position < 24766522)
plot_dataset(rcc_gpnmb, main="Dataset RCC")
vres_gpnmb_rcc <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=gpnmb, dataset2=rcc_gpnmb)


#HGF/crc
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset HGF")
plot_dataset(crc_hgf, main="Dataset CRC")
crc_hgf <- crc_hgf[complete.cases(crc_hgf[ , 21]),]
vres_hgf <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=hgf, dataset2=crc_hgf)

par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(colon_hgf, main="Dataset Colon")
colon_hgf <- colon_hgf[complete.cases(colon_hgf[ , 21]),]
vres_hgf_colon <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=hgf, dataset2=colon_hgf)

par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(rectal_hgf, main="Dataset Rectal")
rectal_hgf <- rectal_hgf[complete.cases(rectal_hgf[ , 21]),]
vres_hgf_rectal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=hgf, dataset2=rectal_hgf)

par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(proximal_hgf, main="Dataset Proximal")
proximal_hgf <- proximal_hgf[complete.cases(proximal_hgf[ , 21]),]
vres_hgf_proximal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=hgf, dataset2=proximal_hgf)

par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(distal_hgf, main="Dataset Distal")
distal_hgf <- distal_hgf[complete.cases(distal_hgf[ , 21]),]
vres_hgf_distal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=hgf, dataset2=distal_hgf)

##rcc analysis
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
#subset chr7 for hgf analysis
rcc_hgf <- rcc %>% filter (chromosome == 7) 
rcc_hgf <- rcc_hgf %>% filter (position > 80229735, position < 83229735)
plot_dataset(rcc_hgf, main="Dataset RCC")
vres_hgf_rcc <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=hgf, dataset2=rcc_hgf)

#RET/crc
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset RET")
plot_dataset(crc_ret, main="Dataset CRC")
crc_ret <- crc_ret[complete.cases(crc_ret[ , 21]),]
vres_ret <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=ret, dataset2=crc_ret)

par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(colon_ret, main="Dataset Colon")
colon_ret <- colon_ret[complete.cases(colon_ret[ , 21]),]
vres_ret_colon <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=ret, dataset2=colon_ret)

par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(rectal_ret, main="Dataset Rectal")
rectal_ret <- rectal_ret[complete.cases(rectal_ret[ , 21]),]
vres_ret_rectal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=ret, dataset2=rectal_ret)

par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(proximal_ret, main="Dataset Proximal")
proximal_ret <- proximal_ret[complete.cases(proximal_ret[ , 21]),]
vres_ret_proximal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=ret, dataset2=proximal_ret)

par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(distal_ret, main="Dataset Distal")
distal_ret <- distal_ret[complete.cases(distal_ret[ , 21]),]
vres_ret_distal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=ret, dataset2=distal_ret)

##rcc analysis
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
#subset chr10 for ret analysis
rcc_ret <- rcc %>% filter (chromosome == 10) 
plot_dataset(rcc_ret, main="Dataset RCC")
vres_ret_rcc <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=ret, dataset2=rcc_ret)

#TLR3/crc
par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset TLR3")
plot_dataset(crc_tlr3, main="Dataset CRC")
crc_tlr3 <- crc_tlr3[complete.cases(crc_tlr3[ , 21]),]
vres_tlr3 <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=tlr3, dataset2=crc_tlr3)

par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(colon_tlr3, main="Dataset Colon")
colon_tlr3 <- colon_tlr3[complete.cases(colon_tlr3[ , 21]),]
vres_tlr3_colon <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=tlr3, dataset2=colon_tlr3)

par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(rectal_tlr3, main="Dataset Rectal")
rectal_tlr3 <- rectal_tlr3[complete.cases(rectal_tlr3[ , 21]),]
vres_tlr3_rectal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=tlr3, dataset2=rectal_tlr3)

par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(proximal_tlr3, main="Dataset Proximal")
proximal_tlr3 <- proximal_tlr3[complete.cases(proximal_tlr3[ , 21]),]
vres_tlr3_proximal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=tlr3, dataset2=proximal_tlr3)

par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(distal_tlr3, main="Dataset Distal")
distal_tlr3 <- distal_tlr3[complete.cases(distal_tlr3[ , 21]),]
vres_tlr3_distal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=tlr3, dataset2=distal_tlr3)

##rcc analysis
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset tlr3")
#subset chr4 for tlr3 analysis
rcc_tlr3 <- rcc %>% filter (chromosome == 4) 
plot_dataset(rcc_tlr3, main="Dataset RCC")
vres_tlr3_rcc <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=tlr3, dataset2=rcc_tlr3)

#WIF1/crc
par(mfrow=c(2,1))
plot_dataset(wif1, main="Dataset WIF1")
plot_dataset(crc_wif1, main="Dataset CRC")
crc_wif1 <- crc_wif1[complete.cases(crc_wif1[ , 21]),]
vres_wif1 <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=wif1, dataset2=crc_wif1)

par(mfrow=c(2,1))
plot_dataset(wif1, main="Dataset wif1")
plot_dataset(colon_wif1, main="Dataset Colon")
colon_wif1 <- colon_wif1[complete.cases(colon_wif1[ , 21]),]
vres_wif1_colon <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=wif1, dataset2=colon_wif1)

par(mfrow=c(2,1))
plot_dataset(wif1, main="Dataset wif1")
plot_dataset(rectal_wif1, main="Dataset Rectal")
rectal_wif1 <- rectal_wif1[complete.cases(rectal_wif1[ , 21]),]
vres_wif1_rectal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=wif1, dataset2=rectal_wif1)

par(mfrow=c(2,1))
plot_dataset(wif1, main="Dataset wif1")
plot_dataset(proximal_wif1, main="Dataset Proximal")
proximal_wif1 <- proximal_wif1[complete.cases(proximal_wif1[ , 21]),]
vres_wif1_proximal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=wif1, dataset2=proximal_wif1)

par(mfrow=c(2,1))
plot_dataset(wif1, main="Dataset wif1")
plot_dataset(distal_wif1, main="Dataset Distal")
distal_wif1 <- distal_wif1[complete.cases(distal_wif1[ , 21]),]
vres_wif1_distal <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=wif1, dataset2=distal_wif1)

##rcc analysis
par(mfrow=c(2,1))
plot_dataset(wif1, main="Dataset wif1")
#subset chr12 for wif1 analysis
rcc_wif1 <- rcc %>% filter (chromosome == 12) 
plot_dataset(rcc_wif1, main="Dataset RCC")
vres_wif1_rcc <- coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,dataset1=wif1, dataset2=rcc_wif1)

##LIVER CANCER using gwasglue package
#esm1
#chrpos_esm1 <- "5:54777867-55185593"
chrpos_esm1 <- "5:44198775-64198775"
#extract twice so that function works
liver_esm1 <- gwasvcf_to_coloc("ieu-b-4953.vcf.gz", "ieu-b-4953.vcf.gz", chrpos_esm1)
#df format
liver_c_esm1 <- as.data.frame(liver_esm1[[1]])
#remove duplicated snps
liver_c_esm1 <- liver_c_esm1 %>% distinct(snp, .keep_all = TRUE)
liver_c_esm1$sdY = 1
liver_c_esm1$position <- liver_c_esm1$pos
#plot
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(liver_c_esm1, main="Dataset Liver Carcinoma")
#run colocalisation analysis
vres_esm1_liver <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,esm1, liver_c_esm1)

#furin
#chrpos_furin <- "15:90668588-91083457"
chrpos_furin <- "15:81425232-101425232"
#extract twice so that function works
liver_furin <- gwasvcf_to_coloc("ieu-b-4953.vcf.gz", "ieu-b-4953.vcf.gz", chrpos_furin)
#df format
liver_c_furin <- as.data.frame(liver_furin[[1]])
#remove duplicated snps
liver_c_furin <- liver_c_furin %>% distinct(snp, .keep_all = TRUE)
liver_c_furin$sdY = 1
liver_c_furin$position <- liver_c_furin$pos
#plot 
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset FURIN")
plot_dataset(liver_c_furin, main="Dataset Liver Carcinoma")
#run colocalisation analysis
vres_furin_liver <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,furin, liver_c_furin)

#gpnmb
#chrpos_gpnmb <- "7:23046746-23475096"
chrpos_gpnmb <- "7:13306141-33306141"
#extract twice so that function works
liver_gpnmb <- gwasvcf_to_coloc("ieu-b-4953.vcf.gz", "ieu-b-4953.vcf.gz", chrpos_gpnmb)
#df format
liver_c_gpnmb <- as.data.frame(liver_gpnmb[[1]])
#remove duplicated snps
liver_c_gpnmb <- liver_c_gpnmb %>% distinct(snp, .keep_all = TRUE)
liver_c_gpnmb$sdY = 1
liver_c_gpnmb$position <- liver_c_gpnmb$pos
#plot
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset GPNMB")
plot_dataset(liver_c_gpnmb, main="Dataset Liver Carcinoma")
#run colocalisation analysis
vres_gpnmb_liver <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,gpnmb, liver_c_gpnmb)

#hgf
#chrpos_hgf <- "7:81499010-81970047"
chrpos_hgf <- "7:71359051-91359051"
#extract twice so that function works
liver_hgf <- gwasvcf_to_coloc("ieu-b-4953.vcf.gz", "ieu-b-4953.vcf.gz", chrpos_hgf)
#df format
liver_c_hgf <- as.data.frame(liver_hgf[[1]])
#remove duplicated snps
liver_c_hgf <- liver_c_hgf %>% distinct(snp, .keep_all = TRUE)
liver_c_hgf$sdY = 1
liver_c_hgf$position <- liver_c_hgf$pos
#plot
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset HGF")
plot_dataset(liver_c_hgf, main="Dataset Liver Carcinoma")
#run colocalisation analysis
vres_hgf_liver <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,hgf, liver_c_hgf)

#ret
#chrpos_ret <- "10:42877069-43327504"
chrpos_ret <- "10:33352894-53352894"
#extract twice so that function works
liver_ret <- gwasvcf_to_coloc("ieu-b-4953.vcf.gz", "ieu-b-4953.vcf.gz", chrpos_ret)
#df format
liver_c_ret <- as.data.frame(liver_ret[[1]])
#remove duplicated snps
liver_c_ret <- liver_c_ret %>% distinct(snp, .keep_all = TRUE)
liver_c_ret$sdY = 1
liver_c_ret$position <- liver_c_ret$pos
#plot
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset RET")
plot_dataset(liver_c_ret, main="Dataset Liver Carcinoma")
#run colocalisation analysis #issue harmonising snps
vres_ret_liver <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,ret, liver_c_ret)

#tlr3
#chrpos_tlr3 <- "4:185869156-186288073"
chrpos_tlr3 <- "4:177004074-197004074"
#extract twice so that function works
liver_tlr3 <- gwasvcf_to_coloc("ieu-b-4953.vcf.gz", "ieu-b-4953.vcf.gz", chrpos_tlr3)
#df format
liver_c_tlr3 <- as.data.frame(liver_tlr3[[1]])
#remove duplicated snps
liver_c_tlr3 <- liver_c_tlr3 %>% distinct(snp, .keep_all = TRUE)
liver_c_tlr3$sdY = 1
liver_c_tlr3$position <- liver_c_tlr3$pos
#plot
par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset TLR3")
plot_dataset(liver_c_tlr3, main="Dataset Liver Carcinoma")
#run colocalisation analysis #issue harmonising snps
vres_tlr3_liver <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,tlr3, liver_c_tlr3)

#wif1
#chrpos_wif1 <- "12:64850626-65321305"
chrpos_wif1 <- "12:55338145-75338145"
#extract twice so that function works
liver_wif1 <- gwasvcf_to_coloc("ieu-b-4953.vcf.gz", "ieu-b-4953.vcf.gz", chrpos_wif1)
#df format
liver_c_wif1 <- as.data.frame(liver_wif1[[1]])
#remove duplicated snps
liver_c_wif1 <- liver_c_wif1 %>% distinct(snp, .keep_all = TRUE)
liver_c_wif1$sdY = 1
liver_c_wif1$position <- liver_c_wif1$pos
#plot
par(mfrow=c(2,1))
plot_dataset(wif1, main="Dataset WIF1")
plot_dataset(liver_c_wif1, main="Dataset Liver Carcinoma")
#run colocalisation analysis #issue harmonising snps
vres_wif1_liver <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,wif1, liver_c_wif1)

##OVARIAN CANCER
#esm1 #oca overall
#extract twice so that function works
ovarian_esm1 <- gwasvcf_to_coloc("ieu-a-1120.vcf.gz", "ieu-a-1120.vcf.gz", chrpos_esm1)
#df format
ovarian_c_esm1 <- as.data.frame(ovarian_esm1[[1]])
#remove duplicated snps
ovarian_c_esm1 <- ovarian_c_esm1 %>% distinct(snp, .keep_all = TRUE)
ovarian_c_esm1$sdY = 1
ovarian_c_esm1$position <- ovarian_c_esm1$pos
#plot
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(ovarian_c_esm1, main="Dataset ovarian cancer")
#run colocalisation analysis
vres_esm1_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,esm1, ovarian_c_esm1)

#esm1 #invasive mucinous ovarian cancer 
#extract twice so that function works
im_ovarian_esm1 <- gwasvcf_to_coloc("ieu-a-1123.vcf.gz", "ieu-a-1123.vcf.gz", chrpos_esm1)
#df format
im_ovarian_c_esm1 <- as.data.frame(im_ovarian_esm1[[1]])
#remove duplicated snps
im_ovarian_c_esm1 <- im_ovarian_c_esm1 %>% distinct(snp, .keep_all = TRUE)
im_ovarian_c_esm1$sdY = 1
im_ovarian_c_esm1$position <- im_ovarian_c_esm1$pos
#plot
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(im_ovarian_c_esm1, main="Dataset invasive mucinous ovarian cancer")
#run colocalisation analysis
vres_esm1_im_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,esm1, im_ovarian_c_esm1)

#esm1 #high grade serous ovarian cancer 
#extract twice so that function works
hg_ovarian_esm1 <- gwasvcf_to_coloc("ieu-a-1121.vcf.gz", "ieu-a-1121.vcf.gz", chrpos_esm1)
#df format
hg_ovarian_c_esm1 <- as.data.frame(hg_ovarian_esm1[[1]])
#remove duplicated snps
hg_ovarian_c_esm1 <- hg_ovarian_c_esm1 %>% distinct(snp, .keep_all = TRUE)
hg_ovarian_c_esm1$sdY = 1
hg_ovarian_c_esm1$position <- hg_ovarian_c_esm1$pos
#plot
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(hg_ovarian_c_esm1, main="Dataset high grade serous ovarian cancer")
#run colocalisation analysis
vres_esm1_hg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,esm1, hg_ovarian_c_esm1)

#esm1 #low grade serous ovarian cancer 
#extract twice so that function works
lg_ovarian_esm1 <- gwasvcf_to_coloc("ieu-a-1122.vcf.gz", "ieu-a-1122.vcf.gz", chrpos_esm1)
#df format
lg_ovarian_c_esm1 <- as.data.frame(lg_ovarian_esm1[[1]])
#remove duplicated snps
lg_ovarian_c_esm1 <- lg_ovarian_c_esm1 %>% distinct(snp, .keep_all = TRUE)
lg_ovarian_c_esm1$sdY = 1
lg_ovarian_c_esm1$position <- lg_ovarian_c_esm1$pos
#plot
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(lg_ovarian_c_esm1, main="Dataset low grade serous ovarian cancer")
#run colocalisation analysis
vres_esm1_lg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,esm1, lg_ovarian_c_esm1)

#esm1 #endometrioid ovarian cancer 
#extract twice so that function works
endo_ovarian_esm1 <- gwasvcf_to_coloc("ieu-a-1125.vcf.gz", "ieu-a-1125.vcf.gz", chrpos_esm1)
#df format
endo_ovarian_c_esm1 <- as.data.frame(endo_ovarian_esm1[[1]])
#remove duplicated snps
endo_ovarian_c_esm1 <- endo_ovarian_c_esm1 %>% distinct(snp, .keep_all = TRUE)
endo_ovarian_c_esm1$sdY = 1
endo_ovarian_c_esm1$position <- endo_ovarian_c_esm1$pos
#plot
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(endo_ovarian_c_esm1, main="Dataset endometrioid ovarian cancer")
#run colocalisation analysis
vres_esm1_endo_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,esm1, endo_ovarian_c_esm1)

#esm1 #clear cell ovarian cancer 
#extract twice so that function works
cc_ovarian_esm1 <- gwasvcf_to_coloc("ieu-a-1124.vcf.gz", "ieu-a-1124.vcf.gz", chrpos_esm1)
#df format
cc_ovarian_c_esm1 <- as.data.frame(cc_ovarian_esm1[[1]])
#remove duplicated snps
cc_ovarian_c_esm1 <- cc_ovarian_c_esm1 %>% distinct(snp, .keep_all = TRUE)
cc_ovarian_c_esm1$sdY = 1
cc_ovarian_c_esm1$position <- cc_ovarian_c_esm1$pos
#plot
par(mfrow=c(2,1))
plot_dataset(esm1, main="Dataset ESM1")
plot_dataset(cc_ovarian_c_esm1, main="Dataset clear cell ovarian cancer")
#run colocalisation analysis
vres_esm1_cc_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,esm1, cc_ovarian_c_esm1)

#furin
#extract twice so that function works
ovarian_furin <- gwasvcf_to_coloc("ieu-a-1120.vcf.gz", "ieu-a-1120.vcf.gz", chrpos_furin)
#df format
ovarian_c_furin <- as.data.frame(ovarian_furin[[1]])
#remove duplicated snps
ovarian_c_furin <- ovarian_c_furin %>% distinct(snp, .keep_all = TRUE)
ovarian_c_furin$sdY = 1
ovarian_c_furin$position <- ovarian_c_furin$pos
#plot - function not working
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset FURIN")
plot_dataset(ovarian_c_furin, main="Dataset ovarian Carcinoma")
#run colocalisation analysis
vres_furin_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,furin, ovarian_c_furin)

#furin #invasive mucinous ovarian cancer 
#extract twice so that function works
im_ovarian_furin <- gwasvcf_to_coloc("ieu-a-1123.vcf.gz", "ieu-a-1123.vcf.gz", chrpos_furin)
#df format
im_ovarian_c_furin <- as.data.frame(im_ovarian_furin[[1]])
#remove duplicated snps
im_ovarian_c_furin <- im_ovarian_c_furin %>% distinct(snp, .keep_all = TRUE)
im_ovarian_c_furin$sdY = 1
im_ovarian_c_furin$position <- im_ovarian_c_furin$pos
#plot
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(im_ovarian_c_furin, main="Dataset invasive mucinous ovarian cancer")
#run colocalisation analysis
vres_furin_im_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,furin, im_ovarian_c_furin)

#furin #high grade serous ovarian cancer 
#extract twice so that function works
hg_ovarian_furin <- gwasvcf_to_coloc("ieu-a-1121.vcf.gz", "ieu-a-1121.vcf.gz", chrpos_furin)
#df format
hg_ovarian_c_furin <- as.data.frame(hg_ovarian_furin[[1]])
#remove duplicated snps
hg_ovarian_c_furin <- hg_ovarian_c_furin %>% distinct(snp, .keep_all = TRUE)
hg_ovarian_c_furin$sdY = 1
hg_ovarian_c_furin$position <- hg_ovarian_c_furin$pos
#plot
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(hg_ovarian_c_furin, main="Dataset high grade serous ovarian cancer")
#run colocalisation analysis
vres_furin_hg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,furin, hg_ovarian_c_furin)

#furin #low grade serous ovarian cancer 
#extract twice so that function works
lg_ovarian_furin <- gwasvcf_to_coloc("ieu-a-1122.vcf.gz", "ieu-a-1122.vcf.gz", chrpos_furin)
#df format
lg_ovarian_c_furin <- as.data.frame(lg_ovarian_furin[[1]])
#remove duplicated snps
lg_ovarian_c_furin <- lg_ovarian_c_furin %>% distinct(snp, .keep_all = TRUE)
lg_ovarian_c_furin$sdY = 1
lg_ovarian_c_furin$position <- lg_ovarian_c_furin$pos
#plot
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(lg_ovarian_c_furin, main="Dataset low grade serous ovarian cancer")
#run colocalisation analysis
vres_furin_lg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,furin, lg_ovarian_c_furin)

#furin #endometrioid ovarian cancer 
#extract twice so that function works
endo_ovarian_furin <- gwasvcf_to_coloc("ieu-a-1125.vcf.gz", "ieu-a-1125.vcf.gz", chrpos_furin)
#df format
endo_ovarian_c_furin <- as.data.frame(endo_ovarian_furin[[1]])
#remove duplicated snps
endo_ovarian_c_furin <- endo_ovarian_c_furin %>% distinct(snp, .keep_all = TRUE)
endo_ovarian_c_furin$sdY = 1
endo_ovarian_c_furin$position <- endo_ovarian_c_furin$pos
#plot
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(endo_ovarian_c_furin, main="Dataset endometrioid ovarian cancer")
#run colocalisation analysis
vres_furin_endo_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,furin, endo_ovarian_c_furin)

#furin #clear cell ovarian cancer 
#extract twice so that function works
cc_ovarian_furin <- gwasvcf_to_coloc("ieu-a-1124.vcf.gz", "ieu-a-1124.vcf.gz", chrpos_furin)
#df format
cc_ovarian_c_furin <- as.data.frame(cc_ovarian_furin[[1]])
#remove duplicated snps
cc_ovarian_c_furin <- cc_ovarian_c_furin %>% distinct(snp, .keep_all = TRUE)
cc_ovarian_c_furin$sdY = 1
cc_ovarian_c_furin$position <- cc_ovarian_c_furin$pos
#plot
par(mfrow=c(2,1))
plot_dataset(furin, main="Dataset furin")
plot_dataset(cc_ovarian_c_furin, main="Dataset clear cell ovarian cancer")
#run colocalisation analysis
vres_furin_cc_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,furin, cc_ovarian_c_furin)

#gpnmb
#extract twice so that function works
ovarian_gpnmb <- gwasvcf_to_coloc("ieu-a-1120.vcf.gz", "ieu-a-1120.vcf.gz", chrpos_gpnmb)
#df format
ovarian_c_gpnmb <- as.data.frame(ovarian_gpnmb[[1]])
#remove duplicated snps
ovarian_c_gpnmb <- ovarian_c_gpnmb %>% distinct(snp, .keep_all = TRUE)
ovarian_c_gpnmb$sdY = 1
ovarian_c_gpnmb$position <- ovarian_c_gpnmb$pos
#plot
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset GPNMB")
plot_dataset(ovarian_c_gpnmb, main="Dataset ovarian Carcinoma")
#run colocalisation analysis
vres_gpnmb_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,gpnmb, ovarian_c_gpnmb)

#gpnmb #invasive mucinous ovarian cancer 
#extract twice so that function works
im_ovarian_gpnmb <- gwasvcf_to_coloc("ieu-a-1123.vcf.gz", "ieu-a-1123.vcf.gz", chrpos_gpnmb)
#df format
im_ovarian_c_gpnmb <- as.data.frame(im_ovarian_gpnmb[[1]])
#remove duplicated snps
im_ovarian_c_gpnmb <- im_ovarian_c_gpnmb %>% distinct(snp, .keep_all = TRUE)
im_ovarian_c_gpnmb$sdY = 1
im_ovarian_c_gpnmb$position <- im_ovarian_c_gpnmb$pos
#plot
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(im_ovarian_c_gpnmb, main="Dataset invasive mucinous ovarian cancer")
#run colocalisation analysis
vres_gpnmb_im_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,gpnmb, im_ovarian_c_gpnmb)

#gpnmb #high grade serous ovarian cancer 
#extract twice so that function works
hg_ovarian_gpnmb <- gwasvcf_to_coloc("ieu-a-1121.vcf.gz", "ieu-a-1121.vcf.gz", chrpos_gpnmb)
#df format
hg_ovarian_c_gpnmb <- as.data.frame(hg_ovarian_gpnmb[[1]])
#remove duplicated snps
hg_ovarian_c_gpnmb <- hg_ovarian_c_gpnmb %>% distinct(snp, .keep_all = TRUE)
hg_ovarian_c_gpnmb$sdY = 1
hg_ovarian_c_gpnmb$position <- hg_ovarian_c_gpnmb$pos
#plot
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(hg_ovarian_c_gpnmb, main="Dataset high grade serous ovarian cancer")
#run colocalisation analysis
vres_gpnmb_hg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,gpnmb, hg_ovarian_c_gpnmb)

#gpnmb #low grade serous ovarian cancer 
#extract twice so that function works
lg_ovarian_gpnmb <- gwasvcf_to_coloc("ieu-a-1122.vcf.gz", "ieu-a-1122.vcf.gz", chrpos_gpnmb)
#df format
lg_ovarian_c_gpnmb <- as.data.frame(lg_ovarian_gpnmb[[1]])
#remove duplicated snps
lg_ovarian_c_gpnmb <- lg_ovarian_c_gpnmb %>% distinct(snp, .keep_all = TRUE)
lg_ovarian_c_gpnmb$sdY = 1
lg_ovarian_c_gpnmb$position <- lg_ovarian_c_gpnmb$pos
#plot
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(lg_ovarian_c_gpnmb, main="Dataset low grade serous ovarian cancer")
#run colocalisation analysis
vres_gpnmb_lg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,gpnmb, lg_ovarian_c_gpnmb)

#gpnmb #endometrioid ovarian cancer 
#extract twice so that function works
endo_ovarian_gpnmb <- gwasvcf_to_coloc("ieu-a-1125.vcf.gz", "ieu-a-1125.vcf.gz", chrpos_gpnmb)
#df format
endo_ovarian_c_gpnmb <- as.data.frame(endo_ovarian_gpnmb[[1]])
#remove duplicated snps
endo_ovarian_c_gpnmb <- endo_ovarian_c_gpnmb %>% distinct(snp, .keep_all = TRUE)
endo_ovarian_c_gpnmb$sdY = 1
endo_ovarian_c_gpnmb$position <- endo_ovarian_c_gpnmb$pos
#plot
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(endo_ovarian_c_gpnmb, main="Dataset endometrioid ovarian cancer")
#run colocalisation analysis
vres_gpnmb_endo_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,gpnmb, endo_ovarian_c_gpnmb)

#gpnmb #clear cell ovarian cancer 
#extract twice so that function works
cc_ovarian_gpnmb <- gwasvcf_to_coloc("ieu-a-1124.vcf.gz", "ieu-a-1124.vcf.gz", chrpos_gpnmb)
#df format
cc_ovarian_c_gpnmb <- as.data.frame(cc_ovarian_gpnmb[[1]])
#remove duplicated snps
cc_ovarian_c_gpnmb <- cc_ovarian_c_gpnmb %>% distinct(snp, .keep_all = TRUE)
cc_ovarian_c_gpnmb$sdY = 1
cc_ovarian_c_gpnmb$position <- cc_ovarian_c_gpnmb$pos
#plot
par(mfrow=c(2,1))
plot_dataset(gpnmb, main="Dataset gpnmb")
plot_dataset(cc_ovarian_c_gpnmb, main="Dataset clear cell ovarian cancer")
#run colocalisation analysis
vres_gpnmb_cc_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,gpnmb, cc_ovarian_c_gpnmb)

#hgf
#extract twice so that function works
ovarian_hgf <- gwasvcf_to_coloc("ieu-a-1120.vcf.gz", "ieu-a-1120.vcf.gz", chrpos_hgf)
#df format
ovarian_c_hgf <- as.data.frame(ovarian_hgf[[1]])
#remove duplicated snps
ovarian_c_hgf <- ovarian_c_hgf %>% distinct(snp, .keep_all = TRUE)
ovarian_c_hgf$sdY = 1
ovarian_c_hgf$position <- ovarian_c_hgf$pos
#plot
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset HGF")
plot_dataset(ovarian_c_hgf, main="Dataset ovarian Carcinoma")
#run colocalisation analysis
vres_hgf_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,hgf, ovarian_c_hgf)

#hgf #invasive mucinous ovarian cancer 
#extract twice so that function works
im_ovarian_hgf <- gwasvcf_to_coloc("ieu-a-1123.vcf.gz", "ieu-a-1123.vcf.gz", chrpos_hgf)
#df format
im_ovarian_c_hgf <- as.data.frame(im_ovarian_hgf[[1]])
#remove duplicated snps
im_ovarian_c_hgf <- im_ovarian_c_hgf %>% distinct(snp, .keep_all = TRUE)
im_ovarian_c_hgf$sdY = 1
im_ovarian_c_hgf$position <- im_ovarian_c_hgf$pos
#plot
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(im_ovarian_c_hgf, main="Dataset invasive mucinous ovarian cancer")
#run colocalisation analysis
vres_hgf_im_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,hgf, im_ovarian_c_hgf)

#hgf #high grade serous ovarian cancer 
#extract twice so that function works
hg_ovarian_hgf <- gwasvcf_to_coloc("ieu-a-1121.vcf.gz", "ieu-a-1121.vcf.gz", chrpos_hgf)
#df format
hg_ovarian_c_hgf <- as.data.frame(hg_ovarian_hgf[[1]])
#remove duplicated snps
hg_ovarian_c_hgf <- hg_ovarian_c_hgf %>% distinct(snp, .keep_all = TRUE)
hg_ovarian_c_hgf$sdY = 1
hg_ovarian_c_hgf$position <- hg_ovarian_c_hgf$pos
#plot
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(hg_ovarian_c_hgf, main="Dataset high grade serous ovarian cancer")
#run colocalisation analysis
vres_hgf_hg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,hgf, hg_ovarian_c_hgf)

#hgf #low grade serous ovarian cancer 
#extract twice so that function works
lg_ovarian_hgf <- gwasvcf_to_coloc("ieu-a-1122.vcf.gz", "ieu-a-1122.vcf.gz", chrpos_hgf)
#df format
lg_ovarian_c_hgf <- as.data.frame(lg_ovarian_hgf[[1]])
#remove duplicated snps
lg_ovarian_c_hgf <- lg_ovarian_c_hgf %>% distinct(snp, .keep_all = TRUE)
lg_ovarian_c_hgf$sdY = 1
lg_ovarian_c_hgf$position <- lg_ovarian_c_hgf$pos
#plot
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(lg_ovarian_c_hgf, main="Dataset low grade serous ovarian cancer")
#run colocalisation analysis
vres_hgf_lg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,hgf, lg_ovarian_c_hgf)

#hgf #endometrioid ovarian cancer 
#extract twice so that function works
endo_ovarian_hgf <- gwasvcf_to_coloc("ieu-a-1125.vcf.gz", "ieu-a-1125.vcf.gz", chrpos_hgf)
#df format
endo_ovarian_c_hgf <- as.data.frame(endo_ovarian_hgf[[1]])
#remove duplicated snps
endo_ovarian_c_hgf <- endo_ovarian_c_hgf %>% distinct(snp, .keep_all = TRUE)
endo_ovarian_c_hgf$sdY = 1
endo_ovarian_c_hgf$position <- endo_ovarian_c_hgf$pos
#plot
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(endo_ovarian_c_hgf, main="Dataset endometrioid ovarian cancer")
#run colocalisation analysis
vres_hgf_endo_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,hgf, endo_ovarian_c_hgf)

#hgf #clear cell ovarian cancer 
#extract twice so that function works
cc_ovarian_hgf <- gwasvcf_to_coloc("ieu-a-1124.vcf.gz", "ieu-a-1124.vcf.gz", chrpos_hgf)
#df format
cc_ovarian_c_hgf <- as.data.frame(cc_ovarian_hgf[[1]])
#remove duplicated snps
cc_ovarian_c_hgf <- cc_ovarian_c_hgf %>% distinct(snp, .keep_all = TRUE)
cc_ovarian_c_hgf$sdY = 1
cc_ovarian_c_hgf$position <- cc_ovarian_c_hgf$pos
#plot
par(mfrow=c(2,1))
plot_dataset(hgf, main="Dataset hgf")
plot_dataset(cc_ovarian_c_hgf, main="Dataset clear cell ovarian cancer")
#run colocalisation analysis
vres_hgf_cc_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,hgf, cc_ovarian_c_hgf)

#ret
#extract twice so that function works
ovarian_ret <- gwasvcf_to_coloc("ieu-a-1120.vcf.gz", "ieu-a-1120.vcf.gz", chrpos_ret)
#df format
ovarian_c_ret <- as.data.frame(ovarian_ret[[1]])
#remove duplicated snps
ovarian_c_ret <- ovarian_c_ret %>% distinct(snp, .keep_all = TRUE)
ovarian_c_ret$sdY = 1
ovarian_c_ret$position <- ovarian_c_ret$pos
#plot
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset RET")
plot_dataset(ovarian_c_ret, main="Dataset ovarian Carcinoma")
#run colocalisation analysis 
vres_ret_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,ret, ovarian_c_ret)

#ret #invasive mucinous ovarian cancer 
#extract twice so that function works
im_ovarian_ret <- gwasvcf_to_coloc("ieu-a-1123.vcf.gz", "ieu-a-1123.vcf.gz", chrpos_ret)
#df format
im_ovarian_c_ret <- as.data.frame(im_ovarian_ret[[1]])
#remove duplicated snps
im_ovarian_c_ret <- im_ovarian_c_ret %>% distinct(snp, .keep_all = TRUE)
im_ovarian_c_ret$sdY = 1
im_ovarian_c_ret$position <- im_ovarian_c_ret$pos
#plot
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(im_ovarian_c_ret, main="Dataset invasive mucinous ovarian cancer")
#run colocalisation analysis
vres_ret_im_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,ret, im_ovarian_c_ret)

#ret #high grade serous ovarian cancer 
#extract twice so that function works
hg_ovarian_ret <- gwasvcf_to_coloc("ieu-a-1121.vcf.gz", "ieu-a-1121.vcf.gz", chrpos_ret)
#df format
hg_ovarian_c_ret <- as.data.frame(hg_ovarian_ret[[1]])
#remove duplicated snps
hg_ovarian_c_ret <- hg_ovarian_c_ret %>% distinct(snp, .keep_all = TRUE)
hg_ovarian_c_ret$sdY = 1
hg_ovarian_c_ret$position <- hg_ovarian_c_ret$pos
#plot
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(hg_ovarian_c_ret, main="Dataset high grade serous ovarian cancer")
#run colocalisation analysis
vres_ret_hg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,ret, hg_ovarian_c_ret)

#ret #low grade serous ovarian cancer 
#extract twice so that function works
lg_ovarian_ret <- gwasvcf_to_coloc("ieu-a-1122.vcf.gz", "ieu-a-1122.vcf.gz", chrpos_ret)
#df format
lg_ovarian_c_ret <- as.data.frame(lg_ovarian_ret[[1]])
#remove duplicated snps
lg_ovarian_c_ret <- lg_ovarian_c_ret %>% distinct(snp, .keep_all = TRUE)
lg_ovarian_c_ret$sdY = 1
lg_ovarian_c_ret$position <- lg_ovarian_c_ret$pos
#plot
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(lg_ovarian_c_ret, main="Dataset low grade serous ovarian cancer")
#run colocalisation analysis
vres_ret_lg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,ret, lg_ovarian_c_ret)

#ret #endometrioid ovarian cancer 
#extract twice so that function works
endo_ovarian_ret <- gwasvcf_to_coloc("ieu-a-1125.vcf.gz", "ieu-a-1125.vcf.gz", chrpos_ret)
#df format
endo_ovarian_c_ret <- as.data.frame(endo_ovarian_ret[[1]])
#remove duplicated snps
endo_ovarian_c_ret <- endo_ovarian_c_ret %>% distinct(snp, .keep_all = TRUE)
endo_ovarian_c_ret$sdY = 1
endo_ovarian_c_ret$position <- endo_ovarian_c_ret$pos
#plot
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(endo_ovarian_c_ret, main="Dataset endometrioid ovarian cancer")
#run colocalisation analysis
vres_ret_endo_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,ret, endo_ovarian_c_ret)

#ret #clear cell ovarian cancer 
#extract twice so that function works
cc_ovarian_ret <- gwasvcf_to_coloc("ieu-a-1124.vcf.gz", "ieu-a-1124.vcf.gz", chrpos_ret)
#df format
cc_ovarian_c_ret <- as.data.frame(cc_ovarian_ret[[1]])
#remove duplicated snps
cc_ovarian_c_ret <- cc_ovarian_c_ret %>% distinct(snp, .keep_all = TRUE)
cc_ovarian_c_ret$sdY = 1
cc_ovarian_c_ret$position <- cc_ovarian_c_ret$pos
#plot
par(mfrow=c(2,1))
plot_dataset(ret, main="Dataset ret")
plot_dataset(cc_ovarian_c_ret, main="Dataset clear cell ovarian cancer")
#run colocalisation analysis
vres_ret_cc_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,ret, cc_ovarian_c_ret)

#tlr3
#extract twice so that function works
ovarian_tlr3 <- gwasvcf_to_coloc("ieu-a-1120.vcf.gz", "ieu-a-1120.vcf.gz", chrpos_tlr3)
#df format
ovarian_c_tlr3 <- as.data.frame(ovarian_tlr3[[1]])
#remove duplicated snps
ovarian_c_tlr3 <- ovarian_c_tlr3 %>% distinct(snp, .keep_all = TRUE)
ovarian_c_tlr3$sdY = 1
ovarian_c_tlr3$position <- ovarian_c_tlr3$pos
#plot
par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset TLR3")
plot_dataset(ovarian_c_tlr3, main="Dataset ovarian Carcinoma")
#run colocalisation analysis 
vres_tlr3_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,tlr3, ovarian_c_tlr3)

#tlr3 #invasive mucinous ovarian cancer 
#extract twice so that function works
im_ovarian_tlr3 <- gwasvcf_to_coloc("ieu-a-1123.vcf.gz", "ieu-a-1123.vcf.gz", chrpos_tlr3)
#df format
im_ovarian_c_tlr3 <- as.data.frame(im_ovarian_tlr3[[1]])
#remove duplicated snps
im_ovarian_c_tlr3 <- im_ovarian_c_tlr3 %>% distinct(snp, .keep_all = TRUE)
im_ovarian_c_tlr3$sdY = 1
im_ovarian_c_tlr3$position <- im_ovarian_c_tlr3$pos
#plot
par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(im_ovarian_c_tlr3, main="Dataset invasive mucinous ovarian cancer")
#run colocalisation analysis
vres_tlr3_im_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,tlr3, im_ovarian_c_tlr3)
print("here1")
#tlr3 #high grade serous ovarian cancer 
#extract twice so that function works
hg_ovarian_tlr3 <- gwasvcf_to_coloc("ieu-a-1121.vcf.gz", "ieu-a-1121.vcf.gz", chrpos_tlr3)
#df format
hg_ovarian_c_tlr3 <- as.data.frame(hg_ovarian_tlr3[[1]])
#remove duplicated snps
hg_ovarian_c_tlr3 <- hg_ovarian_c_tlr3 %>% distinct(snp, .keep_all = TRUE)
hg_ovarian_c_tlr3$sdY = 1
hg_ovarian_c_tlr3$position <- hg_ovarian_c_tlr3$pos
#plot
par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(hg_ovarian_c_tlr3, main="Dataset high grade serous ovarian cancer")
#run colocalisation analysis
vres_tlr3_hg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,tlr3, hg_ovarian_c_tlr3)

#tlr3 #low grade serous ovarian cancer 
#extract twice so that function works
lg_ovarian_tlr3 <- gwasvcf_to_coloc("ieu-a-1122.vcf.gz", "ieu-a-1122.vcf.gz", chrpos_tlr3)
#df format
lg_ovarian_c_tlr3 <- as.data.frame(lg_ovarian_tlr3[[1]])
#remove duplicated snps
lg_ovarian_c_tlr3 <- lg_ovarian_c_tlr3 %>% distinct(snp, .keep_all = TRUE)
lg_ovarian_c_tlr3$sdY = 1
lg_ovarian_c_tlr3$position <- lg_ovarian_c_tlr3$pos
#plot
par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(lg_ovarian_c_tlr3, main="Dataset low grade serous ovarian cancer")
#run colocalisation analysis
vres_tlr3_lg_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,tlr3, lg_ovarian_c_tlr3)

#tlr3 #endometrioid ovarian cancer 
#extract twice so that function works
endo_ovarian_tlr3 <- gwasvcf_to_coloc("ieu-a-1125.vcf.gz", "ieu-a-1125.vcf.gz", chrpos_tlr3)
#df format
endo_ovarian_c_tlr3 <- as.data.frame(endo_ovarian_tlr3[[1]])
#remove duplicated snps
endo_ovarian_c_tlr3 <- endo_ovarian_c_tlr3 %>% distinct(snp, .keep_all = TRUE)
endo_ovarian_c_tlr3$sdY = 1
endo_ovarian_c_tlr3$position <- endo_ovarian_c_tlr3$pos
#plot
par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(endo_ovarian_c_tlr3, main="Dataset endometrioid ovarian cancer")
#run colocalisation analysis
vres_tlr3_endo_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,tlr3, endo_ovarian_c_tlr3)

#tlr3 #clear cell ovarian cancer 
#extract twice so that function works
cc_ovarian_tlr3 <- gwasvcf_to_coloc("ieu-a-1124.vcf.gz", "ieu-a-1124.vcf.gz", chrpos_tlr3)
#df format
cc_ovarian_c_tlr3 <- as.data.frame(cc_ovarian_tlr3[[1]])
#remove duplicated snps
cc_ovarian_c_tlr3 <- cc_ovarian_c_tlr3 %>% distinct(snp, .keep_all = TRUE)
cc_ovarian_c_tlr3$sdY = 1
cc_ovarian_c_tlr3$position <- cc_ovarian_c_tlr3$pos
#plot
par(mfrow=c(2,1))
plot_dataset(tlr3, main="Dataset tlr3")
plot_dataset(cc_ovarian_c_tlr3, main="Dataset clear cell ovarian cancer")
#run colocalisation analysis
vres_tlr3_cc_oca <- coloc::coloc.abf(p1=1E-3,p2=1E-4,p12=1E-5,tlr3, cc_ovarian_c_tlr3)

#BIG LIST OF RESULTS TO EXPORT
dfs <- lapply(ls(pattern="^vres"), function(x) get(x))
summaries <- sapply(dfs, "[[", 1)
colheaders <- ls(pattern = "vres_")
colnames(summaries) <- colheaders
write.table(summaries, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t", file = "~/coloc_results201022.txt")



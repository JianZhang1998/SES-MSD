library(forestplot)
library(haven)
library(TwoSampleMR)
library(MendelianRandomization)
library(MRPRESSO)
library(RadialMR)
library(foreign)
library(rms)
library(tidyverse)
library(readr)
library(RMediation)
library(forestploter)
library(haven)
library(dplyr)


###############edu-bmi two-stage######

ea<- read_delim(file = "E:/07 UKB/23 SES mediator and MSK/00 data/01 education gwas/2018-NG-41.38-educational attainment gwas2-GWAS_EA_excl23andMe.txt",delim="\t")


head(ea)

ea<- ea%>% rename(SNP=MarkerName, chr= CHR,  pos=POS , BETA=Beta,  P =Pval) %>%   # 重命名列
  mutate(EXPOSURE = "education") 

ea1 <- subset(ea, P <= 5e-08)

exp_ea=format_data(ea1,
                   snps = NULL,
                   header = TRUE,
                   snp_col = "SNP",
                   beta_col ="BETA",
                   se_col ="SE",
                   eaf_col = "EAF",
                   effect_allele_col = "A1",
                   other_allele_col = "A2",
                   pval_col = "P")

exposure_dat=clump_data(
  exp_ea,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR")


#mediator:chronotype
chrono<- read_delim(file = "E:/07 UKB/20 sleep and sp/00 data/00 sleep gwas/summary statistics/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt",delim="\t")

head(chrono)

chrono<- chrono%>% rename(chr= CHR,  pos=BP, A1=ALLELE1, A2=ALLELE0, EAF=A1FREQ,P =P_BOLT_LMM) %>%   # 重命名列
  mutate(EXPOSURE = "chronotype")

chrono$HWE_P <- NULL
chrono$INFO <- NULL

chrono1 <- subset(chrono, P <= 5e-08)

exp_chrono=format_data(chrono1,
                       snps = NULL,
                       header = TRUE,
                       snp_col = "SNP",
                       beta_col ="BETA",
                       se_col ="SE",
                       eaf_col = "EAF",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P")

mediator_dat=clump_data(
  exp_chrono,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR")

#mediator：BMI;ieu-b-40  smoking initiation：ieu-b-4877, Cigarettes per Day ieu-b-25 vpa:ebi-a-GCST006098  mvpa:ebi-a-GCST006097

mediator_dat <- extract_instruments(outcomes ="ebi-a-GCST006097", p1 =5e-08, clump = TRUE,p2 =5e-08, r2 = 0.001, kb = 10000)

#去除暴露中介相同snp
# 获取暴露和中介SNP的交集
common_snps <- intersect(exposure_dat$SNP, mediator_dat$SNP)
mediator_dat <- mediator_dat[!mediator_dat$SNP %in% common_snps, ]
exposure_dat <- exposure_dat[!exposure_dat$SNP %in% common_snps, ]

#暴露-中介——效应值A
outcome_ex_med <- extract_outcome_data(exposure_dat$SNP, outcomes = "ieu-b-25") 

data_ex_med <- harmonise_data(exposure_dat, outcome_ex_med)
res_ex_med <- mr(data_ex_med)
res_ex_med <- generate_odds_ratios(res_ex_med)


#中介-结局——效应值B
outcome_med_outcome <- extract_outcome_data(mediator_dat$SNP, outcomes = "ukb-b-4711") 

pain<-read_delim(file = 'E:/07 UKB/23 SES mediator and MSK/00 data/Chronic widespread musculoskeletal pain-2021-33926923-1.txt')

head(pain)

outcome_med_outcome=format_data(pain,
                                type='outcome',
                                snp_col = "SNP",
                                beta_col ="BETA",
                                se_col ="SE",
                                eaf_col = "A1FREQ",
                                effect_allele_col = "A1",
                                other_allele_col = "A2",
                                pval_col = "p")

data_med_outcome <- harmonise_data( mediator_dat, outcome_med_outcome)
res_med_outcome <- mr(data_med_outcome)
res_med_outcome <- generate_odds_ratios(res_med_outcome)

outcomes <- c("ebi-a-GCST007092", "ieu-a-833", "ukb-b-10215","ukb-b-4711","ebi-a-GCST90018797","ukb-b-18596")

###########
Exposures <- data.frame()
Outcomes <- data.frame()
Exposure_Outcome_Results <- data.frame()
Exposure_Mediator_Heterogeneity <- data.frame()
Exposure_Mediator_Horizontal_pleiotropy <- data.frame()
st <- Sys.time()

if (sum(exposure_dat$mr_keep.exposure)!=0){
  for (j in seq_along(outcomes)) {
    ed <- Sys.time()
    if (difftime(ed,st,units = 'secs')>150) {
      print(ed) # 调试用
      Sys.sleep(20)
      st <- Sys.time()
      print(st) # 调试用
    }
    outcome_dat <- TwoSampleMR::extract_outcome_data(snps = mediator_dat$SNP, outcomes = outcomes[j])
    outcome_info <- ao[which(ao$id==outcomes[j]),c('id','year','pmid','population','sample_size','build')]
    if (is.null(dim(outcome_dat))!=TRUE) {
      dat <- harmonise_data(exposure_dat = mediator_dat, outcome_dat = outcome_dat) %>% add_rsq() # 整合数据
      if (sum(mediator_dat$mr_keep.exposure)<2){
        res <- TwoSampleMR::mr(dat, method_list = "mr_wald_ratio") %>% generate_odds_ratios()
      } else if (sum(exposure_dat$mr_keep.exposure)==2) {
        res <- TwoSampleMR::mr(dat = dat, method_list = c('mr_ivw')) %>% generate_odds_ratios() # MR分析
      } else if (sum(exposure_dat$mr_keep.exposure)==3) {
        res <- TwoSampleMR::mr(dat = dat, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median')) %>% generate_odds_ratios() # MR分析
        # res_loo <- TwoSampleMR::mr_leaveoneout(dat)
      } else if (sum(exposure_dat$mr_keep.exposure)>3) {
        res <- TwoSampleMR::mr(dat = dat, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median')) %>% generate_odds_ratios() # MR分析
        # res_loo <- TwoSampleMR::mr_leaveoneout(dat)
        presso_res <- mr_presso(data = dat,BetaOutcome = "beta.outcome",BetaExposure = "beta.exposure",SdOutcome = "se.outcome",
        SdExposure = "se.exposure",OUTLIERtest = TRUE,DISTORTIONtest = FALSE,SignifThreshold = 0.05)
        print(presso_res)
      }
      het <- mr_heterogeneity(dat)
      horizontal_pleiotropy <- mr_pleiotropy_test(dat) # 敏感性分析
      
      # 保存结果
      Exposures <- rbind(Exposures)
      Outcomes <- rbind(Outcomes, outcome_info)
      Exposure_Outcome_Results <- rbind(Exposure_Outcome_Results, res)
      Exposure_Mediator_Heterogeneity <- rbind(Exposure_Mediator_Heterogeneity,het)
      Exposure_Mediator_Horizontal_pleiotropy <- rbind(Exposure_Mediator_Horizontal_pleiotropy, horizontal_pleiotropy)
      
      # 保存后清除缓存结果
      rm(list = 'res','dat','het','horizontal_pleiotropy')
    }
  }
}
# df_Exposure_Mediator <- merge(Exposure_Mediator_Results, Exposure_Mediator_Heterogeneity, by = # c('id.exposure','id.outcome','method','exposure','outcome'),all = TRUE) %>% merge(Exposure_Mediator_Horizontal_pleiotropy, by = # c('id.exposure','id.outcome','exposure','outcome'),all = TRUE)
df_Exposure_Outcome <- Exposure_Outcome_Results
##########
#计算直接效应
outcome_data <- extract_outcome_data(exposure_dat$SNP, outcomes = "ukb-b-18596") 

outcome_med_outcome=format_data(pain,
                                type='outcome',
                                snp_col = "SNP",
                                beta_col ="BETA",
                                se_col ="SE",
                                eaf_col = "A1FREQ",
                                effect_allele_col = "A1",
                                other_allele_col = "A2",
                                pval_col = "p")

data <- harmonise_data( exposure_dat, outcome_med_outcome)
res <- mr(data)
res <- generate_odds_ratios(res)

MV_Res <- rbind(res_ex_med, res_med_outcome,res)

MV_Res <- rbind(res_med_outcome)
write.csv(MV_Res,file="MV_Res_ea_mvpa_speed_twosatge.csv")


######income-chrono-pain/widepain######

exposure_dat<-extract_instruments(outcomes ="ukb-b-7408", p1 =5e-08, clump = TRUE,p2 =5e-08, r2 = 0.001, kb = 10000)

#mediator：BMI; cig:ieu-b-25 vpa:ebi-a-GCST006098

mediator_dat <- extract_instruments(outcomes ="ebi-a-GCST006097", p1 =5e-08, clump = TRUE,p2 =5e-08, r2 = 0.001, kb = 10000)

#去除暴露中介相同snp
# 获取暴露和中介SNP的交集
common_snps <- intersect(exposure_dat$SNP, mediator_dat$SNP)
mediator_dat <- mediator_dat[!mediator_dat$SNP %in% common_snps, ]
exposure_dat <- exposure_dat[!exposure_dat$SNP %in% common_snps, ]

#暴露-中介——效应值A

outcome_ex_med <- extract_outcome_data(exposure_dat$SNP, outcomes = "ieu-b-25") 

data_ex_med <- harmonise_data(exposure_dat, outcome_ex_med)
res_ex_med <- mr(data_ex_med)
res_ex_med <- generate_odds_ratios(res_ex_med)
#中介-结局——效应值B
outcomes <- c("ebi-a-GCST90000025","ukb-b-10215","ukb-b-4711","ukb-b-18596")

outcome_med_outcome <- extract_outcome_data(mediator_dat$SNP, outcomes = "ukb-b-4711") 

outcome_med_outcome=format_data(pain,
                                type='outcome',
                                snp_col = "SNP",
                                beta_col ="BETA",
                                se_col ="SE",
                                eaf_col = "A1FREQ",
                                effect_allele_col = "A1",
                                other_allele_col = "A2",
                                pval_col = "p")

data_med_outcome <- harmonise_data( mediator_dat, outcome_med_outcome  )
res_med_outcome <- mr(data_med_outcome)
res_med_outcome <- generate_odds_ratios(res_med_outcome)


#计算直接效应
outcome_data <- extract_outcome_data(exposure_dat$SNP, outcomes = "ukb-b-4711") 

data <- harmonise_data( exposure_dat, outcome_med_outcome  )
res <- mr(data)
res <- generate_odds_ratios(res)

MV_Res <- rbind(res_ex_med, res_med_outcome)

write.csv(MV_Res,file="MV_Res_income_cig_alm_twosatge.csv")

###############occup-bmi two-stage######

# occupation
occup<- read_delim(file = "E:/07 UKB/23 SES mediator and MSK/00 data/01 education gwas/2022-brain-occupational attainment gwas.tsv",delim="\t")

occup1 <- subset(occup, p_value <= 5e-08)

exp_occup=format_data(occup1,
                      snps = NULL,
                      header = TRUE,
                      snp_col = "variant_id",
                      beta_col ="beta",
                      se_col ="standard_error",
                      eaf_col = NULL,
                      effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele",
                      pval_col = "p_value")

exposure_dat=clump_data(
  exp_occup,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR")

#mediator：BMI;

mediator_dat <- extract_instruments(outcomes ="ieu-b-40", p1 =5e-08, clump = TRUE,p2 =5e-08, r2 = 0.001, kb = 10000)

#去除暴露中介相同snp
# 获取暴露和中介SNP的交集
common_snps <- intersect(exposure_dat$SNP, mediator_dat$SNP)
mediator_dat <- mediator_dat[!mediator_dat$SNP %in% common_snps, ]
exposure_dat <- exposure_dat[!exposure_dat$SNP %in% common_snps, ]

#暴露-中介——效应值A

outcome_ex_med <- extract_outcome_data(exposure_dat$SNP, outcomes = "ieu-b-40") 

data_ex_med <- harmonise_data(exposure_dat, outcome_ex_med)
res_ex_med <- mr(data_ex_med)
res_ex_med <- generate_odds_ratios(res_ex_med)

#中介-结局——效应值B

outcome_med_outcome <- extract_outcome_data(mediator_dat$SNP, outcomes = "ebi-a-GCST90000025") 

data_med_outcome <- harmonise_data(mediator_dat, outcome_med_outcome)
res_med_outcome <- mr(data_med_outcome)
res_med_outcome <- generate_odds_ratios(res_med_outcome)


#计算直接效应
outcome_data <- extract_outcome_data(exposure_dat$SNP, outcomes = "ebi-a-GCST90000025") 

data <- harmonise_data( exposure_dat, outcome_data)
res <- mr(data)
res <- generate_odds_ratios(res)

MV_Res <- rbind(res_ex_med, res_med_outcome,res)

write.csv(MV_Res,file="MV_Res_occup_BMI_alm_twosatge.csv")

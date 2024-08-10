pacman::p_load(tidyverse,
               ggplot2,
               dplyr,
               lmerTest,
               moderndive,
               performance,
               GGally,
               glmmTMB,
               DHARMa,
               stats,
               car,
               redres,
               lmtest,
               nlme)
pacman::p_load(ggcorrplot, inspect, matrix)
pacman::p_load(nlme,
               tidyverse,
               ggplot2,
               lmerTest,
               haven,
               RColorBrewer)

df.raw <- read_tsv("MODS_Rohrlach_EAGER.tsv")
df.raw <- df.raw %>%
  filter(nr_of_input_reads>1e6)
view(df.raw)

#### Process----
df.raw %>% arrange(nr_of_input_reads)

dtf <- df.raw %>%
  mutate(prop_1240k=covered_sn_ps_on_1240k/total_sn_ps_on_1240k,
         percent_endogenous_dna_over_30bp=percent_endogenous_dna_over_30bp/100)

norm <- function(x) {
  return ((x-mean(x))/sd(x))
}

dtf <- dtf %>%
  select(-total_sn_ps_on_1240k,
         -nr_of_input_reads_over_30bp,
         -nr_of_mapped_reads,
         -percent_endogenous_dna,
         -mean_mt_coverage,
         -median_fold_coverage,
         -nr_mt_dna_reads,
         -nr_nuclear_dna_reads,
         -nuclear_contamination_m1_mom,
         -mean_fold_coverage_on_nuclear_genome,
         -nr_of_reads_total)
dtf <- dtf %>%
  select(!contains('genome_covered')) %>%
  select(!contains('damage')) %>%
  select(!contains('error')) %>%
  select(!contains('m2')) %>%
  select(!contains('chromosome')) %>%
  select(!contains('length')) %>%
  select(!contains('relative_coverage')) %>%
  select(!contains('length')) %>%
  select(!contains('length'))

dtf <- dtf %>%
  mutate(SNPs_1240k_N=norm(covered_sn_ps_on_1240k),
         input_reads_N=norm(nr_of_input_reads),
         mapped_reads_N=norm(nr_of_mapped_reads_over_30bp))

dtf <- dtf %>% rename(SNPs_1240k=covered_sn_ps_on_1240k,
                      input_reads=nr_of_input_reads,
                      mapped_reads=nr_of_mapped_reads_over_30bp,
                      endogenous=percent_endogenous_dna_over_30bp,
                      prop_dup=proportion_of_duplicate_reads,
                      mtNuc_ratio=mt_to_nuclear_read_ratio,
                      unique_reads=nr_of_unique_mapped_reads,
                      coverage=mean_fold_coverage,
                      GC=percent_gc_of_unique_reads,
                      SNPs_cont=nr_sn_ps_used_in_contamination_estimation,
                      contam=nuclear_contamination_m1_ml) %>%
  mutate(source=factor(source,levels=c('Twist','1240k')))

df <- dtf %>% group_by(sample)%>%
  filter(n()==2)%>%
  ungroup()%>%
  mutate(sample=factor(sample))

df %>%
  names()%>%
  paste0(c('', paste0(1:11,":"),rep(' ',ncol(dtf)-12)),.)%>%
  paste0(collapse='\n')%>%
  paste0('\n')%>%
  cat()
view(df)

df %>%
  select(-sample)%>%
  select(!contains('_N'))%>%
  select(!contains('prop_1240k'))%>%
  ggpairs(aes(colour = source))


library(lme4)
library(glmmTMB)
library(car)
library(lmtest)
library(dplyr)
library(ggplot2)

#### input reads vs SNPs_1240k----
summary(df$SNPs_1240k)

df %>%
  ggplot(aes(x=input_reads_N,
             y=SNPs_1240k_N,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE, method = "lm")

snps.interaction <- lmer(SNPs_1240k_N ~ input_reads_N*source + (1|sample), data = df)
summary(snps.interaction)

snps.add <- lmer(SNPs_1240k_N ~ input_reads_N + source + (1|sample), data = df)
summary(snps.add)
Anova(snps.interaction) 
waldtest(snps.interaction, snps.add)

snps.input <- lmer(SNPs_1240k_N ~ input_reads_N + (1|sample), data = df)
summary(snps.input)

snps.source <- lmer(SNPs_1240k_N ~ source + (1|sample), data = df)
summary(snps.source)
waldtest(snps.add, snps.input)
waldtest(snps.add, snps.source)
summary(snps.add)

par(mfrow=c(2,2))
plot(snps.add, which = 1)
plot(snps.add, which = 2)
plot(snps.add, which = 3)
plot(snps.add, which = 5)

#### input reads vs mapped_reads----
df %>%
  ggplot(aes(x=input_reads_N,
             y=mapped_reads_N,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

mapped_reads.interaction <- lmer(mapped_reads_N ~ input_reads_N*source + (1|sample), data = df)
summary(mapped_reads.interaction)

mapped_reads.add <- lmer(mapped_reads_N ~ input_reads_N + source + (1|sample), data = df)
summary(mapped_reads.add)
waldtest(mapped_reads.interaction, mapped_reads.add)

mapped_reads.input <- lmer(mapped_reads_N ~ input_reads_N + (1|sample), data = df)
summary(mapped_reads.input)

mapped_reads.source <- lmer(mapped_reads_N ~ source + (1|sample), data = df)
summary(mapped_reads.source)
waldtest(mapped_reads.add, mapped_reads.input)
waldtest(mapped_reads.add, mapped_reads.source)
summary(mapped_reads.add)

par(mfrow=c(2,2))
plot(mapped_reads.add, which = 1)
plot(mapped_reads.add, which = 2)
plot(mapped_reads.add, which = 3)
plot(mapped_reads.add, which = 5)

#### input reads vs endogenous----
df %>%
  ggplot(aes(x=input_reads_N,
             y=endogenous,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

endogenous.interaction <- glmmTMB(endogenous ~ input_reads_N*source + (1|sample), 
                                  data = df, family = beta_family())
summary(endogenous.interaction)

endogenous.add <- glmmTMB(endogenous ~ input_reads_N + source + (1|sample), 
                          data = df, family = beta_family())
summary(endogenous.add)
waldtest(endogenous.interaction, endogenous.add)

endogenous.source <- glmmTMB(endogenous ~ source + (1|sample), 
                             data = df, family = beta_family())
summary(endogenous.source)

endogenous.input <- glmmTMB(endogenous ~ input_reads_N + (1|sample), 
                            data = df, family = beta_family())
summary(endogenous.input)
waldtest(endogenous.add, endogenous.source)
waldtest(endogenous.add, endogenous.input)

par(mfrow=c(2,2))
plot(endogenous.source, which = 1)
plot(endogenous.source, which = 2)
plot(endogenous.source, which = 3)
plot(endogenous.source, which = 5)

#### input reads vs prop_dup----
df %>%
  ggplot(aes(x=input_reads_N,
             y=prop_dup,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

prop_dup.interaction <- glmmTMB(prop_dup ~ input_reads_N*source + (1|sample), 
                                data = df, family = beta_family())
summary(prop_dup.interaction)

prop_dup.add <- glmmTMB(prop_dup ~ input_reads_N + source + (1|sample), 
                        data = df, family = beta_family())
summary(prop_dup.add)
waldtest(prop_dup.interaction, prop_dup.add)

prop_dup.source <- glmmTMB(prop_dup ~ source + (1|sample), 
                           data = df, family = beta_family())
summary(prop_dup.source)

prop_dup.input <- glmmTMB(prop_dup ~ input_reads_N + (1|sample), 
                          data = df, family = beta_family())
summary(prop_dup.input)
waldtest(prop_dup.add, prop_dup.source)
waldtest(prop_dup.add, prop_dup.input)

par(mfrow=c(2,2))
plot(prop_dup.add, which = 1)
plot(prop_dup.add, which = 2)
plot(prop_dup.add, which = 3)
plot(prop_dup.add, which = 5)

#### input reads vs mtNuc_ratio----
df %>%
  ggplot(aes(x=input_reads_N,
             y=mtNuc_ratio,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

mtNuc_ratio.interaction <- lmer(mtNuc_ratio ~ input_reads_N*source + (1|sample), data = df)
summary(mtNuc_ratio.interaction)

mtNuc_ratio.add <- lmer(mtNuc_ratio ~ input_reads_N + source + (1|sample), data = df)
summary(mtNuc_ratio.add)
waldtest(mtNuc_ratio.interaction, mtNuc_ratio.add)

mtNuc_ratio.source <- lmer(mtNuc_ratio ~ source + (1|sample), data = df)
summary(mtNuc_ratio.source)

mtNuc_ratio.input <- lmer(mtNuc_ratio ~ input_reads_N + (1|sample), data = df)
summary(mtNuc_ratio.input)
waldtest(mtNuc_ratio.add, mtNuc_ratio.source)
waldtest(mtNuc_ratio.add, mtNuc_ratio.input)

par(mfrow=c(2,2))
plot(mtNuc_ratio.add, which = 1)
plot(mtNuc_ratio.add, which = 2)
plot(mtNuc_ratio.add, which = 3)
plot(mtNuc_ratio.add, which = 5)

#### input reads vs unique_reads----
df %>%
  ggplot(aes(x=input_reads_N,
             y=unique_reads,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

unique_reads.interaction <- lmer(unique_reads ~ input_reads_N*source + (1|sample), data = df)
summary(unique_reads.interaction)

unique_reads.add <- lmer(unique_reads ~ input_reads_N + source + (1|sample), data = df)
summary(unique_reads.add)
waldtest(unique_reads.interaction, unique_reads.add)

unique_reads.source <- lmer(unique_reads ~ source + (1|sample), data = df)
summary(unique_reads.source)

unique_reads.input <- lmer(unique_reads ~ input_reads_N + (1|sample), data = df)
summary(unique_reads.input)
waldtest(unique_reads.add, unique_reads.source)
waldtest(unique_reads.add, unique_reads.input)

par(mfrow=c(2,2))
plot(unique_reads.add, which = 1)
plot(unique_reads.add, which = 2)
plot(unique_reads.add, which = 3)
plot(unique_reads.add, which = 5)

#### input reads vs coverage----
df %>%
  ggplot(aes(x=input_reads_N,
             y=coverage,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

coverage.interaction <- lmer(coverage ~ input_reads_N*source + (1|sample), data = df)
summary(coverage.interaction)

coverage.add <- lmer(coverage ~ input_reads_N + source + (1|sample), data = df)
summary(coverage.add)
waldtest(coverage.interaction, coverage.add)

coverage.source <- lmer(coverage ~ source + (1|sample), data = df)
summary(coverage.source)

coverage.input <- lmer(coverage ~ input_reads_N + (1|sample), data = df)
summary(coverage.input)
waldtest(coverage.add, coverage.source)
waldtest(coverage.add, coverage.input)

par(mfrow=c(2,2))
plot(coverage.add, which = 1)
plot(coverage.add, which = 2)
plot(coverage.add, which = 3)
plot(coverage.add, which = 5)

#### input reads vs GC linear----
df %>%
  ggplot(aes(x=input_reads_N,
             y=GC,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

gc.interaction <- lmer(GC ~ input_reads_N*source + (1|sample), data = df)
summary(gc.interaction)

gc.add <- lmer(GC ~ input_reads_N + source + (1|sample), data = df)
summary(gc.add)
waldtest(gc.interaction, gc.add)

gc.source <- lmer(GC ~ source + (1|sample), data = df)
summary(gc.source)

gc.input <- lmer(GC ~ input_reads_N + (1|sample), data = df)
summary(gc.input)
waldtest(gc.add, gc.source)
waldtest(gc.add, gc.input)

par(mfrow=c(2,2))
plot(gc.add, which = 1)
plot(gc.add, which = 2)
plot(gc.add, which = 3)
plot(gc.add, which = 5)

####GC log----
gc1.interaction <- glmmTMB(gc ~ input_reads_N*source + (1|sample), 
                           data = df1, family = beta_family())
summary(gc1.interaction)

gc1.add <- glmmTMB(gc ~ input_reads_N + source + (1|sample), 
                   data = df1, family = beta_family())
summary(gc1.add)
waldtest(gc1.interaction, gc1.add)

gc1.source <- glmmTMB(gc ~ source + (1|sample), 
                      data = df1, family = beta_family())
summary(gc1.source)

gc1.input <- glmmTMB(gc ~ input_reads_N + (1|sample), 
                     data = df1, family = beta_family())
summary(gc1.input)
waldtest(gc1.add, gc1.source)
waldtest(gc1.add, gc1.input)

par(mfrow=c(2,2))
plot(gc1.add, which = 1)
plot(gc1.add, which = 2)
plot(gc1.add, which = 3)
plot(gc1.add, which = 5)

#### input reads vs SNPs_cont ----
df %>%
  ggplot(aes(x=input_reads_N,
             y=SNPs_cont,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

snps_cont.interaction <- lmer(SNPs_cont ~ input_reads_N*source + (1|sample), data = df)
summary(snps_cont.interaction)

snps_cont.add <- lmer(SNPs_cont ~ input_reads_N + source + (1|sample), data = df)
summary(snps_cont.add)
waldtest(snps_cont.interaction, snps_cont.add)

par(mfrow=c(2,2))
plot(snps_cont.add, which = 1)
plot(snps_cont.add, which = 2)
plot(snps_cont.add, which = 3)
plot(snps_cont.add, which = 5)

#### input reads vs contam----
df %>%
  ggplot(aes(x=input_reads_N,
             y=contam,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

contam.interaction <- glmmTMB(contam ~ input_reads_N*source + (1|sample), 
                              data = df, family = beta_family())
summary(contam.interaction)

contam.add <- glmmTMB(contam ~ input_reads_N + source + (1|sample), 
                      data = df, family = beta_family())
summary(contam.add)
waldtest(contam.interaction, contam.add)

contam.source <- glmmTMB(contam ~ source + (1|sample), 
                         data = df, family = beta_family())
summary(contam.source)

contam.input <- glmmTMB(contam ~ input_reads_N + (1|sample), 
                        data = df, family = beta_family())
summary(contam.input)
waldtest(contam.add, contam.source)
waldtest(contam.add, contam.input)

par(mfrow=c(2,2))
plot(contam.add, which = 1)
plot(contam.add, which = 2)
plot(contam.add, which = 3)
plot(contam.add, which = 5)

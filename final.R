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

#### input reads vs SNPs_1240k----
snps.interaction <- lm(SNPs_1240k_N ~ input_reads_N*source, data = df)

snps.interaction.random <- lme(fixed = SNPs_1240k_N ~ input_reads_N * source,
                 random = ~ 1 | sample,
                 data = df)

lmer_additive_slope <- lmer(formula = SNPs_1240k_N ~ input_reads_N*source + (1 | sample), 
                            data    = df)

anova(snps.interaction, snps.interaction.random)

snps.add <- lm(SNPs_1240k_N ~ input_reads_N+source, data = df)
summary(snps.add)
Anova(snps.interaction) 
# Up can tell us interaction does not matter
# Below do the same
waldtest(snps.interaction, snps.add)

snps.input <- lm(SNPs_1240k_N ~ input_reads_N, data = df)
summary(snps.input)

snps.source <- lm(SNPs_1240k_N ~ source, data = df)
summary(snps.source)
waldtest(snps.add, snps.input)
# Need source
waldtest(snps.add, snps.source)
# Need input

summary(snps.add)
# SNPs_1240k_N=0.6103+0.3071×input_reads_N−1.1740×source1240k+ϵ
# Twist performance better
par(mfrow=c(2,2))
plot(snps.add, which = 1)
plot(snps.add, which = 2)
plot(snps.add, which = 3)
plot(snps.add, which = 5)
# Outlier 怎么办，如果plot不完美怎么办
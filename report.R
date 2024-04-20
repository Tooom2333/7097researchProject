pacman::p_load(tidyverse,
               ggplot2,
               dplyr,
               lme4,
               lmerTest,
               moderndive,
               performance,
               GGally,
               glmmTMB,
               DHARMa,
               stats,
               car,
               redres,
               lmtest)
library(ggcorrplot, inspect)


df.raw <- read_tsv("MODS_Rohrlach_EAGER.tsv")
df.raw <- df.raw %>%
  filter(nr_of_input_reads>1e6)

view(df.raw)
df.raw

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
view(df)
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

df

df %>%
  select(-sample)%>%
  select(!contains('_N'))%>%
  ggpairs(aes(colour = source))


#### input reads vs SNPs_1240k----
summary(df$SNPs_1240k)

df %>%
  ggplot(aes(x=input_reads_N,
             y=SNPs_1240k_N,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE,
              method = "lm")

snps.interaction <- lm(SNPs_1240k_N ~ input_reads_N*source, data = df)
summary(snps.interaction)

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

#### input reads vs mapped_reads----
df %>%
  ggplot(aes(x=input_reads_N,
             y=mapped_reads_N,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

mapped_reads.interaction <- lm(mapped_reads_N ~ input_reads_N*source, data = df)
summary(mapped_reads.interaction)

mapped_reads.add <- lm(mapped_reads_N ~ input_reads_N+source, data = df)
summary(mapped_reads.add)
# P is high, so they performance same well
waldtest(mapped_reads.interaction, mapped_reads.add)

mapped_reads.input <- lm(mapped_reads_N ~ input_reads_N, data = df)
summary(mapped_reads.add)

mapped_reads.source <- lm(mapped_reads_N ~ source, data = df)
summary(mapped_reads.source)

waldtest(mapped_reads.add, mapped_reads.input)
waldtest(mapped_reads.add, mapped_reads.source)
summary(mapped_reads.add)
par(mfrow=c(2,2))

plot(mapped_reads.add, which = 1)
plot(mapped_reads.add, which = 2)
plot(mapped_reads.add, which = 3)
plot(mapped_reads.add, which = 5)
# Twist better


#### input reads vs endogenous----
df %>%
  ggplot(aes(x=input_reads_N,
             y=endogenous,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

endogenous.interaction <- glm(endogenous ~ input_reads_N*source, 
                              data = df, family = "binomial")
summary(endogenous.interaction)
endogenous.add <- glm(endogenous ~ input_reads_N+source, 
                      data = df, family = "binomial")
summary(endogenous.add)

waldtest(endogenous.interaction, endogenous.add)
summary(endogenous.add)

endogenous.source <- glm(endogenous ~ source, 
                              data = df, family = "binomial")
summary(endogenous.source)
endogenous.input <- glm(endogenous ~ input_reads_N , 
                      data = df, family = "binomial")
summary(endogenous.input)
waldtest(endogenous.add, endogenous.source)

waldtest(endogenous.add, endogenous.input)
summary(endogenous.source)
endogenous.null <- glm(endogenous ~ 1 , 
                        data = df, family = "binomial")
summary(endogenous.null)
waldtest(endogenous.null, endogenous.source)

par(mfrow=c(2,2))
plot(endogenous.source, which = 1)
plot(endogenous.source, which = 2)
plot(endogenous.source, which = 3)
plot(endogenous.source, which = 5)
summary(endogenous.source)

#### input reads vs prop_dup----
df %>%
  ggplot(aes(x=input_reads_N,
             y=prop_dup,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)
# Interaction
prop_dup.interaction <- glm(prop_dup ~ input_reads_N*source, 
                              data = df, family = "binomial")
summary(prop_dup.interaction)

# Add
prop_dup.add <- glm(prop_dup ~ input_reads_N+source, 
                      data = df, family = "binomial")
summary(prop_dup.add)
# Compare
waldtest(prop_dup.interaction, prop_dup.add)
# Winner
summary(prop_dup.add)

# Source
prop_dup.source <- glm(prop_dup ~ source, 
                         data = df, family = "binomial")
summary(prop_dup.source)
# Input
prop_dup.input <- glm(prop_dup ~ input_reads_N , 
                        data = df, family = "binomial")
summary(prop_dup.input)
# Compare
waldtest(prop_dup.add, prop_dup.source)

waldtest(prop_dup.add, prop_dup.input)

prop_dup.null <- glm(prop_dup ~ 1 , 
                      data = df, family = "binomial")
waldtest(prop_dup.source, prop_dup.null)
waldtest(prop_dup.input, prop_dup.null)
waldtest(prop_dup.add, prop_dup.null)
par(mfrow=c(2,2))

plot(prop_dup.null, which = 1)
plot(prop_dup.null, which = 2)
plot(prop_dup.null, which = 3)
plot(prop_dup.null, which = 5)
summary(prop_dup.null)

#### input reads vs mtNuc_ratio----
df %>%
  ggplot(aes(x=input_reads_N,
             y=mtNuc_ratio,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)
# Interaction
summary(df$mtNuc_ratio)
mtNuc_ratio.interaction <- lm(mtNuc_ratio ~ input_reads_N*source, data = df)
summary(mtNuc_ratio.interaction)

# Add
mtNuc_ratio.add <- lm(mtNuc_ratio ~ input_reads_N+source, data = df)
summary(mtNuc_ratio.add)
# Compare
waldtest(mtNuc_ratio.interaction, mtNuc_ratio.add)
# Winner
summary(mtNuc_ratio.add)

# Source
mtNuc_ratio.source <- lm(mtNuc_ratio ~ source, data = df)
summary(mtNuc_ratio.source)
# Input
mtNuc_ratio.input <- lm(mtNuc_ratio ~ input_reads_N, data = df)
summary(mtNuc_ratio.input)
# Compare
waldtest(mtNuc_ratio.add, mtNuc_ratio.source)

waldtest(mtNuc_ratio.add, mtNuc_ratio.input)
summary(mtNuc_ratio.source)

mtNuc_ratio.null <- lm(mtNuc_ratio ~ 1, data = df)
waldtest(mtNuc_ratio.source, mtNuc_ratio.null)
par(mfrow=c(2,2))

plot(mtNuc_ratio.source, which = 1)
plot(mtNuc_ratio.source, which = 2)
plot(mtNuc_ratio.source, which = 3)
plot(mtNuc_ratio.source, which = 5)
summary(mtNuc_ratio.source)

#### input reads vs unique_reads----
df %>%
  ggplot(aes(x=input_reads_N,
             y=unique_reads,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)
# Interaction
unique_reads.interaction <- lm(unique_reads ~ input_reads_N*source, data = df)
summary(unique_reads.interaction)
# Add
unique_reads.add <- lm(unique_reads ~ input_reads_N+source, data = df)
summary(unique_reads.add)
# Compare
waldtest(unique_reads.interaction, unique_reads.add)
# Winner
summary(unique_reads.add)

# Source
unique_reads.source <- lm(unique_reads ~ source, data = df)
summary(unique_reads.source)
# Input
unique_reads.input <- lm(unique_reads ~ input_reads_N, data = df)
summary(unique_reads.input)
# Compare
waldtest(unique_reads.add, unique_reads.source)

waldtest(unique_reads.add, unique_reads.input)
summary(unique_reads.add)
par(mfrow=c(2,2))

plot(unique_reads.add, which = 1)
plot(unique_reads.add, which = 2)
plot(unique_reads.add, which = 3)
plot(unique_reads.add, which = 5)
summary(unique_reads.add)

#### input reads vs coverage----
df %>%
  ggplot(aes(x=input_reads_N,
             y=coverage,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)
# Interaction
coverage.interaction <- lm(coverage ~ input_reads_N*source, data = df)
summary(coverage.interaction)
# Add
coverage.add <- lm(coverage ~ input_reads_N+source, data = df)
summary(coverage.add)
# Compare
waldtest(coverage.interaction, coverage.add)
# Winner
summary(coverage.add)

# Source
coverage.source <- lm(coverage ~ source, data = df)
summary(coverage.source)
# Input
coverage.input <- lm(coverage ~ input_reads_N, data = df)
summary(coverage.input)
# Compare
waldtest(coverage.add, coverage.source)

waldtest(coverage.add, coverage.input)
summary(coverage.add)
par(mfrow=c(2,2))

plot(coverage.add, which = 1)
plot(coverage.add, which = 2)
plot(coverage.add, which = 3)
plot(coverage.add, which = 5)
summary(coverage.add)

#### input reads vs GC linear----
df %>%
  ggplot(aes(x=input_reads_N,
             y=GC,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)
df1 <- df%>%
  mutate(gc=GC/100)

# Interaction
df$GC
gc.interaction <- lm(GC ~ input_reads_N*source, data = df)
summary(gc.interaction)
# Add
gc.add <- lm(GC ~ input_reads_N + source, data = df)
summary(gc.add)
# Compare
waldtest(gc.interaction, gc.add)
# Winner
summary(gc.add)

# Source
gc.source <- lm(GC ~ source, data = df)
summary(gc.source)
# Input
gc.input <- lm(GC ~ input_reads_N, data = df)
summary(gc.input)
# Compare
waldtest(gc.add, gc.source)

waldtest(gc.add, gc.input)
summary(gc.add)
par(mfrow=c(2,2))

plot(gc.add, which = 1)
plot(gc.add, which = 2)
plot(gc.add, which = 3)
plot(gc.add, which = 5)
summary(gc.add)


####GC log----
gc1.interaction <- glm(gc ~ input_reads_N*source, data = df1, family = "binomial")
summary(gc1.interaction)
# Add
gc1.add <- glm(gc ~ input_reads_N + source, data = df1,family = "binomial")
summary(gc1.add)
# Compare
waldtest(gc1.interaction, gc1.add)
# Winner
summary(gc1.add)

# Source
gc1.source <- glm(gc ~ source, data = df1, family = "binomial")
summary(gc1.source)
# Input
gc1.input <- glm(gc ~ input_reads_N, data = df1, family = "binomial")
summary(gc1.input)
# Compare
waldtest(gc1.add, gc1.source)

waldtest(gc1.add, gc1.input)
gc1.null <- glm(gc ~ 1, data = df1, family = "binomial")
summary(gc1.null)
waldtest(gc1.source, gc1.null)
waldtest(gc1.input, gc1.null)
summary(gc1.null)


#### input reads vs SNPs_cont ----
df %>%
  ggplot(aes(x=input_reads_N,
             y=SNPs_cont,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)
# Interaction
snps_cont.interaction <- lm(SNPs_cont ~ input_reads_N*source, data = df)
summary(snps_cont.interaction)
# Add
snps_cont.add <- lm(SNPs_cont ~ input_reads_N + source, data = df)
summary(snps_cont.add)
# Compare
waldtest(snps_cont.interaction, snps_cont.add)
# Winner
summary(snps_cont.interaction)
Anova(snps_cont.interaction)
par(mfrow=c(2,2))

plot(snps_cont.interaction, which = 1)
plot(snps_cont.interaction, which = 2)
plot(snps_cont.interaction, which = 3)
plot(snps_cont.interaction, which = 5)
summary(snps_cont.interaction)

#### input reads vs contam----
df %>%
  ggplot(aes(x=input_reads_N,
             y=contam,
             col=source))+
  theme_bw()+
  geom_point(size = 2)+
  geom_smooth(se = FALSE)

df$contam <- df$contam + 0.1
# Interaction
contam.interaction <- glm(contam ~ input_reads_N*source, 
                            data = df, family = "binomial")
summary(contam.interaction)

# Add
contam.add <- glm(contam ~ input_reads_N+source, 
                    data = df, family = "binomial")
summary(contam.add)
# Compare
waldtest(contam.interaction, contam.add)
# Winner
summary(contam.add)

# Source
contam.source <- glm(contam ~ source, data = df, family = "binomial")
summary(contam.source)
# Input
contam.input <- glm(contam ~ input_reads_N , data = df, family = "binomial")
summary(contam.input)
# Compare
waldtest(contam.add, contam.source)

waldtest(contam.add, contam.input)

contam.null <- glm(contam ~ 1 , data = df, family = "binomial")
waldtest(contam.add, contam.null)
waldtest(contam.input, contam.null)
waldtest(contam.source, contam.null)

summary(contam.null)
par(mfrow=c(2,2))
plot(contam.null, which = 1)
plot(contam.null, which = 2)
plot(contam.null, which = 3)
plot(contam.null, which = 5)

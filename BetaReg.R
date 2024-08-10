pacman::p_load(tidyverse,
               glmmTMB,
               lmtest,
               psych,
               DHARMa,
               emmeans,
               ggplot2,
               betareg)

Data <- read.table(header=TRUE, stringsAsFactors=TRUE, text="
Class            Grade  Pass  Fail
 'Heron Hill'    12     16     4
 'Heron Hill'    12     14     6
 'Heron Hill'    12     15     5
 'Heron Hill'    11     18     2
 'Hawk Haven'    11     10    10
 'Heron Hill'    10     17     3
 'Hawk Haven'    10      9    11
 'Hawk Haven'     9     12     8
 'Cape May'       9      8    12
 'Hawk Haven'     8     10    10
 'Cape May'       8      8    12
 'Hawk Haven'     7     12     8
 'Cape May'       7      4    16
") %>%
  dplyr::mutate(Proportion=Pass/(Pass+Fail),
                Sex=sample(c('XX','XY'),size=n(),replace=T)) %>%
  as_tibble() %>%
  dplyr::arrange(-Proportion)

Data %>%
  ggplot(aes(x=Grade,
             y=Proportion,
             shape=Sex))+
  theme_bw()+
  geom_point(size=3)+
  scale_shape_manual(values=c(21,22))

Data %>%
  ggplot(aes(x=Grade,
             y=Proportion,
             fill=Class,
             shape=Sex))+
  theme_bw()+
  geom_point(size=3)+
  scale_shape_manual(values=c(21,22))
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              aes(col=Sex,group=Sex),
              linetype='dashed',
              se = T)
  
  Data %>%
    ggplot(aes(x=Class,y=Grade,fill=Class))+
    theme_bw()+
    geom_boxplot(outlier.shape=21)
  
  lm(Grade~Class,data=Data) %>%
    summary()

model.beta.int <- betareg(Proportion~Grade*Sex, 
                      data=Data)
model.beta.add <- betareg(Proportion~Grade+Sex, 
                          data=Data)
model.beta.sex <- betareg(Proportion~Sex, 
                          data=Data)
model.beta.grade <- betareg(Proportion~Grade, 
                          data=Data)
model.beta.null <- betareg(Proportion~1, 
                          data=Data)
AIC(model.beta.int,
    model.beta.add,
    model.beta.grade,
    model.beta.sex,
    model.beta.null) %>% 
  dplyr::arrange(AIC)

summary(model.beta.grade)

par(mfrow = c(3, 2))
plot(model.beta.grade, which = 1:6)
par(mfrow = c(1, 1))

# Mixed effects?
model.beta.mef <- glmmTMB(Proportion~1+(1|Class),
                          data=Data, 
                          family=beta_family(link="logit"))
model.beta.mef.sex <- glmmTMB(Proportion~Sex+(1|Class),
                          data=Data, 
                          family=beta_family(link="logit"))
model.beta.mef.grade <- glmmTMB(Proportion~Grade+(1|Class),
                          data=Data, 
                          family=beta_family(link="logit"))
anova(model.beta.mef,
      model.beta.mef.grade)

model.beta.fef <- glmmTMB(Proportion~1-1,
                          data=Data, 
                          family=beta_family(link="logit"))
anova(model.beta.mef,
      model.beta.mef.sex)

summary(model.beta.mef)
anova(model.beta.mef,
      model.beta.fef)

simulateResiduals(model.beta.mef) %>%
  DHARMa::plotQQunif()







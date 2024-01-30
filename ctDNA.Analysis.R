#Analysis of ctDNA Data
library(readxl)
library(ggplot2)
library(dplyr)
library(swimplot)
library(survival)
library(survminer)
library(ggsci)
library(showtext) #Adds sysfonts as dependency
#font_add_google("Roboto Condensed") #Only Once
#showtext_auto() #Only Once

ctDNA.Swim <-  read_excel("ctDNA.xlsx", sheet="Biocartis") %>% 
  mutate_at(c("Patient", "ID", "Result"), as.factor) %>% 
  mutate_at(c("Date"), as.Date) %>% 
  group_by(ID) %>%
  mutate(start.date=min(Date),
         last.date=max(Date),
        days.test=Date-min(Date),
        duration.study=max(Date)-min(Date)) %>%
  arrange(start.date, Date) %>%
  mutate(start.positive=max(ifelse(Date==start.date & Result=="1", 1, 0)),
         last.positive=max(ifelse(Date==last.date & Result=="1", 1, 0)),
         number.positive= sum(as.numeric(as.character(Result)), na.rm = TRUE),
         number.total= sum(!is.na(Result)))

ctDNA.Swim.Index <- read_excel("ctDNA.xlsx", sheet="Index") %>% 
  mutate_at(c("DoB", "Death", "Stage4"), as.Date) %>% 
  mutate_at(c("ID", "Sex", "BRAF", "Stage", "ANA"), as.factor)


# Kaplan Maeier Plots---- =================

ctDNA.Analysis <-   merge(x=ctDNA.Swim.Index, y=ctDNA.Swim[,c("ID",
                                                              "start.date", 
                                                              "duration.study", 
                                                              "start.positive", 
                                                              "last.positive",
                                                              "number.positive",
                                                              "number.total")], 
                          by='ID', all.x = TRUE)
#ctDNA.Analysis <- distinct(ctDNA.Analysis, "ID", "Stage4", "start.date", "number.positive", .keep_all = FALSE)
# ctDNA.Analysis$Death <- as.Date(ctDNA.Analysis$Death) 

ctDNA.Analysis <-  unique(ctDNA.Analysis) %>%
  mutate(duration.followup=ifelse(is.na(Death), Sys.Date()-start.date, Death - start.date),
        death.event=ifelse(is.na(Death), 0, 1),
        age=as.integer(as.numeric(difftime(start.date, DoB, units = "weeks"))/52.25),
        ratio.positive=number.positive/number.total,
        ctDNA.detection.rate=(ratio.positive/duration.followup)*100,
        mets.brain.start1=ifelse(Metastasis==2 & start.positive==1, 2, 
                                  ifelse(Metastasis==2 & start.positive==0, 1, 0)),
        mets.all.start1=ifelse(Metastasis>=1 & start.positive==1, 2,
                               ifelse(Metastasis>=1 & start.positive==0, 1, 0)),
        mets.brain.startend1=ifelse(Metastasis==2 & start.positive==1 & last.positive==1, 1, 
                                    ifelse(Metastasis==2 & !(start.positive==1 & last.positive==1), 0, NA)))
#ctDNA.Analysis$duration.followup <- as.numeric(ctDNA.Analysis$duration.followup)
#View Quantiles
#quantile(survfit(Surv(duration.study, death.event) ~ start.positive, data=ctDNA.Analysis))

#View survival information at each time points
#summary(survfit(Surv(duration.study, death.event) ~ start.positive, data=ctDNA.Analysis))

#View Median Survival Information
#print(survfit(Surv(duration.study, death.event) ~ start.positive, data=ctDNA.Analysis))
# print(survfit(Surv(duration.study, death.event) ~ Metastasis, data=ctDNA.Analysis))

#Log Rank Test
#survdiff(Surv(duration.study, death.event) ~ start.positive, data=ctDNA.Analysis)
# summary(ctDNA.Analysis$duration.followup)
# summary(subset(ctDNA.Analysis, Metastasis == 1 & age >= 18, select=duration.followup))
# summarise(subset(ctDNA.Analysis, Metastasis == 1, select=age))
# 
# print(survfit( Surv(duration.followup, death.event)~Metastasis, data=ctDNA.Analysis))
# 
# df_stat <- ctDNA.Analysis %>% group_by(Metastasis) %>% summarize(
#   count = n(),
#   mean = mean(age, na.rm = TRUE), 
#   sd = sd(age, na.rm = TRUE))
# 
# # T Test
# 
# tTest.data <- subset(ctDNA.Analysis, as.numeric(Metastasis) > 0 & number.total > 1, select=c(Metastasis, ratio.positive)) 
# boxplot(ratio.positive~Metastasis, data = tTest.data)
#   t.test(ratio.positive~Metastasis, data = tTest.data, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
#   t.test(ratio.positive[Metastasis=="1"], ratio.positive[Metastasis=="2"])
  
  
#Median survival time calculation
survfit(Surv(duration.followup, death.event) ~ 1, data = ctDNA.Analysis)
#Compare survival time between groups
survdiff(Surv(duration.followup, death.event) ~ Metastasis, data = subset(ctDNA.Analysis, Metastasis >= 1 & age >= 18)) 


#The chi-square statistic and p-value generated above are for the log rank test of the null
#...hypothesis of no difference between the two survival curve.

#Cox-Proportional Hazard Model
#coxph(Surv(duration.study, death.event) ~ start.positive + age + number.positive, data=ctDNA.Analysis)
#In the above, in place of start.positive, it could be any other or multiple parameters like 
#... Age + Stage + Subtype and so on. Read https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R-Manual/R-Manual21.html

#To remove the variable to calculate using different parameters no need if using different params
#remove(.Survfit) 

#Kaplan Meier of All Patients----=================
.Survfit.All <- survfit(Surv(duration.followup, death.event) ~ 1, data = ctDNA.Analysis)
ggsurvplot(.Survfit.All, palette = "#2E9FDF",
           risk.table = TRUE)



#Kaplan Meier Treatment Settings (Manuscript)----=================
#Duration Follow Up: Considering Alive people to be followed up till current date
.Survfit.Metastasis <- survfit(formula = Surv(duration.followup, death.event) ~ Metastasis,
                                   conf.type ="log", conf.int = 0.95,
                                   type = "kaplan-meier", error = "greenwood",
                                   data=ctDNA.Analysis)

#Duration Study: Considering the duration from first test to last test
# .Survfit.Metastasis <- survfit(formula = Surv(duration.study, death.event) ~ Metastasis,
#                                conf.type ="log", conf.int = 0.95,
#                                type = "kaplan-meier", error = "greenwood", 
#                                data=ctDNA.Analysis)
#.Survfit.Metastasis
ggsurvplot(.Survfit.Metastasis,
           risk.table = TRUE, risk.table.pos = "out", risk.table.col = "black", risk.table.y.text.col = TRUE,
           #tables.theme = theme_cleantable(),
           #ncensor.plot= TRUE, ncensor.plot.height= 0.1,
           font.tickslab = c(12),
           legend.title = "Groups",
           pval = TRUE, pval.method = TRUE,
           #surv.median.line = "hv",
           conf.int = FALSE,
           cumevents = FALSE,
           cumcensor = FALSE,
           break.time.by = 60,
           palette = "cosmic", #lancet, jama, aaas, nejm, npg, jco, cosmic, igv, locuszoom, d3
           legend =c(0.9, 0.4),
           legend.labs = c("Adjuvant", "Mets to other organs", "Mets to the brain"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           title = "All Patients: adjuvant, metastasis to brain & metastasis to other organs",
           ggtheme = theme_classic(base_size=12, base_family = "Roboto Condensed"),
           font.family="Roboto Condensed")

#Kaplan Meier Stage IV Settings (SMR Poster 20211020)----=================
#Duration Follow Up: Considering Alive people to be followed up till current date
ctDNA.Analysis.IV <- filter(ctDNA.Analysis, Metastasis >= 1)
.Survfit.Metastasis.IV <- survfit(formula = Surv(duration.followup, death.event) ~ Metastasis,
                               conf.type ="log", conf.int = 0.95,
                               type = "kaplan-meier", error = "greenwood",
                               data=ctDNA.Analysis.IV)

#Duration Study: Considering the duration from first test to last test
# .Survfit.Metastasis <- survfit(formula = Surv(duration.study, death.event) ~ Metastasis,
#                                conf.type ="log", conf.int = 0.95,
#                                type = "kaplan-meier", error = "greenwood", 
#                                data=ctDNA.Analysis)
#.Survfit.Metastasis
ggsurvplot(.Survfit.Metastasis.IV,
           risk.table = TRUE, risk.table.pos = "out", risk.table.col = "black", risk.table.y.text.col = TRUE,
           #censor.size=0.1, 
           size = 1.2,
           #tables.theme = theme_cleantable(),
           #ncensor.plot= TRUE, ncensor.plot.height= 0.1,
           font.tickslab = c(12),
           legend.title = "Groups",
           pval = TRUE, pval.method = TRUE,
           #surv.median.line = "hv",
           conf.int = FALSE,
           cumevents = FALSE,
           cumcensor = FALSE,
           break.time.by = 60,
           #palette = "cosmic", #lancet, jama, aaas, nejm, npg, jco, cosmic, igv, locuszoom, d3
           palette = c("#51CF66", "#F08C00"),
           legend =c(0.88, 0.32),
           #legend.title = element_text(size = 20, color = "#005BBB", hjust=0.5, vjust = 1),
           legend.labs = c("Mets to other organs", "Mets to the brain"),
           #ggplot.title = element_text(size = 20, color = "#005BBB", hjust=0.5, vjust = 1),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           title = "Stage IV Melanoma; metastasis to other organs vs metastasis to the CNS",
           ggtheme = theme_classic(base_size=12, base_family = "Roboto Condensed"),
           font.family="Roboto Condensed")

#Kaplan Meier All Baseline Positive vs Baseline Negative----=================
.Survfit.Start.Positive <- survfit(formula = Surv(duration.followup, death.event) ~ start.positive,
                                   conf.type ="log", conf.int = 0.95,
                                   type = "kaplan-meier", error = "greenwood", 
                                   data=ctDNA.Analysis)
#.Survfit.Start.Positive

# plot(.Survfit.Start.Positive, col=1:2, lty=1:2, conf.int=FALSE, mark.time=TRUE, xlab = "Days", ylab = "Survival Probability")
# legend("bottomleft", legend = c("Started Negative", "Started Positive"), title = "ctDNA Mutation Detection",
#        col=1:2, lty = 1:2, bty="n")
# title("Kaplan Meier Curve for \nSurvival by Initial Detection")

ggsurvplot(.Survfit.Start.Positive,
           risk.table = TRUE,
           pval = TRUE, pval.method = TRUE,
           surv.median.line = "hv",
           conf.int = FALSE,
           cumevents = FALSE,
           cumcensor = FALSE,
           break.time.by = 60,
           palette = "nejm", #lancet, jama, aaas, nejm, npg, jco, cosmic, igv, locuszoom, d3
           legend =c(0.8,0.8),
           legend.labs = c("Baseline Negative", "Baseline Positive"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           title = "All Patients: Baseline Positive vs Baseline Negative ctDNA",
          ggtheme = theme_classic(base_size=12, base_family = "Roboto Condensed"),
          font.family="Roboto Condensed")

           


#Kaplan Meier Brain Mets Baseline Positive vs Baseline Negative vs Others (manuscript)----================
.Survfit.Mets.Brain.Start1 <- survfit(Surv(duration.followup, death.event) ~ mets.brain.start1, 
                                      conf.type ="log", conf.int = 0.95, 
                                      type = "kaplan-meier", error = "greenwood", data=ctDNA.Analysis)
#.Survfit.Mets.Brain.Start1
ggsurvplot(.Survfit.Mets.Brain.Start1,
           risk.table = TRUE,
           pval = TRUE, pval.method = TRUE,
           #surv.median.line = "hv",
           conf.int = FALSE,
           cumevents = FALSE,
           cumcensor = FALSE,
           break.time.by = 60,
           palette = "jco", #lancet, jama, aaas, nejm, npg, jco, cosmic, igv, locuszoom, d3
           legend =c(0.8,0.65),
           legend.labs = c("No Brain Metastasis", "Brain Mets: Baseline Negative", "Brain Mets: Baseline Positive"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           title = "Patients with Brain Metastasis: Baseline Positive vs Baseline Negative ctDNA",
           ggtheme = theme_classic(base_size=12, base_family = "Roboto Condensed"),
           font.family="Roboto Condensed")


#Kaplan Meier All Mets Baseline Positive vs Baseline Negative----=================
.Survfit.Mets.All.Start1 <- survfit(Surv(duration.followup, death.event) ~ mets.all.start1, 
                                      conf.type ="log", conf.int = 0.95, 
                                      type = "kaplan-meier", error = "greenwood", data=ctDNA.Analysis)
#.Survfit.Mets.All.Start1
ggsurvplot(.Survfit.Mets.All.Start1,
           risk.table = TRUE, risk.table.pos = "out", risk.table.col = "black", risk.table.y.text.col = TRUE,
           #tables.theme = theme_cleantable(),
           #ncensor.plot= TRUE, ncensor.plot.height= 0.1,
           font.tickslab = c(12),
           legend.title = "Groups",
           pval = TRUE, pval.method = TRUE,
           surv.median.line = "hv",
           conf.int = FALSE,
           cumevents = FALSE,
           cumcensor = FALSE,
           break.time.by = 60,
           palette = "jama", #lancet, jama, aaas, nejm, npg, jco, cosmic, igv, locuszoom, d3
           legend =c(0.9, 0.6),
           legend.labs = c("Adjuvant: Baseline Negative", "Metastatic: Baseline Negative", "Metastatic: Baseline Positive"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           title = "Adjuvantly treated Patients with Baseline Negative vs Started Positive vs Started Negative ctDNA",
           ggtheme = theme_classic(base_size=12, base_family = "Roboto Condensed"),
           font.family="Roboto Condensed")

#Kaplan Meier Brain Mets Start and Last Positive----=================
.Survfit.Mets.Brain.StartEnd1 <- survfit(Surv(duration.followup, death.event) ~ mets.brain.startend1, 
                                    conf.type ="log", conf.int = 0.95, 
                                    type = "kaplan-meier", error = "greenwood", data=ctDNA.Analysis)
#.Survfit.Mets.Brain.StartEnd1
ggsurvplot(.Survfit.Mets.Brain.StartEnd1,
           risk.table = TRUE, risk.table.pos = "out", risk.table.col = "black", risk.table.y.text.col = TRUE,
           #tables.theme = theme_cleantable(),
           #ncensor.plot= TRUE, ncensor.plot.height= 0.1,
           font.tickslab = c(12),
           legend.title = "Groups",
           pval = TRUE, pval.method = TRUE,
           #surv.median.line = "hv",
           conf.int = FALSE,
           cumevents = FALSE,
           cumcensor = FALSE,
           break.time.by = 60,
           palette = "jama", #lancet, jama, aaas, nejm, npg, jco, cosmic, igv, locuszoom, d3
           legend =c(0.8, 0.7),
           legend.labs = c("Not Both Baseline and Last Detectable", "Both Baseline and Last Detectable"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           title = "Patients with Brain Metastasis: Start and Last ctDNA Positive vs Not",
           ggtheme = theme_classic(base_size=12, base_family = "Roboto Condensed"),
           font.family="Roboto Condensed")

#Kaplan Meier Ratio Positive----=================
ratio <- as_tibble(ctDNA.Analysis)
ratio.non0 <- filter(ratio, ratio.positive > 0)%>%
  mutate(ratio.group=cut(ratio.positive, breaks=c(0, 0.4, Inf), 
                         labels=c("low.positivity", "high.positivity")))

.Survfit.Ratio.Positive <- survfit(Surv(duration.followup, death.event) ~ ratio.group, 
                                         conf.type ="log", conf.int = 0.95, 
                                         type = "kaplan-meier", error = "greenwood", 
                                         data=ratio.non0)
#.Survfit.Ratio.Positive
ggsurvplot(.Survfit.Ratio.Positive,
           risk.table = TRUE, risk.table.pos = "out", risk.table.col = "black", risk.table.y.text.col = TRUE,
           #tables.theme = theme_cleantable(),
           #ncensor.plot= TRUE, ncensor.plot.height= 0.1,
           font.tickslab = c(12),
           legend.title = "Groups",
           pval = TRUE, pval.method = TRUE,
           #surv.median.line = "hv",
           conf.int = FALSE,
           cumevents = FALSE,
           cumcensor = FALSE,
           break.time.by = 60,
           palette = "igv", #lancet, jama, aaas, nejm, npg, jco, cosmic, igv, locuszoom, d3
           legend =c(0.8, 0.8),
           legend.labs = c("Detection Ratio ≤ 40%", "Detection Ratio > 40%"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           title = "Patients with detectable ctDNA: high vs low detection Ratio",
           ggtheme = theme_classic(base_size=12, base_family = "Roboto Condensed"),
           font.family="Roboto Condensed")


#Kaplan Meier: Rate of Detection (manuscript)----=================
rate <- as_tibble(ctDNA.Analysis)
rate.non0 <- filter(rate, ctDNA.detection.rate > 0)%>%
  mutate(rate.group=cut(ctDNA.detection.rate, breaks=c(0, 0.3, Inf), 
                         labels=c("low.detection.rate", "high.detection.rate")))

.Survfit.Rate.Positive <- survfit(Surv(duration.followup, death.event) ~ rate.group, 
                                   conf.type ="log", conf.int = 0.95, 
                                   type = "kaplan-meier", error = "greenwood", 
                                   data=rate.non0)
#.Survfit.Rate.Positive
ggsurvplot(.Survfit.Rate.Positive,
           risk.table = TRUE, risk.table.pos = "out", risk.table.col = "black", risk.table.y.text.col = TRUE,
           #tables.theme = theme_cleantable(),
           #ncensor.plot= TRUE, ncensor.plot.height= 0.1,
           font.tickslab = c(12),
           legend.title = "Groups",
           pval = TRUE, pval.method = TRUE,
           #surv.median.line = "hv",
           conf.int = FALSE,
           cumevents = FALSE,
           cumcensor = FALSE,
           break.time.by = 60,
           palette = "locuszoom", #lancet, jama, aaas, nejm, npg, jco, cosmic, igv, locuszoom, d3
           legend =c(0.85, 0.8),
           legend.labs = c("Detection Rate ≤ 0.3", "Detection Rate > 0.3"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           #xlim = c(0, 500), #added to increase x limit
           title = "Patients with detectable ctDNA: high vs low rate of detection",
           ggtheme = theme_classic(base_size=12, base_family = "Roboto Condensed"),
           font.family="Roboto Condensed")

#Cox Hazard Ratio Positive (manuscript)----=================
ctDNA.Analysis <- ctDNA.Analysis %>%
  mutate_at(c("Metastasis", "start.positive"), as.factor) %>%
  mutate(ageGroup = ifelse(age >=50, "old","young"))


ctDNA.Analysis$Metastasis <- factor(ctDNA.Analysis$Metastasis, levels = c(0, 1, 2), 
                                    labels = c("Adjuvant", "Metastatic", "Brain Mets"))
ctDNA.Analysis$ageGroup <- factor(ctDNA.Analysis$ageGroup)

survCox <- coxph(Surv(duration.followup, death.event) ~ Sex + BRAF + ageGroup + ANA + start.positive,
                 data = ctDNA.Analysis)
ggforest(survCox, data = ctDNA.Analysis)

#Plotting Rate of Positivity---- 

rate <- as_tibble(ctDNA.Analysis)
rate.non0 <- filter(rate, ctDNA.detection.rate > 0)%>%
  mutate(rate.group=cut(ctDNA.detection.rate, breaks=c(0, 0.4, Inf), 
                        labels=c("low.detection.rate", "high.detection.rate")))
number.patients <- length(unique(rate.non0[["ID"]])) #[ gives a data.frame(list), [[ gives a vector

graph10 = rate.non0 %>%
  ggplot(aes(x=reorder(ID,-as.numeric(ctDNA.detection.rate)), y=ctDNA.detection.rate)) +
  geom_point(size= 2, color="#FF2550", alpha=0.75)+
  #theme_few()+
  theme_minimal() +
  labs(title="ctDNA Detection Rate",
       subtitle="Detection Rate = (# of Positive Tests/Total # of Tests)/Duration of Followup *100",
       x= "Patients",
       y= "ctDNA Detection Rate") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  scale_y_continuous(breaks=seq(0, 2.6, by=0.2))+
  geom_hline(aes(yintercept = median(ctDNA.detection.rate)))+
  scale_x_discrete(labels = 1:number.patients)
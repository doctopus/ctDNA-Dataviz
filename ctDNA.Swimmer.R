#Swimmer Plot and Analysis of ctDNA Data
library(readxl)
library(ggplot2)
library(dplyr)
library(swimplot)
library(survival)
library(survminer)
library(ggsci)
library(showtext) #Adds sysfonts as dependency
#font_add_google("Roboto Condensed")
#showtext_auto()

ctDNA.Swim <-  read_excel("ctDNA.xlsx", sheet="Biocartis") %>% 
  mutate_at(c("Patient", "ID", "Result"), as.factor) %>% 
  mutate_at(c("Date"), as.Date) %>% 
  group_by(ID) %>%
  mutate(start.date=min(Date),
        days.test=Date-min(Date),
        duration.study=max(Date)-min(Date)) %>%
  arrange(start.date, Date) %>%
  mutate(start.positive=max(ifelse(Date==start.date & Result=="1", 1, 0)),
        number.positive= sum(as.numeric(as.character(Result)), na.rm = TRUE))

ctDNA.Swim.Index <- read_excel("ctDNA.xlsx", sheet="Index") %>% 
  mutate_at(c("Death", "Stage4"), as.Date) %>% 
  mutate_at(c("ID"), as.factor)


#Swimmer Plot---- ==============
# ctDNA.Swimmer <-   merge(ctDNA.Swim, ctDNA.Swim.Index[,c(1:5)], by='ID', all.x = TRUE)
# ctDNA.Swimmer <- ctDNA.Swimmer %>%
#   group_by(ID) %>%
#   mutate(days.study=Date-min(Date), days.stage4=Date-Stage4)
#id order could arrange based on numeric values but it does not convert number to factor
#..on the fly, so changing it beforehand.
ctDNA.Analysis$Metastasis <- as.factor(ctDNA.Analysis$Metastasis)

swim <- swimmer_plot(df=ctDNA.Analysis, id="Serial", end="duration.followup",
                        id_order = "Metastasis",
                        col="#0F0F0F", name_fill= "Metastasis", width= .8)
                        
swim.dots <- swim + swimmer_points(df_points=ctDNA.Analysis, id="Serial", time="duration.study",
                                         name_shape=as.factor("mets.all.start1"),
                                         fill="White", col="black")
swim.dots.arrows <- swim.dots + swimmer_arrows(df_arrows = ctDNA.Analysis,
                       id="Serial", arrow_start = "duration.followup",
                       cont = "Sex", type = "open", arrow_positions = c(1, 10))

# swimmer.trial <- read_excel("/Users/i/Documents/dev/play/r/swimmer/Swimmer.Trial.xlsx")
# swimmer.trial.plot <- swimmer_plot(df=swimmer.trial, id="id", end="End_trt")



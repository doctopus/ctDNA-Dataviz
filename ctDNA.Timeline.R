library(ggplot2)
library(gganimate)
library(readxl)
library(dplyr)
library(ggthemes)
library(gifski)
library(ggsci)

ctDNA.Index <- read_excel("ctDNA.xlsx", sheet="Index") %>% 
  mutate_at(c("ID", "Sex", "BRAF", "Metastasis", "Stage", "ANA"), as.factor) %>% 
  mutate_at(c("DoB", "Death", "Stage4", "StageBrain"), as.Date)

ctDNA.Biocartis <-  read_excel("ctDNA.xlsx", sheet="Biocartis")%>% 
  mutate_at(c("Patient", "ID", "Result"), as.factor) %>% 
  mutate_at(c("Date"), as.Date)

#Graph 0 Number of tests over time----
Graph0 = ggplot(ctDNA.Biocartis, aes(format(Date, "%Y-%m"))) +
  geom_bar(stat = "count",
           fill = pal_aaas()(5)[c(1)]) +
  labs(title="ctDNA BRAF Cartridge Usage in Each Month",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Month",
       y= "Number of ctDNA BRAF Tests",
       caption = paste("Total ", nrow(ctDNA.Biocartis), "ctDNA BRAF tests done from",  min(ctDNA.Biocartis$Date), "to",max(ctDNA.Biocartis$Date)))+
  theme_excel_new()+ #theme_fivethirtyeight() +#theme_solarized()+ #theme_economist()+
  theme(title=element_text(color = "#657b83"),
        axis.text.y=element_text(size=15),
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.title.y = element_text(size=15, hjust = 0.5, vjust = 1),
        #axis.title.y = element_blank(),
        plot.title.position = "panel",
        plot.title = element_text(hjust=0.5, vjust = 1),
        plot.subtitle = element_text(hjust=0.5, vjust = 1),
        legend.title =  element_text(inherit.blank = FALSE, size = 10, color = "#000000"))



#Graph 1 Actual Dates in Timeline----
graph1 = ctDNA.Biocartis %>% 
  #ggplot(aes(x=Date, y=reorder(Patient,-patient.num), color=Result, size=result.num)) +
  ggplot(aes(x=Date, y=reorder(Patient,-as.numeric(Patient)), color=Result, size=as.numeric(Result))) +
  geom_line(size= 2, color="#06313C", alpha=0.75) + 
  geom_point(aes(group = seq_along(Date)), # needed, otherwise transition doesn't work
             alpha =0.9, stroke= 1) + 
  theme_economist() +
  labs(title="ctDNA Detection Timeline",
       subtitle="Detection of BRAF mutation in ctDNA",
       x= "Date",
       y= "Patients",
       color="Test Result") +
  guides(size=FALSE,
         color=guide_legend(reverse = TRUE)) +
  theme(title=element_text(color = "#657b83"),
        plot.title.position = "panel",
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "bottom") +
  scale_x_date(date_labels = "%b/%y", date_breaks = "2 month")+
  scale_color_discrete(breaks=c(0, 1),
                     labels=c("Negative", "Positive"))
#  options(repr.plot.width = 4, repr.plot.height = 3)
#  scale_y_reverse() #To reverse the Y axis
#  scale_size(range = c(0, 1))
#graph1
#ggsave("graph1.png", plot = graph1, height = 576, width = 1024, units=c("mm"))
#Graph1 Animation
graph1.animation = graph1 +
  transition_reveal(Date) +
#  view_follow(fixed_y = TRUE) +
#  ease_aes('cubic-in-out') +
#  ease_aes(x = 'bounce-in') +
  enter_grow() + 
  enter_fade() 

animate(graph1.animation, height= 576, width= 1024, fps = 20, duration = 15, start_pause= 10, end_pause = 120)
#anim_save("animation1.gif", animation = last_animation(), path = NULL)

#Graph 2 Days from Start of Study--------------------------------------------------------
#Data Manipulation
ctDNA.Patient <- ctDNA.Biocartis %>%
   group_by(Patient) %>%
   mutate(days=Date-min(Date))
#  summarise(start_date =min(Date))
graph2 =ctDNA.Patient %>%
  ggplot(aes(x=days, y=reorder(Patient,-as.numeric(Patient)), color=Result, size=as.numeric(Result))) +
  geom_line(size= 2, color="#06313C", alpha=0.75) +
  geom_point(aes(group = seq_along(Date)), # needed, otherwise transition doesn't work
             alpha =0.9, stroke= 1) + 
  labs(title="ctDNA Detection Timeline",
       subtitle="Detection of BRAF mutation in ctDNA",
       x= "days",
       y= "Patients",
       color="Test Result") +
  guides(size=FALSE,
         color=guide_legend(reverse = TRUE)) +
  theme_fivethirtyeight() +
  theme(title=element_text(color = "#657b83"),
        plot.title.position = "panel",
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "bottom") +
  scale_color_manual(values=c("#D79A48", "#823939"),
                     breaks=c(0, 1),
                     labels=c("Negative", "Positive"))+
  #scale_x_continuous(limits=c(0,300), breaks = seq(0, 300, by=30))
  scale_x_continuous(breaks=seq(0, max(as.numeric(ctDNA.Patient$days)), by=60))
#  scale_x_date(date_labels = "%b/%y", date_breaks = "2 month")
#graph2
#Graph2 Animation
graph2.animation = graph2 +
  transition_reveal(Date) +
  #  view_follow(fixed_y = TRUE) +
  #  ease_aes('cubic-in-out') +
  #  ease_aes(x = 'bounce-in') +
  enter_grow() + 
  enter_fade() 

#animate(graph2.animation, height= 500, width= 800, fps = 20, duration = 10, end_pause = 20)
animate(graph2.animation, height= 576, width= 1024, fps = 20, duration = 15, start_pause= 10, end_pause = 120)
#anim_save("animation2.gif", animation = last_animation(), path = NULL)



#Graph 3 By ID -Days from Study Start -Mets-----------------------------------------------
#Data Manipulation
ctDNA.ID <- merge(x=ctDNA.Biocartis, y=ctDNA.Index[,c("ID",
                                            "Death", 
                                            "Metastasis")], 
                                            by='ID', all.x = TRUE) %>%
#ctDNA.ID <- merge(ctDNA.Biocartis, ctDNA.Index[,c(1:3)], by='ID', all.x = TRUE) %>%
  group_by(ID) %>%
  mutate(start.date=min(Date),
        days=Date-min(Date),
        brain_mets=if_else(Metastasis==2, 1, 0)) %>%
  arrange(start.date, Date)
#Count number of unique patients; will be used as y-axis
number.patients <- length(unique(ctDNA.ID[["ID"]])) #[ gives a data.frame(list), [[ gives a vector

#Graph3
graph3 =ctDNA.ID %>%
  #ggplot(aes(x=days, y=reorder(ID,desc(start.date)), color=Result, size=ID.result.num)) +
  #ggplot(aes(x=days, y=reorder(ID,desc(start.date)), size=ID.result.num, color=Result)) +
  ggplot(aes(x=days, y=reorder(ID,desc(start.date)), size=Result, color=Result)) +
      geom_line(aes(color=as.factor(brain_mets)), size= 2, alpha=0.75) +
      geom_point(aes(group = seq_along(Date)), # needed, otherwise transition doesn't work
             alpha =0.9, stroke= 1) + 
  labs(title="ctDNA Detection Timeline",
       subtitle="Detection of BRAF mutation in ctDNA",
       x= "Date",
       y= "Patients",
       color="Brain Metastasis",
       size="Test Positivity") +
  guides(size=guide_legend(reverse = TRUE),
         color=guide_legend(reverse = TRUE)) +
  theme_fivethirtyeight() +
  theme(title=element_text(color = "#657b83"),
        plot.title.position = "panel",
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "bottom",
        legend.box = "horizontal") +
  scale_color_manual(values=c("#D79A48", "#823939"),
                     breaks=c(0, 1),
                     labels=c("Negative", "Positive")) +
  scale_size_manual(values = c(0, 2),
                    breaks = c(0, 1),
                    labels=c("Negative", "Positive")) +
  #scale_x_continuous(limits=c(0,500), breaks = seq(0, 500, by=30)) +
  scale_x_continuous(breaks=seq(0, max(as.numeric(ctDNA.ID$days)), by=60))+
  scale_y_discrete(labels = number.patients:1)
#graph3
#Graph3 Animation
graph3.animation = graph3 +
  transition_reveal(Date) +
  enter_grow() + 
  enter_fade() 

#animate(graph3.animation, height= 500, width= 800, fps = 20, duration = 10, end_pause = 60)
animate(graph3.animation, height= 576, width= 1024, fps = 20, duration = 15, start_pause= 10, end_pause = 120)
#anim_save("animation3.gif", animation = last_animation(), path = NULL)

#Graph 4 ID -Chronological from Beginning of Study -Brain -Mets-Adj------------------------------------------

#Data Manipulation

ctDNA.ID <- merge(x=ctDNA.Biocartis, y=ctDNA.Index[,c("ID",
                                            "Death", 
                                            "Metastasis")], 
                  by='ID', all.x = TRUE) %>%
#ctDNA.ID <- merge(ctDNA.Biocartis, ctDNA.Index[,c(1:3)], by='ID', all.x = TRUE) %>%
  group_by(ID) %>%
  mutate(start.date=min(Date),
        days=Date-min(Date),
        duration=max(Date)-min(Date)) %>%
  arrange(start.date, Date) %>%
  mutate(start.1= ifelse(Date==start.date & Result=="1", 1, 0),
        start.positive=max(start.1),
        number.positive= sum(as.numeric(as.character(Result)), na.rm = TRUE))
#Count number of unique patients; will be used as y-axis
number.patients <- length(unique(ctDNA.ID[["ID"]])) #[ gives a data.frame(list), [[ gives a vector

#ctDNA.ID <-   merge(ctDNA.ID, ctDNA.Index[,c(1:2)], by='ID', all.x = TRUE)
#ctDNA.ID$Date <- as.Date(as.character(ctDNA.ID$Date), format="%Y-%m-%d")

#ID.result.num <- as.factor(ctDNA.ID$Result)
#ID.mets.num <- as.factor(ctDNA.ID$Metastasis)
#ID.end.date <- ctDNA.ID %>% group_by(ID) %>% as.Date(max(Date), format="%Y-%m-%d")
#ID.end.date <- ctDNA.ID %>% group_by(ID) %>% max(ctDNA.ID$Date)
#ID.start.date <- ctDNA.ID %>% group_by(ID) %>% as.Date.numeric(min(ctDNA.ID$Date))

#ID.length.days <- ctDNA.ID %>% group_by(ID) %>% ID.end.date - ID.start.date

#ctDNA.ID$Brain_Mets <- as.factor(ctDNA.ID$Brain_Mets)
#ctDNA.ID.Plot <- arrange(ctDNA.ID, start.date, Date) 
#ctDNA.HiddenRows <- ctDNA.ID[!complete.cases("start.date","days"), ] #To see which rows have NA
#Graph4
graph4 =ctDNA.ID %>%
  ggplot(aes(x=days, y=reorder(ID,desc(start.date)), size=Result, color=Result)) +
  geom_line(aes(color=Metastasis), size= 2, alpha=0.75) +
  geom_point(aes(group = seq_along(Date)), # needed, otherwise transition doesn't work
             alpha =0.9, stroke= 1, 
             color = pal_jco()(5)[c(1)])+
  labs(title="ctDNA Detection Timeline",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Days",
       y= "Patients",
       color="Metastasis",
       size="ctDNA Detection",
       caption = paste("ctDNA Data from",  min(ctDNA.ID$Date), "to",max(ctDNA.ID$Date))) +
  guides(size=guide_legend(reverse = TRUE),
       color=guide_legend(reverse = FALSE)) +
  #theme_fivethirtyeight() +
  #theme_economist()+
  #theme_gdocs() +
  #theme_stata()+
  #theme_solarized()+
  #theme_calc()+
  theme_excel_new()+
  theme(title=element_text(color = "#657b83"),
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.title.y=element_blank(),
        plot.title.position = "panel",
        plot.title = element_text(hjust=0.5, vjust = 1),
        plot.subtitle = element_text(hjust=0.5, vjust = 1),
        legend.position = "bottom",
        #legend.position = c(0.3, 0),
        legend.title =  element_text(inherit.blank = FALSE, size = 10, color = "#000000"),
        legend.box = "horizontal") +
  scale_color_manual(values = pal_jama("default")(5)[c(2,4,3)], labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  #scale_color_manual(values=c("#425903", "#A5512B", "#632526"),
  #scale_color_jama(breaks=c(0, 1, 2),labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  #scale_color_manual(values=c("#36A6E5", "#D76A38", "#A23E52"),breaks=c(0, 1, 2),
  #                   labels=c("Adjuvant", "Metastatic", "Brain Mets")) +
  scale_size_manual(values = c(0, 3),
                    breaks = c(0, 1),
                    labels=c("Negative", "Positive")) +
  #scale_x_continuous(limits=c(0,490), breaks = seq(0, 490, by=90)) +
  scale_x_continuous(breaks=seq(0, max(as.numeric(ctDNA.ID$duration)), by=60))+
  scale_y_discrete(labels = number.patients:1) #+
  #facet_grid(start.positive ~ .)
  #facet_grid(Metastasis ~ ., scales="free")
#graph4
#Graph4 Animation
graph4.animation = graph4 +
  transition_reveal(Date) +
  enter_grow() + 
  enter_fade() 

#animate(graph4.animation, height= 500, width= 800, fps = 20, duration = 10, end_pause = 60)
animate(graph4.animation, height= 700, width= 935, fps = 40, duration = 20, start_pause= 10, end_pause = 220)
#anim_save("animation4.gif", animation = last_animation(), path = NULL)


#Graph 5 ID -Days from Beginning of Study -Mets Grouped------------------------------------------

#Data Manipulation

ctDNA.ID <- merge(x=ctDNA.Biocartis, y=ctDNA.Index[,c("ID",
                                            "Death", 
                                            "Metastasis")], 
                  by='ID', all.x = TRUE) %>%
  group_by(ID) %>%
  mutate(start.date=min(Date),
         days=Date-min(Date),
         duration=max(Date)-min(Date)) %>%
  arrange(start.date, Date) %>%
  mutate(start.1= ifelse(Date==start.date & Result=="1", 1, 0),
         start.positive=max(start.1),
         number.positive= sum(as.numeric(as.character(Result)), na.rm = TRUE))
  #Count number of unique patients; will be used as y-axis
  number.patients <- length(unique(ctDNA.ID[["ID"]])) #[ gives a data.frame(list), [[ gives a vector

graph5 =ctDNA.ID %>%
  #ggplot(aes(x=days, y=reorder(ID,desc(start.date)), size=Result, color=Result)) +
  ggplot(aes(x=days, y=reorder(ID, -as.numeric(Metastasis), mean), size=Result, color=Result)) +
  geom_line(aes(color=Metastasis), size= 2, alpha=0.75) +
  geom_point(aes(group = seq_along(Date)), # needed, otherwise transition doesn't work
      alpha =0.9, stroke= 1,
      #color = c("forestgreen")) + 
      #color = pal_jco()(1))+
      color = pal_jco()(5)[c(1)])+
  labs(title="ctDNA BRAF Detection Timeline",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Days",
       y= "Patients",
       color="Metastasis:",
       size="ctDNA Detection:",
       caption = paste("ctDNA Data from",  min(ctDNA.ID$Date), "to",max(ctDNA.ID$Date))) +
  guides(size=guide_legend(reverse = TRUE),
  color=guide_legend(reverse = FALSE)) +
  theme_excel_new()+ #theme_fivethirtyeight() +#theme_solarized()+ #theme_economist()+
  theme(title=element_text(color = "#657b83"),
      axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      axis.title.y = element_blank(),
      plot.title.position = "panel",
      plot.title = element_text(hjust=0.5, vjust = 1),
      plot.subtitle = element_text(hjust=0.5, vjust = 1),
      legend.position = "bottom",
      #legend.position = c(0.3, 0),
      legend.title =  element_text(inherit.blank = FALSE, size = 10, color = "#000000"),
      legend.box = "horizontal") +
  scale_color_manual(values = pal_jama("default")(5)[c(2,4,3)], labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  #scale_color_manual(values=c("#425903", "#A5512B", "#632526"),

  #scale_color_jama(breaks=c(0, 1, 2), labels=c("Adjuvant", "Metastatic", "Brain Mets"))+

  #scale_color_manual(values=c("#36A6E5", "#D76A38", "#A23E52"),breaks=c(0, 1, 2),
      #labels=c("Adjuvant", "Metastatic", "Brain Mets")) +
  scale_size_manual(values = c(0, 3),
                    breaks = c(0, 1),
                    labels=c("Negative", "Positive")) +
  #scale_x_continuous(limits=c(0,490), breaks = seq(0, 490, by=90)) +
  scale_x_continuous(breaks=seq(0, max(as.numeric(ctDNA.ID$duration)), by=60))+
  scale_y_discrete(labels = number.patients:1) #+
  #facet_grid(start.positive ~ ., scales = "free", space="free")#+
  #facet_grid(Metastasis ~ ., scales="free")

#Graph5 Animation
  # graph5.animation = graph5 +
  #   transition_reveal(Date) +
  #   enter_grow() +
  #   enter_fade()
    
#animate(graph5.animation, height= 500, width= 800, fps = 20, duration = 10, end_pause = 60)
#animate(graph5.animation, height= 576, width= 1024, fps = 20, duration = 15, start_pause= 10, end_pause = 120)
#anim_save("animation5.gif", animation = last_animation(), path = NULL)





#Graph 6 Unfaceted categorization Start Negative vs Positive===============================
ctDNA.ID <- merge(x=ctDNA.Biocartis, y=ctDNA.Index[,c("ID",
                                            "Sex",
                                            "Death", 
                                            "Metastasis")], 
                  by='ID', all.x = TRUE) %>%
  group_by(ID) %>%
  mutate(start.date=min(Date),
         days=Date-min(Date),
         duration=max(Date)-min(Date)) %>%
  #arrange(start.date, Date) %>%
  mutate(start.1= ifelse(Date==start.date & Result=="1", 1, 0),
         start.positive=max(start.1),
         number.positive= sum(as.numeric(as.character(Result)), na.rm = TRUE)) %>%
  #ungroup() %>%
  group_by(start.positive) %>%
  arrange(Metastasis, .by_group = TRUE)
#Count number of unique patients; will be used as y-axis
number.patients <- length(unique(ctDNA.ID[["ID"]])) #[ gives a data.frame(list), [[ gives a vector

y.order = ctDNA.ID %>% 
  group_by(start.positive, duration, ID) %>% #Adding duration here makes the IDs ordered per it
  summarise(mets.order=mean(as.numeric(as.character(Metastasis)))) %>% 
  #summarise(duration.order=mean(as.numeric(duration))) %>% 
  ungroup() %>%
  arrange(desc(start.positive), desc(mets.order)) %>% 
  #arrange(desc(.data[["start.positive"]]), desc(.data[["mets.order"]]), .data[["duration"]]) %>% 
  #arrange(desc(start.positive), mets.order, duration.order, .by_group = TRUE) %>% 
  pull(ID)
ctDNA.ID$ID = factor(ctDNA.ID$ID, levels=y.order)

graph6 =ctDNA.ID %>%
  arrange(start.positive, Metastasis) %>%
  #ggplot(aes(x=days, y=reorder(ID,desc(start.date)), size=Result, color=Result)) +
  #ggplot(aes(x=days, y=reorder(ID, -as.numeric(Metastasis), mean), size=Result, color=Result)) +
  ggplot(aes(x=days, y=ID, size=Result, color=Result)) +
  geom_line(aes(color=Metastasis), size= 2, alpha=0.75) +
  geom_point(aes(group = seq_along(Date)), # needed, otherwise transition doesn't work
             alpha =0.9, stroke= 1,
             #color = c("forestgreen")) + 
             #color = pal_jco()(1))+
             color = pal_jco()(5)[c(1)])+
  labs(title="ctDNA BRAF Detection Timeline",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Days",
       y= "Patients",
       color="Metastasis:",
       size="ctDNA Detection:",
       caption = paste("ctDNA Data from",  min(ctDNA.ID$Date), "to",max(ctDNA.ID$Date))) +
  guides(size=guide_legend(reverse = TRUE),
         color=guide_legend(reverse = FALSE)) +
  theme_excel_new()+ #theme_fivethirtyeight() +#theme_solarized()+ #theme_economist()+
  theme(title=element_text(color = "#657b83"),
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.title.y = element_blank(),
        plot.title.position = "panel",
        plot.title = element_text(hjust=0.5, vjust = 1),
        plot.subtitle = element_text(hjust=0.5, vjust = 1),
        legend.position = "bottom",
        #legend.position = c(0.3, 0),
        legend.title =  element_text(inherit.blank = FALSE, size = 10, color = "#000000"),
        legend.box = "horizontal") +
  scale_color_manual(values = pal_jama("default")(5)[c(2,4,3)], labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  #scale_color_manual(values=c("#425903", "#A5512B", "#632526"),
  
  #scale_color_jama(breaks=c(0, 1, 2), labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  
  #scale_color_manual(values=c("#36A6E5", "#D76A38", "#A23E52"),breaks=c(0, 1, 2),
  #labels=c("Adjuvant", "Metastatic", "Brain Mets")) +
  scale_size_manual(values = c(0, 3),
                    breaks = c(0, 1),
                    labels=c("Negative", "Positive")) +
  #scale_x_continuous(limits=c(0,490), breaks = seq(0, 490, by=90)) +
  scale_x_continuous(breaks=seq(0, max(as.numeric(ctDNA.ID$duration)), by=60))+
  scale_y_discrete(labels = number.patients:1) #+
#facet_grid(start.positive ~ ., scales = "free")#+
#facet_grid(Metastasis ~ ., scales="free")


#Graph 7 Metastatic categorization Ordered============================
#Data Manipulation

ctDNA.ID <- merge(x=ctDNA.Biocartis, y=ctDNA.Index[,c("ID",
                                                      "Death", 
                                                      "Metastasis")], 
                  by='ID', all.x = TRUE) %>%
  group_by(ID) %>%
  mutate(start.date=min(Date),
         days=Date-min(Date),
         duration=max(Date)-min(Date)) %>%
  arrange(start.date, Date) %>%
  mutate(start.1= ifelse(Date==start.date & Result=="1", 1, 0),
         start.positive=max(start.1),
         number.positive= sum(as.numeric(as.character(Result)), na.rm = TRUE))
#Count number of unique patients; will be used as y-axis
number.patients <- length(unique(ctDNA.ID[["ID"]])) #[ gives a data.frame(list), [[ gives a vector

length.order = ctDNA.ID %>% 
  group_by(Metastasis, duration, ID) %>% #Adding duration here makes the IDs ordered per it
  summarise(mets.order=mean(as.numeric(as.character(Metastasis)))) %>% 
  #summarise(duration.order=mean(as.numeric(duration))) %>% 
  #ungroup() %>%
  arrange(desc(mets.order), duration) %>% 
  #arrange(desc(.data[["start.positive"]]), desc(.data[["mets.order"]]), .data[["duration"]]) %>% 
  #arrange(desc(start.positive), mets.order, duration.order, .by_group = TRUE) %>% 
  pull(ID)
ctDNA.ID$ID = factor(ctDNA.ID$ID, levels=length.order)

graph7 =ctDNA.ID %>%
  #ggplot(aes(x=days, y=reorder(ID,desc(start.date)), size=Result, color=Result)) +
  ggplot(aes(x=days, y=ID, size=Result, color=Result)) +
  geom_line(aes(color=Metastasis), size= 2, alpha=0.75) +
  geom_point(aes(group = seq_along(Date)), # needed, otherwise transition doesn't work
             alpha =0.9, stroke= 1,
             #color = c("forestgreen")) + 
             #color = pal_jco()(1))+
             color = pal_jco()(5)[c(1)])+
  labs(title="ctDNA BRAF Detection Timeline",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Days",
       y= "Patients",
       color="Metastasis:",
       size="ctDNA Detection:",
       caption = paste("ctDNA Data from",  min(ctDNA.ID$Date), "to",max(ctDNA.ID$Date))) +
  guides(size=guide_legend(reverse = TRUE),
         color=guide_legend(reverse = FALSE)) +
  theme_excel_new()+ #theme_fivethirtyeight() +#theme_solarized()+ #theme_economist()+
  theme(title=element_text(color = "#657b83"),
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.title.y = element_blank(),
        plot.title.position = "panel",
        plot.title = element_text(hjust=0.5, vjust = 1),
        plot.subtitle = element_text(hjust=0.5, vjust = 1),
        legend.position = "bottom",
        #legend.position = c(0.3, 0),
        legend.title =  element_text(inherit.blank = FALSE, size = 10, color = "#000000"),
        legend.box = "horizontal") +
  scale_color_manual(values = pal_jama("default")(5)[c(2,4,3)], labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  #scale_color_manual(values=c("#425903", "#A5512B", "#632526"),
  
  #scale_color_jama(breaks=c(0, 1, 2), labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  
  #scale_color_manual(values=c("#36A6E5", "#D76A38", "#A23E52"),breaks=c(0, 1, 2),
  #labels=c("Adjuvant", "Metastatic", "Brain Mets")) +
  scale_size_manual(values = c(0, 3),
                    breaks = c(0, 1),
                    labels=c("Negative", "Positive")) +
  #scale_x_continuous(limits=c(0,490), breaks = seq(0, 490, by=90)) +
  scale_x_continuous(breaks=seq(0, max(as.numeric(ctDNA.ID$duration)), by=60))+
  scale_y_discrete(labels = number.patients:1)

#Graph 8 Timeline with followup time----
ctDNA.Timeline <- merge(x=ctDNA.Biocartis, 
                        y=ctDNA.Index[, c("ID",
                                          "Death", 
                                          "Metastasis",
                                          "Stage4",
                                          "StageBrain")], 
                  by='ID', all.x = TRUE) %>%
  group_by(ID) %>%
  mutate(start.date=min(Date),
         days=Date-min(Date),
         duration=max(Date)-min(Date),
         death.duration=Death-start.date,
         followup=ifelse(!is.na(Death), Death-start.date, Sys.Date()-start.date),
         start.1= ifelse(Date==start.date & Result=="1", 1, 0),
         start.positive=max(start.1),
         number.positive= sum(as.numeric(as.character(Result)), na.rm = TRUE)) %>% 
  mutate_at(c("days", "duration", "death.duration", "followup"), as.numeric) %>% 
  arrange(start.date, Date) #To order the table from older to newer recruitment

#Count number of unique patients; will be used as y-axis
number.patients <- length(unique(ctDNA.Timeline[["ID"]])) #[ gives a data.frame(list), [[ gives a vector

#Custom Order the IDs
length.order = ctDNA.Timeline %>% 
  group_by(Metastasis, duration, followup, ID) %>% #Adding duration here makes the IDs ordered per it
  #group_by(Metastasis, duration, ID) %>% #Adding duration here makes the IDs ordered per it
  summarise(mets.order=mean(as.numeric(as.character(Metastasis)))) %>% 
  #summarise(duration.order=mean(as.numeric(duration))) %>% 
  #ungroup() %>%
  arrange(desc(mets.order), duration, followup) %>% 
  #arrange(desc(.data[["start.positive"]]), desc(.data[["mets.order"]]), .data[["duration"]]) %>% 
  #arrange(desc(start.positive), mets.order, duration.order, .by_group = TRUE) %>% 
  pull(ID)
ctDNA.Timeline$ID = factor(ctDNA.Timeline$ID, levels=length.order)


#Graph
graph8 =ctDNA.Timeline %>% ggplot()+
  #ggplot(aes(x=days, y=reorder(ID,desc(start.date)), size=Result, color=Result)) +
  #ggplot(aes(x=days, y=ID, size=Result, color=Result)) +
  geom_line(aes(x=days, y=ID, color=Metastasis), size= 2, alpha=0.75) +
  geom_linerange(aes(xmin=duration, xmax=followup, y=ID, color=Metastasis), size= 0.3, alpha =0.8) +
  geom_point(aes(x=days, y=ID, size=Result, group = seq_along(Date)), # needed, otherwise transition doesn't work
             alpha =0.9, stroke= 1,
             #color = c("forestgreen")) + 
             #color = pal_jco()(1))+
             color = pal_jco()(5)[c(1)])+
  #geom_linerange(aes(xmin=duration, xmax=followup, y=ID, color=Metastasis), size= 0.3, alpha =0.8) +
  geom_point(aes(x=death.duration, y=ID), shape = 15, color = "red", fill= "#234168") +
  labs(title="ctDNA BRAF Detection Timeline",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Days",
       y= "Patients",
       color="Metastasis:",
       size="ctDNA Detection:",
       caption = paste("ctDNA Data from",  min(ctDNA.Timeline$Date), "to",max(ctDNA.Timeline$Date))) +
  guides(size=guide_legend(reverse = TRUE),
         color=guide_legend(reverse = FALSE)) +
  theme_excel_new()+ #theme_fivethirtyeight() +#theme_solarized()+ #theme_economist()+
  theme(title=element_text(color = "#657b83"),
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.title.y = element_blank(),
        plot.title.position = "panel",
        plot.title = element_text(hjust=0.5, vjust = 1),
        plot.subtitle = element_text(hjust=0.5, vjust = 1),
        legend.position = "bottom",
        #legend.position = c(0.3, 0),
        legend.title =  element_text(inherit.blank = FALSE, size = 10, color = "#000000"),
        legend.box = "horizontal",
        #panel.background = element_rect(fill = '#0E2439', color = 'purple'), #plot area color
        #plot.background = element_rect(fill = "#0E2439"), # Background color of the plot Darktheme
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y=element_blank()) +
  #scale_color_manual(values = pal_jama("default")(5)[c(2,4,3)], labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  scale_color_manual(values=c("#9ac026", "#E2A26E", "#35C0F0"), labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  #scale_color_jama(breaks=c(0, 1, 2), labels=c("Adjuvant", "Metastatic", "Brain Mets"))+
  #scale_color_manual(values=c("#36A6E5", "#D76A38", "#A23E52"),breaks=c(0, 1, 2),labels=c("Adjuvant", "Metastatic", "Brain Mets")) +
  #scale_color_jama()+
  scale_size_manual(values = c(0, 3),
                    breaks = c(0, 1),
                    labels=c("Negative", "Positive")) +
  #scale_x_continuous(limits=c(0,490), breaks = seq(0, 490, by=90)) +
  scale_x_continuous(breaks=seq(0, max(as.numeric(ctDNA.Timeline$duration)), by=60)) #+
  #scale_y_discrete(labels = number.patients:1)
#E7AD7B Brown
#9ac026 lemon green
#73D393 Torques green
#0E2439 Indiehacker color
#scales::show_col(pal_nejm("default")(10))
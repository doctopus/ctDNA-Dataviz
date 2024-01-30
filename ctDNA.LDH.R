#LDH Plot
library(readxl)
library(dplyr)
library(ggthemes)
library(ggplot2)
library(ggbeeswarm)
library(ggsci)
library(plotly)

ctDNA.LDH <-  read_excel("ctDNA.xlsx", sheet="LDH")%>% 
  mutate_at(c("ID"), as.factor) %>% 
  mutate_at(c("LDH"), as.numeric)

ctDNA.LDH$Date <- as.Date(as.numeric(ctDNA.LDH$Date), origin = "1899-12-30")

ctDNA.Biocartis <-  read_excel("ctDNA.xlsx", sheet="Biocartis")%>% 
  mutate_at(c("Patient", "ID", "Result"), as.factor) %>% 
  mutate_at(c("Date"), as.Date)

ctDNA.Index <- read_excel("ctDNA.xlsx", sheet="Index") %>% 
  mutate_at(c("ID", "Metastasis"), as.factor) %>% 
  mutate_at(c("Death", "Stage4"), as.Date)

#Beeswarm Plot Data Prep----
#Count number of unique patients; will be used as y-axis
number.patients <- length(unique(ctDNA.LDH[["ID"]])) #[ gives a data.frame(list), [[ gives a vector


ctDNA.LDH.ordered <- merge(x=ctDNA.LDH, y=ctDNA.Index[,c("ID",
                                          "Death", 
                                          "Metastasis")], 
                  by='ID', all.x = TRUE)

mets.order =  ctDNA.LDH.order  %>%
  group_by(Metastasis, ID) %>%
  summarise(mets.order=mean(as.numeric(as.character(Metastasis)))) %>%
  arrange(desc(mets.order), ID) %>% 
  # mutate(start.date=min(Date),
  #        days=Date-min(Date),
  #        duration=max(Date)-min(Date)) %>%
  pull(ID)

ctDNA.LDH.ordered$ID = factor(ctDNA.LDH.ordered$ID, levels=mets.order)

#Rename Mets
ctDNA.LDH.ordered$Metastasis <- factor(ctDNA.LDH.ordered$Metastasis,
                                               labels = c("Adjuvant",
                                                          "Metastasis to Other Organs", 
                                                          "Metastasis to Brain"))

#Beeswarm Plot -----
ggplot(ctDNA.LDH.ordered, aes(x = LDH, y = ID, color=Metastasis)) +
  geom_beeswarm(groupOnX = FALSE, priority = "density")+
  geom_vline(xintercept = 135, linetype = "longdash") +
  geom_vline(xintercept = 225, linetype = "longdash") +
  #coord_flip() ++
  theme_minimal()+
  scale_x_continuous(breaks=seq(0, max(ctDNA.LDH.ordered$LDH, na.rm=T), by=500))+
  scale_y_discrete(labels = number.patients:1)+

labs(title="LDH levels of all time among treatment groups",
     #subtitle="Detection of BRAF mutation in ctDNA",
     x= "LDH",
     y= "Patients",
     color="Metastasis:") +
theme(title=element_text(color = "#657b83"),
      axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      axis.title.y = element_blank(),
      plot.title.position = "panel",
      plot.title = element_text(hjust=0.5, vjust = 1),
      plot.subtitle = element_text(hjust=0.5, vjust = 1),
      legend.position = "bottom",
      #legend.position = c(0.3, 0),
      legend.title =  element_text(inherit.blank = FALSE, size = 10, color = "#000000"),
      legend.box = "horizontal")


##Line Chart of Timeline Data Prep-----

#Combining Biocartis + ID to get Start Date
ctDNA.Biocartis.Index <- merge(x=ctDNA.Biocartis, y=ctDNA.Index[,c("ID",
                                                                   "Death", 
                                                                   "Metastasis",
                                                                   "Stage4")], 
                         by='ID', all.x = TRUE) %>%
  group_by(ID) %>%
  mutate(start.date=min(Date))

ctDNA.BI.unique <- unique(ctDNA.Biocartis.Index[,c("ID",
                                                   "Death", 
                                                   "Metastasis",
                                                   "Stage4",
                                                   "start.date")]) 

#====All LDH data with reference to start of Biocartis
ctDNA.LDH.Biocartis.Index <- merge(x=ctDNA.LDH, y=ctDNA.BI.unique[,c("ID",
                                                                     "Death", 
                                                                     "Metastasis",
                                                                     "Stage4",
                                                                     "start.date")], 
                                   by='ID', all.x = TRUE) %>%

  mutate(biocartis.start.to.ldh= difftime(Date, start.date,units = "weeks"),
         biocartis.start.to.death=difftime(Death, start.date, units = "weeks"),
         biocartis.start.to.stage4=difftime(Stage4, start.date, units = "weeks"))

#Rename Facets
ctDNA.LDH.Biocartis.Index$Metastasis <- factor(ctDNA.LDH.Biocartis.Index$Metastasis,
                                                               labels = c("Adjuvant",
                                                                          "Metastasis to Other Organs", 
                                                                          "Metastasis to Brain"))




#==LDH data after Stage IV wrt start of Biocartis
ctDNA.LDH.Biocartis.Index.After.Stage4 <- filter(ctDNA.LDH.Biocartis.Index, biocartis.start.to.stage4 < biocartis.start.to.ldh)

#Rename Facets
ctDNA.LDH.Biocartis.Index.After.Stage4$Metastasis <- factor(ctDNA.LDH.Biocartis.Index.After.Stage4$Metastasis,
                                                               labels = c("Adjuvant",
                                                                          "Metastasis to Other Organs", 
                                                                          "Metastasis to Brain"))




#==LDH data only after the start of study
ctDNA.LDH.Biocartis.Index.After.Biocartis <- filter(ctDNA.LDH.Biocartis.Index, biocartis.start.to.ldh > 0)

#Rename Facets
ctDNA.LDH.Biocartis.Index.After.Biocartis$Metastasis <- factor(ctDNA.LDH.Biocartis.Index.After.Biocartis$Metastasis,
                                                               labels = c("Adjuvant",
                                                                          "Metastasis to Other Organs", 
                                                                          "Metastasis to Brain"))
# Graph of All LDH----
graph.ldh <- ctDNA.LDH.Biocartis.Index %>%
    ggplot(aes(x=biocartis.start.to.ldh, y=LDH)) +
    geom_line(aes(color= ID), size= 0.5, alpha=0.75) +
    
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 135, ymax = 225, fill = "#257FFF", alpha = .2, color = NA)+
    theme_minimal()+
    theme(legend.position = "none") +
    scale_y_continuous(trans='log10')+
    scale_x_continuous()+
  
    facet_grid(Metastasis ~ .)+
    #facet_grid(Metastasis ~ ., scales = "free")+
  
  labs(title="All-time LDH Levels in Treatment Groups",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Weeks",
       y= "LDH in Log Scale",
       color="Metastasis:") +  
  theme(title=element_text(color = "#657b83"),
      axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      plot.title = element_text(hjust=0.5, vjust = 1),
      strip.text = element_text(color = "#657b83")
      )
#Make Graph Interactive
#ggplotly(graph.ldh)

# Graph After Stage4----
graph.ldh.stage4 <- ctDNA.LDH.Biocartis.Index.After.Stage4 %>%
  ggplot(aes(x=biocartis.start.to.ldh, y=LDH)) +
  geom_line(aes(x=biocartis.start.to.ldh, y=LDH, color=ID), size= 0.7, alpha= 0.75)+
  
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 135, ymax = 225, fill = "#257FFF", alpha = .2, color = NA)+
  theme_minimal()+
  theme(legend.position = "none") +
  scale_y_continuous(trans='log10')+
  scale_x_continuous()+
  
  facet_grid(Metastasis ~ .)+
  #facet_grid(Metastasis ~ ., scales = "free")
  
  labs(title="LDH Levels Since Stage4 Diagnosis in Treatment Groups",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Weeks",
       y= "LDH in Log Scale",
       color="Metastasis:") +  
  theme(title=element_text(color = "#657b83"),
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        plot.title = element_text(hjust=0.5, vjust = 1),
        strip.text = element_text(color = "#657b83")
  )

#Make Graph Interactive
#ggplotly(graph.ldh)


# Graph After Biocartis----
graph.ldh.biocartis <- ctDNA.LDH.Biocartis.Index.After.Biocartis %>%
  ##ctDNA.LDH.Biocartis.Index %>%
  ggplot(aes(x=biocartis.start.to.ldh, y=LDH)) +
  geom_line(aes(color= ID), size= 0.7, alpha=0.75) +
  
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 135, ymax = 225, fill = "#257FFF", alpha = .2, color = NA)+
  theme_minimal()+
  theme(legend.position = "none") +
  scale_y_continuous(trans='log10')+
  scale_x_continuous()+
  
  facet_grid(Metastasis ~ .)+
  #facet_grid(Metastasis ~ ., scales = "free")
  
  labs(title="LDH Levels During the Study in Treatment Groups",
       #subtitle="Detection of BRAF mutation in ctDNA",
       x= "Weeks",
       y= "LDH in Log Scale",
       color="Metastasis:") +  
  theme(title=element_text(color = "#657b83"),
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        plot.title = element_text(hjust=0.5, vjust = 1),
        strip.text = element_text(color = "#657b83")
  )

#Make Graph Interactive
#ggplotly(graph.ldh)


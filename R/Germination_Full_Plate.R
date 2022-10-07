library(readxl)
library(ggplot2)
library(ggpubr)

t = seq(0,120,1)

plate = read_xlsx("Downloads/20220707_germinationassay.xlsx",sheet=2) %>% 
  mutate(Time=t) %>% gather(key=Well,value = OD,-Time) %>% 
  mutate(Row = substr(Well,1,1), Column = substr(Well,2,3),
         Condition = ifelse(Column %in% seq(1:4),"Control","10mM Alanine"))

strains = c("PY79","JLG370","JLG413","JLG3140","JLG2641","JLG3162","JLG2457","Blank")
rows = unique(plate$Row)

platelist = list()

for (i in 1:8){
  
  plate2 = plate %>% filter(Row == rows[i]) %>% mutate(Strain = strains[i])
  
  platelist[[i]] = plate2
  
  rm(plate2)
}

plate_final = bind_rows(platelist) %>% select(Time,OD,Condition,Strain)

plate_summ = plate_final %>% group_by(Time,Condition,Strain) %>% 
  summarise(MeanOD = mean(OD), SD=sd(OD))

control = plate_summ %>% filter(Condition == "Control")
alanine = plate_summ %>% filter(Condition != "Control")

cplot = ggplot(control, aes(x=Time,y=MeanOD,color=Strain)) + theme_bw() +
  geom_point() + geom_line() + ylab("OD 580") + xlab("Time (min)") +
  ggtitle("Control") + ylim(0.1,0.4) + xlim(0,30)

aplot = ggplot(alanine, aes(x=Time,y=MeanOD,color=Strain)) + theme_bw() +
  geom_point() + geom_line() + ylab("OD 580") + xlab("Time (min)") +
  ggtitle("10mM Alanine")  + ylim(0.1,0.4) + xlim(0,30)

ggarrange(cplot,aplot,nrow = 1,common.legend = TRUE,legend = "bottom")


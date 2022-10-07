###### Package installation and Library loading ######
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","readxl","ggforce","ggpubr")

###### Data loading and preparation #####

# The function %in% allows to filter for values on a vector
# Negating it allows us to filter for values that are not in the vector
'%ni%' <- Negate("%in%")

# Read the excel file that contains the measurements
snake = read_xlsx("Desktop/ECB/Snakes_cell_shape/Measures.xlsx", sheet=1) 

# Calculate Surface Area (SA), Volume (V) and Surface to Volume ratio (AR) of the cells

snake = snake %>% # %>% is called the "pipeline" function 
  # This allows to link multiple functions of a data frame in one line of code
  mutate(# mutate creates a new variable and allows for operations with other columns
    SA=pi*Snake_Length*Snake_Width_Mean,
    V = ((pi*Snake_Length*Snake_Width_Mean^2)/4)-(pi*Snake_Width_Mean^3)/12,
    AR = SA/V)

# All strains with the JLG names as used on the microscopy files

# We know that the strains will have at most 7 characters (JLG+max.4 numbers)
# With a simple loop we can filter create a vector of just strain names

# First, create an empty vector where the strain names will be stored
Strains= c()

# Create a position variable that will allow us to store data in the vector
# This will increase in 1 every time we add a new name to the vector
k=1

# We start by analyzing the unique names of our Species column (which contains species and strains)
# If we go one by one (i.e every name) it would take longer and we would need to compress to the unique ones anyway
for (i in 1:length(unique(snake$Species))) {
# Store the name in anew variable
  Species = unique(snake$Species)[[i]]
  # The "if" statement will proceed only when the condition is fulfilled
  # In this case, if the number of characters is lower or equal to 7
  if (nchar(Species) <= 7) {
  # This will add the strain name to the vector
  Strains[k] = Species   
  # And increase k so as to not overwrite the existing name
  k = k+1    
  }
}
rm(k,Species)

PBPdel = snake %>% filter(Species %in% Strains & !is.na(AR) & Type != "ClassB" & Snake_Length<6) %>% 
  mutate(Genotype = fct_reorder(Genotype, Snake_Width_Mean, .fun = 'mean')) %>% rename(Strain=Species)

snake = snake %>% filter(Species %ni% Strains & !is.na(AR)) %>% 
  mutate(Species = fct_reorder(Species, Snake_Width_Mean, .fun = 'mean'))

sumtable = snake %>% group_by(Species) %>% 
  summarise(N=n(),MeanL=mean(Snake_Length),SDL=sd(Snake_Length),
            Meanw=mean(Snake_Width_Mean),SDw=sd(Snake_Width_Mean),
            MeanSA=mean(SA),SDSA=sd(SA),
            MeanV=mean(V),SDV=sd(V),
            MeanAR=mean(AR),SDAR=sd(AR)) 

###### Species Violin Plots #####
SAplot = ggplot(snake,aes(x=Species,y=SA,fill=Species)) + theme_bw() +
  xlab("Species") + ylab("Surface Area (µm^2)") +
  geom_violin() + ylim(c(0,30)) + xlab("") +
  scale_fill_brewer(type="qual",palette="Dark2") +
  stat_summary(fun=median, geom="crossbar", size=1,
               color="black",position = position_dodge(width = 0.9)) + 
  scale_color_brewer(type="qual",palette="Dark2") + 
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x = element_blank(),legend.text = element_text(size=50),
        legend.title = element_text(size=50))

Vplot = ggplot(snake,aes(x=Species,y=V,fill=Species)) + theme_bw() +
  xlab("Species") + ylab("Volume (µm^3)") +
  geom_violin() + ylim(c(0,10)) + xlab("") +
  scale_fill_brewer(type="qual",palette="Dark2") +
  stat_summary(fun=median, geom="crossbar", size=1,
               color="black",position = position_dodge(width = 0.9)) + 
  scale_color_brewer(type="qual",palette="Dark2") + 
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x = element_blank(),legend.text = element_text(size=50),
        legend.title = element_text(size=50))

lplot=ggplot(snake,aes(x=Species,y=Snake_Length,fill=Species)) + theme_bw() +
  xlab("Species") + ylab("Cell length (µm)") +
  geom_violin() + ylim(c(0,7.5)) + xlab("") +
  scale_fill_brewer(type="qual",palette="Dark2") +
  stat_summary(fun=median, geom="crossbar", size=1,
               color="black",position = position_dodge(width = 0.9)) + 
  scale_color_brewer(type="qual",palette="Dark2") + 
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x = element_blank(),legend.text = element_text(size=50),
        legend.title = element_text(size=50))

wplot=ggplot(snake,aes(x=Species,y=Snake_Width_Mean,fill=Species)) + theme_bw() +
  xlab("Species") + ylab("Cell width (µm)") +
  geom_violin() + ylim(c(0,2)) + xlab("") +
  scale_fill_brewer(type="qual",palette="Dark2") +
  stat_summary(fun=median, geom="crossbar", size=1,
               color="black",position = position_dodge(width = 0.9)) + 
  scale_color_brewer(type="qual",palette="Dark2") + 
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x = element_blank(),legend.text = element_text(size=50),
        legend.title = element_text(size=50))

ggarrange(lplot,wplot,SAplot,Vplot,common.legend = TRUE,legend="bottom")
#Save as a 4000 width x 1600 height file

###### Species Scatter plots ######
wSAplot = ggplot(snake,aes(x=Snake_Width_Mean,y=SA)) + stat_smooth(inherit.aes = F, data=snake,aes(x=Snake_Width_Mean,y=SA),color="black",lty=2,se=F,size=5) +
  geom_jitter(aes(color=Species),size=3) + theme_bw() +
  xlab("Cell width (µm)") + ylab("Surface Area (µm^2)") + 
  stat_smooth(method = "lm",formula= y~exp(-x),aes(color=Species),size=5,se=F) +
  geom_point(inherit.aes = F, data = sumtable,aes(x=Meanw,y=MeanSA),color="black",size=3)  +
  scale_color_brewer(type="qual",palette="Dark2") + ylim(c(0,30)) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))

wVplot = ggplot(snake,aes(x=Snake_Width_Mean,y=V)) + stat_smooth(inherit.aes = F, data=snake,aes(x=Snake_Width_Mean,y=V),color="black",lty=2,se=F,size=5) +
  geom_jitter(aes(color=Species),size=3) + theme_bw() +
  xlab("Cell width (µm)") + ylab("Volume (µm^3)") + 
  stat_smooth(method = "lm",formula= y~exp(-x),aes(color=Species),size=5,se=F) +
  geom_point(inherit.aes = F, data = sumtable,aes(x=Meanw,y=MeanAR),color="black",size=3)  +
  scale_color_brewer(type="qual",palette="Dark2") + ylim(c(0,10))+ 
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))

wARplot = ggplot(snake,aes(x=Snake_Width_Mean,y=AR)) + stat_smooth(inherit.aes = F, data=snake,aes(x=Snake_Width_Mean,y=AR),color="black",lty=2,se=F,size=5,formula = y~I(x^-1)) +
  geom_jitter(aes(color=Species),size=3) + theme_bw() +
  xlab("Cell width (µm)") + ylab("Surface to Volume ratio (1/µm)") + 
  stat_smooth(method = "lm",formula = y~I(x^-1),aes(color=Species),size=5,se=F) +
  geom_point(inherit.aes = F, data = sumtable,aes(x=Meanw,y=MeanAR),color="black",size=7)  +
  scale_color_brewer(type="qual",palette="Dark2") + ylim(c(2.5,9.5)) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))

lSAplot = ggplot(snake,aes(x=Snake_Length,y=SA)) + stat_smooth(inherit.aes = F, data=snake,aes(x=Snake_Length,y=SA),color="black",lty=2,method="gam",se=F,size=5) +
  geom_jitter(aes(color=Species),size=3) + theme_bw() +
  xlab("Cell length (µm)") + ylab("Surface Area (µm^2)") + ylim(c(0,30)) + 
  stat_smooth(method = "lm",formula=y ~ exp(x),aes(color=Species),size=5,se=F) +
  geom_point(inherit.aes = F, data = sumtable,aes(x=MeanL,y=MeanAR),color="black",size=3)  +
  scale_color_brewer(type="qual",palette="Dark2") + xlim(c(0,7.5)) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))

lVplot = ggplot(snake,aes(x=Snake_Length,y=V)) + stat_smooth(inherit.aes = F, data=snake,aes(x=Snake_Length,y=V),color="black",lty=2,method="gam",se=F,size=5) +
  geom_jitter(aes(color=Species),size=3) + theme_bw() +
  xlab("Cell length (µm)") + ylab("Volume (µm^2)") + ylim(c(0,10))+ 
  stat_smooth(method = "lm",formula=y ~ exp(x),aes(color=Species),size=5,se=F) +
  geom_point(inherit.aes = F, data = sumtable,aes(x=MeanL,y=MeanAR),color="black",size=3)  +
  scale_color_brewer(type="qual",palette="Dark2") + xlim(c(0,7.5)) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))

lARplot = ggplot(snake,aes(x=Snake_Length,y=AR)) + stat_smooth(inherit.aes = F, data=snake,aes(x=Snake_Length,y=AR),color="black",lty=2,method="gam",se=F,size=5,formula = y~I(x^-1)) +
  geom_jitter(aes(color=Species),size=3) + theme_bw() +
  xlab("Cell length (µm)") + ylab("Surface to Volume ratio (1/µm)") + 
  stat_smooth(method = "lm",formula = y~I(x^-1),aes(color=Species),size=5,se=F) +
  geom_point(inherit.aes = F, data = sumtable,aes(x=MeanL,y=MeanAR),color="black",size=7)  +
  scale_color_brewer(type="qual",palette="Dark2") + xlim(c(0,7.5)) + ylim(c(2.5,9.5)) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))

# All metrics
ggarrange(wSAplot, lSAplot, wVplot, lVplot,wARplot,lARplot,
          nrow=3, ncol=2,common.legend = TRUE, legend="bottom")

# Surface to volume ratio
ggarrange(wARplot,lARplot,nrow=1,ncol=2,common.legend = TRUE,legend = "bottom")

###### PBPs violin plots #####
SAplot = ggplot(PBPdel,aes(x=Genotype,y=SA,fill=Genotype)) + theme_bw() +
  xlab("Genotype") + ylab("Surface Area (µm^2)") +
  geom_violin() + ylim(c(0,11)) + xlab("") +
  scale_fill_brewer(type="diverging",palette="Spectral") +
  stat_summary(fun=median, geom="crossbar", size=1,
               color="black",position = position_dodge(width = 0.9)) + 
  scale_color_brewer(type="qual",palette="Dark2") + 
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x = element_blank(),legend.text = element_text(size=50),
        legend.title = element_text(size=50))

Vplot = ggplot(PBPdel,aes(x=Genotype,y=V,fill=Genotype)) + theme_bw() +
  xlab("Genotype") + ylab("Volume (µm^3)") +
  geom_violin() + ylim(c(0,2.3)) + xlab("") +
  scale_fill_brewer(type="diverging",palette="Spectral") +
  stat_summary(fun=median, geom="crossbar", size=1,
               color="black",position = position_dodge(width = 0.9)) + 
  scale_color_brewer(type="qual",palette="Dark2") + 
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x = element_blank(),legend.text = element_text(size=50),
        legend.title = element_text(size=50))

lplot=ggplot(PBPdel,aes(x=Genotype,y=Snake_Length,fill=Genotype)) + theme_bw() +
  xlab("Genotype") + ylab("Cell length (µm)") +
  geom_violin() + ylim(c(1,4.1)) + xlab("") +
  scale_fill_brewer(type="diverging",palette="Spectral") +
  stat_summary(fun=median, geom="crossbar", size=1,
               color="black",position = position_dodge(width = 0.9)) + 
  scale_color_brewer(type="qual",palette="Dark2") + 
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x = element_blank(),legend.text = element_text(size=50),
        legend.title = element_text(size=50))

wplot=ggplot(PBPdel,aes(x=Genotype,y=Snake_Width_Mean,fill=Genotype)) + theme_bw() +
  xlab("Genotype") + ylab("Cell width (µm)") +
  geom_violin() + ylim(c(0.3,1.1)) + xlab("") +
  scale_fill_brewer(type="diverging",palette="Spectral") +
  stat_summary(fun=median, geom="crossbar", size=1,
               color="black",position = position_dodge(width = 0.9)) + 
  scale_color_brewer(type="qual",palette="Dark2") + 
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x = element_blank(),legend.text = element_text(size=50),
        legend.title = element_text(size=50))

ggarrange(lplot,wplot,SAplot,Vplot,common.legend = TRUE,legend="bottom")
#Save as a 4000 width x 1600 height file

###### PBPS scatter plots ####
wSAplot = ggplot(PBPdel,aes(x=Snake_Width_Mean,y=SA)) + stat_smooth(inherit.aes = F, data=PBPdel,aes(x=Snake_Width_Mean,y=SA),color="black",lty=2,se=F,size=5) +
  geom_jitter(aes(color=Genotype),size=5) + theme_bw() +
  xlab("Cell width (µm)") + ylab("Surface Area (µm^2)")+ 
  stat_smooth(method = "lm",formula= y~exp(-x),aes(color=Genotype),size=5) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50)) + 
  scale_color_brewer(type="diverging",palette="Spectral") +
  geom_mark_ellipse(aes(fill = PonA_del), alpha = 0.15) + xlim(c(0.3,1.1)) + ylim(c(1,11))

wVplot = ggplot(PBPdel,aes(x=Snake_Width_Mean,y=V)) + stat_smooth(inherit.aes = F, data=PBPdel,aes(x=Snake_Width_Mean,y=V),color="black",lty=2,se=F,size=5) +
  geom_jitter(aes(color=Genotype),size=5) + theme_bw() +
  xlab("Cell width (µm)") + ylab("Volume (µm^3)")+ 
  stat_smooth(method = "lm",formula= y~exp(-x),aes(color=Genotype),size=5) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50)) + 
  scale_color_brewer(type="diverging",palette="Spectral") +
  geom_mark_ellipse(aes(fill = PonA_del), alpha = 0.15) + xlim(c(0.3,1.1)) + ylim(c(0,2.3))

wARplot = ggplot(PBPdel,aes(x=Snake_Width_Mean,y=AR)) + stat_smooth(inherit.aes = F, data=PBPdel,aes(x=Snake_Width_Mean,y=AR),color="black",lty=2,se=F,size=5) +
  geom_jitter(aes(color=Genotype),size=5) + theme_bw() +
  xlab("Cell width (µm)") + ylab("Surface to Volume ratio (1/µm)")+ 
  stat_smooth(method = "lm",formula= y~exp(-x),aes(color=Genotype),size=5) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50)) + 
  scale_color_brewer(type="diverging",palette="Spectral") +
  geom_mark_ellipse(aes(fill = PonA_del), alpha = 0.15) + xlim(c(0.3,1.1)) + ylim(c(3.5,12)) 


lSAplot = ggplot(PBPdel,aes(x=Snake_Length,y=SA)) + stat_smooth(inherit.aes = F, data=PBPdel,aes(x=Snake_Length,y=SA),color="black",lty=2,se=F,size=5) +
  geom_jitter(aes(color=Genotype),size=5) + theme_bw() +
  xlab("Cell length (µm)") + ylab("Surface Area (µm^2)")+ 
  stat_smooth(method = "lm",formula= y~exp(-x),aes(color=Genotype),size=5) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50)) + 
  scale_color_brewer(type="diverging",palette="Spectral") +
  geom_mark_ellipse(aes(fill = PonA_del), alpha = 0.15) + xlim(c(1,6)) + ylim(c(1,11))

lVplot = ggplot(PBPdel,aes(x=Snake_Length,y=V)) + stat_smooth(inherit.aes = F, data=PBPdel,aes(x=Snake_Length,y=V),color="black",lty=2,se=F,size=5) +
  geom_jitter(aes(color=Genotype),size=5) + theme_bw() +
  xlab("Cell length (µm)") + ylab("Volume (µm^3)")+ 
  stat_smooth(method = "lm",formula= y~exp(-x),aes(color=Genotype),size=5) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50)) + 
  scale_color_brewer(type="diverging",palette="Spectral") +
  geom_mark_ellipse(aes(fill = PonA_del), alpha = 0.15) + xlim(c(1,6)) + ylim(c(0,2.3))

lARplot = ggplot(PBPdel,aes(x=Snake_Length,y=AR)) + stat_smooth(inherit.aes = F, data=PBPdel,aes(x=Snake_Length,y=AR),color="black",lty=2,se=F,size=5) +
  geom_jitter(aes(color=Genotype),size=5) + theme_bw() +
  xlab("Cell length (µm)") + ylab("Surface to Volume ratio (1/µm)")+ 
  stat_smooth(method = "lm",formula= y~exp(-x),aes(color=Genotype),size=5) +
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50)) + 
  scale_color_brewer(type="diverging",palette="Spectral") +
  geom_mark_ellipse(aes(fill = PonA_del), alpha = 0.15) + xlim(c(1,6)) + ylim(c(3.5,12))

# All metrics
ggarrange(wSAplot, lSAplot, wVplot, lVplot,wARplot,lARplot,
          nrow=3, ncol=2,common.legend = TRUE, legend="bottom")

# Surface to volume ratio
ggarrange(wARplot,lARplot,nrow=1,ncol=2,common.legend = TRUE,legend = "bottom")

#Figure in report
ggarrange(wplot,wARplot,nrow=1,ncol=2,common.legend = TRUE,legend="bottom")

###### Fit the width and S/V ratio data ######

# We assume an inverse relationship between width and S/V ratio
x = snake$Snake_Width_Mean
y = snake$AR

fit = lm(y ~ I(x^-1), snake)

# We can then see the relationship between width and S/V
gamma = fit$coefficients[[2]]
# Considering our function, this means that S/V = 4.49/w
# Next, we calculate the standard error 
k=length(fit$coefficients)-1 #Subtract one to ignore intercept
SSE=sum(fit$residuals**2)
n=length(fit$residuals)
SE = sqrt(SSE/(n-(1+k))) #Residual Standard Error

# Round standard error and value of fit to the last significant digit
for (i in 1:10){
  SE2 = round(SE,i)
  if (SE2 != 0){
    break
  }
}
# As this is the number were our error is, we cannot have more decimals than the error 
gamma = round(gamma,i)

# Our new formula is S/V = γ/w; γ = 4.5 ± 0.2

# We can then plot the data and add our new fit and the standard errors
ggplot(snake,aes(x=Snake_Width_Mean,y=AR)) + 
  geom_jitter(aes(color=Species),size=3) + theme_bw() +
  xlab("Cell width (µm)") + ylab("Surface to Volume ratio (1/µm)") + 
  stat_smooth(method = "lm",formula = y~I(x^-1),aes(color=Species),size=3,se=F) +
  scale_color_brewer(type="qual",palette="Dark2") + ylim(c(2,9.5)) +
  geom_line(aes(y=gamma/Snake_Width_Mean),size=2,lty=2) +
  geom_line(aes(y=((gamma-SE2)/Snake_Width_Mean)),size=2,lty=3) +
  geom_line(aes(y=((gamma+SE2)/Snake_Width_Mean)),size=2,lty=3) #+
  theme(axis.text=element_text(size=55),axis.title = element_text(size=50),
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))
  


###### Fit the width and S/V ratio data for each species ######
  
  fitData = data.frame(Species = unique(snake$Species), gamma=1:4,SSE=1:4)
  colorvec = RColorBrewer::brewer.pal(length(unique(snake$Species)), name = "Dark2")
  specPlot = list()
  
for (j in 1:4){
    
    spec = as.vector(unique(snake$Species))
    bact = snake %>% filter(Species == spec[j])  
    bactsum = sumtable %>% filter(Species == spec[j]) 
    
    # We assume an inverse relationship between width and S/V ratio
    x = bact$Snake_Width_Mean
    y = bact$AR
    
    fit = lm(y ~ I(x^-1), bact)
    
    # We can then see the relationship between width and S/V
    gamma = fit$coefficients[[2]]
    # Considering our function, this means that S/V = 4.49/w
    # Next, we calculate the standard error 
    k=length(fit$coefficients)-1 #Subtract one to ignore intercept
    SSE=sum(fit$residuals**2)
    n=length(fit$residuals)
    SE = sqrt(SSE/(n-(1+k))) #Residual Standard Error
    
    # Round standard error and value of fit to the last significant digit
    for (i in 1:10){
      SE2 = round(SE,i)
      if (SE2 != 0){
        break
      }
    }
    # As this is the number were our error is, we cannot have more decimals than the error 
    gamma = round(gamma,i)
    
    fitData$gamma[j] = gamma
    fitData$SSE[j] = SE2
    
    FitFormula = paste("S/V = γ/w; ","γ = ",gamma," ± ",SE2,sep="")
    
    p = ggplot(bact,aes(x=Snake_Width_Mean,y=AR,color=Species)) + theme_bw() + 
      geom_jitter(color=colorvec[j],size=4) + ggtitle(paste(spec[j])) +
      stat_smooth(method="lm",formula = y~I(x^-1),color=colorvec[j],se=F,size=4) +
      geom_point(inherit.aes = F, data = bactsum,aes(x=Meanw,y=MeanAR),color="black",size=3)  +
      xlab("Cell width (µm)") + ylab("Surface to Volume ratio (1/µm)") +
      geom_text(x=max(x)-0.1,y=max(y)-0.1,label=FitFormula,color="black",size=30) +
      theme(axis.text=element_text(size=60),axis.title = element_text(size=60),
            title = element_text(size=80))
    
    specPlot[[j]] = p
    
}
  
ggarrange(plotlist = specPlot, nrow=2,ncol = 2)

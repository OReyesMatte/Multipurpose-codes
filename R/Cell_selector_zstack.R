# This code processes cell imaging data in order to find the best focal plane for each single cell
# It consists in a series of nested loops to properly analyze each single image taken
## The data is obtained from analyzed images using the Fiji plugins JFilament and MicrobeJ
## The necessary data from the cells is the length (L), width (w), x and y centroids

# The data is 'nested' in the following levels:
# Medium -> Time point -> Condition -> Image -> Slice
# Keep these in order to properly process the data
# New experiments can add extra variables (i.e. genotype, multiple conditions) 
# If this happens, add loops in a new code as necessary

###### Package installation and Library loading #####
# pacman: installing and loading necessary packages
# tidyverse: data wrangling and plotting tools
# readxl: reading excel files
# modeest: calculator of data mode (most frequent value)
# ggpubr: extra plotting tools

if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","readxl","modeest",'ggpubr')

###### Custom functions ######

# We define a custom function to enumerate quantiles of different parameters
# It requires specifying a data frame and a column that has a numerical value
# Alternatively, you can use a specific probability and remove previously defined quantiles
quantilizer = function(data,column,probability=0.2,QuantRemover=FALSE){
  
  if (QuantRemover==TRUE){
    
    data = data %>% select(-contains('Quantile',ignore.case = TRUE))
  }
  
  a = quantile(column,probs=seq(0,1,probability),type = 7,names = FALSE)
  
  b = column
  
  dfs = list()
  
  for (i in 2:length(a)){
    if (i == 2){
      df = data %>% filter(b >= a[i-1] & b <= a[i]) %>% 
        mutate(Quantile = as.factor(paste('Quantile ',i-1,sep = ''))) 
      
      dfs[[i]] = df}
    else {
      df = data %>% filter(b > a[i-1] & b <= a[i]) %>% 
        mutate(Quantile = paste('Quantile ',i-1,sep = '')) 
      
      dfs[[i]] = df  
    }
  }
  
  df = bind_rows(dfs)
  
}

###### Data loading and preparation #####

snake_data = read_xlsx('~/Desktop/Bsubtilis_growthcurve_cellshape.xlsx', # Read the file
                       sheet=4 # Select the sheet where the data is stored
                       ) %>% 
  select(-c(Cell,Experiment,ID) # Discard columns without useful information
         ) %>% filter(w > 0.3 # Do a prior filter
                      )  %>% 
  mutate(temp_col = str_split(Name, stringr::fixed("_"),n=4)
         # Name has all the necessary info to sort the data
         # Take the info and create new columns
         ) %>%
  mutate(Medium = map_chr(temp_col,1), # Growth medium
         Time = map_chr(temp_col,2), # Time point
         Image = map_chr(temp_col,3), # Image number
         Condition = 'None' # We have no conditions, so we create an 'empty' column
         ) %>% 
  select(-c(temp_col,Name) # Discard columns without useful information
         ) %>% 
  mutate(Time = as.numeric(substr(Time,2,3)) # Take the actual numbers of the time points
         ) %>% # Transform to hour values
  mutate(Time = ifelse(Time <=20, Time*15/60, (15*20+(Time-20)*30)/60))

# Create an empty list to store processed data and plots of each image 

dataList = list()
plotL = list()

# Create two increasing indexes to not overwrite the data
h=1
i=1

###### Data processing ####

# Start the loop according to the different media available, thus splitting by this condition
for (a in 1:length(unique(snake_data$Medium))){

  # Define a new vector of the different media
  Media = unique(snake_data$Medium)
  # Use the positions of the vector to filter by each media in each iteration
  snakeM = snake_data %>% filter(Medium == Media[a])  

  # Define a new vector of the different conditions  
  Conditions = unique(snakeM$Condition)
  
  # Start the next loop to split by the different conditions
  for (b in 1:length(Conditions)){
  
    # Use the positions of the vector to filter by each condition in each iteration
    snakeC = snakeM %>% filter(Condition == Conditions[b]) 
    
    # Define a new vector of the different time points
    Times = unique(snakeM$Time)
    
    # Start the next loop to split by the different time points
    for (c in 1:length(Times)){
      
      # Use the positions of the vector to filter by each time point in each iteration
      snakeT = snakeC %>% filter(Time == Times[c])

      # Define a new vector of the different images  
      Images = unique(snakeT$Image)
      
      # Start the next loop to split by the different images taken
      for (d in 1:length(Images)){
  
        # Use the positions of the vector to filter by each image in each iteration
        snakeI = snakeT %>% filter(Image == Images[d])
  
        # Now the processing starts!
        
        # Calculate the number of cells in each slice and filter by the most frequent value (mode)
        cellNum = snakeI %>% group_by(Slice) %>% summarise(Ncell = n()) %>% 
        filter(Ncell %in% mlv(Ncell))

        # In case there are two possible values for the high frequency, maximize number of cells
        if (length(unique(cellNum$Ncell)) > 1){
        cellNum = cellNum %>% filter(Ncell == max(Ncell))}
        
        # Filter the slices to be analyzed and create a new column of the cell ID
        snakeS = snakeI %>% filter(Slice %in% cellNum$Slice) %>% mutate(Cell=NA)

        # Use the first slice to number the cells from 1 to total
        snakeS$Cell[1:unique(cellNum$Ncell)] = c(1:unique(cellNum$Ncell))
        
        # Define a new vector of the different slices  
        Slices = unique(snakeS$Slice)

        # An easy way to track a single cell in each image is by the distance
        # The closest (considering the x and y centroids) cell in the following slice should be the same cell
        # We create an empty vector to store the measured distances
        distvec = c()
        
        # Start the cell annotation
        for (e in 1:(length(Slices)-1)){
          
          # Make two data frames: one for the current slice and one for the next slice
          sl1 = snakeS %>% filter(Slice == Slices[e])
          sl2 = snakeS %>% filter(Slice == Slices[e+1])
          
          # Start the distance calculation, going from each cell in the first slice...
          for (f in 1:mlv(cellNum$Ncell)){
            
            # ...to compare it to each available cell in the next slice
            for (g in 1:mlv(cellNum$Ncell)){
              
              # Calculate the euclidean distance (shortest possible distance)
              # It follows the formula of Pythagoras' theorem (a^2 + b^2 = c^2)
              # 'a' is the difference in x coordinates
              # 'b' is the difference in y coordinates
              eucdist = sqrt((sl1$x[f]-sl2$x[g])^2 + (sl1$y[f]-sl2$y[g])^2)
              
              # Store the distance values in a new vector
              distvec[g] = eucdist
              
              # At the end of this loop, all the distances have been calculated 
              # for one cell in the current slice to all cells in the next slice
              # Match according to the index
              n = match(min(distvec),distvec)
            }
            # Anotate the indexes in the next slice
            sl2$Cell[f] = n
            } 
          # Anotate the indexes in the data of the slice
          snakeS$Cell[((mlv(cellNum$Ncell)*e)+1):((mlv(cellNum$Ncell)*e)+mlv(cellNum$Ncell))] = sl2$Cell
          
          # Remove the distance value and the data frames
          rm(eucdist,sl1,sl2)
          }
        # Remove the distance vector
        rm(distvec)
        
        # Select each single cell based on the maximum width
        cellw = aggregate(w ~ Cell, snakeS, function(x) max(x))
        
        # Create a line plot of the image data, with the slices as x axis and width as y axis
        #p = ggplot(snakeS,aes(x=as.integer(Slice),y=w,color=as.factor(Cell))) + theme_bw() +
        #  geom_point(size=3) + geom_line(size=2) +
        #  # Axis names
        #  xlab("Slice") + ylab('Cell width (µm)') + 
        #  # Plot title
        #  ggtitle(paste(Media[a],Conditions[b],Times[c],Images[d],sep = ", ")) +
        #  # Remove legends
        #  theme(legend.position = 'none') + 
        #  # Define range of the y axis
        #  ylim(0.3,1.2)
        
      #  # Store in the plot list and increase the index to prevent overwriting
      #  plotL[[h]] = p
      #  h = h+1
        
        # Filter on the data frame of the slice
        snakeS = snakeS %>% filter(w %in% cellw$w)
        
        # Store the processed data frame in the list
        dataList[[i]] = snakeS
        
        # Increase the index to not overwrite
        i = i+1
        
        # Remove data frames that will be overwritten
        rm(snakeS,cellw)
        }
      # Remove data frames that will be overwritten
      rm(cellNum,snakeT,snakeI)
    }
  }
  
  # When the code reaches the end:
  if (a == length(Media)){
    
    # Join all the data in a new data frame 
    # Calculate cell parameters based on length and width (Ojkic et al. 2019)
    snakes = bind_rows(dataList) %>% select(-c(Slice,Image)) %>% 
      mutate(SA = pi*L*w, # Surface area
             V = SA*w/4 - pi*(w^3)/12, # Volume
             SVR = SA/V, # Surface to volume ratio
             AR = L/w)  # Aspect ratio 
      
    
    # Remove all things that will not be used
    rm(Conditions,Images,Media,Slices,Times,snakeM,snakeC,snake_data,dataList,
     a,b,c,d,e,f,g,h,i,n,p)
  }
  }

###### Time-series Boxplots #####
# Construct boxplots to show the behavior of cell shape in each medium

# Create a plot list 
boxplots = list()

# In this loop, Length, Width, Volume and S/V ratio will be plotted as a function of time
# It will be split by media so it is easier to look at the data


for (a in 1:length(unique(snakes$Medium))){

Media = unique(snakes$Medium)  
  
colorval = c('dark blue', 'dark red')    

snakes2 = snakes %>% filter(Medium == Media[a])
  
p1 = ggplot(snakes2,aes(x=as.factor(Time),y=w,fill=Medium)) + theme_bw() + geom_boxplot() +
  xlab('') + ylab('Cell width (µm)') +
  stat_smooth(aes(group=Medium,color=Medium),size=1.5,se=F) +
  scale_fill_manual(values = colorval[a]) + ggtitle(Media[a]) +
  scale_color_manual(values = colorval[a]) + ylim(0.6,1) #+  
  #theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
  #      axis.text.x=element_text(size=55,angle = 90),legend.text = element_blank(),
  #      plot.title = element_text(size=55))

p2 = ggplot(snakes2,aes(x=as.factor(Time),y=L,fill=Medium)) + theme_bw() + geom_boxplot() +
  xlab('') + ylab('Cell length (µm)') +
  stat_smooth(aes(group=Medium,color=Medium),size=1.5,se=F) +
  scale_fill_manual(values = colorval[a]) + ylim(0,10) +
  scale_color_manual(values = colorval[a]) + ggtitle(Media[a]) #+  
#theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
#      axis.text.x=element_text(size=55,angle = 90),legend.text = element_blank(),
#      plot.title = element_text(size=55))

p3 = ggplot(snakes2,aes(x=as.factor(Time),y=V,fill=Medium)) + theme_bw() + geom_boxplot() +
  xlab('') + ylab('Volume (µm^3)') +
  stat_smooth(aes(group=Medium,color=Medium),size=1.5,se=F) +
  scale_fill_manual(values = colorval[a]) + ylim(0,4) +
  scale_color_manual(values = colorval[a]) + ggtitle(Media[a]) #+  
#theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
#      axis.text.x=element_text(size=55,angle = 90),legend.text = element_blank(),
#      plot.title = element_text(size=55))

p4 =ggplot(snakes2,aes(x=as.factor(Time),y=SVR,fill=Medium)) + theme_bw() + geom_boxplot() +
  xlab('') + ylab('Surface to Volume ratio (1/µm)') +
  stat_smooth(aes(group=Medium,color=Medium),size=1.5,se=F) +
  scale_fill_manual(values = colorval[a]) + ylim(4,10) +
  scale_color_manual(values = colorval[a]) + ggtitle(Media[a]) #+  
#theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
#      axis.text.x=element_text(size=55,angle = 90),legend.text = element_blank(),
#      plot.title = element_text(size=55))

boxplots[[(a*4-3)]] = p1 
boxplots[[(a*4-2)]] = p2 
boxplots[[(a*4-1)]] = p3 
boxplots[[(a*4)]] = p4

if(a == length(unique(snakes$Medium))){
  
  rm(Media,colorval,p1,p2,p3,p4,snakes2)
}
}

ggarrange(plotlist = boxplots,ncol =length(boxplots)/a ,nrow = a,legend = 'none')

###### Line Plots: time vs cell shape #####

# Construct line plots to better compare between the media
plots = list()

sumtable = snakes %>% group_by(Time,Medium) %>% 
  summarise(MeanL=mean(L),SDL=sd(L),
            Meanw=mean(w),SDw=sd(w),
            MeanV=mean(V),SDV=sd(V),
            MeanSVR=mean(SVR),SDSVR=sd(SVR),
            MeanAR=mean(AR),SDAR=sd(AR))

startValsL= sumtable %>% filter(Time < 0 & Medium =='Lb')
startValsS = sumtable %>% filter(Time < 0 & Medium =='S750')

breakvec = c(seq(-1,8,1),24)

p1 = ggplot(snakes,aes(x=Time,y=w,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=5) + 
  geom_point(inherit.aes = FALSE,data=sumtable,
             aes(x=Time,y=Meanw,color=Medium),size=10) +
  xlab('Time (hours)') + ylab('Cell Width (µm)') +
  geom_errorbar(inherit.aes = FALSE,data=sumtable,
                aes(x=Time,y=Meanw,color=Medium,
                    ymax=Meanw+SDw,ymin=Meanw-SDw),position = position_dodge(0.05),size=5) + 
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +
  geom_hline(yintercept = startValsL$Meanw,size=1,color='dark red',alpha=0.5) +
  geom_hline(yintercept = startValsS$Meanw,size=1,color='dark blue',lty=2,alpha=0.5) +  
  theme(axis.text.y=element_text(size=105),axis.title = element_text(size=100),
        axis.text.x=element_text(size=105),legend.text = element_blank(),
        plot.title = element_text(size=105))


p2 = ggplot(snakes,aes(x=Time,y=L,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=5) + 
  geom_point(inherit.aes = FALSE,data=sumtable,
             aes(x=Time,y=MeanL,color=Medium),size=10) +
  xlab('Time (hours)') + ylab('Cell Length (µm)') +
  geom_errorbar(inherit.aes = FALSE,data=sumtable,
                aes(x=Time,y=MeanL,color=Medium,
                    ymax=MeanL+SDL,ymin=MeanL-SDL),position = position_dodge(0.05),size=5) + 
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +
  geom_hline(yintercept = startValsL$MeanL,size=1,color='dark red',alpha=0.5) +
  geom_hline(yintercept = startValsS$MeanL,size=1,color='dark blue',lty=2,alpha=0.5) +  
  theme(axis.text.y=element_text(size=105),axis.title = element_text(size=100),
        axis.text.x=element_text(size=105),legend.text = element_blank(),
        plot.title = element_text(size=105))

p3 = ggplot(snakes,aes(x=Time,y=V,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=5) + 
  geom_point(inherit.aes = FALSE,data=sumtable,
             aes(x=Time,y=MeanV,color=Medium),size=10) +
  xlab('Time (hours)') + ylab('Volume (µm^3)') +
  geom_errorbar(inherit.aes = FALSE,data=sumtable,
                aes(x=Time,y=MeanV,color=Medium,
                    ymax=MeanV+SDV,ymin=MeanV-SDV),position = position_dodge(0.05),size=5) + 
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +
  geom_hline(yintercept = startValsL$MeanV,size=1,color='dark red',alpha=0.5) +
  geom_hline(yintercept = startValsS$MeanV,size=1,color='dark blue',lty=2,alpha=0.5)+  
  theme(axis.text.y=element_text(size=105),axis.title = element_text(size=100),
        axis.text.x=element_text(size=105),legend.text = element_blank(),
        plot.title = element_text(size=105))

p4 = ggplot(snakes,aes(x=Time,y=SVR,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=5) + 
  geom_point(inherit.aes = FALSE,data=sumtable,
             aes(x=Time,y=MeanSVR,color=Medium),size=10) +
  xlab('Time (hours)') + ylab('Surface to Volume ratio (1/µm)') +
  geom_errorbar(inherit.aes = FALSE,data=sumtable,
                aes(x=Time,y=MeanSVR,color=Medium,
                    ymax=MeanSVR+SDSVR,ymin=MeanSVR-SDSVR),position = position_dodge(0.05),size=5) + 
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +
  geom_hline(yintercept = startValsL$MeanSVR,size=1,color='dark red',alpha=0.5) +
  geom_hline(yintercept = startValsS$MeanSVR,size=1,color='dark blue',lty=2,alpha=0.5)+  
  theme(axis.text.y=element_text(size=105),axis.title = element_text(size=100),
        axis.text.x=element_text(size=105),legend.text = element_blank(),
        plot.title = element_text(size=105)) 

plots[[1]]=p1
plots[[2]]=p2
plots[[3]]=p3
plots[[4]]=p4

ggarrange(plotlist = plots,ncol = 2,nrow=2,legend='none')

###### Growth curve data ######

growthcurve = read_xlsx('~/Desktop/Bsubtilis_growthcurve_cellshape.xlsx',sheet = 1) %>% 
  mutate(Time = ifelse(Time==28,24,ifelse(Time <=20, Time*15/60, (15*20+(Time-20)*30)/60)))%>% 
  group_by(Time,Medium) %>% summarise(OD=mean(OD600),SD=sd(OD600)) %>% filter(!is.na(OD))

ggplot(growthcurve,aes(x=Time,y=log(OD),color=Medium)) + theme_bw() + xlim(0,8) +
  geom_point(size=10) + geom_line(size=5) + xlab('Time (hours)') +
  scale_color_manual(values = c('dark red','dark blue')) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55),legend.position = 'none')

kvec = c()
mlist = list()

for (j in 1:length(unique(growthcurve$Medium))){
  
  Media = unique(growthcurve$Medium)
  medium = growthcurve %>% filter(Medium == Media[j])
  
  a = length(unique(medium$Time))-1
  
  for (i in 1:a){
    
    t = unique(medium$Time)
    OD = medium$OD
    
    kvec[i] = ((OD[i+1] - OD[i])/(t[i+1] - t[i]))/OD[i]
    
    
  }
  medium = medium %>% filter(Time %in% t[1:a]) 
  
  medium$k = kvec
  
  mlist[[j]] = medium
}

mlistk = list()

for (j in 1:length(mlist)){
  
  medium = sumtable %>% filter(Medium == Media[j] & Time %in% mlist[[j]]$Time) %>% 
    left_join(mlist[[j]]) 
  
  mlistk[[j]] = medium
}

sumtable_k = bind_rows(mlistk)

###### Line plots :including the growth rate ####
plots_k = list()

snakes_k = snakes %>% filter(Time %in% sumtable_k$Time) %>% left_join(sumtable_k) %>% 
  select(Medium,Time,k,L,w,V,SVR)

startValsL= sumtable_k %>% filter(Time == 0 & Medium =='Lb')
startValsS = sumtable_k %>% filter(Time == 0 & Medium =='S750')

breakvec = seq(0,8,1)

p1_k = ggplot(snakes_k,aes(x=Time,y=w,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=5) + 
  geom_point(inherit.aes = FALSE,data=sumtable_k,
             aes(x=Time,y=Meanw,color=Medium,size=10,alpha=k)) +
  xlab('Time (hours)') + ylab('Cell Width (µm)') +
  geom_errorbar(inherit.aes = FALSE,data=sumtable_k,
                aes(x=Time,y=Meanw,color=Medium,
                    ymax=Meanw+SDw,ymin=Meanw-SDw,alpha=k),position = position_dodge(0.05),size=5) + 
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) + theme(legend.key.size = unit(2, 'cm'))
geom_hline(yintercept = startValsL$Meanw,size=1,color='dark red',alpha=0.5) +
  geom_hline(yintercept = startValsS$Meanw,size=1,color='dark blue',lty=2,alpha=0.5) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55))


p2_k = ggplot(snakes_k,aes(x=Time,y=L,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=5) + 
  geom_point(inherit.aes = FALSE,data=sumtable_k,
             aes(x=Time,y=MeanL,color=Medium,size=10,alpha=k)) +
  xlab('Time (hours)') + ylab('Cell Length (µm)') +
  geom_errorbar(inherit.aes = FALSE,data=sumtable_k,
                aes(x=Time,y=MeanL,color=Medium,
                    ymax=MeanL+SDL,ymin=MeanL-SDL,alpha=k),position = position_dodge(0.05),size=5) + 
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +
  geom_hline(yintercept = startValsL$MeanL,size=1,color='dark red',alpha=0.5) +
  geom_hline(yintercept = startValsS$MeanL,size=1,color='dark blue',lty=2,alpha=0.5) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55))

p3_k = ggplot(snakes_k,aes(x=Time,y=V,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=5) + 
  geom_point(inherit.aes = FALSE,data=sumtable_k,
             aes(x=Time,y=MeanV,color=Medium,size=10,alpha=k)) +
  xlab('Time (hours)') + ylab('Volume (µm^3)') +
  geom_errorbar(inherit.aes = FALSE,data=sumtable_k,
                aes(x=Time,y=MeanV,color=Medium,
                    ymax=MeanV+SDV,ymin=MeanV-SDV,alpha=k),position = position_dodge(0.05),size=5) + 
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +
  geom_hline(yintercept = startValsL$MeanV,size=1,color='dark red',alpha=0.5) +
  geom_hline(yintercept = startValsS$MeanV,size=1,color='dark blue',lty=2,alpha=0.5) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55))

p4_k = ggplot(snakes_k,aes(x=Time,y=SVR,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=5) + 
  geom_point(inherit.aes = FALSE,data=sumtable_k,
             aes(x=Time,y=MeanSVR,color=Medium,size=10,alpha=k)) +
  xlab('Time (hours)') + ylab('Surface to Volume ratio (1/µm)') +
  geom_errorbar(inherit.aes = FALSE,data=sumtable_k,
                aes(x=Time,y=MeanSVR,color=Medium,
                    ymax=MeanSVR+SDSVR,ymin=MeanSVR-SDSVR,alpha=k),position = position_dodge(0.05),size=5) + 
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +
  geom_hline(yintercept = startValsL$MeanSVR,size=1,color='dark red',alpha=0.5) +
  geom_hline(yintercept = startValsS$MeanSVR,size=1,color='dark blue',lty=2,alpha=0.5) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55))

plots_k[[1]]=p1_k
plots_k[[2]]=p2_k
plots_k[[3]]=p3_k
plots_k[[4]]=p4_k

ggarrange(plotlist = plots_k,ncol = 2,nrow=2,legend='none')#,common.legend = TRUE)

###### Scatter plots: growth rate vs cell shape #####

plots_ks = list()

p1_ks = ggplot(snakes_k,aes(x=k,y=w,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=3) + geom_jitter(aes(alpha=Time),size=5) +
  xlab('Growth rate (1/hours)') + ylab('Cell Width (µm)') +
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55))

p2_ks = ggplot(snakes_k,aes(x=k,y=L,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=3) + geom_jitter(aes(alpha=Time),size=5) +
  xlab('Growth rate (1/hours)') + ylab('Cell Length (µm)') +
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55))

p3_ks = ggplot(snakes_k,aes(x=k,y=V,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=3) + geom_jitter(aes(alpha=Time),size=5) +
  xlab('Growth rate (1/hours)') + ylab('Volume (µm^3)') +
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55))

p4_ks = ggplot(snakes_k,aes(x=k,y=SVR,color=Medium)) + theme_bw() +
  stat_smooth(aes(lty=Medium),se=F,size=3) + geom_jitter(aes(alpha=Time),size=5) +
  xlab('Growth rate (1/hours)') + ylab('Surface to Volume ratio (1/µm)') +
  scale_color_manual(values=c('dark red','dark blue')) + 
  scale_x_continuous(breaks=breakvec)  +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55))

plots_ks[[1]]=p1_ks
plots_ks[[2]]=p2_ks
plots_ks[[3]]=p3_ks
plots_ks[[4]]=p4_ks

ggarrange(plotlist = plots_ks,ncol = 2,nrow=2,legend='none')#,common.legend = TRUE)


###### Binning data by cell shape ####
snakes = quantilizer(data=snakes,column=c(snakes$L),QuantRemover = TRUE) %>% 
  rename(Length.Quantile = Quantile)

snakes = quantilizer(data=snakes,column=snakes$w) %>% 
  rename(Width.Quantile = Quantile) %>% 
  mutate(Length.Quantile = fct_reorder(Length.Quantile,L,.fun = 'min'))

p = ggplot(snakes,(aes(x=L,y=w))) + theme_bw() + geom_jitter(aes(shape=Medium,color=Width.Quantile),size=2.5) +
  stat_smooth(aes(lty=Medium),se=F,color='black',method='lm',size=2) + 
  xlab('Length') + ylab('Width') + ylim(min(snakes$w),max(snakes$w)) +  
  theme(axis.text.y=element_text(size=55),axis.title = element_text(size=50),
        axis.text.x=element_text(size=55)) + 
  scale_color_brewer(type='seq',palette = "YlOrRd")

a = quantile(snakes$L,probs=seq(0,1,0.2),type = 7,names = FALSE)

ggplot(snakes,aes(x=L)) + theme_bw() + geom_histogram(color='purple',fill='white',size=5,binwidth = 0.3) +
  xlab('Cell length (µm)') + ylab('Counts') +
  theme(axis.text.y=element_text(size=105),axis.title = element_text(size=100),
        axis.text.x=element_text(size=105))

ggplot(snakes,aes(x=L)) + theme_bw() + geom_histogram(aes(fill=Length.Quantile),color='black',size=5,binwidth = 0.3) +
  xlab('Cell length (µm)') + ylab('Counts') +
  scale_fill_brewer(type = 'seq',palette = 'Purples') +
  theme(axis.text.y=element_text(size=105), axis.title = element_text(size=100),
        axis.text.x=element_text(size=105),legend.position = 'none')

facet(p, scales = "free_x",facet.by = c("Medium","Length.Quantile"),
      short.panel.labs = FALSE,nrow=2)

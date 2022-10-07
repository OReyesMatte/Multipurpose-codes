###### Package installation and Library loading ######
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse","readxl","car")

# Create a "list" object to store averaged data plates 
plates = list()

# Average replicates within a plate. In the first 'for' function change the last number to the number of sheets in the excel file 

for (j in 1:3){

# Read data file. The sheets will be read in order
# Change the path accordingly to your own files
growth <- read_xlsx("~/Desktop/ECB/Growth assays/Growth_selected_promoters.xlsx",sheet=j)

#Create a time vector of minutes. The vector created will have the length of the total kinetic 
t = seq(0,(nrow(growth)-1)*2,2)

# Change the Time column on the data sheet with the new time vector
growth$Time = t

# In R is easier to work with less columns even if you have a lot of rows
# The 'gather' function takes a "wide" data frame and transforms the existing columns into a new category

growth2 = gather(growth, # The data frame to be transformed 
                 key = Promoter, # The name of the new column that will have the previous column names 
                 value = OD, # The name of the column that will have the measures 
                 -Time) # Which column to leave intact (optional)

for (i in 1:length(growth2$Promoter)){ #Select only the promoter name (i.e. Pcv1)
  # This depends on the length of your Promoter/Genotype/Category and it's optional
  
  growth2$Promoter[i] = substr(growth2$Promoter[i],1,5)
  # 1 and 5 define the limits of the string to be considered for the name of the promoter
}

# In my case, I changed the original wt name. Again, this is optional
# Make sure the name you provide is as exactly as you have it so the function will work
growth2$Promoter = car::recode(growth2$Promoter,"'wt...'='WT'")

# Average the values for each category
# %>% is a pipeline argument, which allows you to chain functions into one line as long as it's all related to the same data object

growth_summ = growth2 %>% group_by(Time,Promoter) %>% # Columns that will be used to cuantify values
  # Grouping by these two means that you'll get a mean value of each category for every time point
  summarise(Mean = mean(OD),SD=sd(OD))# Calculate means and standard deviations

# Store the summary table in the list object
plates[[j]] = growth_summ 

}

# Now, we'll average between plates to calculate the total average curve
growth_summ = bind_rows(plates) %>% 
  # 'bind_rows' joins data frames that have the same columns and stack the rows
  ungroup %>% # Ungroup just in case, then repeat the averaging process
  group_by(Time,Promoter) %>% 
  summarise(StrainMean = mean(Mean),StrainSD=sd(Mean))

# As we want the OD values minus the blank OD, we separate these values to then substract them
blanks_mean = subset(growth_summ, Promoter == "Blank") 
growth_data = subset(growth_summ, Promoter != "Blank")

# Now we plot! 'ggplo' is the general plotting start using the ggplot2 library (included in tidyverse)
# For these plots new arguments are chained using a "+"
ggplot(growth_data, # Where the data to be plotted is  
       aes(#'aes' is used to define x and y values, as well as coloring and filling categories
         x=Time/60, # Define the values for x axis (in this case, time in hours)
         y=StrainMean-blanks_mean$StrainMean, # Define the values for y axis
         # In this case, it's going to be the data of the non-blank curves with the substraction of the man blank OD
         color=Promoter # Based on which category the data will be colored
         )) + geom_point() + # Plot data as points
  geom_line() + # Plot data as a line
  theme_minimal() + # Have a white background with a light grey grid
  scale_color_brewer(type="qual",palette="Paired") + # Define a color palette 
  xlim(c(0,660/60)) + # Define the x-axis limits. Make sure the numbers are a subset of the ones used in 'aes'
  ylim(c(0,0.5)) + # Define the y-axis limits. Make sure the numbers are a subset of the ones used in 'aes'
  xlab("Time (h)") + # Re-label the x-axis
  ylab("OD600") + # Re-label the y-axis
  ggtitle("Growth curve assay") # Add a plot title

# While this may give us information in the behaviour of the populations, we need 
# to calculate the growth rate (in exponential) to determine which grows better/faster

# Subset the data until the point where bacteria are finishing exponential growth
# Look at the previous plot and define the time accordingly
times =  subset(growth_data, Time <= 240)

# Define a vector with your categories
promoters = unique(growth_data$Promoter)

# Define a vector of exponential time points
tdiv = unique(times$Time[-1])

# Create a list to store the calculated growth rates
Rates = list()

# In this loop  we'll create new data frames for storing calculated growth rates
for (j in 1:length(promoters)){
  # Create the new data frame, defining the different columns and the lenght
  RateData = data.frame(Time=tdiv[-1], Promoter=1:(length(tdiv)-1),k=1:(length(tdiv)-1))  
  
  # Subset the category of interest
  promData = subset(times, Promoter == promoters[j])
  
  # Create an empty vector for storing growth rate values
  growthRate = c()
  
  # Calculate the growth rate as the exact derivative of the population growth
  for (i in 1:length(promData$Time)-1){
    
    growthRate[i] = ((promData$StrainMean[i+1] - 
                        promData$StrainMean[i])/(promData$Time[i+1]-promData$Time[i]))/(promData$StrainMean[i])
    
  }
  # Add the calculated growth rates to the data frames and store in the list
  RateData = RateData %>% mutate(Promoter = promoters[j], k=growthRate)
  Rates[[j]] = RateData
}

# Join growth rate calculations and filter out negative, 0 and NaN values
RateData = bind_rows(Rates) %>% filter(k>0 & !is.nan(k))

# Plot boxplots of the growth rate values
# The structure is similar as in the other plot, just changing "color" for "fill both in aes and in the palette function
ggplot(RateData, aes(x=Promoter,y=k*60,fill=Promoter)) + geom_boxplot() +  
  theme_minimal() + scale_fill_brewer(type="qual",palette="Paired")  + 
  xlab("Promoter") + ylab("Growth rate (1/h)") + ggtitle("Calculation of growth rates")

# Normality test
normtest = shapiro.test(RateData$k*60)

# Homoscedasticty test
leveneTest(k*60~Promoter, RateData)

# ANOVA or kurskal-wallis test depending on results of normality test
if (normtest$p.value >= 0.05){
  summary(aov(k*60~Promoter, RateData))
} else {
kruskal.test(k*60~Promoter, RateData)
}

# Store mean growth rates to compare them
meanRate = RateData %>% group_by(Promoter) %>% summarise(Meank = mean(k*60))




# Multipurpose codes
Different codes that I use for crunching, wrangling and analysing (and sometimes generating) different types of biological data.

The different folders guide to the specific languages the codes are written in. 

## R codes

Written mostly using R studio, traditionally I use mostly [Tidyverse](https://www.tidyverse.org/) for the data wrangling and also loading

### Cell_selector_z-stack.R

The objective of this code is to "compress" image analysis data into a more manageable table that makes it easier for plotting and representation. It 'tracks' cells based on their position across different focal points of a Z-stack image.
The input data is a table of cell measurements obtained with the [MicrobeJ](https://www.microbej.com/) plugin for FIJI. For this to be used, it is necessary to have at least measurements of cell length and width, as well as the x and y coordinates of the cell centroids.


## MATLAB codes

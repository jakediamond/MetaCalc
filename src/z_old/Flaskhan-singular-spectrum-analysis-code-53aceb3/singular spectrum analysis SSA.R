# Author: Hares Khan 2019, (hkhan.ch@gmail.com)
# This is the R code used in the study "Khan H., Laas A., Marcé R., Obrador B., "Major effects of alkalinity on the relationship between metabolism and dissolved inorganic carbon dynamics in lakes"
# The code was used in the mentioned study to isolate the diel metabolic signal (seasonal pattern of 24 hours) from High Frequency Measurements time series of dissolved oxygen and dissolved inorganic carbon.
# It can be useful to extract seasonal patterns, trends or noise signals from a time series using Singular Spectrum Analysis (SSA).

# This code can be used freely and redistributed.
# This code comes with no warranty

### IMPORTANT! FOR SSA, your datafile must not contain any missing data. Also, your time-series must always have the same time interval!

# Install and download required packages
install.packages("Rssa")
library(Rssa)
library(xts)
library(lattice)  #For multipanel graphs

# Set your working directory in the folder where you keep your datafile
# Load your datafile
yourdata <- read.csv("datafile.csv", header = TRUE, sep = ",", na.strings = "NA")
yourdata <- df


yourdata <- filter(yourdata, year > 1994) %>%
  select(datetime, SpC) %>%
  imputeTS::na_kalman()

# Put your date-time column into POSIXct format
yourdata$pr_dt <- as.POSIXct(yourdata$datetime, format = "%m/%d/%Y %H:%M")
# Create a zooreg object for your variable
yourvariable<-yourdata$SpC
yourvariable<- zooreg(yourvariable, order.by=as.POSIXct(yourdata$pr_dt, "%m/%d/%Y %H:%M"))
# summary(yourvariable)
# Decomposition stage: this decomposes the time-series into its fundamental Eigenvectors of oscillations, trends, and noise
E <- ssa(yourvariable)
summary(E)

# Diagnosing the decomposed time-series to find out how to reconstruct it.
plot(E) # Eigenvalues
plot(E, type = "vectors", main="yourdata Eigenvectors for yourvariable") ### This shows the most important Eigenvectors. Here you can identify Eigenvectors responsible for trends, seasonality, and noise. 
plot(E, type = "paired", main="yourdata pairs of Eigenvectors for yourvariable") ### Pairs of eigenvectors helps to identify combinations of Eigenvectors responsible for seasonality such as diel signals. If the paired Eigenvectors create a circular pattern, then they are responsible for a seasonal pattern and should be combined in the reconstruction phase. 
plot(wcor(E)) # w-correlation matrix plot. This shows correlations between the Eigenvectors. This also helps to identify which Eigenvectors should be combined together in the reconstruction phase. The first Eigenvectors that are strongly correlated together should be combined in the reconstruction phase. 

# Reconstruction stage
# Here you reconstruct the time series by combining Eigenvectors according to your previous diagnosis using the plots.
recon_yourvariable <- reconstruct(E, groups = list(c(2,3), c(4,5))) # The numbers need to be changed to group the Eigenvectors that are strongly correlated in the previous diagnosis stage. In this example, the reconstruction will group Eigenvectors 1 and 2, and Eigenvectors 3 and 4. 

# Calculate the residuals
res_yourvariable <- residuals(recon_yourvariable)

# Here you visualize your reconstructed series
plot(recon_yourvariable, ylab=recon_yourvariable$name, main="yourdata reconstructed series for yourvariable" )

#plot reconstructed series (trends, seasonal patterns)
plot(yourvariable)
pattern1<-recon_yourvariable$F1 
pattern2<-recon_yourvariable$F2 

# You can also visualize all 50 Eigenvectors at once. This shows you how SSA decomposes your signal.To do so, run the two following lines:
recon_yourvariable <- reconstruct(E)
plot(recon_yourvariable, ylab=recon_yourvariable$name)


################### Multivariate Singular Spectrum Analysis (MSSA)
# This is a Multivariat SSA that allows you to isolate common trends or oscillations between to variables.
your2variables <- yourdata[, c("SpC", "NEP_mmolO2m3")]
H <- ssa(your2variables, kind = "mssa")

# same diagnosis procedure as for the simple SSA
plot(H) # Eigenvalues
plot(H, type = "vectors", main="Eigenvectors") # Eigenvectors
plot(H, type = "paired", main="paired Eigenvectors") # Pairs of eigenvectors
plot(wcor(E)) # w-correlation matrix plot

# Reconstruction phase, similar as for SSA
recon_your2variables <- reconstruct(H, groups = list(c(2,3), c(4,5))) #change numbers as needed according to the previous diagnosis stage.
plot(recon_your2variables, ylab=recon_your2variables$name, main="yourdata MSSA")
# To visualize correlations between the reconstructed series of your two variables

plot(recon_your2variables$F1[1:24])

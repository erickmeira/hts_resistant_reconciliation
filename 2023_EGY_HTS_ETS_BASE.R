#  Meira et al. (2023) - A novel reconciliation approach for hierarchical electricity consumption forecasting based on resistant regression----
# Code: Reconciling ETS base forecasts of electricity consumption
# Data: Electricity Consumption across the divisions of the Brazilian National Interlinked System

# Load the required libraries----
library(hts)
library(readr)
library(MASS)
library(dplyr)
library(plyr)
library(tidyr)
library(quantreg)
library(nnls)

# Evaluation setup----
h = 4 # Number of forecasting steps (forecast lead time)

# Train set 
train.beg.m = 1; train.beg.y = 2004 # Month and year of the first train set observation
train.end.m = 7; train.end.y = 2021 # Month and year of the last train set observation

# Test set
test.beg.m = 8; test.beg.y = 2021 # Month and year of the first test set observation
test.end.m = 11; test.end.y = 2021 # Month and year of the last test set observation


# Reading data and adjusting its format to time series to serve as input to the hts function
data.0 <-  read.csv("HTS_Brazilian_Power_System.csv")
data.1 <- ts(data.0,start=c(2004,1), end=c(test.end.y,test.end.m), frequency=12)

# Generating the hts object
data.2 <- hts(data.1,characters=c(3,2))

# Extracting the matrix which defines the hierarchy, i.e., the linear dependencies between the time series
X <- smatrix(data.2)

# Defining test and training sets in the hts object
train <- window(data.2, start=c(train.beg.y,train.beg.m), end = c(train.end.y ,train.end.m))
test <- window(data.2, start=c(test.beg.y,test.beg.m), end = c(test.end.y,test.end.m))

# Extracting the sets of time series that are present at each level of the hierarchy----
# Train set
train.0 <- aggts(train, levels=0) # Top level (total national electricity consumption)
train.1 <- aggts(train, levels=1) # Intermediate level (Consumption across classes - commercial, residential...)
train.2 <- aggts(train, levels=2) # Bottom level (Consumption across the distinct geographic regions within classes)

# Test set
test.0 <- aggts(test,levels=0)
test.1 <- aggts(test,levels=1)
test.2 <- aggts(test,levels=2)

# 1. Base forecasts - Generating the independent (base) forecasts for each time series in the hierarchy using ETS as forecasting method---- 

# h-steps ahead base forecasts for the time series at the top level
# (total national electricity consumption)
allf.0 <- forecast(ets(train.0[,1]), h = h)$mean 
# Note: if the column 'train.0[,1]' is not specified, the forecast function does not understand the input as a time series

# Accuracy metrics for the ETS base forecasts for the time series at the top level
acc.0 <- accuracy(allf.0,test.0)

# h-steps ahead base forecasts for the four time series at the intermediate level
# (commercial, residential, industrial and others)
allf.1 <- matrix(NA, nrow = h, ncol = ncol(train.1))
acc.1<-matrix(NA, nrow = ncol(train.1), ncol = length(acc.0) )
for(j in 1:ncol(train.1)){
  allf.1[,j] <- forecast(ets(train.1[,j]), h = h)$mean
  acc.1[j,]<-accuracy(allf.1[,j],test.1[,j])
}

# h-steps ahead base forecasts for the twenty time series at the bottom level
# (Commercial North, Commercial Northeast,..., Commercial South, Residential North, ..., Others South)
allf.2 <- matrix(NA, nrow = h, ncol = ncol(train.2))
acc.2<-matrix(NA, nrow = ncol(train.2), ncol = length(acc.0))
for(j in 1:ncol(train.2)){
  allf.2[,j] <- forecast(ets(train.2[,j]), h = h)$mean
  acc.2[j,]<- accuracy(allf.2[,j],test.2[,j])
}

# Column binding all h-steps ahead base forecasts into a single table
allf <- cbind(allf.0, allf.1, allf.2)

# Compiling in a single table the accuracy measures (ME, RMSE, MAE, MAPE and MPE) for the ETS base forecasts
acc.base <- rbind( acc.0, acc.1, acc.2)[,c(1,2,3,5,4)]

# 2. Reconciled forecasts - Generating reconciled forecasts using benchmarks and the proposed LAD resistant approach----

# 2.1 Benchmark methods - Reconciling the h-steps ahead ETS base forecasts using BU, OLS, MinT(S), TDGSF and TDGSA----
fcastsbu     <- forecast(train, h = h, method = "bu", fmethod="ets")  # Bottom-up
fcastsols    <- forecast(train, h = h, weights="ols", fmethod="ets")  # OLS
fcastsmint   <- forecast(train, h = h, weights="mint", covariance = "shr", fmethod="ets") # Mint-Shrink (MinT-S)
fcaststdgsf  <- forecast(train, h = h, method = "tdgsf", fmethod="ets") # TDGSF
fcaststdgsa  <- forecast(train, h = h, method = "tdgsa", fmethod="ets") # TDGSA

# Note that the hts function generates reconciled forecasts for the benchmark methods using the forecast function once only
# For the other methods (the robust - HUBER 1 and HUBER 2 - methods and the proposed resistant reconciliation approach,
# we must first generate the base forecasts, store the results and then use specific functions to reconcile such forecasts

# Computing the accuracy metrics for the forecasts reconciled using benchmark methods
acc.bu     <- t(accuracy(fcastsbu, test))[,1:5]   
acc.ols    <- t(accuracy(fcastsols, test))[,1:5]
acc.mint   <- t(accuracy(fcastsmint, test))[,1:5] 
acc.tdgsf  <- t(accuracy(fcaststdgsf, test))[,1:5] 
acc.tdgsa  <- t(accuracy(fcaststdgsa, test))[,1:5]

# 2.2 Regression methods - Reconciling the h-steps ahead ETS base forecasts using the HUBER 1, HUBER 2 and LAD approaches----
# Defining the forecasting matrices for the Regression (HUBER 1, HUBER 2 and LAD) approaches 
fcasthuber1 <- fcasthuber2 <- fcastlad <- matrix(NA, nrow = h, ncol = ncol(allf))

for(j in 1:nrow(allf)){
  suppressWarnings(fcasthuber1[j,]<-t(X%*%rlm(X, allf[j,], psi = psi.huber,maxit=1000,scale.est="MAD")$coefficients))
  suppressWarnings(fcasthuber2[j,]<-t(X%*%rlm(X, allf[j,], psi = psi.huber,maxit=1000,scale.est="proposal 2")$coefficients))
  suppressWarnings(fcastlad[j,]<-t(X%*%rq(allf[j,]~X-1, tau = 0.5)$coefficients))
}

# Accuracy metrics for the forecasts reconciled using regression methods
# Defining the accuracy matrices
acc.huber1 <- acc.huber2 <- acc.lad <- matrix(NA, nrow =nrow(acc.base), ncol = ncol(acc.base))

# Computing the accuracy metrics for the regression-based reconciled forecasts
for(j in 1:ncol(allf)){
  acc.huber1[j,] <- accuracy(ts(fcasthuber1[,j],
                                start=c(test.beg.y,test.beg.m),
                                end = c(test.end.y,test.end.m),frequency=12),
                             ts(allts(test)[,j],
                                start=c(test.beg.y,test.beg.m),
                                end = c(test.end.y,test.end.m),frequency=12))[,c(1,2,3,5,4)]
  
  acc.huber2[j,] <- accuracy(ts(fcasthuber2[,j],
                                start=c(test.beg.y,test.beg.m),
                                end = c(test.end.y,test.end.m),frequency=12),
                             ts(allts(test)[,j],
                                start=c(test.beg.y,test.beg.m),
                                end = c(test.end.y,test.end.m),frequency=12))[,c(1,2,3,5,4)]
  
  acc.lad[j,] <- accuracy(ts(fcastlad[,j],
                             start=c(test.beg.y,test.beg.m),
                             end = c(test.end.y,test.end.m),frequency=12),
                          ts(allts(test)[,j],
                             start=c(test.beg.y,test.beg.m),
                             end = c(test.end.y,test.end.m),frequency=12))[,c(1,2,3,5,4)]
  }

# 3. Evaluation of the results-----
# Obtaining the relative accuracy metrics according to hierarchical aggregation

# Creating a function to compute the geometric mean
geomean <- function(x){exp(mean(log(x)))}

# Applying the geometric mean function for each method
grand.mean.results <- rbind(apply(acc.huber1[,2:4] /acc.base[,2:4],2,geomean), 
                   apply(acc.huber2[,2:4] /acc.base[,2:4],2,geomean),
                   apply(acc.lad[,2:4]/acc.base[,2:4],2,geomean),
                   apply(acc.bu[,2:4]/acc.base[,2:4],2,geomean),
                   apply(acc.ols[,2:4]/acc.base[,2:4],2,geomean),
                   apply(acc.mint[,2:4]/acc.base[,2:4],2,geomean),
                   apply(acc.tdgsf[,2:4]/acc.base[,2:4],2,geomean),
                   apply(acc.tdgsa[,2:4]/acc.base[,2:4],2,geomean))
  
rownames(grand.mean.results) <- c("huber1", "huber2", "lad",
                                  "bu", "ols", "mint", "tdgsf", "tdgsa")

# EGY Paper Table 2 Results-----
paper.table2 <- round(grand.mean.results[, c(2,1,3)], 3)

# View results
View(paper.table2)
# Note that the results will be those for the right panel of Table 2 if h was defined as 4 steps-ahead
# If h was defined as 3 (h = 3), the results will be those of the left panel of Table 2

# EGY Paper Table 4 Results-----

# Top level average relative metrics
one<-1
top.results <- rbind(acc.huber1[one,2:4] /acc.base[one,2:4], 
                             acc.huber2[one,2:4] /acc.base[one,2:4],
                             acc.lad[one,2:4]/acc.base[one,2:4],
                             acc.bu[one,2:4]/acc.base[one,2:4],
                             acc.ols[one,2:4]/acc.base[one,2:4],
                             acc.mint[one,2:4]/acc.base[one,2:4],
                             acc.tdgsf[one,2:4]/acc.base[one,2:4],
                             acc.tdgsa[one,2:4]/acc.base[one,2:4])

rownames(top.results) <- c("huber1", "huber2", "lad",
                                   "bu", "ols", "mint", "tdgsf", "tdgsa")

# Intermediate level average relative metrics
one<-1
inter<-seq(1, ncol(allf.1),1)+one

intermediate.results <- rbind(apply(acc.huber1[inter,2:4] /acc.base[inter,2:4],2,geomean), 
                              apply(acc.huber2[inter,2:4] /acc.base[inter,2:4],2,geomean),
                              apply(acc.lad[inter,2:4]/acc.base[inter,2:4],2,geomean),
                              apply(acc.bu[inter,2:4]/acc.base[inter,2:4],2,geomean),
                              apply(acc.ols[inter,2:4]/acc.base[inter,2:4],2,geomean),
                              apply(acc.mint[inter,2:4]/acc.base[inter,2:4],2,geomean),
                              apply(acc.tdgsf[inter,2:4]/acc.base[inter,2:4],2,geomean),
                              apply(acc.tdgsa[inter,2:4]/acc.base[inter,2:4],2,geomean))

rownames(intermediate.results) <- c("huber1", "huber2", "lad",
                                    "bu", "ols", "mint", "tdgsf", "tdgsa")    

# Bottom level average relative metrics

one<-1
inter<-seq(1, ncol(allf.2),1)+one+ncol(allf.1)

bottom.results <- rbind(apply(acc.huber1[inter,2:4] /acc.base[inter,2:4],2,geomean), 
                        apply(acc.huber2[inter,2:4] /acc.base[inter,2:4],2,geomean),
                        apply(acc.lad[inter,2:4]/acc.base[inter,2:4],2,geomean),
                        apply(acc.bu[inter,2:4]/acc.base[inter,2:4],2,geomean),
                        apply(acc.ols[inter,2:4]/acc.base[inter,2:4],2,geomean),  
                        apply(acc.mint[inter,2:4]/acc.base[inter,2:4],2,geomean),
                        apply(acc.tdgsf[inter,2:4]/acc.base[inter,2:4],2,geomean),
                        apply(acc.tdgsa[inter,2:4]/acc.base[inter,2:4],2,geomean))

rownames(bottom.results) <- c("huber1", "huber2", "lad",
                              "bu", "ols", "mint", "tdgsf", "tdgsa")

paper.table4.part1 <- round(top.results[, c(2,1,3)], 3)
paper.table4.part2 <- round(intermediate.results[, c(2,1,3)], 3)
paper.table4.part3 <- round(bottom.results[, c(2,1,3)], 3)

View(paper.table4.part1)
View(paper.table4.part2)
View(paper.table4.part3)
# Note that the results will be those for the right panel of Table 4 if h was defined as 4 steps-ahead
# If h was defined as 3 (h = 3), the results will be those of the left panel of Table 4
library(lavaan)
library(tidyverse)

model <- '
#residual covariance
x ~~ 0.64*y

#factor
x =~ NA*x1 + 0.8*x1 + 0.8*x2
y =~ NA*y1 + 0.8*y1 + 0.8*y2

#residual variance
x1 ~~ 0.36*x1
x2 ~~ 0.36*x2
y1 ~~ 0.36*y1
y2 ~~ 0.36*y2
'
# simulate data
data20 <- simulateData(model = model, model.type = "sem",
            sample.nobs = 20, seed = 1000); data20

data50 <- simulateData(model = model, model.type = "sem",
                       sample.nobs = 50, seed = 1000)

data100 <- simulateData(model = model, model.type = "sem",
                       sample.nobs = 100, seed = 1000)


c_names = c("x1","x2","y1","y2")
prc_missing = 0.15



z1 <- data20 %>%
  add_column(id = 1:20) %>%
  gather(var, value, -id) %>%    # reshape data
  mutate(r = runif(nrow(.)),   # simulate a random number from 0 to 1 for each row
         value = ifelse(var %in% c_names & r <= prc_missing, NA, value)) %>%  # if it's one of the variables you specified and the random number is less than your threshold update to NA
  select(-r) %>%                 # remove random number
  spread(var, value)             # reshape back to original format

##### multiverse ####
library(ggplot2)
library(Rmisc)
library(caret)
library(lm.beta)
library(reghelper)
library(glmnet)

annlist <- c(1) #number of iteration
#deplist <- c("A")
all.data.multiverses <- list()
all.cor <- list()
avg.cor <- list()

for (i in 1:length(annlist)){
  #rm(list = setdiff(ls(), c("annlist", "all.cor", "all.data.multiverses", "i")))
  #ann[a] <- annlist[a]

    model <- '
    #residual covariance
    x ~~ 0.64*y
    
    #factor
    x =~ NA*x1 + 0.8*x1 + 0.8*x2
    y =~ NA*y1 + 0.8*y1 + 0.8*y2
    
    #residual variance
    x1 ~~ 0.36*x1
    x2 ~~ 0.36*x2
    y1 ~~ 0.36*y1
    y2 ~~ 0.36*y2
    '

  no.it <- 500
  no.x <- 2
  no.y <- 2
  no.z <- 2
  
  data.multiverse <- array(list(), dim = c(no.it, no.x, no.y, no.z))
  cor.multiverse <- array(0, dim = c(no.it, no.x, no.y, no.z))
  
  for (a in 1:no.it){
  data <- simulateData(model = model, model.type = "sem",
                         sample.nobs = 100)
  data.proc <- data
  
  data.proc.init <- data.proc
  
  for (j in 1:no.x) {
    for (k in 1:no.y) {
      for (l in 1:no.z) {

        data.proc <- data.proc.init
        
        if (j == 1) {data.proc$x <- data$x1} else {data.proc$x <- data$x2}
        if (k == 1) {data.proc$y <- data$y1} else {data.proc$y <- data$y2}
        if (l == 1) {data.proc <- data.proc} 
        else {data.proc <- data.proc[(data.proc$x < 1.5),]} # listwise deletion
        
        data.multiverse[[a, j, k, l]] = data.proc
        }
      }
    }
  }
  ### analysis ###
  for (a in 1:no.it) {
    for (j in 1:no.x) {
      for (k in 1:no.y) {
         for (l in 1:no.z) {

        # data.proc$x <- factor(data.proc$x)
        # data.proc$y <- factor(data.proc$y)
        
        #an = coef(summary(lm(y ~ x, data.multiverse[[j, k, l]])))
          dt <- data.multiverse[[a, j, k, l]] %>%
            select(x, y)
          an = cor(dt$y, dt$x)
          cor.multiverse[a, j, k, l] <- an
        
        #se.multiverse[j, k, l] <- an[2,2]
          avg.cor[[a]] = mean(cor.multiverse[a,,,])
        }
      }     
    } 
   
  }
  
  all.data.multiverses[[i]] <- data.multiverse
  all.cor[[i]] <- cor.multiverse
}


# plot corr
barplot(table(round(unlist(avg.cor), digits = 2)), xlab = "Correlation", ylab = "Freq", main = "Correlation: 500 iterations (N = 100)")


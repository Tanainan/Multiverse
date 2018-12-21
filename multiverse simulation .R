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

annlist <- c(1)
#deplist <- c("A")
all.data.multiverses <- list()
all.se <- list()

for (iii in 1:length(annlist)) {
  #rm(list = setdiff(ls(), c("annlist", "deplist", "all.se", "all.data.multiverses", "iii")))
  ann <- annlist[iii]
  
  if (ann == 1) {
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
    data20 <- simulateData(model = model, model.type = "sem",
                           sample.nobs = 20)
  }
  
  
  
  no.x <- 2
  no.y <- 2
  no.z <- 2
  
  data.multiverse <- array(list(), dim = c(no.x, no.y, no.z))
  se.multiverse <- array(0, dim = c(no.x, no.y, no.z))
  
  data.proc <- data20
  #data.proc <- data50
  #data.proc <- data100
  
  data.proc.init <- data.proc
  
  for (j in 1:no.x) {
    for (k in 1:no.y) {
      for (l in 1:no.z) {
        
        data.proc <- data.proc.init
        
        if (j == 1) {data.proc$x <- data20$x1} else {data.proc$x <- data20$x2}
        if (k == 1) {data.proc$y <- data20$y1} else {data.proc$y <- data20$y2}
        if (l == 1) {data.proc <- data.proc} 
        else {data.proc <- data.proc[(data.proc$x < 1.5),]} # listwise deletion
        
        data.multiverse[[j, k, l]] = data.proc
      }
    }
  }
  
  ### analysis ###
  for (j in 1:no.x) {
    for (k in 1:no.y) {
      for (l in 1:no.z) {
        data.proc$x <- factor(data.proc$x)
        data.proc$y <- factor(data.proc$y)
        
        an = coef(summary(lm(y ~ x, data.multiverse[[j, k, l]])))
        
        se.multiverse[j, k, l] <- an[2,2]
        
      }
    }     
  }     
  all.data.multiverses[[iii]] <- data.multiverse
  all.se[[iii]] <- se.multiverse
}

# average SE
mean(unlist(all.se))

# plot SE
graphnames <- c("Standard Error for N = 20")

hists <- list()
pv <- list()
ylabs=c("Frequency")
xlabs=c("SE")
for (iii in 1:length(annlist)) local({
  ann <- annlist[iii]
  se <- all.se[[ann]]
  if (ann == 1) {
    cat1 <- rep(c(1:15), 8)
    cat2 <- rep(1:8, each = 15)
  } else {
    cat1 <- rep(c(1:15), 14)
    cat2 <- rep(1:14, each = 15)
  }
  df <- data.frame(category1 = cat1, category2 = cat2, value = (as.vector(se[!is.na(se)])))
  df[["sign"]] = ifelse(df[["value"]] <= 0.05, "significant", "nonsignificant")
  pv[[ann]]=df$value
  hists[[ann]] <<- qplot(pv[[ann]], geom = "histogram", binwidth = 0.01) + xlim(0,1) + geom_histogram(colour = "black", fill = "white", binwidth = 0.01) + 
    xlab(xlabs[[ann]]) + ylab(ylabs[[ann]]) + #geom_vline(xintercept = 0.05, colour = "red", linetype = "longdash") + 
    ggtitle(graphnames[ann]) + theme(plot.title = element_text(lineheight = 0.8, face = "bold")) + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16))
  #windows(8, 5)
  print(hists[[ann]])
  rm(se)
  rm(df)
})


###### n = 50 ##############################################################
annlist <- c(1)
#deplist <- c("A")
all.data.multiverses50 <- list()
all.se50 <- list()

for (iii in 1:length(annlist)) {
  #rm(list = setdiff(ls(), c("annlist", "deplist", "all.se", "all.data.multiverses", "iii")))
  ann50 <- annlist[iii]
  
  if (ann50 == 1) {
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
    data50 <- simulateData(model = model, model.type = "sem",
                           sample.nobs = 50, seed = 1000)
  }
  
  
  
  no.x <- 2
  no.y <- 2
  no.z <- 2
  
  data.multiverse50 <- array(list(), dim = c(no.x, no.y, no.z))
  se.multiverse50 <- array(0, dim = c(no.x, no.y, no.z))
  
  data.proc50 <- data50
  #data.proc <- data50
  #data.proc <- data100
  
  data.proc.init50 <- data.proc50
  
  for (j in 1:no.x) {
    for (k in 1:no.y) {
      for (l in 1:no.z) {
        
        data.proc50 <- data.proc.init50
        
        if (j == 1) {data.proc50$x <- data50$x1} else {data.proc50$x <- data50$x2}
        if (k == 1) {data.proc50$y <- data50$y1} else {data.proc50$y <- data50$y2}
        if (l == 1) {data.proc50 <- data.proc50} 
        else {data.proc50 <- data.proc50[(data.proc50$x < 1.5),]} # listwise deletion
        
        data.multiverse50[[j, k, l]] = data.proc50
      }
    }
  }
  
  ### analysis ###
  for (j in 1:no.x) {
    for (k in 1:no.y) {
      for (l in 1:no.z) {
        data.proc50$x <- factor(data.proc50$x)
        data.proc50$y <- factor(data.proc50$y)
        
        an50 = coef(summary(lm(y ~ x, data.multiverse50[[j, k, l]])))
        
        se.multiverse50[j, k, l] <- an50[2,2]
        
      }
    }     
  }     
  all.data.multiverses50[[iii]] <- data.multiverse50
  all.se50[[iii]] <- se.multiverse50
}

# average SE
mean(unlist(all.se50))

# plot SE
graphnames <- c("Standard Error for N = 50")

hists <- list()
pv <- list()
ylabs=c("Frequency")
xlabs=c("SE")
for (iii in 1:length(annlist)) local({
  ann50 <- annlist[iii]
  se50 <- all.se50[[ann50]]
  if (ann50 == 1) {
    cat1 <- rep(c(1:15), 8)
    cat2 <- rep(1:8, each = 15)
  } else {
    cat1 <- rep(c(1:15), 14)
    cat2 <- rep(1:14, each = 15)
  }
  df <- data.frame(category1 = cat1, category2 = cat2, value = (as.vector(se50[!is.na(se50)])))
  df[["sign"]] = ifelse(df[["value"]] <= 0.05, "significant", "nonsignificant")
  pv[[ann50]]=df$value
  hists[[ann50]] <<- qplot(pv[[ann50]], geom = "histogram", binwidth = 0.01) + xlim(0,1) + geom_histogram(colour = "black", fill = "white", binwidth = 0.01) + 
    xlab(xlabs[[ann50]]) + ylab(ylabs[[ann50]]) + #geom_vline(xintercept = 0.05, colour = "red", linetype = "longdash") + 
    ggtitle(graphnames[ann50]) + theme(plot.title = element_text(lineheight = 0.8, face = "bold")) + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16))
  #windows(8, 5)
  print(hists[[ann50]])
  rm(se50)
  rm(df)
})



###### n = 100 ##############################################################
annlist <- c(1)
#deplist <- c("A")
all.data.multiverses100 <- list()
all.se100 <- list()

for (iii in 1:length(annlist)) {
  #rm(list = setdiff(ls(), c("annlist", "deplist", "all.se", "all.data.multiverses", "iii")))
  ann100 <- annlist[iii]
  
  if (ann100 == 1) {
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
    data100 <- simulateData(model = model, model.type = "sem",
                           sample.nobs = 100, seed = 1000)
  }
  
  
  
  no.x <- 2
  no.y <- 2
  no.z <- 2
  
  data.multiverse100 <- array(list(), dim = c(no.x, no.y, no.z))
  se.multiverse100 <- array(0, dim = c(no.x, no.y, no.z))
  

  data.proc100 <- data100
  
  data.proc.init100 <- data.proc100
  
  for (j in 1:no.x) {
    for (k in 1:no.y) {
      for (l in 1:no.z) {
        
        data.proc100 <- data.proc.init100
        
        if (j == 1) {data.proc100$x <- data100$x1} else {data.proc100$x <- data100$x2}
        if (k == 1) {data.proc100$y <- data100$y1} else {data.proc100$y <- data100$y2}
        if (l == 1) {data.proc100 <- data.proc100} 
        else {data.proc100 <- data.proc100[(data.proc100$x < 1.5),]} # listwise deletion
        
        data.multiverse100[[j, k, l]] = data.proc100
      }
    }
  }
  
  ### analysis ###
  for (j in 1:no.x) {
    for (k in 1:no.y) {
      for (l in 1:no.z) {
        data.proc100$x <- factor(data.proc100$x)
        data.proc100$y <- factor(data.proc100$y)
        
        an100 = coef(summary(lm(y ~ x, data.multiverse100[[j, k, l]])))
        
        se.multiverse100[j, k, l] <- an100[2,2]
        
      }
    }     
  }     
  all.data.multiverses100[[iii]] <- data.multiverse100
  all.se100[[iii]] <- se.multiverse100
}

# average SE
mean(unlist(all.se100))

# plot SE
graphnames <- c("Standard Error for N = 100")

hists <- list()
pv <- list()
ylabs=c("Frequency")
xlabs=c("SE")
for (iii in 1:length(annlist)) local({
  ann100 <- annlist[iii]
  se100 <- all.se100[[ann100]]
  if (ann100 == 1) {
    cat1 <- rep(c(1:15), 8)
    cat2 <- rep(1:8, each = 15)
  } else {
    cat1 <- rep(c(1:15), 14)
    cat2 <- rep(1:14, each = 15)
  }
  df <- data.frame(category1 = cat1, category2 = cat2, value = (as.vector(se100[!is.na(se100)])))
  df[["sign"]] = ifelse(df[["value"]] <= 0.05, "significant", "nonsignificant")
  pv[[ann100]]=df$value
  hists[[ann100]] <<- qplot(pv[[ann100]], geom = "histogram", binwidth = 0.01) + xlim(0,1) + geom_histogram(colour = "black", fill = "white", binwidth = 0.01) + 
    xlab(xlabs[[ann100]]) + ylab(ylabs[[ann100]]) + #geom_vline(xintercept = 0.05, colour = "red", linetype = "longdash") + 
    ggtitle(graphnames[ann100]) + theme(plot.title = element_text(lineheight = 0.8, face = "bold")) + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16))
  #windows(8, 5)
  print(hists[[ann100]])
  rm(se100)
  rm(df)
})

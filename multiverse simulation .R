##### multiverse ####
library(ggplot2)
library(Rmisc)
library(caret)
library(lm.beta)
library(reghelper)
library(glmnet)
library(DescTools)
library(lavaan)
library(tidyverse)
library(graphics)

#deplist <- c("A") # only need when want to test with different dependent variables (see Steeegen et al. 2016)
# For example, dependent variable A, B, C >> deplist <- c("A", "B", "C")
annlist <- c(1) # analysis 1
# create empty lists for multiverse outputs
all.data.multiverses <- list()
all.cor <- list() # for correlation outputs
all.z <- list() # for fisher transformation outputs
all.p <- list() # for all p-value of outputs
avg.cor <- list() # for mean of all correlation outputs in each iteration
sd.cor <- list() # for sd of all correlation outputs in each iteration
max.cor <- list() # for max of all correlation outputs in each iteration
signi.cor <- list() # for only significant correlation outputs in each iteration
mm.cor <- list() # for mean of multiverse means

for (i in 1:length(annlist)){ # since annlist = 1 >> compute 1 time for the whole process

  # simulate data from a model with given conditions
    model <- '
    #residual covariance
    x ~~ 0.64*y
    
    #factor
    x =~ 0.8*x1 + 0.8*x2
    y =~ 0.8*y1 + 0.8*y2
    
    #residual variance
    x1 ~~ 0.36*x1
    x2 ~~ 0.36*x2
    y1 ~~ 0.36*y1
    y2 ~~ 0.36*y2
    '

  no.it <- 10 # numbers of iteration
  no.a <- 4 # x2, (x1 + y2), (x2 + y1), (x2 + y2) as x2
  no.z <- 2 # complete and missing data
  no.x <- 2 # 2 ways to analyze a variable x, namely, x1 and x2
  no.y <- 2 # y1 and y2


  
  # indicates that there will be 4 dimensions for each matrix
  data.multiverse <- array(list(), dim = c(no.it, no.a, no.z, no.x, no.y))
  cor.multiverse <- array(0, dim = c(no.it, no.a, no.z, no.y, no.y))
  cor.signi <- array(0, dim = c(no.it, no.a, no.z, no.x, no.y))
  z.multiverse <- array(0, dim = c(no.it, no.a, no.z, no.x, no.y))
  p.multiverse <- array(0, dim = c(no.it, no.a, no.z, no.x, no.y))
  # (I'm not sure what is the difference between list() and 0. They might have something to do with the data source)
  
  # start simulating the data based on the model above and numbers of iterations
  for (a in 1:no.it){ # a is the order of iteration
  data <- simulateData(model = model, model.type = "sem",
                         sample.nobs = 100) # sample size
  data.proc <- data
  
  data.proc.init <- data.proc # this line is a MUST
  
  # prepare the analysis with all possible combinations of variables x and y and missing value or not
  for (j in 1:no.a) {
    data.proc <- data.proc.init # this line is also a MUST and needs to be right after the first for loop
                                # when there is a condition for that loop before going to the next for loop (e.g., "preparing the data")
                                # for example, in this case, I needed to first alter the variable x2 to something else
                                # before running no.z, no.x, no.y
                                # In Steegen et al., this line is after for loops when dealing with variable categorizations 
                                # because they didn't have to prepare the data overall
                                
    if (j == 1) {data.proc$x2 <- data.proc$x2}
    else if (j == 2) {data.proc$x2 <- (data.proc$x1 + data.proc$y2)}
    else if (j == 3) {data.proc$x2 <- (data.proc$x2 + data.proc$y1)}
    else if (j == 4) {data.proc$x2 <- (data.proc$x2 + data.proc$y2)}
    
  for (k in 1:no.z) { # need to create a completed vs missing data set first

        
        
        # we can have a complete data set vs not missing at random
        # c_names = c("x1","x2","y1","y2")
        # prc_missing = 0.15 # let's say 15% missing 
        # 
        # z1 <- data %>%
        #   add_column(id = 1:100) %>% # must change to sample size 1:N
        #   gather(var, value, -id) %>%    # reshape data
        #   mutate(r = runif(nrow(.)),   # simulate a random number from 0 to 1 for each row
        #          value = ifelse(var %in% c_names & r <= prc_missing, NA, value)) %>%  # if it's one of the variables you specified and the random number is less than your threshold update to NA
        #   select(-r) %>%                 # remove random number
        #   spread(var, value)             # reshape back to original format
        # 
        # if (k == 1) {data.proc <- data.proc} 
        # else {data.proc <- z1}
        
        # or we can try the listwise deletion for outliers in x1 and x2
        # MUST "#" the lines above from c_names to else {data.proc <- z1}
        if (k == 1) {data.proc <- data.proc}
        else {data.proc <- data.proc[(data.proc$x1 < 1.5 & data.proc$x2 < 1.5),]}
        # listwise deletion for those that are considered to be outliers 
        # in this case, I picked a number of 1.5 (for no reason)


        for (l in 1:no.x) { # then deal with variables x and y
          for (m in 1:no.y) {
            
            
        if (l == 1) {data.proc$x <- data.proc$x1} else {data.proc$x <- data.proc$x2}
        # this would indicate that we can use x1 (if k == 1) or x2 (if k ==2) as x in the real analysis
        if (m == 1) {data.proc$y <- data.proc$y1} else {data.proc$y <- data.proc$y2}
        
        # Here's what I did for data & variable preparation:
            # 1. determined what should be x2 (no.a)
            # 2. added a condition of missing values (random or not random) (no.z)
            # 3. determined between x1 and x2 which to be chosen as x in the analysis
            # 4. determined between y1 and y2 which to be chosen as y in the analysis
            
        data.multiverse[[a, j, k, l, m]] = data.proc # save the new data set into matrix
          }
        }
      }
    }
  }
  
  ### analysis ###
  for (a in 1:no.it){
    for (j in 1:no.a){
      for (k in 1:no.z){
        for (l in 1:no.x){
          for (m in 1:no.y){
        
          # since there are 4 dimensions, I only select x and y columns to find the correlation between them
          dt <- data.multiverse[[a, j, k, l, m]] %>%
            select(x, y)
          
          # test for correlations between x and y
          an = cor(dt$y, dt$x)
          cor.multiverse[a, j, k, l, m] <- an # save the output as a matrix
          
          zs <- FisherZ(an) # Fisher transformation
          z.multiverse[a, j, k, l, m] <- zs
          
          # with cor.test, we can get the other info, such as p-value
          test <- cor.test(dt$y, dt$x)
          if (test$p.value < .05) # only select the correlations that are significant, otherwise save as NA in the matrix
          {cor.signi[a, j, k, l, m] <- test$estimate} else {cor.signi[a, j, k, l, m] <- NA}
        
          # all p-values of correlations
          p.multiverse[a, j, k, l, m] <- test$p.value
          
          # for each iteration a, compute mean, sd, max, and mean of significant correlations and save them into the list we created at the beginning
          avg.cor[[a]] = mean(cor.multiverse[a,,,,]) 
          sd.cor[[a]] = sd(cor.multiverse[a,,,,])
          max.cor[[a]] = max(cor.multiverse[a,,,,])
          signi.cor[[a]] = mean(cor.signi[a,,,,], na.rm = T) # mean of only significant correlations in iteration a, excluding non-significant
          mm.cor = mean(cor.multiverse[,,,,])
          # these must be inside the same for loops of the correlation analysis, otherwise, the output be wrong
        
          }
        }
      }     
    } 
  }
  
  # save the matrix results into the list
  all.data.multiverses[[i]] <- data.multiverse 
  all.cor[[i]] <- cor.multiverse
  all.z[[i]] <- z.multiverse
  all.p[[i]] <- p.multiverse
}
# the orders of dimensions are flexible.

# plot multiverse distributions for each iteration
par(mfrow = c(2,2)) # make 4 distributions to be on the same page
for (i in 1:10){ # 10 iterations
  barplot(table(round(unlist(cor.multiverse[i,,,,]), digits = 2)), xlab = "Correlation", ylab = "Freq", main = "10 iterations (N = 100)")
}

# plot mean of each iteration
par(mfrow = c(1,1))
barplot(table(round(unlist(avg.cor), digits = 2)), xlab = "Mean", ylab = "Freq", main = "10 iterations (N = 10)")

# SD
barplot(table(round(unlist(sd.cor), digits = 2)), xlab = "SD", ylab = "Freq", main = "10 iterations (N = 10)")

# max
barplot(table(round(unlist(max.cor), digits = 2)), xlab = "Max", ylab = "Freq", main = "10 iterations (N = 10)")

# mean of all significant multiverse
barplot(table(round(unlist(signi.cor), digits = 2)), xlab = "Mean of Correlations of Significant Multiverses", ylab = "Freq", main = "10 iterations (N = 10)")

# mean of multiverse means ????
barplot(table(round(unlist(cor.multiverse), digits = 2)), xlab = "Mean of Multiverse Means", ylab = "Freq", main = "N = 100 for each iteration")


# grid for p-value
grids <- list()
for (iii in c(1)){ 
  ann <- annlist[iii]
  p <- all.p[[ann]]
  p.grid <- array(0,dim=c(no.it, no.z, no.x, no.y, no.a))  # change the dimensions of the p multiverse for visualization purposes
  for (jj in 1:3){
    for (jjj in 1:3){
      p.grid[, , , , jjj] <- p[, jjj, , , ]
    }
  }
  cat1 <- rep(c(1:20), 16) # 20 columns (a*j) and 16 rows (k*l*m)
  cat2 <- rep(1:16, each = 20)
  df <- data.frame(category1 = cat1, category2 = cat2, value = (as.vector
                                                                (p.grid[!is.na(p.grid)])))
  df[["sign"]] = ifelse(df[["value"]] <= 0.05, "significant", "nonsignificant")
  grids[[ann]] <- ggplot(df, aes(x = category1, y = category2, fill = sign)) +
    geom_tile(colour = "black") +
    geom_text(label = round((df$value), 2), size = 3, colour = "black") + 
    # for columns 
    # draw no.z (complete vs missing) branches vertical
    geom_segment(aes(x = 5.5, y = -1.7, xend = 5.5, yend = -0.3)) + 
    geom_segment(aes(x = 15.5, y = -1.7, xend = 15.5, yend = -0.3)) + 
    # draw no.z (complete vs missing) branches horizontal
    geom_segment(aes(x = 5.5, y = -1.7, xend = 15.5, yend = -1.7)) + 
    # draw no.it (iteration) branches vertical
    geom_segment(aes(x = 1, y = -0.3, xend = 1, yend = 0.5)) + 
    geom_segment(aes(x = 2, y = -0.3, xend = 2, yend = 0.5)) + 
    geom_segment(aes(x = 3, y = -0.3, xend = 3, yend = 0.5)) + 
    geom_segment(aes(x = 4, y = -0.3, xend = 4, yend = 0.5)) + 
    geom_segment(aes(x = 5, y = -0.3, xend = 5, yend = 0.5)) + 
    geom_segment(aes(x = 6, y = -0.3, xend = 6, yend = 0.5)) + 
    geom_segment(aes(x = 7, y = -0.3, xend = 7, yend = 0.5)) + 
    geom_segment(aes(x = 8, y = -0.3, xend = 8, yend = 0.5)) + 
    geom_segment(aes(x = 9, y = -0.3, xend = 9, yend = 0.5)) + 
    geom_segment(aes(x = 10, y = -0.3, xend = 10, yend = 0.5)) + 
    geom_segment(aes(x = 11, y = -0.3, xend = 11, yend = 0.5)) + 
    geom_segment(aes(x = 12, y = -0.3, xend = 12, yend = 0.5)) + 
    geom_segment(aes(x = 13, y = -0.3, xend = 13, yend = 0.5)) + 
    geom_segment(aes(x = 14, y = -0.3, xend = 14, yend = 0.5)) + 
    geom_segment(aes(x = 15, y = -0.3, xend = 15, yend = 0.5)) + 
    geom_segment(aes(x = 16, y = -0.3, xend = 16, yend = 0.5)) + 
    geom_segment(aes(x = 17, y = -0.3, xend = 17, yend = 0.5)) + 
    geom_segment(aes(x = 18, y = -0.3, xend = 18, yend = 0.5)) + 
    geom_segment(aes(x = 19, y = -0.3, xend = 19, yend = 0.5)) + 
    geom_segment(aes(x = 20, y = -0.3, xend = 20, yend = 0.5)) + 
    # draw no.it (iteration) branches horizontal
    geom_segment(aes(x = 1, y = -0.3, xend = 10, yend = -0.3)) + 
    geom_segment(aes(x = 11, y = -0.3, xend = 20, yend = -0.3)) + 
    
    # for rows (from the outiest lines)
    # draw no.a (x2 criteria) branches horizontal
    geom_segment(aes(x = 23.5, y = 2.5, xend = 25.5, yend = 2.5)) + 
    geom_segment(aes(x = 23.5, y = 6.5, xend = 25.5, yend = 6.5)) + 
    geom_segment(aes(x = 23.5, y = 10.5, xend = 25.5, yend = 10.5)) + 
    geom_segment(aes(x = 23.5, y = 14.5, xend = 25.5, yend = 14.5)) + 
    # draw no.a (x2 criteria) branches vertical
    geom_segment(aes(x = 25.5, y = 2.5, xend = 25.5, yend = 14.5)) + 
    # draw no.y (y) branches horizontal
    geom_segment(aes(x = 21.5, y = 1.5, xend = 23.5, yend = 1.5)) + 
    geom_segment(aes(x = 21.5, y = 3.5, xend = 23.5, yend = 3.5)) + 
    geom_segment(aes(x = 21.5, y = 5.5, xend = 23.5, yend = 5.5)) + 
    geom_segment(aes(x = 21.5, y = 7.5, xend = 23.5, yend = 7.5)) + 
    geom_segment(aes(x = 21.5, y = 9.5, xend = 23.5, yend = 9.5)) + 
    geom_segment(aes(x = 21.5, y = 11.5, xend = 23.5, yend = 11.5)) + 
    geom_segment(aes(x = 21.5, y = 13.5, xend = 23.5, yend = 13.5)) + 
    geom_segment(aes(x = 21.5, y = 15.5, xend = 23.5, yend = 15.5)) + 
    # draw no.y (y) branches vertical
    geom_segment(aes(x = 23.5, y = 1.5, xend = 23.5, yend = 3.5)) + 
    geom_segment(aes(x = 23.5, y = 5.5, xend = 23.5, yend = 7.5)) + 
    geom_segment(aes(x = 23.5, y = 9.5, xend = 23.5, yend = 11.5)) + 
    geom_segment(aes(x = 23.5, y = 13.5, xend = 23.5, yend = 15.5)) + 
    # draw no.x (x) branches horizontal
    geom_segment(aes(x = 20.5, y = 1, xend = 21.5, yend = 1)) + 
    geom_segment(aes(x = 20.5, y = 2, xend = 21.5, yend = 2)) + 
    geom_segment(aes(x = 20.5, y = 3, xend = 21.5, yend = 3)) + 
    geom_segment(aes(x = 20.5, y = 4, xend = 21.5, yend = 4)) + 
    geom_segment(aes(x = 20.5, y = 5, xend = 21.5, yend = 5)) + 
    geom_segment(aes(x = 20.5, y = 6, xend = 21.5, yend = 6)) + 
    geom_segment(aes(x = 20.5, y = 7, xend = 21.5, yend = 7)) + 
    geom_segment(aes(x = 20.5, y = 8, xend = 21.5, yend = 8)) + 
    geom_segment(aes(x = 20.5, y = 9, xend = 21.5, yend = 9)) + 
    geom_segment(aes(x = 20.5, y = 10, xend = 21.5, yend = 10)) + 
    geom_segment(aes(x = 20.5, y = 11, xend = 21.5, yend = 11)) + 
    geom_segment(aes(x = 20.5, y = 12, xend = 21.5, yend = 12)) + 
    geom_segment(aes(x = 20.5, y = 13, xend = 21.5, yend = 13)) + 
    geom_segment(aes(x = 20.5, y = 14, xend = 21.5, yend = 14)) + 
    geom_segment(aes(x = 20.5, y = 15, xend = 21.5, yend = 15)) + 
    geom_segment(aes(x = 20.5, y = 16, xend = 21.5, yend = 16)) + 
    # draw no.x (x) branches vertical
    geom_segment(aes(x = 21.5, y = 1, xend = 21.5, yend = 2)) + 
    geom_segment(aes(x = 21.5, y = 3, xend = 21.5, yend = 4)) + 
    geom_segment(aes(x = 21.5, y = 5, xend = 21.5, yend = 6)) + 
    geom_segment(aes(x = 21.5, y = 7, xend = 21.5, yend = 8)) + 
    geom_segment(aes(x = 21.5, y = 9, xend = 21.5, yend = 10)) + 
    geom_segment(aes(x = 21.5, y = 11, xend = 21.5, yend = 12)) + 
    geom_segment(aes(x = 21.5, y = 13, xend = 21.5, yend = 14)) + 
    geom_segment(aes(x = 21.5, y = 15, xend = 21.5, yend = 16)) + 
    annotate("text", x = c(5.5, 15.5), y = -2.2, label = c("z1", "z2")) + 
    annotate("text", x = 1:20, y = -0.8, label = rep(c("it1", "it2", "it3", 
                                                       "it4", "it5", "it6", "it7",
                                                       "it8", "it9", "it10"), 2)) + 
    annotate("text", x = 24.5, y = c(1.7, 5.7, 9.7, 13.7), label = c("a1", "a2", 
                                                        "a3", "a4")) + 
    annotate("text", x = 22.5, y = c(0.8, 2.8), label = c("y1", "y2")) + 
    annotate("text", x = 22.5, y = c(4.8, 6.8), label = c("y1", "y2")) + 
    annotate("text", x = 22.5, y = c(8.8, 10.8), label = c("y1", "y2")) + 
    annotate("text", x = 22.5, y = c(12.8, 14.8), label = c("y1", "y2")) + 
    
    annotate("text", x = 21, y = c(0.7, 1.7, 2.7, 3.7, 4.7, 5.7, 6.7, 7.7, 8.7, 
                                   9.7, 10.7, 11.7, 12.7, 13.7, 14.7, 15.7), label = rep(c("x1", "x2"), 8)) + 
    scale_fill_manual(values = c(significant = "grey", nonsignificant = "white")) + 
    scale_x_discrete(expand = c(0, 0)) + scale_y_reverse() + #ggtitle(graphnames[ann]) + 
    theme(plot.title = element_text(lineheight = 0.8, face = "bold")) + 
    theme(panel.grid.minor = element_blank()) + theme(panel.grid.major = element_blank()) + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank()) + theme(panel.background = element_rect(fill = "transparent")) + 
    theme(legend.position = "none") + theme() + xlab("") + ylab("")  
  # windows(30, 20)
  #windows(10, 7)
  
  print(grids[ann])
  rm(df)
  rm(p)
  rm(p.grid)
}




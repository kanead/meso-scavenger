#' Plot equilibrium data for extension of scavenger community DE model 
#' Adam Kane
#' 21/07/2023

#' housekeeping
rm(list=ls())
graphics.off()

#' load packages
library(tidyverse)
library(reshape2)
library(patchwork)
library(wesanderson)

#' load data
mydata <- read_csv("increase_n_f.csv", col_names = TRUE)

#' reshape your data
mydata <- melt(mydata, id.vars=c("Fraction"))
mydata <- mydata %>% mutate(fraction = Fraction)

#' plot mammals first
levels(as.factor(mydata$variable))
mydata1 <-
  mydata %>% filter(
      variable == "Jackals" |
      variable == "Hyenas" |
      variable == "Lions" 
  )

p1 <- ggplot(data = mydata1, aes(x = fraction , y = value, color = variable, fill = "species")) + geom_line(size=1.5)  + ylim(c(0,0.7)) +
  xlab(expression("fraction"~n[f])) + ylab(expression(paste("density", " (individuals/", km^2,")")))  + theme_bw() + 
  labs(color="group") +  #xlim(c(0,15000)) +
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73")) + ggtitle("B")
p1 


#' plot vultures 
mydata2 <-
  mydata %>% filter(
     variable == "Vultures")

p2 <- ggplot(data = mydata2, aes(x = fraction, y = value, color = variable)) + geom_line(size = 1.5) +
  xlab(expression("fraction"~n[f])) + ylab(expression(paste("density", " (individuals/", km^2,")"))) + theme_bw() + 
  labs(color="group") + #xlim(c(0,15000)) +
  scale_color_manual(values=c("#D55E00")) + ggtitle("A")
p2



#' plot carrion
mydata3 <-
  mydata %>% filter(
      variable == "Carrion 1"
  )

p3 <- ggplot(data = mydata3, aes(x = fraction, y = value, color = variable)) + geom_line(size = 1.5) +
  xlab(expression("fraction"~n[f])) +  ylab(expression(density (kg/km^2))) + theme_bw() + 
  labs(color="group") + ylim(c(0,0.005)) +
  scale_color_manual(labels = "Carrion C",values=c("#0072B2")) + ggtitle("C")
p3




#' plot carrion again
mydata4 <-
  mydata %>% filter(
      variable == "Carrion 2"
  )

p4 <- ggplot(data = mydata4, aes(x = fraction, y = value, color = variable)) + geom_line(size = 1.5) +
  xlab(expression("fraction"~n[f])) +  ylab(expression(density (kg/km^2))) + theme_bw() + 
  labs(color="group") + #xlim(c(0,15000)) +
  scale_color_manual(labels = "Carrion N",values=c("#0072B2")) + ggtitle("D")
p4

# build a panel of the plots 
(p2 | p1)/(p3 | p4)



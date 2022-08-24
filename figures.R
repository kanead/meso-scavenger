#' Plot data for scavenger community DE model 
#' Adam Kane
#' 16/04/2021

#' housekeeping
rm(list=ls())
graphics.off()

#' load packages
library(tidyverse)
library(reshape2)
library(patchwork)
library(wesanderson)

#' load data
mydata <- read_csv("maintextscenario.csv", col_names = TRUE)
# options 
# maintextscenario.csv
# halvedCin.csv 
# halvedKs.csv

#' reshape your data
mydata <- melt(mydata, id.vars=c("Time"))

#' could plot years instead of days on x axis
mydata <- mydata %>% mutate(years = Time/365)

#' plot mammals first
levels(as.factor(mydata$variable))
mydata1 <-
  mydata %>% filter(
      variable == "Jackals" |
      variable == "Hyenas" |
      variable == "Lions" 
  )

p1 <- ggplot(data = mydata1, aes(x = years , y = value, color = variable, fill = "species")) + geom_line(size=1.5)  + ylim(c(0,0.7)) +
  xlab("years") + ylab(expression(paste("density", " (individuals/", km^2,")")))  + theme_bw() + 
  labs(color="group") +  #xlim(c(0,15000)) +
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73")) + ggtitle("B")
p1 

#' years example
p1a <- ggplot(data = mydata1, aes(x = years, y = value, color = variable, fill = "species")) + geom_line(size=1.5)  + ylim(c(0,0.7)) +
  xlab("years") + ylab(expression(paste("density", " (individuals/", km^2,")")))  + theme_bw() + 
  labs(color="group") + # xlim(c(0,15000)) +
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73")) + ggtitle("A")
p1a 

#' plot vultures 
mydata2 <-
  mydata %>% filter(
     variable == "Vultures")

p2 <- ggplot(data = mydata2, aes(x = years, y = value, color = variable)) + geom_line(size = 1.5) +
  xlab("years") + ylab(expression(paste("density", " (individuals/", km^2,")"))) + theme_bw() + 
  labs(color="group") + #xlim(c(0,15000)) +
  scale_color_manual(values=c("#D55E00")) + ggtitle("A")
p2


#' plot carcass removal rates
mydata3 <-
  mydata %>% filter(
    variable == "JackalRemovalRate" |
      variable == "HyenaRemovalRate" |
      variable == "LionRemovalRate" |
      variable == "VultureRemovalRate" |
      variable == "DecayRate"
      
  )

p3 <- ggplot(data = mydata3, aes(x = years, y = value, color = variable)) + geom_line(size = 1.5) +
  xlab("years") + ylab(expression(paste("carcass removal (kg", " /", km^2,"/day)"))) + theme_bw() + 
  labs(color="group") + #xlim(c(0,15000)) 
  ggtitle("D") +
  scale_color_manual(labels = c("Jackals", "Hyenas", "Lions", "Vultures", "Other"), values=c("#E69F00","#56B4E9","#009E73","#D55E00","#000000")) 
p3



#' plot carrion 
mydata4 <-
  mydata %>% filter(
      variable == "Log10Carrion"
  )

p4 <- ggplot(data = mydata4, aes(x = years, y = value, color = variable)) + geom_line(size = 1.5) +
  xlab("years") +  ylab(expression(log[10]^{""}~(kg/km^2))) + theme_bw() + 
  labs(color="group") + #xlim(c(0,15000)) +
  scale_color_manual(labels = "Carrion",values=c("#0072B2")) + ggtitle("C")
p4

# build a panel of the plots 
(p2 | p1)/(p4 | p3)



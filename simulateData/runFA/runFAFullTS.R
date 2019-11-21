## AFGR 2019 Nov
## This code will be used to load simulated matlab EEG data
## run a FFT on the data - across subjects
## and then run a EFA on these data

## Load library(s)
install_load("R.matlab", "tidyverse","reshape2", "psych")

## Declare statics and set wd
setwd("~/Documents/dFAonEEG/simulateData/runFA/")

## Read the data
## This will be done in a for loop
in.data <- system("ls ../calhounScripts/ss2_create_sim/", intern = T)
out.values <- NULL
for(subj in in.data[1:3]){
  subj <- paste("../calhounScripts/ss2_create_sim/", subj, sep='')
  tmp.in <- readMat(subj)
  out.values <- apply(tmp.in$sources.sim, c(1,2), mean)
}

## Create a visualization for the sampiling procedure
# I will need to create a variable with the onsets and length of each sampiling window
# here is the matlab code to do so: onsets = round(1:(dur-olap/100*dur):(size(sources_sim,2)-dur)); %% considers overlap
# 0.0001    0.0188    0.0376    0.0563    0.0751
dur <- 750
olap <- 75
length_ts <- dim(tmp.in$sources.sim)[2]
onsets <-  floor(seq(1, length_ts, by=(dur-olap/100*dur)))
window <- onsets + 750
## Now organize these in a data frame for plotting
plot.data <- data.frame(start.points=onsets,
                        end.points=window)
plot.data$color <- 1:dim(plot.data)[1]
out.plot <- ggplot(plot.data[1:3,], aes(xmin=start.points, xmax=end.points, ymin=0, ymax=1, color=factor(color), fill=factor(color))) +
  geom_rect(alpha=.3) +
  theme_bw() +
  theme(legend.position = "none")

plot.data <- reshape2::melt(data = plot.data, direction = 'wide', idvar = 'color',measure.vars=c('start.points','end.points'))
out.plot.2 <- ggplot(data=plot.data, aes(x=value, y=color, group=color, color=factor(color))) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Epoch") +
  xlab("Time (ms)")

index <- which(plot.data$color %in% 1:10)
out.plot.2.1 <- out.plot.2 <- ggplot(data=plot.data[index,], aes(x=value, y=color, group=color, color=factor(color))) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Epoch") +
  xlab("Time (ms)")

## Now create a corllation matrix for the first epoch
subj <- paste("../calhounScripts/ss2_create_sim/", in.data[4], sep='')
tmp.in <- readMat(subj)
plot(tmp.in$outdata[1,,1], tmp.in$outdata[2,,1])
pc.loadings.one <- principal(r = as.matrix(tmp.in$outdata[,,1]), nfactors = 20)
cor.vals <- cor(tmp.in$outdata[,,1])
cor.vals <- melt(cor.vals)
out.plot.3 <- ggplot(data = cor.vals, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
## AFGR 2019 Nov
## This code will be used to load simulated matlab EEG data
## run a FFT on the data - across subjects
## and then run a EFA on these data

## Load library(s)
install_load("R.matlab", "tidyverse","reshape2", "psych", "greta","gganimate", "plotly")

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
  theme(legend.position = "none",text = element_text(size=20)) +
  ylab("Epoch") +
  xlab("Time (ms)")

index <- which(plot.data$color %in% 1:10)
out.plot.2.1 <- out.plot.2 <- ggplot(data=plot.data[index,], aes(x=value, y=color, group=color, color=factor(color))) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none",text = element_text(size=20)) +
  ylab("Epoch") +
  xlab("Time (ms)")

## Now create a corllation matrix for the first epoch
subj <- paste("../calhounScripts/ss2_create_sim/", in.data[4], sep='')
tmp.in <- readMat(subj)
out.fa <- fa()
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
                                   size = 12, hjust = 1),legend.position='none')+
  coord_fixed()

## Now plot the correlation matrix across time in a animated fashion
## In order to do this we are going to iterate through every epoch and create a corellation matrix
## and then collapse it into one object
dataForAnim <- NULL
for(t in 1:dim(tmp.in$outdata)[3]){
  cor.vals <- cor(tmp.in$outdata[,,t])
  cor.vals <- melt(cor.vals)
  cor.vals$epoch <- t
  dataForAnim <- rbind(dataForAnim, cor.vals)
}
out.plot.4 <- ggplot(data=dataForAnim, aes(x=Var1, y=Var2, fill=value, frame=epoch)) +
  geom_tile() +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),legend.position='none')+
  coord_fixed()
p <- ggplotly(out.plot.4)

plot(tmp.in$outdata[1,,1], tmp.in$outdata[2,,200])
pc.loadings.one <- principal(r = as.matrix(tmp.in$outdata[,,1]), nfactors = 4)
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
                                   size = 12, hjust = 1),legend.position='none') +
  coord_fixed()

## Now perform FA on the spatiospectral data
# first we are going to have to reformat the input data
my_matrix<-data.frame()
for(j in 1:dim(tmp.in$outdata)[3]){
  my_matrix <- rbind(my_matrix, tmp.in$outdata[,,j])
}
## Now we can fun the FA on the three d data collapsed across time
fa.res <- fa(my_matrix, nfactors = 4, rotate = "promax")

## Now I need to look at factor scores across time
score.mat <- matrix(NA, dim(my_matrix)[1], 4)
dim(score.mat) <- c(60,4,379)
for(q in 1:dim(score.mat)[3]){
  score.vals <- predict(fa.res, tmp.in$outdata[,,q])
  score.mat[,,q] <- score.vals
}

## Now try to produce the heat map from figure 2
# I am going to "unmix" the sources based on the FA weights
# in order to do this I am going to multiply the resultant scores by the individual weights from each components/factor
vals.one <- score.mat[,1,1] %*% t(fa.res$weights[,1])
tmp <- tmp.in$outdata
for(q in 1:dim(score.mat)[3]){
  tmp[,,q] <- score.mat[,1,q]%*% t(fa.res$weights[,1])
}
vals.one <- apply(tmp, c(1,2), mean)
tmp <- tmp.in$outdata
for(q in 1:dim(score.mat)[3]){
  tmp[,,q] <- score.mat[,2,q]%*% t(fa.res$weights[,2])
}
vals.two <- apply(tmp, c(1,2), mean)
vals.three <- score.mat[,3,1] %*% t(fa.res$weights[,3])
tmp <- tmp.in$outdata
for(q in 1:dim(score.mat)[3]){
  tmp[,,q] <- score.mat[,3,q]%*% t(fa.res$weights[,3])
}
vals.three <- apply(tmp, c(1,2), mean)
tmp <- tmp.in$outdata
for(q in 1:dim(score.mat)[3]){
  tmp[,,q] <- score.mat[,4,q]%*% t(fa.res$weights[,4])
}
vals.four <- apply(tmp, c(1,2), mean)

# This is essentially what figure 2 is plotting - acorss averaged across all epochs
# Now plot these bad boys
to.plot.one <- melt(vals.one)
to.plot.two <- melt(vals.two)
to.plot.three <- melt(vals.three)
to.plot.four <- melt(vals.four)
plot.one <- ggplot(to.plot.one, aes(x=Var1, y=Var2, fill=range01(value))) +
  geom_tile() +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue",  high = "red", mid = "white", 
                       midpoint=.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),legend.position='none') +
  coord_fixed()
plot.two <- ggplot(to.plot.two, aes(x=Var1, y=Var2, fill=range01(value))) +
  geom_tile() +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", mid = "red", high = "white", 
                       midpoint = .5, limit = c(range(to.plot.two$value)[1],range(to.plot.two$value)[2]), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),legend.position='none') +
  coord_fixed()
plot.three <- ggplot(to.plot.three, aes(x=Var1, y=Var2, fill=range01(value))) +
  geom_tile() +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", mid = "red", high = "white", 
                       midpoint = .5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),legend.position='none') +
  coord_fixed()
plot.four <- ggplot(to.plot.four, aes(x=Var1, y=Var2, fill=range01(value))) +
  geom_tile() +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", mid = "red", high = "white", 
                       midpoint = .5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),legend.position='none') +
  coord_fixed()

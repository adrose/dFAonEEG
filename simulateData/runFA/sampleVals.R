## AFGR 2019

##Usage##
## This script will be used to sample from a given distribution with known skew and kurtosis
## Possible sampiling distributions include: noraml & logistic
## However these will expand to include: %AFGR fill this out as needed%

# Turn off warnings
options(warn=-1)

###################################################################
# Load required libraries
###################################################################
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(psych)))


###################################################################
# Read in options
###################################################################
option_list = list(
  make_option(c("-d", "--dist"), action="store", default=NA, type='character',
              help="A type of distribution. Acceptable input: normal; logistic"),
  make_option(c("-m", "--meanVal"), action="store", default=NA, type='numeric',
              help="Mean value for the distribution."),
  make_option(c("-s", "--sigmaVal"), action="store", default=NA, type='numeric',
              help="Sigma value for the distribution."),
  make_option(c("-k", "--kurtosisVal"), action="store", default=NA, type='numeric',
              help="Kurtosis value for the distribution."),
  make_option(c("-S", "--skewVal"), action="store", default=NA, type='numeric',
              help="Skewness value for the distribution."),
  make_option(c("-n", "--sampleSize"), action="store", default=NA, type='numeric',
              help="output sample size."),
  make_option(c("-o", "--outputFile"), action="store", default=NA, type='character',
              help="CSV location for output.")
)
opt = parse_args(OptionParser(option_list=option_list))

if(is.na(opt$dist)){
  opt$dist <- 'normal'
}
opt$dist <- tolower(opt$dist) 
if(is.na(opt$dist)) {
  opt$dist <- "normal"
} 
dist <- opt$dist
if (is.na(opt$meanVal)) {
  meanVal <- 0
} else{
  meanVal <- opt$meanVal
}
if (is.na(opt$sigmaVal)) {
  sigmaVal <- 1
} else {
  sigmaVal <- opt$sigmaVal
}
if (is.na(opt$kurtosisVal)) {
  kurtosisVal <- 4
} else{
  kurtosisVal <- opt$kurtosisVal
}
if (is.na(opt$skewVal)) {
  skewVal <- 1.5
} else{
  skewVal <- opt$skewVal
}
if (is.na(opt$sampleSize)){
  cat('User did not specify sample size.\n')
  cat('Use sampleVals.R -h for an expanded usage menu.\n')
  quit()
}
sampleSize <- opt$sampleSize
out <- opt$outputFile

###################################################################
# Sample from normal distribution here
###################################################################
if(dist == "normal"){
  suppressMessages(suppressWarnings(library(PearsonDS)))
  moments <- c(mean = meanVal,variance = sigmaVal,skewness = skewVal, kurtosis = kurtosisVal)
  toReturn <- rpearson(sampleSize, moments = moments)
}

###################################################################
# Sample from logistic distribution here
###################################################################
if(dist == "logistic"){
  suppressMessages(suppressWarnings(library(genlogis)))
  toReturn <- genlogis::rgenlog_sk(sampleSize, mu=meanVal, skew = skewVal, a = sqrt(2/pi), b = 0.5, p = 2)
}

###################################################################
# Write the output
###################################################################
write.csv(toReturn, out, quote = F, row.names = F)
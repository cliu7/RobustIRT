% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robust.R
\name{choose.tuco}
\alias{choose.tuco}
\title{Plot histogram of residuals along plot of weight (dependent on TuCo) vs residuals}
\usage{
choose.tuco(r, H = NULL, B = NULL)
}
\arguments{
\item{r}{A vector of residuals}

\item{H}{Huber tuning parameter}

\item{B}{Bisquare tuning parameter}
}
\value{
Histogram plot of residuals beneath a graph of the weight functions vs. the residuals
}
\description{
This function plots a histogram of residuals along the graph of the weighting function (dependent on the tuning parameter) as a function of the residual
}
\details{
The goal of this plot is to visualize the proportion of residuals that are downweighted based on the tuning parameter and allow the researcher to choose a tuning parameter that suits their data well.
}
\examples{
########## Choose.tuco example - Unidimensional IRT ######
n=40
thetas<-matrix(seq(0,2, by=.05), ncol=1)#generate real thetas
a<-matrix(runif(n, .5, 1.5), ncol=1) #set item slope
d<-rnorm(n)
#Generate probabilities of correct response as if working at a suboptimal level of -1
theta.drop<-.5
chng.pt<-0.6
probs<-rbind(apply(thetas, 1, function(x) probs.gen(x, matrix(a[1:(chng.pt*n)], ncol=1), d[1:(chng.pt*n)])$P), apply(thetas-theta.drop, 1, function(x) probs.gen(x, matrix(a[(chng.pt*n+1):n,], ncol=1), d[(chng.pt*n+1):n])$P))
dat<-apply(probs, c(1, 2), function(x) rbinom(1, 1, x))
example<-theta.est(dat, a, d, iter=30, cutoff=.01, init.val=rep(0,ncol(a)), weight.type="equal", tuning.par=NULL)
theta_plots(dat, a, d=d, iter=30, cutoff=0.01, H=1, B=4, type='MIRT')
choose.tuco(r=matrix(na.omit(example$residual), ncol=1), B=4)

########## Choose.tuco example - MIRT#######
library(mirt)
data(SAT12)
SAT12[SAT12 == 8] <- NA #set 8 as a missing value

#correct answer key
key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
scoredSAT12 <- key2binary(SAT12, key)
specific <- c(2, 3, 2, 3, 3, 2, 1, 2, 1, 1, 1, 3, 1, 3, 1, 2, 1, 1, 3, 3, 1, 1, 3, 1, 3, 3, 1, 3, 2, 3, 1,2) #which factor each item loads on
b_mod1 <- mirt(scoredSAT12, specific)
ipars<-matrix(unlist(coef(b_mod1))[1:(32*6)], nrow = length(key), byrow=T) #item parameters
a <- ipars[,1:3]
d<- ipars[,4]

########### Choose.tuco example - GRM ######
n=40
nthresh<-4
thetas<-seq(-2,2.1, by=.1)#generate real thetas

a<-runif(n, .90, 2.15) #set item slope
b<- matrix(runif(n*nthresh, -2.5,2.5), nrow = n, ncol =nthresh)
b<-t(apply(b, 1, sort)) #category threshold parameters
probs<-probs.gen.grm(thetas, a, b)
dat<-data.gen(probs$P)
abdat<-dat
chng.pt<-.6 #random guessing for latter 30\% of the exam
abdat[(chng.pt*n+1):n, ]<-sample(c(1:(nthresh+1)), length(thetas)*(n-chng.pt*n), replace = T)
mle<-theta.est.grm(dat, a, b, iter=30, cutoff=0.01, init.val=0, weight.type="equal")
choose.tuco(matrix(mle$residual), H=.1, B=.8)

######## theta_plots example - MIRT ##########
data(SAT12)
SAT12[SAT12 == 8] <- NA #set 8 as a missing value

#correct answer key
key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
scoredSAT12 <- key2binary(SAT12, key)
specific <- c(2, 3, 2, 3, 3, 2, 1, 2, 1, 1, 1, 3, 1, 3, 1, 2, 1, 1, 3, 3, 1, 1, 3, 1, 3, 3, 1, 3, 2, 3, 1,2) #which factor each item loads on
b_mod1 <- mirt(scoredSAT12, specific)
ipars<-matrix(unlist(coef(b_mod1))[1:(32*6)], nrow = length(key), byrow=T) #item parameters

### Set Parameters
a <- ipars[,1:3]
d<- ipars[,4]
dat<-scoredSAT12[!is.na(rowSums(scoredSAT12)),] # removing vectors with missing data
colnames(dat)<-NULL
dat<-scoredSAT12[!is.na(rowSums(scoredSAT12)),] # removing vectors with missing data
colnames(dat)<-NULL

out<-theta.est(t(dat), a, d, iter=30, cutoff=.01, weight.type="equal")
choose.tuco(matrix(out$residual[,,2]), H=1, B=4)
}

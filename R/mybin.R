#' Binomial Experiment
#'
#' This function simulates binomial experiments by sampling n Bernoulli
#' trials with probability p of success, repeated "iter" times.
#' It also produces a barplot of the emperical proportions and returns
#' the table.
#'
#' @param iter Integer. Number of iterations (Simulated experiments).
#' @param n Integer. Number of trials in each experiment.
#' @param p Numerical. Probability of success in each trail.
#'
#' @returns A named table of empirical proportions for 0 to n successes.
#' @export
#'
#' @examples
#' mybin(iter = 100, n = 10, p = 0.5)
#' mybin(iter = 1000, n = 20, p = 0.3)
mybin=function(iter=100,n=10, p=0.5){
  # make a matrix to hold the samples
  #initially filled with NA's
  sam.mat=matrix(NA,nrow=n,ncol=iter, byrow=TRUE)
  #Make a vector to hold the number of successes in each trial
  succ=c()
  for( i in 1:iter){
    #Fill each column with a new sample
    sam.mat[,i]=sample(c(1,0),n,replace=TRUE, prob=c(p,1-p))
    #Calculate a statistic from the sample (this case it is the sum)
    succ[i]=sum(sam.mat[,i])
  }
  #Make a table of successes
  succ.tab=table(factor(succ,levels=0:n))
  #Make a barplot of the proportions
  barplot(succ.tab/(iter), col=rainbow(n+1), main="Binomial simulation", xlab="Number of successes")
  succ.tab/iter
}

#'threshold tuning
#'
#'tunes the thrsehold for the FLOC algorithm to reach a desired false alarm
#'probability
#'
#'@param NJ integer; number of data points in a jump bin
#'@param NK integer; number of data points in a kink bin
#'@param tau integer; number of observations so that P(hat tau < tau) <= eta
#'@param k integer; number of observations that are available as historical data
#'@param eta numerical between 0 ane 1, probability such that
#'P(hat tau < tau) <= eta
#'@param arl integer; desired average run length, tau has to be NULL to tune for
#'the average run length
#'@param r integer; number of simulated maximal test statistics
#'@param data list with at least k + tau rows and r columns;
#'data to tune the threshold on, if NULL standard normal data is generated
#'
#'@return list with:
#'jump - threshold for only jump detector
#'kink - threshold for only kink detector
#'both_jump - threshold for jump part if running both detectors
#'both_kink -threshold for kink part if running both detectors
#'
#'
#'@export

t_tuning <- function(NJ, NK, k,tau = NULL, eta = NULL, arl = NULL, r = 1000,
                     data = NULL){

  #tune for the average run length in the case, that tau is set NULL
  if(is.null(tau)){
    tau <- arl
    eta <- 1 - (1 / exp(1))
  }

  if(is.null(data)){
    #generates r sets of data and calculates the maximal test statistics
    max_Tstat <- replicate(r, max_Tstat_calc(NJ, NK, tau, k), simplify = FALSE)
  }
  else{
    r <- ncol(data)
    max_Tstat <- sapply(1:r,function(i) max_Tstat_calc(NJ, NK, tau,
                                              k, data[,i]), simplify = FALSE)
  }


  #sort the maximal test statistic
  max_Tstat <- do.call(rbind,max_Tstat)
  max_Tstat_sorted <- apply(max_Tstat, 2,sort)

  #the sought quantile
  q <- floor((1-eta) * r) + 1

  #setting the jump and kink threshold as the quantile
  rhoj <- max_Tstat_sorted[q,"jump"]
  rhok <- max_Tstat_sorted[q,"kink"]

  #looking for the correct combined quantile
  while(!max_quantile(q, r, eta, max_Tstat, max_Tstat_sorted)){
    q <- q + 1
  }

  #setting the combined threshold as the quantile
  rhobj <- max_Tstat_sorted[q,"jump"]
  rhobk <- max_Tstat_sorted[q,"kink"]

  thresholds = list(jump =  unname(rhoj), kink = unname(rhok), both_jump = unname(rhobj), both_kink = unname(rhobk))

  return(thresholds)
}

#'maximal Test statistic calculate
#'
#'Internal function that generates data and calculates the maximal test
#'statistic
#'
#'@param NJ integer; number of data points in a jump bin
#'@param NK integer; number of data points in a kink bin
#'@param tau integer; number of observations so that P(hat tau < tau) <= eta
#'@param k integer; number of observations that are available as historical data
#'@param datavec numerical vector; data to calculate the maximal Test statistic
#'on
#'
#'@return numerical vector; the maximal jump test statistic and the maximal kink
#'test statistic
#'
#'@importFrom stats rnorm
#'@importFrom utils tail
#'
#'@keywords internal

max_Tstat_calc <- function(NJ, NK, tau, k, datavec = NULL){

  if(is.null(datavec)){
  #generate historical data if no data is given
  data <- rnorm(k)
  }
  else{
  data <- head(datavec, k)
  }

  #estimate prechange distribution
  x <- (1-k):0
  f <- lm(data ~ x)


  if(is.null(datavec)){
    #generate data if no data vector is given
    data <- c(tail(data, 3 * max(NJ, NK)), rnorm(tau))
  }
  else{
    data <- c(tail(data, 3 * max(NJ, NK)), tail(datavec, -k ))
  }

  data <- data - unname(predict(f,list(x = (1 - 3 * max(NJ, NK)):tau)))

  #initialize bins
  SJ <- S_init(data[1:(3*NJ)], NJ)
  SK <- S_init(data[1:(3*NK)], NK)
  W <- W_init(data[1:(3*NK)], NK)

  #initialize vector of test statistics
  Tstat_jump <- c()
  Tstat_kink <- c()

  #iterates to calculate the test statistics
  for(i in 1:tau){

    #number of observations in last bin
    rJ <- (i -1) %% NJ + 1
    rK <- (i -1) %% NK + 1

    #shift the bins if a new bin is started
    if(rJ == 1){
      SJ <- c(SJ[-1],0)
    }

    if(rK == 1){
      SK <- c(SK[-1],0)
      W <- c(W[-1],0)
    }

    #update bins
    SJ[3] <- SJ[3] + data[i + 3 * NJ]
    SK[3] <- SK[3] + data[i + 3 * NK]
    W[3] <- W[3] + rK * data[i + 3 * NK]

    #calculate new test statistics
    Tstat_jump <- c(Tstat_jump,Tstat_j_calc(SJ, NJ, rJ))
    Tstat_kink <- c(Tstat_kink, Tstat_k_calc(SK, W, NK, rK))
  }

  #take maximal test statistic
  max_Tstat = c(max(abs(Tstat_jump)), max(abs(Tstat_kink)))
  names(max_Tstat) = c("jump", "kink")

  return(max_Tstat)
}

#'maximal test statistic quantile
#'
#'Internal function that takes a quantile of the jump and kink and checks if the
#'desired quantile for both test statistics combined is reached.
#'
#'@param q integer; quantile of the singular test statistics
#'@param r integer; number of maximal test statistics
#'@param eta numerical between 0 and 1; desired probability for a false alarm
#'@param max_Tstat numercial list; maximal test statistics
#'@param max_Tstat_sorted numerical list; maximal test statistics sorted
#'
#'@return boolean; true if the desired quantile is reached
#'
#'@keywords internal

max_quantile <- function(q, r, eta,  max_Tstat, max_Tstat_sorted){

  #counts test statistics that are greater than  the threshold at quantile q
  detect <- apply(max_Tstat, 1, function(x) any(x >= max_Tstat_sorted[q,]))
  detect <- sum(detect)

  #returns true if the count is below the desired quantile
  return(detect < r * eta)

}

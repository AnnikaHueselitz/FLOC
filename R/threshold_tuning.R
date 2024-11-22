#'threshold tuning
#'
#'tunes the thrsehold for the FLOC algorithm to reach a desired false alarm
#'probability
#'
#'@param N integer; number of data points in a bin
#'@param tau integer; number of observations so that P(hat tau < tau) <= eta
#'@param k integer; number of observations that are available as historical data
#'@param eta numerical between 0 ane 1, probability such that
#'P(hat tau < tau) <= eta
#'@param r integer; number of simulated maximal test statistics
#'
#'@return list with:
#'rhoj - threshold for only jump
#'rhok - threshold for only kink
#'rhob - thrshold for jump and kink
#'
#'@export

t_tuning <- function(N, tau, k, eta = 0.5, r = 1000){

  #generates r sets of data and calculates the maximal test statistics
  m_Tstat <- replicate(r, m_Tstat_calc(N, tau, k), simplify = FALSE)

  #sort the maximal test statistic
  m_Tstat <- do.call(rbind,m_Tstat)
  m_Tstat_sorted <- apply(m_Tstat, 2,sort)

  #the sought quantile
  q <- floor((1-eta) * r) + 1

  #setting the jump and kink threshold as the quantile
  rhoj <- m_Tstat_sorted[q,"jump"]
  rhok <- m_Tstat_sorted[q,"kink"]

  #looking for the correct combined quantile
  while(!m_quantile(q, r, eta, m_Tstat, m_Tstat_sorted)){
    q <- q + 1
  }

  #setting the combined threshold as the quantile
  rhob <- c(m_Tstat_sorted[q,"jump"], m_Tstat_sorted[q,"kink"])

  thresholds = list(jump =  unname(rhoj), kink = unname(rhok), both = rhob)

  return(thresholds)
}

#'maximal Test statistic calculate
#'
#'Internal function that generates data and calculates the maximal test
#'statistic
#'
#'@param N integer; number of data points in a bin
#'@param tau integer; number of observations so that P(hat tau < tau) <= eta
#'@param k integer; number of observations that are available as historical data
#'
#'@return numerical vector; the maximal jump test statistic and the maximal kink
#'test statistic
#'
#'@importFrom stats rnorm
#'@importFrom utils tail
#'
#'@keywords internal

m_Tstat_calc <- function(N, tau, k){

  #generate historical data
  data <- rnorm(k)

  #estimate prechange distribution
  x <- (1-k):0
  f <- lm(data ~ x)

  #generate new data
  data <- c(tail(data, 3 * N),rnorm(tau))
  data <- data - unname(predict(f,list(x = (1 - 3 * N):tau)))

  #initialize bins
  S <- S_init(data[1:(3*N)], N)
  W <- W_init(data[1:(3*N)], N)

  #initialize vector of test statistics
  Tstat_jump <- c()
  Tstat_kink <- c()

  #iterates to calculate the test statistics
  for(i in 1:tau){

    #number of observations in last bin
    r <- (i -1) %% N + 1

    #shift the bins if a new bin is started
    if(r == 1){
      S <- c(S[-1],0)
      W <- c(W[-1],0)
    }

    #update bins
    S[3] <- S[3] + data[i + 3 * N]
    W[3] <- W[3] + r * data[i + 3 * N]

    #calculate new test statistics
    Tstat_jump <- c(Tstat_jump,Tstat_j_calc(S, N, r))
    Tstat_kink <- c(Tstat_kink, Tstat_k_calc(S, W, N, r))
  }

  #take maximal test statistic
  m_Tstat = c(max(abs(Tstat_jump)), max(abs(Tstat_kink)))
  names(m_Tstat) = c("jump", "kink")

  return(m_Tstat)
}

#'maximal test statistic quantile
#'
#'Internal function that takes a quantile of the jump and kink and checks if the
#'desired quantile for both test statistics combined is reached.
#'
#'@param q integer; quantile of the singular test statistics
#'@param r integer; number of maximal test statistics
#'@param eta numerical between 0 and 1; desired probability for a false alarm
#'@param m_Tstat numercial list; maximal test statistics
#'@param m_Tstat_sorted numerical list; maximal test statistics sorted
#'
#'@return boolean; true if the desired quantile is reached
#'
#'@keywords internal

m_quantile <- function(q, r, eta,  m_Tstat, m_Tstat_sorted){

  #counts test statistics that are greater than  the threshold at quantile q
  detect <- apply(m_Tstat, 1, function(x) any(x >= m_Tstat_sorted[q,]))
  detect <- sum(detect)

  #returns true if the count is below the desired quantile
  return(detect < r * eta)

}

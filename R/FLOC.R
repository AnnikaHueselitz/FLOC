#'FLOC initialize
#'
#'Initializes the FLOC, Fast Limited-memory Optimal Change detector.
#'
#'@param data numerical vector; historical data that is available to estimate 
#'the prechange distribution from.
#'@param N integer; number of data points in a bin
#'@param rhoj numerical; threshold for the jump part of the algorithm.
#'@param rhok numerical; threshold for the kink part of the algorithm.
#'
#'@return list with
#'detect_jump - Boolean; True if a detection was made by the jump test statistic
#' in this iteration
#'detect_kink - Boolean; True if a detection was made by the kink test statistic
#' in this iteration
#'iteration - integer; number of the iteration
#'S - numerical vector; Sums of the bins
#'W - numerical vector; weighted sums of the bins
#'Tstat - numerical vector; the current test statistic
#'f - objedt of class "lm"; estimate of the pre change distribution
#'N - integer; number of data points in a bin
#'rhoj - numerical; threshold for the jump part of the algorithm.
#'rhok - numerical; threshold for the kink part of the algorithm.
#'
#'@export
#'
#'@importFrom stats lm predict

FLOC_init <- function(data, N, rhoj, rhok){
  
  #indices of the historical data
  x <- (1-length(data)):0
  
  #prediction of the prechange distribution
  f <- lm(data ~ x)
  data <- data - unname(predict(f,list(x)))
  
  #initialize the bins
  S <- S_init(data, N)
  W <- W_init(data, N)
  
  #initialize the test statistic
  Tstat <- c(Tstat_j_calc(S,N,N),Tstat_k_calc(S,W,N,N))
  
  return(list(detect_jump = FALSE,
              detect_kink = FALSE,
              iteration = 0,
              S = S,
              W = W,
              Tstat = Tstat,
              f = f,
              N = N,
              rhoj = rhoj,
              rhok = rhok
              ))
}

#'FLOC iteration
#'
#'Iterates the FLOC, Fast Limited-memory Optimal Change detector, for a new 
#'observation and detectcs if a change happened.
#'
#'@param data numerical; newest observation 
#'@param detector list; as returned by FLOC_init or FLOC_iter
#'
#'@return list with
#'detect_jump - Boolean; True if a detection was made by the jump test statistic
#' in this iteration
#'detect_kink - Boolean; True if a detection was made by the kink test statistic
#' in this iteration
#'iteration - integer; number of the iteration
#'S - numerical vector; Sums of the bins
#'W - numerical vector; weighted sums of the bins
#'Tstat - numerical vector; the current test statistic
#'f - objedt of class "lm"; estimate of the pre change distribution
#'N - integer; number of data points in a bin
#'rhoj - numerical; threshold for the jump part of the algorithm.
#'rhok - numerical; threshold for the kink part of the algorithm.
#'
#'@export
#'
#'@importFrom stats predict

FLOC_iter <- function(data,detector){
  
  #needed variables from detector
  iteration = detector$iteration
  S = detector$S
  W = detector$W
  f = detector$f
  N = detector$N
  rhoj = detector$rhoj
  rhok = detector$rhok
  
  #count iteration
  iteration <- iteration + 1
  
  #number of observations in last bin
  r <- (iteration -1) %% N + 1
  
  #shift the bins if a new bin is started
  if(r == 1){
    S <- c(S[-1],0)
    W <- c(W[-1],0)
  }
  
  #substract prechange prediction
  data <- data - unname(predict(f,list(x = iteration)))
  
  #add observation to bins
  S[3] <- S[3] + data
  W[3] <- W[3] + r * data
  
  #calculate test statistics
  Tstat_jump <- Tstat_j_calc(S, N, r)
  Tstat_kink <- Tstat_k_calc(S, W, N, r)
  
  #do detections
  detect_jump <- Tstat_jump >= rhoj
  detect_kink <- Tstat_kink >= rhok
  
  return(list(detect_jump = detect_jump,
              detect_kink = detect_kink,
              iteration = iteration,
              S = S,
              W = W,
              Tstat = c(Tstat_jump,Tstat_kink),
              f = f,
              N = N,
              rhoj = rhoj,
              rhok = rhok
  ))
}

#'S initialize
#'
#'Internal Function that initialize the sum bins.
#'
#'@param data numerical vector; historical data to initialize the bins from
#'@param N integer; number of observations in a bin
#'
#'@return S - numerical vector; the sum bins
#'
#'@keywords internal

S_init <- function(data, N){
  
  #only the last 3*N data points are needed
  data <- tail(data, 3 * N)
  
  #sort data in three bins of length N
  S1 <- sum(data[1:N])
  S2 <- sum(data[(N+1):(2*N)])
  S3 <- sum(data[(2 * N + 1):(3 * N)])
  
  return(c(S1,S2,S3))
}

#'W initialize
#'
#'Internal Function that initialize the weighted sum bins.
#'
#'@param data numerical vector; historical data to initialize the bins from
#'@param N integer; number of observations in a bin
#'
#'@return W - numerical vector; the weighted sum bins
#'
#'@keywords internal

W_init <- function(data, N){
  
  #sort data in three bins of length N
  data <- tail(data, 3 * N)
  
  #sort data in three bins of length N and weight it accordingly
  W1 <- sum(data[1:N] * (1:N))
  W2 <- sum(data[(N+1):(2*N)] * (1:N) )
  W3 <- sum(data[(2 * N + 1):(3 * N)] * (1:N))
  
  return(c(W1,W2,W3))
}

#'Test statistic jump calculate
#'
#'Calculates the jump test statistic.
#'
#'@param S numerical vector; Sums of the bins
#'@param N integer; number of observations in a bin
#'@param r integer; number of observations in the last bin
#'
#'@return numerical; the jump test statistic
#'
#'@keywords internal

Tstat_j_calc <- function(S, N, r){
  
  #number of observations in all bins
  M <- 2 * N + r
  
  #calculate the test statistic
  Tstat <- sum(S) / M
  
  return(abs(Tstat))
}

#'Test statistic kink calculate
#'
#'Calculates the kink test statistic.
#'
#'@param S numerical vector; Sums of the bins
#'@param W numerical vector; Weighted sums of the bins
#'@param N integer; number of observations in a bin
#'@param r integer; number of observations in the last bin
#'
#'@return numerical; the kink test statistic
#'
#'@keywords internal


Tstat_k_calc <- function(S, W, N, r){
  
  #number of observations in all bins
  M <- 2 * N + r
  
  #factor for weighting the test statistic
  d <- sum((1:M) * (1:M))
  
  #calculate the test statistic
  Tstat <- sum(W) + 2 * N * S[3] + N * S[2]
  Tstat <- Tstat/d
  
  return(abs(Tstat))
}
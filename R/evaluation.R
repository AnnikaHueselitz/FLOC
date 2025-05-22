#'estimate false alarm probability
#'
#'Estimates the probability of FLOC that a false alarm occurs before tau observations
#'for given bin size and thresholds
#'
#'@param NJ integer; number of data points in a jump bin
#'@param NK integer; number of data points in a kink bin
#'@param rhoj numerical; threshold for the jump part of the algorithm.
#'@param rhok numerical; threshold for the kink part of the algorithm.
#'@param tau integer; number of observation to check for a false alarm
#'@param k integer; number of observations that are available as historical data
#'@param m integer; = 200; number of repetitions
#'
#'@return numerical between 0 and 1; estimate of the false alarm probability
#'
#'@export

est_faprob <- function(NJ, NK, rhoj, rhok, tau, k, m = 200){

  #simulates m times data and checks for a false alarm
  fa <- replicate(m, false_alarm(NJ, NK, rhoj, rhok, tau, k))

  return(sum(fa)/m)
}

#'#false alarm probability
#'
#'Internal function that checks if a false alarm occurs before tau observations
#'
#'@param NJ integer; number of data points in a jump bin
#'@param NK integer; number of data points in a kink bin
#'@param rhoj numerical; threshold for the jump part of the algorithm.
#'@param rhok numerical; threshold for the kink part of the algorithm.
#'@param tau integer; number of observation to check for a false alarm
#'@param k integer; number of observations that are available as historical data
#'
#'@return boolean; true if a false alarm was occurs
#'
#'@importFrom stats rnorm
#'
#'@keywords internal

false_alarm <- function(NJ, NK, rhoj, rhok, tau, k){

  #generate historic data and initialize the detector
  data <- rnorm(k)
  detector <- FLOC_init(data, NJ, NK, rhoj, rhok)


  while(detector$iteration <= tau){

    #iterate the detector with a new data point
    data <- rnorm(1)
    detector <- FLOC_iter(data, detector)

    if(detector$detect_jump | detector$detect_kink){

      #a false detection was made
      return(TRUE)

    }

  }

  #No false alarm in tau observations
  return(FALSE)
}

#'estimate expected detection delay
#'
#'Estimates the expected detection delay for given bin size, thresholds and
#'change sizes
#'
#'@param NJ integer; number of data points in a jump bin
#'@param NK integer; number of data points in a kink bin
#'@param rhoj numerical; threshold for the jump part of the algorithm.
#'@param rhok numerical; threshold for the kink part of the algorithm.
#'@param k integer; number of observations that are available as historical data
#'@param jump numerical; size of the jump, size of change in the intercept at
#'the changepoint
#'@param kink numerical; size of the kink, size of change in the slope
#'@param m integer; = 200; number of repetitions
#'
#'@return integer; estimate of the expected detection delay rounded to the next
#'integer
#'
#'@export

est_detdel <- function(NJ, NK, rhoj, rhok, k, jump, kink, m = 200){

  #generates data m times and reports the detection delay
  detdel <- replicate(m, detection_delay(NJ, NK, rhoj, rhok,  k, jump, kink))

  return(round(mean(detdel)))
}

#'detection delay
#'
#'Internal function that generates data and reports the detection delay of FLOC
#'
#'@param NJ integer; number of data points in a jump bin
#'@param NK integer; number of data points in a kink bin
#'@param rhoj numerical; threshold for the jump part of the algorithm.
#'@param rhok numerical; threshold for the kink part of the algorithm.
#'@param k integer; number of observations that are available as historical data
#'@param jump numerical; size of the jump, size of change in the intercept at
#'the changepoint
#'@param kink numerical; size of the kink, size of change in the slope
#'@param max_n integer; number of observations after which the algoithm stops
#'
#'@return integer; detection delay
#'
#'@importFrom stats rnorm
#'
#'@keywords internal

detection_delay <- function(NJ, NK, rhoj, rhok, k, jump, kink, max_n = 1e7){

  #generate historical data under the null and initialize detector
  data <- rnorm(k)
  detector <- FLOC_init(data, NJ, NK, rhoj, rhok)

  while(detector$iteration <= max_n){

    #iterate the detector with observation after change
    data <- rnorm(1) + jump + (detector$iteration + 1) * kink
    detector <- FLOC_iter(data, detector)

    if(detector$detect_jump | detector$detect_kink){

      #return the iteration of the detection
      return(detector$iteration)

    }
  }

  #No change detected in max_n observations after the change
  return(Inf)
}

#'segmented linear function
#'
#'A function that returns a segmented linear function at times x
#'
#'@param x points to evaluate the function at; vector of integers
#'@param tau changepoint; integer
#'@param alpha_minus intercept at tau of first segment; numeric
#'@param beta_minus slope on first segment; numeric
#'@param alpha_plus intercept at tau of second segment; numeric
#'@param beta_plus slope on second segment; numeric
#'
#'@return vector that describes a segmented linear function
#'
#'@export
seg_lin_fun <- function(x,
                        tau,
                        alpha_minus,
                        beta_minus,
                        alpha_plus,
                        beta_plus) {

  #first segment
  f1 <- (x[x <= tau] - tau) * beta_minus + alpha_minus

  #second segment
  f2 <- (x[x > tau] - tau) * beta_plus + alpha_plus

  return(c(f1,f2))
}



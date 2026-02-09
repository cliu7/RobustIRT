#' Convert category threshold probabilities (\eqn{P^*}) to the probabilities of responding in each category
#'
#' Calculate \eqn{P}, the probabilities of responding in each category, from the \eqn{P^*} threshold values using the graded response model (GRM; Samejima, 1969).
#' The probability that a subject responds in or above a category \eqn{k} for item \eqn{j} is \eqn{P^*_{jk}(\theta) = \frac{1}{1+ e^{-a_j (\theta-b_{jk})}}},
#' for \eqn{K} categories and \eqn{K-1} threshold parameters  (\eqn{b_{j,1}, ..., b_{j,K-1}}), where \eqn{b_{j,k}} separates response category \eqn{k} and \eqn{k+1} (\eqn{k=1,..K-1}) (Embretson & Reise, 2000).
#' \eqn{a_j} is the item discrimination parameter.  The probability of endorsing exactly category \eqn{k} is \eqn{P_{jk}(\theta) = P^*_{j,k}(\theta) - P^*_{j,k+1}(\theta),} where \eqn{P^*_{j1}(\theta) \equiv 1.0} and \eqn{P^*_{jK}(\theta) \equiv 0.0.}
#' @references Embretson, S. E., & Reise, S. P. (2000). \emph{Item response theory for psychologists.} Mahwah, N.J: L. Erlbaum Associates.
#' @references Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. \emph{Psychometrika Monograph Supplement, 34} (4, Pt. 2), 100–100.
#' @param Pstar A \eqn{J \times K-1 \times N} array of \eqn{P^*} threshold probability values, for \eqn{K} categories, \eqn{J} items, and \eqn{N} subjects
#' @return The probabilities \eqn{P} of responding in each category
#' @examples
#' One-subject case
#' Pstar <- matrix(c(0.85, 0.50, 0.20,0.70, 0.40, 0.10,0.90, 0.60, 0.30), nrow = 3, byrow = TRUE)
#' pstar_to_p(Pstar)
#' 
#' Multi-subject case
#' J <- 2   # items
#' K <- 4   # categories (0–3)
#' N <- 3   # persons
#' Simulate P* values
#' Pstar <- array(runif(J * (K - 1) * N), dim = c(J, K - 1, N))
#' Pstar <- aperm(apply(Pstar, c(1, 3), sort, decreasing = TRUE), c(2, 1, 3))
#' pstar_to_p(Pstar)

pstar_to_p<-function(Pstar){
  stopifnot(all(Pstar >= 0 & Pstar <= 1))
  
  if(length(dim(Pstar))==3){
    # If there is more than one subject
    
    stopifnot(all(apply(Pstar, c(1, 3), function(x) all(diff(x) <= 0))))
    
    J <- dim(Pstar)[1] # number of items
    thresh <- dim(Pstar)[2] # number of thresholds
    N <- dim(Pstar)[3] # number of subjects
    K <- thresh + 1 # number of categories
    
    # Initialize array for probabilities
    P.in <- array(NA, dim = c(J, K + 1, N))
    
    # Bind P(X>1)=1.0 and P(X>K+1)=0.0
    P.in[, 1, ] <- 1
    P.in[, 2:K, ] <- Pstar
    P.in[, K + 1, ] <- 0
    
    # Difference along category dimension
    P <- P.in[, 1:K, ] - P.in[, 2:(K + 1), ]
    
  }else{ # If there is only one subject
    
    stopifnot(is.matrix(Pstar))
    stopifnot(all(apply(Pstar, 1, function(x) all(diff(x) <= 0))))
    
    J <- nrow(Pstar) # number of items
    thresh <- dim(Pstar)[2] # number of thresholds
    K <- thresh + 1 # number of categories
    
    # Initialize array for probabilities
    P.in <- cbind(1, Pstar, 0)
    
    # Difference along category dimension
    P <- P.in[, 1:K] - P.in[, 2:(K + 1)]
  }
  
  # Return the matrix/array of category response probabilities
  return(P)
}

#' Item Response Probability

#' Computes item response probabilities for select IRT models (1PL, 2PL, MIRT, GRM, and MGRM), given ability and item paratemeters.
#' by constructing the appropriate linear predictors and applying the logistic function. Returns item response probabilities for dichotomous data or item category response probabilities for polytomous data.
#' @references 
#' @references 
#' @param theta A numeric vector or matrix of latent trait values. 
#' @param ipars A matrix of item parameters. See examples for how to structure the columns of the matrix based on the model utilitized.
#' @param model A character string specifying which IRT model to use.
#' @param D A positive scaling constant used in the logistic function. Defaults to 1.7.
#' @return For model accomodating dichotomous data ("1PL", "2PL", "MIRT"), returns an \eqn{N \times J} matrix of response probabilities P(X = 1)
#' @return For models accomodating polytomous data (GRM, MGRM), returns a list with: {pstar}: an array of cumulative probabilities \eqn{P^*(X \geq k)} 
#'         and {P}: an array of category probabilities \eqn{P(X = k)}
#' @examples

#' 1PL case
#' N <- 4 # subjects
#' J <- 2 # items
#' theta <- rnorm(N) # generate ability values
#' ipars <- rnorm(J) # generate item difficulties
#' item.prob(theta, "1PL", ipars)

#' 2PL case
#' N <- 5 # subjects
#' J <- 3 # items
#' theta <- rnorm(N)
#' ipars <- cbind(a = runif(J, 0, 1), # generate item discrimination
#'                b = rnorm(J)) # generate item difficulty
#' item.prob(theta, "2PL", ipars) 
 
#' MIRT case
#' N <- 7 # subjects
#' J <- 2 # items
#' L <- 4 # dimensions
#' theta <- matrix(rnorm(N * L), ncol = L) # N x L ability matrix
#' ipars <- cbind(matrix(runif(J * L, 0, 1), ncol = L), d = rnorm(J)) # generate slopes + intercept
#' item.prob(theta, "MIRT", ipars)

#' GRM case
#' N <- 5 # subjects
#' J <- 3 # items
#' K <- 4 # categories
#' theta <- rnorm(N)
#' a <- runif(J, 0, 1) # J discriminations
#' b <- matrix(sort(rnorm(J * (K - 1))), nrow = J) # J x (K-1) thresholds
#' ipars <- cbind(a, b) # Combine into a J x K matrix
#' item.prob(theta, "GRM", ipars)
 
#' MGRM case
#' N <- 5 # subjects
#' J <- 3 # items
#' l <- 2 # dimensions
#' K <- 4 # categories
#' theta <- matrix(rnorm(N * L), ncol = L) 
#' a <- matrix(runif(J * L, 0, 1), ncol = L) # slopes 
#' b <- matrix(rnorm(J * (K - 1)), nrow = J) # thresholds 
#' ipars <- cbind(a, b)
#' item.prob(theta, "MGRM", ipars)

#' Can handle 1PL, 2PL, MIRT, GRM, MGRM
item.prob<-function(theta, model, ipars, D=1.7){
  model<-toupper(model)
  
  # Ensure theta is a matrix
  if(is.vector(theta)){
    Theta <- matrix(theta, ncol = 1)
  }else{
    Theta <- as.matrix(theta)
  }
  
  #extract model parameters
  N <- nrow(Theta) # number of subjects
  L <- ncol(Theta) # number of dimensions
  J<-nrow(ipars) # number of items
  
  # Logistic function with scaling constant D = 1.7
  invlogit <-function(x) 1/(1+exp(-D*x))
  
  
  # Compute linear predictor (ex) depending on the model
  
  # 1PL predictor: theta - b
  # sapply loops over each person's theta value (Theta[,1]) 
  # For each x = Theta[n,1], compute x - ipars (vector of item difficulties) 
  # sapply returns J x N, then t() makes it N x J
   if(model=="1PL"){
    ex<-t(sapply(Theta[,1], function(x) x-ipars))
  }
  
  # 2PL predictor: a * theta - b
  if(model=="2PL" | model=="Rasch"){
    ex<-t(sapply(Theta[,1], function(x) ipars[,1]*(x-ipars[,2]))) # ipars[,1] = a_j (discrimination), ipars[,2] = b_j (difficulty)
  }
  
  # MIRT predictor: a·theta + d
  if(model=="MIRT"){
    a<-ipars[,1:L] # a: J x L matrix of discriminations
    d<-ipars[,L+1] # d: J-length vector of intercepts
    ex<-(t(a%*%t(Theta))+d)
  }
  
  # GRM predictor: a * (theta - b_k)
  if(model == "GRM"){
    a<-as.vector(ipars[,1]) # a: J-length vector of discriminations
    b<-ipars[,-1] # b: J x K matrix of thresholds (one column per threshold)
    thresh<-ncol(b) # number of thresholds K

    # for each Theta[n,1] = x, compute J x K matrix a_j * (x - b_{j,k}) 
    ex <- vapply(Theta[,1], function(x) {a*(x-b)}, matrix(0, nrow = J, ncol = thresh))
  }
  
  # MGRM predictor: (a·theta) - (sum(a_j) * b_jk)
  if(model=="MGRM"){
    
    a<- ipars[, 1:L] # a: J x L matrix (first L columns: discrimination parameters)
    b<- ipars[, (L+1):ncol(ipars)] # b: J x K matrix (category threshold parameters)
    thresh<-ncol(b) # number of thresholds K
    
    # int: J x N matrix of a·theta, then transposed to J x N
    int <- t(Theta %*% t(a)) 
    a.sum <- rowSums(a) 
    
    # for each k, compute J x N matrix: int - a.sum * b[,k]
    ex <- vapply(1:thresh, function(x) {int - a.sum*b[,x]}, matrix(0, nrow = J, ncol = N))
    ex <- aperm(ex, c(1, 3, 2)) # reorders dimensions from J x N x K to J x K x N
  }
  
  # Apply logistic function and return probabilities for dichotomous models
  if(model %in% c("1PL", "2PL", "MIRT")){
    return(P=invlogit(ex)) #Returns N x J matrix of P(X = 1)
  }
  
  # Apply logistic function and return probabilities for polytomous models
  if(model %in% c("GRM", "MGRM")){
    pstar<-invlogit(ex) # cumulative probabilities P*(X >= k)
    return(list(pstar=pstar, P=pstar_to_p(pstar))) #converts cumulative probabilities to category probabilities P(X = k)
  }
  
}


#' Bisquare Weighting Function
#'
#' Calculate Tukey's bisquare weight (Mosteller & Tukey, 1977) given a residual and bisquare tuning parameter
#' @param r A residual that measures the inconsistency of a response from the subject's assumed response model, on one item. Residuals that are NA are given a weight of 0.
#' @param B Bisquare tuning parameter. Larger values lead to less downweighting
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' @return Bisquare weight value
bisquare<-function(r, B){
  w<-ifelse(is.nan(r), 0, 
            ifelse(abs(r) <= B, (1-(r/B)^2)^2, 0))
  
  return(w)
}

#' Huber Weighting Function
#'
#' Calculate the Huber weight (Huber, 1981) given a residual and Huber tuning parameter
#' @param r A residual that measures the inconsistency of a response from the subject's assumed response model, on one item. Residuals that are NA are given a weight of 0.
#' @param H Huber tuning parameter. Larger values lead to less downweighting.
#' @references Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250
#' @return Huber weight value
huber<-function(r, H){
  w<-ifelse(is.nan(r), 0, 
            ifelse(abs(r) <= H, 1, H/abs(r)))
  return(w)
}

#' Data generation for dichotomous and Likert outcomes based on response probabilities
#' 
#' Bernoulli sampling is conducted for dichotomous items, while multinomial sampling is conducted for polytomous items.
#' @param P Probabilities of a correct response (in the case of dichotomous outcomes) or an array of item category response probabilities. For polytomous data, the number of columns should be equal to the number of response categories. Array structure of P with three dimensions assumes polytomous data generation.
#' @param anchor Value of the lowest category value. Typical values are 0 or 1; default is 0.
#' @param polytomous Is polytomous data desired? Must be specified when P is a matrix of category response probabilities for polytomous data. 
#' 
dat.gen<-function(P, anchor = 0, polytomous = FALSE, seed=NULL){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(!is.null(dim(P)[3])){
    # If dealing with array of polytomous category probabilities
    out <- t(apply(P, 1, function(p) sample(1:length(p), size = 1, prob = p))) - (1-anchor)
    
  }else if(polytomous==T){
    # If dealing with matrix of polytomous category probabilities for one subject
    out <- t(apply(P, c(1, 3), function(p) sample(1:length(p), size = 1, prob = p))) - (1-anchor)
    
  }else if(polytomous==F){
    # If dealing with matrix or vector of item success probabilities on dichotomous items
    U<-matrix(runif(length(P)), ncol = ncol(P), nrow = nrow(P))
    out <- ifelse(P>U, 1, 0)
  }
  return(out)
}

#' Standard error function for 2PL, MIRT, GRM, MGRM
#' 
#' Computes per person: information SE (expected Fisher), asymptotic SE (robust expected, incorporating weights), sandwich SE, Bayesian SE & bayesian sandwich SE (2PL & GRM only for Bayesian)
standard.errors<-function(theta, ipars, dat, model, D=1.7, weight.type = "equal", tuning.par = NULL, bayes = NULL){
  
  # Function to determine weights:
  weight.func <- function(r, weight.type2="equal", tuning.par2=NULL){
    if(weight.type2 == "equal") return(rep(1, length(r)))
    if(weight.type2 == "Huber") return(huber(r, tuning.par))
    if(weight.type2 == "bisquare") return(bisquare(r, tuning.par))
    if(length(weight.type2) == length(r)) return(weight.type2)
    stop("Invalid weight.type")
  }
  
  P<-item.prob(theta, model, ipars, D=D)
  r<-(dat-P)/sqrt(P*(1-P))
  
}
#' Calculate standard errors of ability estimates for MIRT data
#' 
#' Calculate the standard errors of ability estimates using the Fisher Information matrix for multidimensional dichotomous data using the MIRT model
#' @param theta Vector of latent traits (abilities) for an individual
#' @param d Vector of difficulty parameters for the items
#' @param a Matrix of discrimination parameters for the items (rows are items, columns are dimensions)
#' @param D Scaling constant. Default is 1.7 to scale the item parameters to align with a normal ogive model. 1.0 may be used alternatively.)
#' @return Vector of standard errors
#' @export
#'
#' @examples
#' #Define the parameters
#' theta <- c(0.5, -1.2) # Example vector of latent traits
#' d <- c(0.0, -0.5, 1.0) # Example vector of difficulty parameters for 3 items
#' a <- matrix(c(1.0, 0.5,
#'               0.7, 1.2,
#'               1.5, 1.0), nrow = 3, byrow = TRUE) # Example matrix of discrimination parameters for 3 items
#' dat <- c(1,1,0)
#' D <- 1.7 # Scaling constant

#' #Calculate the standard errors
#' std.errs <- std.err.dichotomous(theta, d, a, D)
#' print(std.errs)

std.err.dichotomous<- function(theta, d, a, dat, D = 1.7, residual = "standardized", weight.type="equal", tuning.par=NULL){
  
  # Number of items and dimensions
  n <- nrow(a)
  dim <- ncol(a)
  
  # Calculate probabilities
  probs<-item.prob(theta, "2PL", cbind(a, d))
  
  # Calculate expected value of response
  expected.value<-probs
  
  # Calculate standardized residual # check that this is true
  if(residual == "standardized"){
    resid<- (dat-expected.value)/sqrt(probs*(1-probs))
  }else if(residual=="information"){
    resid<- a%*%matrix(theta)+d
  }else{
    stop(paste(residual, "is not a valid residual type."))
  }
  
  weighting.term <- NULL
  if (weight.type == "bisquare") { #Bisquare weighting
    weighting.term <- bisquare(resid, tuning.par)
  } else if (weight.type == "Huber") { # Huber weighting
    weighting.term <- huber(resid, tuning.par)
  } else if (weight.type == "equal") { # Regular Maximum likelihood estimation
    weighting.term <- rep(1, n)
  } else{ # User-specified weights
    weighting.term <- weight.type
  }
  if (is.null(weighting.term)) {
    stop("Cannot determine the weighting function.")
  }
  
  # Initialize  matrices
  info.mat2.num<-info.mat2.den<-info.mat <- A<-B<-matrix(0, nrow = dim, ncol = dim)
  
  # Loop over each item
  for (j in 1:n) {
    
    Hess<-matrix(ncol = dim, nrow = dim) # Initialize Hessian matrix
    D1<-matrix(nrow=dim) # Initialize First derivative matrix
    
    # Calculate the first derivative (D1) and the Hessian matrix
    for(k in 1:dim){
      D1[k] <- D*weighting.term[j]*a[j,k]*(dat[j]-probs[j,]) # First derivative
      for(l in 1:dim){
        Hess[k,l] <- ((D)^2)*weighting.term[j]*a[j,k]*a[j,l]*(1-probs[j,])*probs[j,] # Hessian matrix of 2nd derivatives
      }
    }
    
    # Sandwich SE: Update item contributions to A and B matrices
    B<-B+ D1 %*% t(D1)
    A<-A-Hess
    
    # Calculate unweighted item information
    I_j <- D^2*probs[j,] *(1-probs[j,]) *(a[j,] %*% t(a[j,])) # need to figure out how to add the weighting term into this
    
    # Weighted SE: Update item contribution to Fisher Information (Also denominator for Magis SE)
    info.mat <- info.mat + weighting.term[j] *I_j
    
    # Magis SE: Update numerator
    info.mat2.num<-info.mat2.num + weighting.term[j]^2 * I_j
  }
  
  # Check if Information matrix is singular
  if(any(is.na(info.mat))|any(is.nan(info.mat))){
    standard_errors <- magis_ase<- NA
  }else if(is.singular.matrix(info.mat)){ 
    standard_errors <- magis_ase<- "singular"
  }else {
    # Calculate the inverse of the Fisher Information matrix
    inv_info_mat <- if(length(info.mat)>1) solve(info.mat) else 1/info.mat
    
    # Weighted SE: square roots of diagonal elements of inverse Fisher Info matrix
    standard_errors <- sqrt(diag(inv_info_mat))
    
    magis_ase<- sqrt(diag(solve(info.mat)%*%info.mat2.num%*%solve(info.mat)))
  }
  
  # Check if A matrix is singular
  if(any(is.na(A))|any(is.nan(A))){
    sandwich_se<- NA
  }else if(is.singular.matrix(A)){ 
    sandwich_se<- "singular"
  }else {
    # Calculate sandwich SE
    sandwich_se <- sqrt(diag(solve(A)%*%B%*%solve(A)))
  }
  
  return(list(standard_errors=standard_errors, ase_magis = magis_ase, sandwich_se=sandwich_se))
}

#' Calculate standard errors of ability estimates for GRM data
#' 
#' 
#' Calculate the standard errors of ability estimates using the Fisher Information matrix or Huber-White sandwich standard errors for Likert-type response data using the GRM model
#' @param theta Vector of latent traits (abilities) for an individual
#' @param dat Vector of observed item responses for the individual, of length $J$ for $J$ items
#' @param b Matrix of threshold parameters for the items, with number of rows corresponding to the number of items and number of columns corresponding to the number of thresholds
#' @param a Vector of discrimination parameters for the items of length $J$
#' @param weight.type Type of weighting function to be used: "equal", "Huber", or "bisquare"
#' @param tuning.par Tuning parameter to be used for the Huber or bisquare weight function
#' @param D Scaling constant. Default is 1.7 to scale the item parameters to align with a normal ogive model. 1.0 may be used alternatively.)
#' @return Standard Error Standard error of $\theta$ based on the Fisher information (expected information) matrix
#' @return Sandwich SE Huber-White sandwich estimator of the standard error
#' @export
#'
#' @examples
std.err.poly<- function(theta, dat, b, a, weight.type, tuning.par, D = 1.7){
  
  # Takes 1 subject and is unidimensional (theta is a scalar)
  
  # Number of items
  j<-length(a) 
  # Number of thresholds between categories
  nthresh<-length(b)/j
  
  # Calculate probabilities
  probs<-item.prob(theta, "GRM", cbind(a, b))
  
  # Initialize the Fisher Information (scalar)
  #    Fisher Information (Dodd, De Ayala, and Koch, 1995)
  I_test<-B<-A<-0
  
  pstars<-cbind(rep(1, j), probs$pstar[,,1], rep(0, j))
  # above lowest threshold gets probability 1; above highest threshold gets probability 0
  
  expected.value<-probs$P[,,1]%*%matrix(c(1:(nthresh+1))) 
  # Calculate standardized residual
  residual<-(dat-expected.value)/sqrt(rowSums(apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*probs$P[,,1]))  
  
  # Loop over each item
  for (i in 1:j) {
    
    # Discrimination vector for the j-th item
    a_j <- a[i]
    
    # Compute item response weights based on specified weight function (bisquare, Huber, equal)
    weighting.term <- NULL
    if (weight.type == "bisquare") {
      weighting.term <- bisquare(residual[i], tuning.par)
    } else if (weight.type == "Huber") {
      weighting.term <- huber(residual[i], tuning.par)
    } else {
      weighting.term <- 1
    }
    # Check if weighting term is determined
    if (is.null(weighting.term)) {
      stop("Cannot determine the weighting function.")
    }
    
    D2<-D1<-D2E<-0
    
    # Looping over each category
    for(k in 1:(nthresh+1)){
      ps0<-pstars[i,k]
      qs0<-1-ps0
      ps1<-pstars[i,(k+1)]
      qs1<-1-ps1
      P<-probs$P[i,k,1]
      
      #Add contribution of kth category to the information
      # First and second derivatives of the log-likelihood
      # Observed Information
      D1<-D1+ D*a_j*weighting.term*(ps0*qs0-ps1*qs1)/P
      D2<-D2+ D^2*a_j^2*weighting.term*( (ps0*qs0*(qs0-ps0)-ps1*qs1*(qs1-ps1))/P - (ps0*qs0-ps1*qs1)^2/P^2 )
      
      # Expected Information (D2*P)
      D2E<-D2E+ D^2*a_j^2*weighting.term*( (ps0*qs0*(qs0-ps0)-ps1*qs1*(qs1-ps1)) - (ps0*qs0-ps1*qs1)^2/P )
      
    }
    
    # Calculate the information of the jth item (Dodd, DeAyala, & Koch, 1995; Embretson & Reise 2001)
    I_j <- -D2E
    
    # Add the jth contribution to the B matrix
    B<- B + D1^2
    
    # Add the jth contribution to the A matrix
    A<- A - D2
    
    # Add the contribution to the test information
    I_test <- I_test + I_j
  }
  
  # The standard errors are the square roots of the inverse Fisher Information
  standard_error <- sqrt(I_test^(-1))
  
  sandwich_se <- sqrt(A^(-1)*B*A^(-1))
  
  return(list("Standard Error" = standard_error, "Sandwich SE" = sandwich_se))
}


#' Ability Estimation Function Using Robust Estimation - Dichotomous Data
#'
#' Calculate robust ability estimates using the MIRT item response function with the given weight function, fixed item parameters, and item responses
#' @param dat A \eqn{J \times N} matrix of dichotomously-coded data where 1 represents an endorsed response and 0 represents a non-endorsed response. \emph{J} is the test length and \emph{N} is the number of subjects.
#' @param a A \eqn{L \times J} matrix containing fixed item slope parameters for \emph{J} items and \emph{L} dimensions
#' @param d A vector of length \emph{J} containing fixed intercept parameters for J items. If the model is unidimensional, this vector should contain difficulty parameters corresponding to \emph{b} in the unidimensional 2PL model.
#' @param iter Maximum number of iterations for the Newton-Raphson method. Default is 30.
#' @param cutoff Threshold value to terminate the iteration when the likelihood changes below this value, which means that the estimation is converged. Default is 0.01.
#' @param init.val A vector of length \emph{L} containing initial latent trait values for the maximum likelihood estimation. The default is 0 for each of the \emph{L} dimensions.
#' @param weight.category The weighting strategy to use: "equal", "bisquare", or "Huber". Default is "equal", which is equally weighted as in standard maximum likelihood estimation.
#' @param tuning.par The tuning parameter for "bisquare" or "Huber" weighting functions. Greater tuning parameters result in less downweighting.
#' @param residual Type of residual if using bisquare or Huber weight functions. Default is "information" while "standardized" can alternatively be used.
#' @details The goal of robust estimation is to downweigh potentially aberrant responses to lessen their impact on the estimation of \eqn{\boldsymbol{\theta}}. Robust estimates resist the harmful effects of response disturbances and tend to be less biased estimates of true ability than maximum likelihood estimates.
#'               Using the multidimensional IRT model for dichotomous data (McKinley & Reckase, 1983), the probability of a correct response (e.g., "1") on item \emph{j} is given by the formula \eqn{P(X_{ij}=1 | \boldsymbol{\theta}_i) =  \frac{1}{1+ e^{-1.7(\boldsymbol{a}_j\boldsymbol\theta_i+d_j)}},} where \eqn{\boldsymbol{a}_j} is a \eqn{L \times J} matrix of item slope parameters for \emph{J} items and \emph{L} dimensions. The intercept parameter is \eqn{d_j}.  
#`               The 2PL IRT model for unidimensional data subsumes the MIRT model and can be used for unidimensional robust estimation (when \emph{L}=1), using the formula \eqn{P(X_{ij} =1 | \theta_i) = \frac{1}{1+ e^{-1.7a_{j}(\theta_i-b_{j})}}} where \eqn{a_j} is the same as the multidimensional case where \emph{L=1}. The parameter \eqn{b_j} differs from MIRT and can be interpreted as the item difficulty.
#`               The contribution of item \emph{j} to the overall log-likelihood for one subject is weighted with a weight \eqn{\omega(r_j)} as a function of a residual \eqn{r_{ij}} for the item \emph{j} and person \emph{i}:
#`               \deqn{\sum_j^J \omega(r_{ij}) \frac{\partial}{\partial\boldsymbol\theta_i} \ln L(\boldsymbol\theta_j;x_{ij}) = 0 }
#`               The residual, which measures the inconsistency of a response from the subject's assumed response model, is \deqn{r_{ij} = \textbf a_j\boldsymbol\theta_i + d_j } in the multidimensional case.
#`               Two types of weight functions are used: Tukey's bisquare weighting function (Mosteller & Tukey, 1977)
#`               \deqn{\omega(r_{ij})=\begin{cases}[1-(r_{ij}/B)^2]^2, & \text{if} |r_{ij}|\leq B.\\0, & \text{if} |r_{ij}|>B.\end{cases}}
#`               and the Huber weighting function (Huber, 1981)
#`               \deqn{\omega(r_{ij})=\begin{cases}1, & \text{if} |r_{ij}|\leq H.\\H/|r_{ij}|, & \text{if} |r_{ij}|>H.\end{cases}}
#`               Both functions are effective in estimating more accurate scores with aberrant data, although the bisquare weight function may lead to nonconvergence when using data containing a high proportion of incorrect responses (Schuster & Yuan, 2011).
#' @references Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250
#' @references McKinley, R. L., & Reckase, M. D. (1983, August). \emph{An Extension of the Two-Parameter Logistic Model to the Multidimensional Latent Space} (Research No. ONR83-2). Iowa City, IA: American College Testing Program.
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' @references Schuster, C., & Yuan, K.-H. (2011). Robust Estimation of Latent Ability in Item Response Models. \emph{Journal of Educational and Behavioral Statistics}, 36(6), 720–735. https://doi.org/10.3102/1076998610396890
#' @return theta A \eqn{N \times L} matrix of ability estimates for \emph{N} subjects and \emph{L} dimensions
#' @return standard.errors A \eqn{N \times L} matrix of the standard errors of ability estimates for \emph{N} subjects and \emph{L} dimensions, calculated by the square root of the reciprocal of the Fisher information. NAs replace nonconverging values. 
#' @return convergence \eqn{N \times L} matrix containing indicators of convergence for \emph{N} subjects: a “0” indicates the value converged; a “1” indicates the maximum likelihood estimation did not converge to any value; and “Singular” indicates the Hessian matrix was singular and could not be used to continue the maximum likelihood estimation.
#' @return theta.progression A \eqn{L \times p \times N} array for \emph{L} dimensions, \emph{p} number of iterations supplied to the input, and \emph{N} subjects. Each column provides the updated theta estimate at each iteration of the Newton-Raphson algorithm until the change in log-likelihood for that subject reaches the cutoff value, infinite values (nonconverged), or encounters a singular matrix error.
#' @return residual A \eqn{J \times N} matrix with residuals corresponding to the ability estimate for \emph{N} subjects respective to the \emph{J} test items
#' @export
#'
#' @examples
#' # Test length
#' n <- 30
#'
#' ## Number of iterations of Newton's method
#' iters<-50
#'
#' ## Number of dimensions
#' dim <- 3
#'
#' ## Correlation between one person's thetas
#' cor <- .6
#'
#' ## Covariance matrix for generating thetas
#' sigma <- matrix(cor, ncol = dim, nrow = dim)
#' diag(sigma) <- 1
#' s <- rep(1, dim)
#' sigma <- s*sigma*s
#'
#' ## Generate real thetas
#' thetas <- matrix(mvrnorm(n = 10, rep(0,dim), Sigma = sigma), ncol = dim)
#'
#' ## Generate slope parameters
#' a <- matrix(runif(n*dim, 1.5, 3), nrow = n, ncol = dim)
#' ## Generate intercept parameters
#' d <- matrix(rnorm(n), ncol = 1)
#' ## Calculate probabilities
#' probs <- item.prob(thetas, "MIRT", cbind(a,d))
#' ## Generate data from probabilities
#' dat <- apply(probs, c(1,2), function(x) rbinom(1,1,x))
#'
#' ## Estimate thetas
#' theta.est(dat, a, d, iters, init.val=matrix(rep(0,dim)), weight.type="Huber", tuning.par=1)

theta.est<-function(dat, a, d, D=1.7, iter=30, cutoff=.01, init.val=rep(0,ncol(a)), weight.type="equal", tuning.par=NULL, residual = "information"){
  # Check if the turning parameter is given when the weight.type is not "equal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The turning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
    if (is.null(residual)) {
      stop(paste("The residual cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  dim<-ncol(a) #number of dimensions
  n<-ncol(dat) #number of items
  N<-nrow(dat) #number of subjects
  
  # Initialize matrices to store results
  theta.out<-standard.errs<-robust.se<-magis.ase<-matrix(nrow = N, ncol = dim) # Estimated thetas, standard errors and residuals
  residual.out<-matrix(nrow = N, ncol = n)
  convergence<-matrix(0, nrow=N, ncol = dim) # Convergence status
  theta.progression<-array(NA, dim = c(dim, iter, N)) # Progression of thetas over iterations
  
  for (i in 1:N){ # Loop over subjects if more than one
    # Initialize theta values
    if(length(init.val)>dim){
      theta<-matrix(init.val[i,], nrow = ncol(a), byrow = T)
    }else{
      theta<-matrix(init.val, nrow = ncol(a), byrow = T)
    }
    P0<-0 # Starting probability for calculating likelihood for convergence
    for(j in 1:iter){ # Loop over max number of iterations for Newton-Raphson Algorithm
      if(dim==1){
        theta<-as.vector(theta)
        ex<-a*(theta-d) # Calculate exponent for probability function if unidimensional
      }else{
        ex<-a%*%theta+d # Calculate exponent for multiple dimensions
      } 
      P<-1/(1+exp(-D*ex)) # Calculate probability
      
      Hess<-matrix(ncol = dim, nrow = dim) # Initialize Hessian matrix
      D1<-matrix(nrow=dim) # Initialize First derivative matrix
      
      # Calculate residual based on specified type
      if(!is.null(residual)){
        if(residual=="information"){
          resid.r<-ex
        }else if(residual=="standardized"){
          resid.r<-(dat[i,]-P)/sqrt(P*(1-P))
        }
      }
      
      
      # Calculate weight based on specified type
      weighting.term <- NULL
      if (weight.type == "bisquare") { #Bisquare weighting
        weighting.term <- bisquare(resid.r, tuning.par)
      } else if (weight.type == "Huber") { # Huber weighting
        weighting.term <- huber(resid.r, tuning.par)
      } else { # Regular Maximum likelihood estimation
        weighting.term <- 1
      }
      if (is.null(weighting.term)) {
        stop("Cannot determine the weighting function.")
      }
      # Caluclae the first derivative (D1) and the Hessian matrix
      for(k in 1:dim){
        D1[k] <- D*sum(weighting.term*a[,k]*(dat[i,]-P)) # First derivative
        for(l in 1:dim){
          Hess[k,l] <- (-(D)^2)*sum(weighting.term*a[,k]*a[,l]*(1-P)*(P)) # Hessian matrix of 2nd derivatives
        }
      }
      # Check if Hessian matrix is singular
      
      if(any(is.na(Hess))|any(is.nan(Hess))){
        convergence[i,]<-rep(NA, dim)
        theta<-NA
        break
      }else if(is.singular.matrix(Hess)){ 
        convergence[i,]<-rep("Singular", dim)
        for(h in 1:dim){ # Bound theta estimates between -3 and 3
          if(theta[h]>3){
            theta[h]<-NA
          }else if(theta[h]<(-3)){
            theta[h]<-NA
          }
        }
        break
      }
      # Check for nonconverging theta updates
      if(any(is.na(theta-solve(Hess)%*%D1))){
        theta.out[i,]<-theta<-theta-solve(Hess)%*%D1
        convergence[i,]<-as.numeric(is.na(theta.out[i,]))
        break # Break if nonconverging
      }
      # Update theta estimates 
      theta<-theta-solve(Hess)%*%D1 
      
      # Compare log-likelihood with log-likelihood from previous iteration for convergence criterion
      log_like<-sum(log(P))-sum(log(P0)) 
      if(!is.nan(log_like) & abs(log_like)<cutoff){
        break
      }
      P0<-P # Update previous probability
      theta.progression[,j,i]<-theta # Store theta from this iteration
    }
    # If theta never converged to tolerance, set theta to NA
    if(j==iter){
      theta.out[i,]<-theta<-rep(NA, dim)
      convergence[i,]<-rep(1, dim)
    }
    # Bound theta estimates between -3 and 3
    theta<-ifelse(abs(theta)>3, sign(theta)*3, theta)
    
    # Store final theta estimates
    theta.out[i,]<-theta
    residual.out[i,]<-resid.r
    # subset the hessian if one or more theta is na
    if(all(is.na(theta))){
      standard.errs[i,]<-robust.se[i,]<-magis.ase[i,]<-NA
    }else{
      std.err<-std.err.dichotomous(as.matrix(theta[!is.na(theta)]), d, as.matrix(a[,!is.na(theta)]), dat[i,], D, residual, weight.type, tuning.par)
      standard.errs[i,!is.na(theta)]<-std.err$standard_errors
      robust.se[i,!is.na(theta)]<-std.err$sandwich_se
      #magis.ase[i,!is.na(theta)]<-std.err$magis_ase
    }
    
  }
  return(list(theta = theta.out, standard.errors=standard.errs, robust.se=robust.se, magis.se = magis.ase, convergence = convergence, theta.progression = theta.progression, residual=residual.out))
}

#' Ability Estimation Function Using Robust Estimation (GRM)
#'
#' Calculate robust ability estimates using the GRM item response function with the given weight function, fixed item parameters, and item responses
#' @param dat A \eqn{J \times N} matrix of polytomously-scored data (e.g., Likert-type) for \emph{J} items and \emph{N} subjects.
#' @param a Vector of slope parameters for \emph{J} items
#' @param b A \eqn{J \times (K-1)} matrix of category threshold parameters for \emph{K} categories
#' @param iter Max number of iterations. Default is 100
#' @param cutoff Threshold value to terminate the iteration when the likelihood changes below this value, which means that the estimation is converged.
#' @param init.val Vector of initial latent trait for the maximum likelihood estimation for \emph{N} subjects. If a single value is provided, that initial value will be used for all subjects. Default is 0.
#' @param weight.category The weighting strategy to use: "equal", "bisquare" and "Huber". Default is "equal", which is equally weighted as in standard maximum likelihood estimation.
#' @param tuning.par The tuning parameter for "bisquare" or "Huber" weighting functions. Greater tuning parameters result in less downweighting in robust estimation.
#' @details The goal of robust estimation is to downweigh potentially aberrant responses to lessen their impact on the estimation of \eqn{\theta_i}. Robust estimates resist the harmful effects of response disturbances and tend to be less biased estimates of true ability than maximum likelihood estimates.
#'               Under the graded response model (GRM; Samejima, 1969), the probability that a subject responds in or above a category \emph{k} for item \emph{j} is \eqn{P^*_{jk}(\theta_i) = \frac{1}{1+ e^{-a_j (\theta_i-b_{jk})}}}  (Embretson & Reise, 2000). \eqn{a_j} is the item discrimination parameter. There are \emph{K} categories and \eqn{K-1} threshold parameters (\eqn{b_{j,1}, ..., b_{j,K-1}}), where \eqn{b_{j,k}} separates response category \eqn{k} and \eqn{k+1} (\eqn{k=1,..K-1}).
#'               The probability of endorsing exactly category \eqn{k} is therefore: \eqn{P_{jk}(\theta_i) = P^*_{j,k}(\theta_i) - P^*_{j,k+1}(\theta_i),} where \eqn{P^*_{j1}(\theta_i) \equiv 1.0} and \eqn{P^*_{jK}(\theta_i) \equiv 0.0.}
#'               The contribution of item \emph{j} to the overall log-likelihood for one subject is weighted with a weight \eqn{\omega(r_{ij})} as a function of a residual \eqn{r_{ij}} for the item:
#'               \deqn{\sum^J_{j=1} \omega(r_{ij}) \sum^K_{k=1} u_{jk}\text{log}P_{jk} = 0 }
#'               \eqn{u_{jk}} is an indicator function: \deqn{u_{jk} = \begin{cases}
#'                                                            1 & \text{if } X_{ij} = k; \\
#'                                                            0 & \text{otherwise}.
#'                                                            \end{cases} }
#'               The residual, which measures the inconsistency of a response from the subject's assumed response model, is \deqn{r_{ij} = \frac{1}{\sigma_{X_{ij}}}\left[X_{ij} - E(X_{ij}|\hat{\theta}_i)\right]} for the GRM.
#'               The difference in fit is determined between the observed response \eqn{X_{ij}} and expected score \eqn{E(X_{ij}|\hat{\theta}_i) = \sum_{k=1}^KkP_{jk}(\hat{\theta}_i)}, and scaled by the variance \eqn{\sigma_{X_{ij}}^2 = \sum_{k=1}^K (X_{ijk}-E[X_{ij}|\hat{\theta}_i])^2P_{jk}(\hat{\theta}_i).}
#'               Two types of weight functions are used: Tukey's bisquare weighting function (Mosteller & Tukey, 1977)
#'                 \deqn{\omega(r_{ij})=\begin{cases}[1-(r_{ij}/B)^2]^2, & \text{if} |r_{ij}|\leq B.\\0, & \text{if} |r_{ij}|>B.\end{cases}}
#'               and the Huber weighting function (Huber, 1981)
#'                 \deqn{\omega(r_{ij})=\begin{cases}1, & \text{if} |r_{ij}|\leq H.\\H/|r_{ij}|, & \text{if} |r_{ij}|>H.\end{cases}}
#'               Both functions are effective in estimating more accurate scores with aberrant data, although the bisquare weight function may lead to nonconvergence when using data containing a high proportion of incorrect responses (Schuster & Yuan, 2011).
#' @references Embretson, S. E., & Reise, S. P. (2000). \emph{Item response theory for psychologists.} Mahwah, N.J: L. Erlbaum Associates.
#' @references Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' @references Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. \emph{Psychometrika Monograph Supplement, 34} (4, Pt. 2), 100–100.
#' @references Schuster, C., & Yuan, K.-H. (2011). Robust Estimation of Latent Ability in Item Response Models. \emph{Journal of Educational and Behavioral Statistics}, 36(6), 720–735. https://doi.org/10.3102/1076998610396890
#' @return theta Ability estimates for \emph{N} subjects. NAs replace values that did not converge to any value. Estimates that converged to values less than -3.0 were replaced with -3.0, while estimates that converged to values greater than 3.0 were replaced with 3.0.
#' @return convergence Indicators of convergence for \emph{N} subjects: a “0” indicates the value converged, while a “1” indicates the maximum likelihood estimation did not converge to any value.
#' @return standard.error Standard errors of the theta estimates for \emph{N} subjects, given by the square root of the reciprocal of the Fisher information. NAs replace nonconverging values. 
#' @return theta.progression A matrix with rows corresponding to each subject and columns corresponding to the number of iterations supplied to the input. Each column provides the updated theta estimate at each iteration of the Newton-Raphson algorithm until the change in log-likelihood for that subject reaches the cutoff value or the value is nonconverged (reaches infinite values).
#' @return residual A \eqn{J \times N \times p} array containing residuals corresponding to the ability estimate for \emph{N} subjects respective to the \emph{J} test items at each iteration until convergence within maximum \emph{p} iterations, nonconvergence, or singular matrix is reached.
#' @export
#' @examples
#' # Test Length
#' n<-30
#' 
#' # Number of thresholds (5-point Likert scale)
#' nthresh<-4
#' 
#' # Number of iterations of Newton's method
#' iter <- 15
#' 
#' # Set critical value for convergence criteria
#' crit.val<-0.01
#' 
#' # Set real thetas - 5 subjects
#' thetas<-c(-2,-1,0,1,2)
#' 
#' # Set item slope
#' a<-runif(n, .90, 2.15)
#' 
#' # Set threshold parameters
#' b<-t(apply(matrix(runif(n*4, -2.5,2.5), nrow = n, ncol =4), 1, sort))
#' 
#' # Calculate Probabilities
#' probs<-item.prob(thetas, "GRM", cbind(a, b))
#' 
#' # Generate Likert data
#' dat<-data.gen(probs$P)
#' 
#' # Make the data aberrant by reverse coding 20% items
#' ab.prop<-0.2
#' index<-sample(c(1:n), ab.prop*n)
#' ab.dat<-dat
#' ab.dat[index, ]<-apply(matrix(dat[index,]), c(1,2), function(x) return(nthresh+2-x))
#'
#' 
#' # Calculate MLE (Non-robust)
#' mle<-theta.est.grm(ab.dat, a, b, iter, crit.val, init.val=0, weight.type="equal")
#' 
#' # Use MLE as starting value, or 0 if NA
#' start.val<-apply(mle$theta, c(1,2), function(x) ifelse(is.na(x), 0, x))
#' 
#' # Calculate bisquare- and Huber-weighted robust estimates
#' b.est<- theta.est.grm(ab.dat, a, b, iter, crit.val, init.val=start.val, weight.type="bisquare", tuning.par=4)
#' h.est<-theta.est.grm(ab.dat, a, b, iter, crit.val, init.val=start.val, weight.type="Huber", tuning.par=1)
#' 
#' # Compare robust ability estimates with MLE
#' b.est$theta
#' h.est$theta
#' mle$theta
#' 
theta.est.grm<-function(dat, a, b, iter=30, cutoff=0.01, init.val=0, weight.type="equal", tuning.par=NULL){
  # Check if the turning parameter is given when the weight.type is not "equal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The turning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  # Number of subjects
  l<-ncol(dat) 
  # Test length
  n<-nrow(dat) 
  # Number of threshold parameters
  nthresh<-ncol(b) 
  
  # Array to store threshold values of interest for each response
  dat.b<-array(dim = list(n, l, 2))
  
  # Arrays/matrices to store final values of theta estimates, standard errors, nonconvergence rates, and the residuals and thetas at each iteration
  theta.est2<- standard.error<- sandwich.sd<-matrix(data=NA, nrow=l)
  convergence<-matrix(0, nrow=l)
  theta.progression<-matrix(NA, nrow = l, ncol = iter)
  residual<-matrix(data=NA, nrow = n, ncol = l)
  
  # Select threshold values (b) based on response categories
  # b_{jk} and b_{j,k+1} are used to calculate the probability of the response k on item j
  for (i in 1:l){
    for(j in 1:n){
      if(dat[j,i]==1){
        # For responses in lowest category
        dat.b[j,i,1]<- -100000 # very small b_{jk} will give probability 1
        dat.b[j,i,2]<-b[j,1]
      }else if(dat[j,i]>nthresh){
        # For responses in highest category
        dat.b[j,i,1]<-b[j,dat[j,i]-1]
        dat.b[j,i,2]<-100000 # very large b_{j,k+1} will give probability 0
      }else{
        # For responses in middle categories
        dat.b[j,i,1]<-b[j,dat[j,i]-1]
        dat.b[j,i,2]<-b[j,dat[j,i]]
      }
    }
  }
  
  # Loop to estimate theta for each subject 
  for(i in 1:l){
    # Initialize theta value
    theta<-ifelse(length(init.val)>1, init.val[i], init.val)
    P0<-0
    
    # Iterative loop for maximum likelihood estimation of theta
    for (k in 1:iter){
      
      # Compute item response probability for the response category and next category
      exponent_0<-a*(theta-dat.b[,i,1])
      exponent_1<-a*(theta-dat.b[,i,2])
      # Calculate P*_k and P*_{k+1}
      ps0<-1/(1+exp(-1.7*exponent_0))
      ps1<-1/(1+exp(-1.7*exponent_1))
      qs0<-1-ps0
      qs1<-1-ps1
      # Compute P_k, the probability of response k
      P<-ps0-ps1
      # Calculate expected response
      pstar<-1/(1+exp(-1.7*(c(a)*(theta-b))))
      probs<-pstar_to_p(pstar)
      expected.value<-probs%*%matrix(c(1:(nthresh+1))) 
      # Calculate standardized residual
      residual[,i]<-(dat[,i]-expected.value)/sqrt(rowSums(apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*probs))  
      # Compute item response weights based on specified weight function (bisquare, Huber, equal)
      weighting.term <- NULL
      if (weight.type == "bisquare") {
        weighting.term <- bisquare(residual[,i], tuning.par)
      } else if (weight.type == "Huber") {
        weighting.term <- huber(residual[,i], tuning.par)
      } else {
        weighting.term <- 1
      }
      # Check if weighting term is determined
      if (is.null(weighting.term)) {
        stop("Cannot determine the weighting function.")
      }
      
      # First and second derivatives of the log-likelihood
      D1<-sum(1.7*a*weighting.term*(ps0*qs0-ps1*qs1)/P) 
      D2<-sum(1.7^2*a^2*weighting.term*( (ps0*qs0*(qs0-ps0)-ps1*qs1*(qs1-ps1))/P - (ps0*qs0-ps1*qs1)^2/P^2 )) 
      
      # Check for NAs & record nonconvergence
      if(is.na(theta-D1/D2)){
        theta.est2[i]<-theta<-NA
        convergence[i,1]<-1
        break
      }
      
      # Update and store theta for this iteration
      theta<-theta.progression[i,k]<-theta-D1/D2
      
      # Stop Newton-Raphson method if log-likelihood difference is converged / less than cutoff
      log_like<-sum(log(P))-sum(log(P0))
      if(abs(log_like)<cutoff){
        break
      }
      
      # Update probability for loglikelihood comparison
      P0<-P
    }
    
    # Store final theta estimate for subject
    theta.est2[i]<-theta
    
    # Compute standard errors via test information
    if(is.na(theta)){
      standard.error[i]<-sandwich.sd[i]<-NA
    } else{
      se.s<-std.err.poly(theta, dat[,i], b, a, weight.type, tuning.par, D = 1.7)
      standard.error[i]<-se.s$"Standard Error"
      sandwich.sd[i]<-se.s$"Sandwich SE"
    }
    
    # Handle cases where theta didn't converge withing number of iterations
    if((k==iter&(abs(log_like)>cutoff))){
      theta.est2[i]<-standard.error[i]<-sandwich.sd[i]<-NA
      convergence[i,1]<-1
    }else if(!is.na(theta) & abs(theta)>3){ 
      # Replace thetas that converged outside [-3, 3]
      theta<-sign(theta)*3
      se.s<-std.err.poly(theta, dat[,i], b, a, weight.type, tuning.par, D = 1.7)
      standard.error[i]<-se.s$"Standard Error"
      sandwich.sd[i]<-se.s$"Sandwich SE"
      
      theta.est2[i]<-theta
    }
  }
  return(list(theta = theta.est2, convergence=convergence, standard.error = standard.error, robust.se = sandwich.sd, theta.progression = theta.progression, residual=residual))
}

#' Ability Estimation Function Using Robust Estimation (MGRM)
#'
#' Calculate robust ability estimates using the MGRM item response function with the given weight function, fixed item parameters, and item responses
#' @param dat A \eqn{N \times J} matrix of polytomously-scored data (e.g., Likert-type) for \emph{N} subjects and \emph{J} items. Indexing begins at 0.
#' @param a A \eqn{J \times L} matrix of fixed slope parameters for \emph{J} items and \emph{L} dimensions
#' @param b A \eqn{J \times (K-1)} matrix of category threshold parameters for \emph{J} items and \emph{K} categories
#' @param D Constant to scale the normal ogive model. 
#' @param weight.category The weighting strategy to use: "equal", "bisquare" and "Huber". Default is "equal", which is equally weighted as in standard maximum likelihood estimation.
#' @param tuning.par The tuning parameter for "bisquare" or "Huber" weighting functions. Greater tuning parameters result in less downweighting in robust estimation.
#' @param init.val Initial latent traits for the estimation. Accepts a \eqn{N \times L} matrix for \emph{N} subjects and \emph{L} dimensions, a vector of length \emph{L} to be used for all subjects, or a single value be used for all subjects and dimensions. Default is 0.
#' @param iter Max number of iterations for the Newton-Rapshon algorithm. Default is 30
#' @param cutoff Threshold value to terminate the Newton-Rapshon algorithm when the likelihood change is below this value, to deem the estimation converged. Default is 0.01.

theta.est.mgrm<-function(dat, a, b, D=1.7, weight.type="equal", tuning.par=NULL, init.val=0, iter=30, cutoff=0.01){
  
  # Check if the turning parameter is given when the weight.type is not "equal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The turning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  
  # Number of subjects
  N<-nrow(dat) 
  # Test length
  J<-ncol(dat) 
  # Number of threshold parameters
  nthresh<-ncol(b) 
  # Number of dimensions
  dim<-ncol(a)
  
  # Initialize Matrices
  theta.out<-matrix(ncol=dim, nrow=N)
  residual.out<-matrix(nrow=N, ncol=J)
  convergence<-matrix(0, nrow=N, ncol = dim) # Convergence status
  theta.progression<-array(NA, dim = c(iter, dim, N)) # Progression of thetas over interations
  
  for(i in 1:N){ # loop over each subject
    
    # Initialize theta
    if(is.matrix(init.val)){
      if(nrow(init.val>1)){
        theta<-as.matrix(init.val[i,])
      } else{
        theta<-as.matrix(init.val)
      }
    }else{
      if(length(init.val)==dim){
        theta<-as.matrix(init.val)
      }else{
        theta<-as.matrix(rep(init.val, dim))
      }
    }
    xi<-dat[i,]
      
    P0<-0
    for(p in 1:iter){
      probs<-item.prob(t(as.matrix(theta)), "MGRM", cbind(a, b))
      
      P<-probs$P[,,1]
      if(J==1){ # if test length only 1
        pstars.all<-cbind(rep(1,J), matrix(probs$pstar[,,1], ncol=nthresh), rep(0, J))
        
        # Calculate standardized residual
        exp.x<-P%*%matrix(c(1:(nthresh+1)))
        sd.x<-sqrt(sum(sapply(matrix(1:(nthresh+1)),  function(x) x-exp.x)^2*P))  
        
      } else {
        
        pstars.all<-cbind(rep(1,J), probs$pstar[,,1], rep(0, J))
        
        # Calculate standardized residual
        exp.x<-P%*%matrix(c(1:(nthresh+1)))
        sd.x<-sqrt(rowSums(apply(matrix(1:(nthresh+1)), 1, function(x) x-exp.x)^2*P))  
      }
      resid.r<-(xi-exp.x)/sd.x
        
        # Calculate weight based on specified type
        weighting.term <- NULL
        if (weight.type == "bisquare") { #Bisquare weighting
          weighting.term <- bisquare(resid.r, tuning.par)
        } else if (weight.type == "Huber") { # Huber weighting
          weighting.term <- huber(resid.r, tuning.par)
        } else { # Regular Maximum likelihood estimation
          weighting.term <- 1
        }
        if (is.null(weighting.term)) {
          stop("Cannot determine the weighting function.")
        }
        
        # Extract probabilities based on responses
        ps0<-pstars.all[cbind(1:J,xi)]
        ps1<-pstars.all[cbind(1:J,xi+1)]
        qs0<-1-ps0
        qs1<-1-ps1
        P1<-P[cbind(1:J, xi)]
        
        # Initialize matrices for 1st and 2nd derivatives
        D1<-matrix(nrow=dim)
        Hess<-matrix(ncol=dim, nrow=dim)
        
        # Calculate 1st and 2nd derivatives
        for(h in 1:dim){
          D1[h,]<-sum(D*a[,h]*weighting.term*(ps0*qs0-ps1*qs1)/P1)
          for(l in 1:dim){
            Hess[h,l]<- sum( D^2*a[,h]*a[,l]*weighting.term*( (ps0*qs0*(1-2*ps0)-ps1*qs1*(1-2*ps1))/P1
                                                              -(ps0*qs0-ps1*qs1)^2/P1^2 ) )
          }
        }
        if(any(is.na(Hess))|any(is.nan(Hess))){
          theta.out[i,]<-theta<-NA
          convergence[i,]<-1
          break
        }
        # Check if Hessian matrix is singular
        if(is.singular.matrix(Hess)){ 
          #convergence[i,]<-rep("Singular", dim)
          Hess <- Hess + diag(1e-6, dim) # add regularization term to hessian
          if(is.singular.matrix(Hess)){
            theta.out[i,]<-NA
            convergence[i,]<-"singular"
            break
          }
        }
        
        theta<-theta-solve(Hess)%*%D1
        
        # Stop Newton-Raphson method if log-likelihood difference is converged / less than cutoff
        log_like<-sum(log(P1))-sum(log(P0))
        if(is.nan(log_like)|abs(log_like)<cutoff){
          break
        }
        
        # Update probability for log-likelihood comparison
        P0<-P1
        theta.progression[p,,i]<-theta # Store theta from this iteration
      
    }
    theta<-ifelse(abs(theta)>3, sign(theta)*3, theta)
    theta.progression[p,,i]<-theta.out[i,]<-theta
    if(!is.null(resid.r)){
      residual.out[i,]<-resid.r
    }
  }
  return(list(thetas=theta.out, residual=residual.out, theta.progression=theta.progression))
}

#' Standard error for MGRM ability estimation

std.err.mgrm<- function(theta, d, a, dat, D = 1.7, weight.type="equal", tuning.par=NULL){
  #takes dat with values of k=0,1,..,nthresh
  # Number of items and dimensions
  J <- nrow(a)
  dim <- ncol(a)
  nthresh<-ncol(d)
  
  # Calculate probabilities
  probs<-item.prob(theta, "MGRM", cbind(a, d))
  
  P1<-probs$P[,,1]
  pstars.all<-cbind(rep(1,J), probs$pstar[,,1], rep(0, J))
  
  exp.x<-P1%*%matrix(c(0:(nthresh)))
  sd.x<-sqrt(rowSums(apply(matrix(0:(nthresh)), 1, function(x) x-exp.x)^2*P1))  
  resid<-(dat-exp.x)/sd.x
  
  # Initialize the Fisher Information matrix
  info.mat <- A<-B<-matrix(0, nrow = dim, ncol = dim)
  
  # Loop over each item
  for (i in 1:J) {
    
    # Discrimination vector for the j-th item
    a_j <- a[i,]
    
    # Compute item response weights based on specified weight function (bisquare, Huber, equal)
    weighting.term <- NULL
    if (weight.type == "bisquare") {
      weighting.term <- bisquare(resid[i], tuning.par)
    } else if (weight.type == "Huber") {
      weighting.term <- huber(resid[i], tuning.par)
    } else {
      weighting.term <- 1
    }
    # Check if weighting term is determined
    if (is.null(weighting.term)) {
      stop("Cannot determine the weighting function.")
    }
    
    Hess<-D2E<-matrix(0,ncol = dim, nrow = dim) # Initialize Hessian matrix
    D1<-D12<-matrix(0,nrow=dim) # Initialize First derivative matrix
    
    # Looping over each category
    for(k in 1:(nthresh+1)){
      ps0<-pstars.all[i,k]
      qs0<-1-ps0
      ps1<-pstars.all[i,(k+1)]
      qs1<-1-ps1
      P<-P1[i,k]
      
      #Add contribution of kth category to the information
      # First and second derivatives of the log-likelihood
      # Calculate 1st and 2nd derivatives & information
      for(h in 1:dim){
        # Observed Information
        D1[h,]<-D1[h,]+D*a_j[h]*weighting.term*(ps0*qs0-ps1*qs1)/P
        D12[h,]<-D12[h,]+D1[h,]^2/P
        for(l in 1:dim){
          # Expected Information (D2*P)
          D2E[h,l]<-D2E[h,l]+ D^2*a_j[h]*a_j[l]*weighting.term*( (ps0*qs0*(qs0-ps0)-ps1*qs1*(qs1-ps1)) - (ps0*qs0-ps1*qs1)^2/P )
          
          Hess[h,l]<-Hess[h,l]+ D^2*a_j[h]*a_j[l]*weighting.term*( (ps0*qs0*(1-2*ps0)-ps1*qs1*(1-2*ps1))/P
                                                                   -(ps0*qs0-ps1*qs1)^2/P^2 ) 
        }
      }
    }
    
    # Calculate the information of the jth item (Dodd, DeAyala, & Koch, 1995; Embretson & Reise 2001)
    I_j <- -D2E
    
    # Add the jth contribution to the B matrix
    B<- B + D1 %*% t(D1)
    
    # Add the jth contribution to the A matrix
    A<- A - Hess
    
    # Add the contribution to the test information
    info.mat <- info.mat + I_j
  }
  
  
  # Check if Info  matrix is singular
  if(any(is.na(info.mat))|any(is.nan(info.mat))){
    standard_errors <- sandwich_se<- NA
  }else if(is.singular.matrix(info.mat)){ 
    standard_errors <- sandwich_se<- "singular"
  }else{
    # Calculate the inverse of the Fisher Information matrix
    inv_info_mat <- if(length(info.mat)>1) solve(info.mat) else 1/info.mat
    
    # The standard errors are the square roots of the diagonal elements of the inverse Fisher Information matrix
    standard_errors <- sqrt(diag(inv_info_mat))
    
    sandwich_se <- sqrt(diag(solve(A)%*%B%*%solve(A)))
  }
  
  return(list(standard_errors=standard_errors, sandwich_se=sandwich_se))
}

#' Modified Standardized Residual
#'
#' Calculate the modified standardized residual (MSR) used to quantify misfit between a response and the IRT model (Yu & Cheng, 2019): 
#' \deqn{ \frac{y_j - E(Y_j | \hat{\theta})}{P(y_j|\hat{\theta})}}
#' given data and parameters that follow the multidimensional graded response model (MGRM; Muraki & Carlson, 1995). Handles unidimensional and dichotomous data as well.
#' @param dat \eqn{N \times J} matrix of data, beginning indexing at 1 for \eqn{N} respondents and \eqn{J} items
#' @param theta \eqn{N \times L} matrix of \eqn{\theta}s, for \eqn{L} dimensions
#' @param a \eqn{J \times L} matrix of discrimination parameters
#' @param b \eqn{J \times (K-1)} matrix of category threshold parameters, for \eqn{K} Likert-type categories
#' @references  Yu, X & Cheng, Y. A change-point analysis procedure based on weighted residuals to detect back random responding. \emph{Psychological Methods} (Oct. 2019), pp. 658–674. DOI: 10.1037/met0000212.
#' @references Muraki, E. & Carlson, J. E. Full-Information Factor Analysis for Polytomous Item Responses. \emph{Applied Psychological Measurement} (Mar. 1995), pp. 73–90. DOI: 10.1177/014662169501900109.
#' @return \eqn{N \times J} matrix of residuals

msr<-function(dat, theta, a, b){
  N<-nrow(dat) # respondents
  J<-ncol(dat) # items
  nthresh<-ncol(b) # number of category thresholds
  residual<-matrix(nrow=N, ncol=J)
  for(i in 1:N){
    # Calculate expected response
    probs<-item.prob(theta[i,], "MGRM", cbind(a, b))$P[,,1]
    expected.value<-probs%*%matrix(c(1:(nthresh+1))) 
    # Extract the probabilities corresponding to the response
    if(J==1){
      probs.i <- probs[ dat[i,]]
    } else {
      probs.i <- probs[cbind(1:J, dat[i,])]
    }
    # Calculate MSR residual
    residual[i,]<-(dat[i,]-expected.value)/probs.i
    
  }
  
  return(residual)
}

#' Change-Point Analysis for Back Random Responding
#'
#'  Uses the modified standardized residual (MSR) to identify shifts in response behavior from normal responding to random responding. Details: This change-point analysis, as outlined in Yu and Cheng, 2019, is designed to detect respondents who exhibits back-random responding (BRR; Clark et al., 2003), which is defined as switching from a normal response behavior to random responding for the remaining test items. For each possible change point, the MSR, 
#'  #' \deqn{ \frac{y_j - E(Y_j | \hat{\theta})}{P(y_j|\hat{\theta})}}
#'   of items beforehe difference in the average absolute weighted residual (ABWR) of items before the change point is taken from that of items after the change point. A subject who exhibits BRR is expected to have large absolute residuals, indicating greater misfit, after the change point, as compared to smaller residuals before the change point. An \eqn{R_j} greater than the user-specified critical value indicates that the subject is likely exhibiting BRR. 
#' given data and parameters that follow the multidimensional graded response model (MGRM; Muraki & Carlson, 1995). Handles unidimensional and dichotomous data as well.
#' @param dat \eqn{N \times J} matrix of data, beginning indexing at 1 for \eqn{N} respondents and \eqn{J} items
#' @param a \eqn{J \times L} matrix of discrimination parameters
#' @param b \eqn{J \times (K-1)} matrix of category threshold parameters, for \eqn{K} Likert-type categories
#' @param crit.val Critical value to separate normal responders from those exhibiting BRR; an \eqn{R_j} greater than this value is classified as exhibiting BRR
#' @references  Yu, X & Cheng, Y. A change-point analysis procedure based on weighted residuals to detect back random responding. \emph{Psychological Methods} (Oct. 2019), pp. 658–674. DOI: 10.1037/met0000212.
#' @references Clark, M. E., Gironda, R. J., & Young, R. W. (2003). Detection of back random responding: Effectiveness of MMPI-2 and Personality Assessment Inventory validity indices. \emph{Psychological Assessment}, 15 (2), 223–234. DOI: 10.1037/1040-3590.15.2.223
#' @return flagged vector of binary indicators: 1 if classified as aberrant, 0 if not.
#' @return change.point the item number with the largest \eqn{R_j} for each respondent. For response vectors flagged, BRR begins at the next item.  
#' @return max.residual the largest \eqn{R_j} for each respondent
#' @export
#' @examples
#' data(BFI)
#' BFI<-BFI[,grep("t1_bfi_N", colnames(BFI))]
#' 
#' # Estimate model parameters using ltm
#' mod<-grm(BFI) 
#' ipars<-matrix(unlist(mod$coefficients), nrow = ncol(BFI), byrow = T)
#' a<-ipars[,5]
#' b<-ipars[,1:4]
#' colnames(BFI)<-rownames(BFI)<-NULL
#' 
#' # reparameterize b to align with MGRM 
#' cpa.brr(BFI, a, b/a, crit.val=80)


cpa.brr<-function(dat, a, b, crit.val=75){
  N<-nrow(dat) # respondents
  J<-ncol(dat) # items
  
  # Initialize matrix of residuals
  out.rj<-matrix(NA, N, J-1)
  for(j in 1:(J-1)){
    # Handle cases where only one item is used
    if(j==1){b1<-matrix(b[1:j,], nrow=1)} else {b1<-b[1:j,]}
    if(j==(J-1)){b2<-matrix(b[(j+1):J,], nrow=1)} else {b2<-b[(j+1):J,]}
    
    # Calculate theta with data up to item j 
    thetas<-theta.est.mgrm(as.matrix(dat[,1:j]), as.matrix(a[1:j]), b1)$thetas 
    # Calculate residuals before and after item j
    resids.1<-msr(as.matrix(dat[,1:j]), thetas, as.matrix(a[1:j]), b1)
    resids.2<-msr(as.matrix(dat[,(j+1):J]), thetas, as.matrix(a[(j+1):J]), b2)
    # Calculate difference in mean
    out.rj[, j-1]<-rowMeans(abs(resids.2))-rowMeans(abs(resids.1))
  }
  
  max.r<-unlist(apply(out.rj, 1, which.max), use.names=F)
  row_max <- apply(out.rj, 1, function(x) max(x, na.rm=T))
  flagged<-rep(0, N)
  flagged[which(row_max>crit.val)]<-1
  return(list("flagged"=flagged, "change.point"=max.r, "max.residual"=row_max))
}


#' Plot histogram of residuals along plot of weight (dependent on TuCo) vs residuals
#'
#' Plot a histogram of residuals along the graph of the weighting function (dependent on the tuning parameter) as a function of the residual
#' @param r A vector of residuals
#' @param H Huber tuning parameter
#' @param B Bisquare tuning parameter
#' @details The goal of this plot is to visualize the proportion of residuals that are downweighted based on the tuning parameter and allow the researcher to choose a tuning parameter that suits their data well.
#'               For a set of residuals with larger variance, a larger tuning parameter should be used.
#'               Generally, the tail end of the weighting function should approach the tail end of the distribution of residuals.
#'               To increase the downweighting applied in estimation, use a smaller tuning parameter. To decrease the amount of downweighting, use a greater tuning parameter.
#'               The function will plot the histogram of residuals below (1) the Huber weight curve (Huber, 1981) if \emph{H} is supplied to the function, (2) Tukey's bisquare weight curve (Mosteller & Tukey, 1977) if \emph{B} is supplied, or (3) both the Huber and bisquare weight curves if both tuning parameters are supplied.
#'               If \emph{H} is supplied, vertical lines will be displayed at \emph{H} and \emph{-H} to highlight the amount of data that is downweighted (a residual greater than \emph{|H|}) versus not downweighted.
#'               If no tuning parameter is supplied, just the histogram of residuals is generated.
#' @references Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' @return Histogram plot of residuals beneath a graph of the weight functions vs. the residuals
#' @export
#' @examples
#' ## Unidimensional IRT Example
#' n=40
#' # Generate real thetas
#' thetas<-matrix(seq(0,2, by=.05), ncol=1)
#' # Set item slope and difficulty
#' a<-matrix(runif(n, .5, 1.5), ncol=1) 
#' b<-rnorm(n)
#'
#' # Introduce response disturbances: working at a suboptimal level (theta minus 1 standard deviation), for last 40% of items
#' theta.drop<-1
#' chng.pt<-0.6
#' probs<-cbind(item.prob(thetas, "2PL", cbind(a[1:(chng.pt*n)], b[1:(chng.pt*n)])), 
#'              item.prob(thetas-theta.drop, "2PL", cbind(a[(chng.pt*n+1):n], b[(chng.pt*n+1):n])))
#' dat<-apply(probs, c(1, 2), function(x) rbinom(1, 1, x))
#' 
#' Estimate thetas
#' example<-theta.est(dat, a, d, iter=30, cutoff=.01, init.val=rep(0,ncol(a)), weight.type="equal", tuning.par=NULL)
#' choose.tuco(r=matrix(na.omit(example$residual), ncol=1), B=4)
#'
#' ## GRM example
#' n=40
#' nthresh<-4
#' # Generate real thetas
#' thetas<-seq(-2,2.1, by=.1)
#'
#' # Set item slope
#' a<-runif(n, .90, 2.15) 
#' # Set category threshold parameters
#' b<- matrix(runif(n*nthresh, -2.5,2.5), nrow = n, ncol =nthresh)
#' b<-t(apply(b, 1, sort)) 
#' 
#' # Calculate response probabilities and generate data
#' probs<-item.prob(thetas, "GRM", cbind(a, b))
#' dat<-data.gen(probs$P)
#' 
#' Introduce response disturbance: random guessing for latter 40% of the exam
#' abdat<-dat
#' chng.pt<-.6 
#' abdat[(chng.pt*n+1):n, ]<-sample(c(1:(nthresh+1)), length(thetas)*(n-chng.pt*n), replace = T)
#' Calculate ability estimates and residuals
#' mle<-theta.est.grm(dat, a, b, iter=30, cutoff=0.01, init.val=0, weight.type="equal")
#' choose.tuco(matrix(mle$residual), H=.1, B=.8)
#'
#' ## MIRT Example
#' data(SAT12)
#' SAT12[SAT12 == 8] <- NA #set 8 as a missing value
#'
#' # Correct answer key
#' key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
#' scoredSAT12 <- key2binary(SAT12, key)
#' specific <- c(2, 3, 2, 3, 3, 2, 1, 2, 1, 1, 1, 3, 1, 3, 1, 2, 1, 1, 3, 3, 1, 1, 3, 1, 3, 3, 1, 3, 2, 3, 1,2) #which factor each item loads on
#' b_mod1 <- mirt(scoredSAT12, specific)
#' ipars<-matrix(unlist(coef(b_mod1))[1:(32*6)], nrow = length(key), byrow=T) #item parameters
#'
#' ## Set Parameters
#' a <- ipars[,1:3]
#' d<- ipars[,4]
#' # Remove vectors with missing data
#' dat<-scoredSAT12[!is.na(rowSums(scoredSAT12)),] 
#' colnames(dat)<-NULL
#'
#' # Calculate theta estimates and residuals
#' out<-theta.est(t(dat), a, d, iter=30, cutoff=.01, weight.type="equal")
#' choose.tuco(matrix(out$residual[,,2]), H=1, B=4)

choose.tuco<-function(r, H=NULL, B=NULL){
  # r is a vector of residuals
  
  residuals<-data.frame(Residual =r)
  hist.out<-ggplot(residuals, aes(x=Residual))+geom_histogram(aes(y = ..density..), bins=50)+ ylab("Density")
  if(!is.null(H) & !is.null(B)){
    hist.out<-hist.out+ geom_vline(xintercept = -H, linetype="dashed", color = "grey")+ 
      geom_vline(xintercept = H, linetype="dashed", color = "grey")
    weight.out<-ggplot()+stat_function(fun=function(x) huber(x, H), aes(colour = "Huber"))+
      stat_function(fun=function(x) bisquare(x, B),  aes(colour = "Bisquare"))+ 
      geom_vline(xintercept = -H, linetype="dashed", color = "grey")+ 
      geom_vline(xintercept = H, linetype="dashed", color = "grey")+
      xlim(min(r, na.rm =T), max(r, na.rm =T)) + ylab("Weight")+
      scale_color_manual(name = "Function", breaks=c('Bisquare', 'Huber'), values=c('Bisquare'="darkcyan", 'Huber'='firebrick')) +
      theme(legend.position = c(.9, .74))+
      ggtitle("Weights Applied in Estimation")
    return(do.call(ggarrange, c(list(weight.out+xlab(NULL), hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  }else if(is.null(H) & !is.null(B)){
    weight.out<-ggplot()+stat_function(fun=function(x) bisquare(x, B), aes(colour = "Bisquare"))+
      xlim(min(r, na.rm =T), max(r, na.rm =T)) + ylab("Weight")+
      scale_color_manual(name = "Function", breaks=c('Bisquare'), values=c('Bisquare'="darkcyan")) +
      theme(legend.position = c(.9, .74))+
      ggtitle("Weights Applied in Estimation")
    return(do.call(ggarrange, c(list(weight.out+xlab(NULL), hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  }else if(is.null(B) & !is.null(H)){
    hist.out<-hist.out+ geom_vline(xintercept = -H, linetype="dashed", color = "grey")+ 
      geom_vline(xintercept = H, linetype="dashed", color = "grey")
    weight.out<-ggplot()+stat_function(fun=function(x) huber(x, H), aes(colour = "Huber"))+
      xlim(min(r, na.rm =T), max(r, na.rm =T)) + ylab("Weight")+
      geom_vline(xintercept = -H, linetype="dashed", color = "grey")+ 
      geom_vline(xintercept = H, linetype="dashed", color = "grey") + 
      scale_color_manual(name = "Function", breaks=c('Huber'), values=c('Huber'='firebrick')) +
      theme(legend.position = c(.9, .74))+
      ggtitle("Weights Applied in Estimation")
    return(do.call(ggarrange, c(list(weight.out+xlab(NULL), hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  }else{
    return(hist.out+ggtitle("Histogram of Residuals"))
  }
}

#' Plot to compare robust estimates with MLE
#'
#' Generate a scatterplot of robust estimates versus the maximum likelihood estimate (MLE)
#' @param dat \eqn{J \times N} matrix of response data for \emph{J} items and \emph{N} subjects
#' @param a \eqn{J \times L} matrix of slope parameters for \emph{J} items and \emph{L} dimensions (\emph{L=1} if using the GRM or unidimensional 2PL model)
#' @param b If type = “GRM”, an \eqn{J \times (K-1)} matrix of intercept parameters
#' @param d If type = “MIRT”, a vector of discrimination parameters for \emph{J} items
#' @param iter Maximum number of iterations. Default is 30
#' @param cutoff Threshold value to terminate the iteration when the likelihood changes below this value, which means that the estimation is converged. Default is 0.01.
#' @param H Huber tuning parameter
#' @param B Bisquare tuning parameter
#' @param same.plot.dim If TRUE and type = “MIRT”, estimates across all \emph{L} dimensions will be plotted on the same graph. If FALSE (default) and type = “MIRT”, one plot per dimension will be generated.
#' @param same.plot If TRUE (default) and both \emph{H} and \emph{B} are supplied, generates both the Huber and bisquare plots in the same image frame. If FALSE, the Huber and bisquare plots are generated on separate images.
#' @param type Type of data: "Dichotomous" for dichotomous data (multidimensional or unidimensional) or "GRM" for Likert-type data
#' @details When the data is not disturbed, robust estimates should not differ greatly from the maximum likelihood estimate (MLE).
#'                                       By plotting the robust estimates against the MLE, the user can identify possible aberrant trends, if the robust estimates are far from the MLE, as indicated by the distance from the \eqn{y=x} identity line.
#'                                       Larger discrepancies between the point plotted for a subject and the identity line suggest there may be some disturbance in this subject’s data that the robust estimation may be correcting.
#'                                       At least one tuning parameter \emph{H} or \emph{B} must be supplied to the function; if both are supplied, the function will return separate plots for both weighting systems.

#' @return ‘Summary Statistics (Huber)’ If \emph{H} is supplied, a dataframe where each row provides a subject’s ID, MLE, the Huber-weighted robust estimate, the minimum distance between the point on the plot and the identity line, and their response vector. The subjects are organized by greatest to least distance.
#' @return ‘Summary Statistics (Bisquare)’ If \emph{B} is supplied, a dataframe where each row provides a subject’s ID, MLE, the bisquare-weighted robust estimate, the minimum distance between the point on the plot and the identity line, and their response vector. The subjects are organized by greatest to least distance.
#' @return Plots If same.plot = TRUE and both \emph{H} and \emph{B} are supplied, each robust ability estimate is plotted against the MLE; graphs for each of the Huber- and bisquare-weighted estimates are generated separately but on the same image frame. The identity line \eqn{y=x} is plotted as a reference line.
#' @return `Huber Plot` If same.plot = FALSE or \emph{B} is not supplied, each Huber-weighted robust ability estimate is plotted against the MLE with the identity line \eqn{y=x} as reference.
#' @return `Bisquare Plot` If same.plot = FALSE or \emph{H} is not supplied, each bisquare-weighted robust ability estimate is plotted against the MLE with the identity line \eqn{y=x} as reference.
#' @export
#' @examples
#' ## Test length
#' n <- 30
#'
#' ## Number of iterations of newton's method
#' iter <- 15
#'
#' ## Number of thresholds
#' nthresh <- 4
#'
#' ## Generate real thetas
#' thetas <- seq(-2, 2, by=.1)
#'
#' ## Generate item slope
#' a <- runif(n, .90, 2.15)
#'
#' ## Generate category threshold parameters
#' b <- matrix(runif(n*nthresh, -2.5,2.5), nrow = n, ncol =nthresh)
#' b <- t(apply(b, 1, sort))
#'
#' ## Calculate probabilities
#' probs <- item.prob(thetas, "GRM", cbind(a, b))
#'
#' ## Generate input data from probabilities
#' abdat <- data.gen(probs$P)
#'
#' ## Introduce aberrant responses: random guessing for latter 40% of the exam
#' chng.pt <- .6
#' abdat[(chng.pt*n+1):n, ] <- sample(c(1:(nthresh+1)), length(thetas)*(n-chng.pt*n), replace = T)
#'
#' ## Plot the GRM
#' out<-theta_plots(abdat, a, b=b, iter=30, cutoff=0.01, H=.1, B=1, same.plot = F, type="GRM")
#' ## Check bisquare plot
#' out$`Bisquare Plot`
#' ## Check Huber summary
#' out$`Summary Statistics (Huber)`
theta_plots<-function(dat, a, d=NULL, b=NULL, iter=30, cutoff=0.01, H=NULL, B=NULL, same.plot.dim = F, same.plot = T, type){
  if(type != "Dichotomous" & type != "GRM"){
    return(print("Please enter a valid type of model (e.g., 'Dichotomous' or 'GRM')."))
  }else if(type == "Dichotomous"){
    if(is.null(d)){ return(print("Please enter a vector of intercept values (d)."))}
    h.plots<-b.plots<-list()
    dim<-ncol(a)
    if(!same.plot.dim){
      theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "equal")$theta
      dim<-ncol(theta_estimate) #number of dimensions
      n<-nrow(theta_estimate) #number of subjects
      
      if(!is.null(H) & !is.null(B)){
        huber_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "Huber", tuning.par = H)$theta
        bisquare_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "bisquare", tuning.par = B)$theta
        pnt.h<-matrix(apply(cbind(matrix(theta_estimate), matrix(huber_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        pnt.b<-matrix(apply(cbind(matrix(theta_estimate), matrix(bisquare_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.h = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)
        Distance.b = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)
        stats.h<-data.frame(Dis=apply(Distance.h, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        stats.b<-data.frame(Dis=apply(Distance.b, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        h.names<-b.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.h<-cbind(stats.h, MLE = theta_estimate[,i],
                         Huber = huber_theta_estimate[,i],
                         Distance =Distance.h[,i])
          stats.b<-cbind(stats.b, MLE = theta_estimate[,i],
                         Bisquare = bisquare_theta_estimate[,i],
                         Distance =Distance.b[,i])
          h.names<-c(h.names, paste0("MLE", i), paste0("Huber", i), paste0("Distance", i))
          b.names<-c(b.names, paste0("MLE", i), paste0("Bisquare", i), paste0("Distance", i))
        }
        colnames(stats.h)<-h.names
        colnames(stats.b)<-b.names
        sum.stats.h<-cbind(stats.h, t(dat))%>%arrange(desc(Dis))
        sum.stats.b<-cbind(stats.b, t(dat))%>%arrange(desc(Dis))
        
        for(i in 1:dim){
          #message(i)
          h.plots[[i]] <- local({
            i <- i
            huberplot<- ggplot(mapping = aes (x = theta_estimate[,i], y = huber_theta_estimate[,i]))+ geom_abline(color = "red", slope = 1) +
              geom_point() + xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) +
              ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE, Dimension" ~ .(i) ))
          })
          b.plots[[i]] <- local({
            i <- i
            bisquareplot<- ggplot(mapping = aes (x = theta_estimate[,i], y = bisquare_theta_estimate[,i]))+ geom_abline(color = "red", slope = 1) +
              geom_point()+ xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) +
              ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE, Dimension" ~ .(i) ))
          })
        }
        return(list("Summary Statistics (Huber)" = sum.stats.h[,-1], "Summary Statistics (Bisquare)" = sum.stats.b[,-1], "Huber Plot" = do.call(ggarrange, c(h.plots, ncol = 1, nrow = dim, common.legend = T)), "Bisquare Plot" = do.call(ggarrange, c(b.plots, ncol = 1, nrow = dim, common.legend = T))))
      }else if(is.null(H) & !is.null(B)){
        bisquare_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "bisquare", tuning.par = B)$theta
        
        pnt.b<-matrix(apply(cbind(matrix(theta_estimate), matrix(bisquare_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.b = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)
        stats.b<-data.frame(Dis=apply(Distance.b, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        b.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.b<-cbind(stats.b, MLE = theta_estimate[,i],
                         Bisquare = bisquare_theta_estimate[,i],
                         Distance =Distance.b[,i])
          b.names<-c(b.names, paste0("MLE", i), paste0("Bisquare", i), paste0("Distance", i))
        }
        colnames(stats.b)<-b.names
        sum.stats.b<-cbind(stats.b, t(dat))%>%arrange(desc(Dis))
        
        for(i in 1:dim){
          b.plots[[i]] <- local({
            i <- i
            bisquareplot<- ggplot(mapping = aes (x = theta_estimate[,i], y = bisquare_theta_estimate[,i]))+ geom_abline(color = "red", slope = 1) +
              geom_point()+ xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) +
              ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE, Dimension" ~ .(i) ))
          })
        }
        return(list("Summary Statistics" = sum.stats.b[,-1], "Bisquare Plots" =do.call(ggarrange, c(b.plots, ncol = 1, nrow = dim, common.legend = T))))
      }else if(!is.null(H) & is.null(B)){
        huber_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "Huber", tuning.par = H)$theta
        
        pnt.h<-matrix(apply(cbind(matrix(theta_estimate), matrix(huber_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.h = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)
        stats.h<-data.frame(Dis=apply(Distance.h, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        h.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.h<-cbind(stats.h, MLE = theta_estimate[,i],
                         Huber = huber_theta_estimate[,i],
                         Distance =Distance.h[,i])
          h.names<-c(h.names, paste0("MLE", i), paste0("Huber", i), paste0("Distance", i))
        }
        colnames(stats.h)<-h.names
        sum.stats.h<-cbind(stats.h, t(dat))%>%arrange(desc(Dis))
        
        for(i in 1:dim){
          h.plots[[i]] <- local({
            i <- i
            huberplot<- ggplot(mapping = aes (x = theta_estimate[,i], y = huber_theta_estimate[,i]))+ geom_abline(color = "red", slope = 1) +
              geom_point() + xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) +
              ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE, Dimension" ~ .(i) ))
          })
          
        }
        return(list("Summary Statistics (Huber)" = sum.stats.h[,-1], "Huber Plots" = do.call(ggarrange, c(h.plots, ncol = 1, nrow = dim, common.legend = T))))
      }else{ return(print("A valid tuning parameter is needed."))}
    }else{ # if not same plot
      theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "equal")$theta
      dim<-ncol(theta_estimate) #number of dimensions
      n<-nrow(theta_estimate) #number of subjects
      
      if(!is.null(H) & !is.null(B)){
        huber_theta_estimate<-theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "Huber", tuning.par = H)$theta
        bisquare_theta_estimate<-theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "bisquare", tuning.par = B)$theta
        
        pnt.h<-matrix(apply(cbind(matrix(theta_estimate), matrix(huber_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        pnt.b<-matrix(apply(cbind(matrix(theta_estimate), matrix(bisquare_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.h = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)
        Distance.b = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)
        stats.h<-data.frame(Dis=apply(Distance.h, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        stats.b<-data.frame(Dis=apply(Distance.b, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        h.names<-b.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.h<-cbind(stats.h, MLE = theta_estimate[,i],
                         Huber = huber_theta_estimate[,i],
                         Distance =Distance.h[,i])
          stats.b<-cbind(stats.b, MLE = theta_estimate[,i],
                         Bisquare = bisquare_theta_estimate[,i],
                         Distance =Distance.b[,i])
          h.names<-c(h.names, paste0("MLE", i), paste0("Huber", i), paste0("Distance", i))
          b.names<-c(b.names, paste0("MLE", i), paste0("Bisquare", i), paste0("Distance", i))
        }
        colnames(stats.h)<-h.names
        colnames(stats.b)<-b.names
        sum.stats.h<-cbind(stats.h, t(dat))%>%arrange(desc(Dis))
        sum.stats.b<-cbind(stats.b, t(dat))%>%arrange(desc(Dis))
        
        dat<-data.frame(MLE=matrix(theta_estimate),
                        Huber=matrix(huber_theta_estimate),
                        Bisquare=matrix(bisquare_theta_estimate),
                        Dimension=as.factor(rep(1:dim, each=n)))
        
        huberplot<- ggplot(data = dat, aes(x=MLE, y=Huber)) + geom_abline(color = "red", slope = 1) +
          geom_point(aes(color=Dimension)) + xlab(bquote(hat(theta)[MLE])) +
          ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) + ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE"))
        
        bisquareplot<- ggplot(data = dat, aes(x=MLE, y=Bisquare)) + geom_abline(color = "red", slope = 1) +
          geom_point(aes(color=Dimension)) + xlab(bquote(hat(theta)[MLE])) +
          ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) + ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE"))
        
        return(list("Summary Statistics (Huber)" = sum.stats.h[,-1], "Summary Statistics (Bisquare)" = sum.stats.b, "Huber Plots" = huberplot, "Bisquare Plots" = bisquareplot))
      }else if(is.null(H) & !is.null(B)){
        bisquare_theta_estimate<-theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "bisquare", tuning.par = B)$theta
        
        pnt.b<-matrix(apply(cbind(matrix(theta_estimate), matrix(bisquare_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.b <- sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)
        stats.b<-data.frame(Dis=apply(Distance.b, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        b.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.b<-cbind(stats.b, MLE = theta_estimate[,i],
                         Bisquare = bisquare_theta_estimate[,i],
                         Distance =Distance.b[,i])
          b.names<-c(b.names, paste0("MLE", i), paste0("Bisquare", i), paste0("Distance", i))
        }
        colnames(stats.b)<-b.names
        sum.stats.b<-cbind(stats.b, t(dat))%>%arrange(desc(Dis))
        
        
        dat<-data.frame(MLE=matrix(theta_estimate),
                        Bisquare=matrix(bisquare_theta_estimate),
                        Dimension=as.factor(rep(1:dim, each=n)))
        bisquareplot<- ggplot(data = dat, aes(x=MLE, y=Bisquare)) + geom_abline(color = "red", slope = 1) +
          geom_point(aes(color=Dimension)) + xlab(bquote(hat(theta)[MLE])) +
          ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) + ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE"))
        return(list("Sumary Statistics (Bisquare)" = sum.stats.b[,-1], "Bisquare Plot" = bisquareplot))
      }else if(!is.null(H) & is.null(B)){
        huber_theta_estimate<-theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "Huber", tuning.par = H)$theta
        
        pnt.h<-matrix(apply(cbind(matrix(theta_estimate), matrix(huber_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.h <- sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)
        stats.h<-data.frame(Dis=apply(Distance.h, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        h.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.h<-cbind(stats.h, MLE = theta_estimate[,i],
                         Huber = huber_theta_estimate[,i],
                         Distance =Distance.h[,i])
          h.names<-c(h.names, paste0("MLE", i), paste0("Huber", i), paste0("Distance", i))
        }
        colnames(stats.h)<-h.names
        sum.stats.h<-cbind(stats.h, t(dat))%>%arrange(desc(Dis))
        
        dat<-data.frame(MLE=matrix(theta_estimate),
                        Huber=matrix(huber_theta_estimate),
                        Dimension=as.factor(rep(1:dim, each=n)))
        
        huberplot<- ggplot(data = dat, aes(x=MLE, y=Huber)) + geom_abline(color = "red", slope = 1) +
          geom_point(aes(color=Dimension)) + xlab(bquote(hat(theta)[MLE])) +
          ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) + ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE"))
        
        return(list("Summary Statistics (Huber)" = sum.stats.h[,-1], "Huber Plots" = huberplot))
      }else{return(print("A valid tuning parameter is needed."))}
    }
  }else if(type == "GRM"){
    
    theta_estimate<-theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="equal")$theta
    if(!is.null(H)& !is.null(B)){
      huber_theta_estimate<-theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="Huber", tuning.par=H)$theta
      bisquare_theta_estimate<-theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="bisquare", tuning.par=B)$theta
      
      h.plot<- ggplot(mapping = aes (x = theta_estimate, y = huber_theta_estimate))+ geom_abline(color = "red", slope = 1) +
        geom_point() + xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) +
        ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE" ))
      
      b.plot<- ggplot(mapping = aes (x = theta_estimate, y = bisquare_theta_estimate))+ geom_abline(color = "red", slope = 1) +
        geom_point()+ xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) +
        ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE" ))
      
      pnt.h<-apply(cbind(theta_estimate, huber_theta_estimate), 1, function(x) sum(x)/2)
      pnt.b<-apply(cbind(theta_estimate, bisquare_theta_estimate), 1, function(x) sum(x)/2)
      sum.stats.h<-cbind(data.frame(ID = 1:nrow(theta_estimate),
                                    MLE = theta_estimate,
                                    Huber = huber_theta_estimate,
                                    Distance = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)), t(dat))%>%arrange(desc(Distance))
      sum.stats.b<-cbind(data.frame(ID = 1:nrow(theta_estimate),
                                    MLE = theta_estimate,
                                    Bisquare = bisquare_theta_estimate,
                                    Distance = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)), t(dat)) %>%arrange(desc(Distance))
      if(same.plot){ #allows user to have sperate plots or see both Huber and Bisquare at same time
        return(list("Summary Statistics (Huber)" = sum.stats.h, "Summary Statistics (Bisquare)" = sum.stats.b, "Plots" = do.call(ggarrange, c(list(h.plot, b.plot), ncol = 1, nrow = 2))))
      }else{
        return(list("Summary Statistics (Huber)" = sum.stats.h, "Summary Statistics (Bisquare)" = sum.stats.b, "Huber Plot" = h.plot, "Bisquare Plot" =b.plot))
      }
    }else if(!is.null(H)& is.null(B)){
      huber_theta_estimate<-theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="Huber", tuning.par=H)$theta
      
      h.plot<- ggplot(mapping = aes (x = theta_estimate, y = huber_theta_estimate))+ geom_abline(color = "red", slope = 1) +
        geom_point() + xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) +
        ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE" ))
      
      pnt.h<-apply(cbind(theta_estimate, huber_theta_estimate), 1, function(x) sum(x)/2)
      sum.stats.h<-cbind(data.frame(ID = 1:nrow(theta_estimate),
                                    MLE = theta_estimate,
                                    Huber = huber_theta_estimate,
                                    Distance = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)), t(dat))%>%arrange(desc(Distance))
      
      return(list("Summary Statistics (Huber)" = sum.stats.h,  "Huber Plot" = h.plot))
    }else if(is.null(H)& !is.null(B)){
      bisquare_theta_estimate<-theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="Huber", tuning.par=H)$theta
      
      b.plot<- ggplot(mapping = aes (x = theta_estimate, y = bisquare_theta_estimate))+ geom_abline(color = "red", slope = 1) +
        geom_point()+ xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) +
        ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE" ))
      
      pnt.b<-apply(cbind(theta_estimate, bisquare_theta_estimate), 1, function(x) sum(x)/2)
      sum.stats.b<-cbind(data.frame(ID = 1:nrow(theta_estimate),
                                    MLE = theta_estimate,
                                    Bisquare = bisquare_theta_estimate,
                                    Distance = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)), t(dat)) %>%arrange(desc(Distance))
      return(list("Summary Statistics (Bisquare)" = sum.stats.b, "Bisquare Plot" = b.plot))
      
    }else{return(print("A valid tuning parameter is needed."))}
    
  }
}

# An overall function for computing standard errors for multiple models
standard.errors <- function(a, thetas, dat, d=NULL, b=NULL, residual = "standardized", weight.function = "equal", tuning.par=NULL, model="2PL",  D=1.7){
  # a = J × K matrix of item discrimination parameters
  # thetas = vector of K theta values
  # dat = vector of J responses
  # d = Jx1 matrix of item difficulty parameters for MIRT model; 
  #     JxK matrix of category threshold parameters for MGRM model
  # b = vector of J item difficulty parameters for 2PL model;
  #     JxK matrix of category threshold parameters for GRM model
  # residual = information, standardized, MSR, DWMSR, DWSTZ
  # weight.function = Huber, biweight, equal
  # tuning.par = tuning parameter for weight function
  # model = MIRT, 2PL, GRM, MGRM
  # D = scaling constant
  
  if(model=="MIRT"){
    
    # Item Response probability
    P<-item.prob(matrix(thetas), "MIRT", cbind(a, d))
    
    # Calculate the residual for each item
    if(residual=="information"){
      # Information residual
      r<-a%*%thetas+d  
    }else if(residual=="standardized"){
      # Standardized residual
      r<-(dat-P)/sqrt(P*(1-P))  
    }else if(residual=="MSR"){
      # Modified standardized residual
      r<-(dat-P)/P
    }else if(residual=="DWMSR"){
      # Information residual and MSR
      r<-cbind(a%*%thetas+d, (dat-P)/P)
    }else if(residual=="DWSTZ"){
      # Information residual and STZ
      r<-cbind(a%*%thetas+d, (dat-P)/sqrt(P*(1-P)))
    }else{
      # If the residual is not identifiable
      stop("Cannot determine the residual type.")
    }
    
    # Calculate the weight based on the residual
    if(residual=="DWMSR"|residual=="DWSTZ"){
      # need to add the two weights for the dual methods
      if(weight.function =="Huber"){
        w<-huber(r[,1], tuning.par)+huber(r[,2], tuning.par)
      }else if(weight.function =="bisquare"){
        w<-bisquare(r[,1], tuning.par)+bisquare(r[,2], tuning.par)
      }else if(weight.function =="equal"){
        w<-1
      }
    }else{
      if(weight.function =="Huber"){
        w<-huber(r, tuning.par)
      }else if(weight.function =="bisquare"){
        w<-bisquare(r, tuning.par)
      }else if(weight.function =="equal"){
        w<-1
      }
    }
    
    # Weighted contribution to information (A matrix)
    wpq   <- c(w*P*(1 - P))
    # Weighted contribution to V matrix in expected SE
    w2pq  <- c(w^2*P*(1 - P))
    
    # KxK matrix A for expected and observed information (Hessian)
    A <- D^2 * t(a *wpq)%*%a 
    # KxK matrix V for expected information
    V <- D^2 * t(a*w2pq) %*%a   
    
    # First derivative of log likelihood
    D1 <- a * c(D*(w*(dat - P)))
    # K x K B Matrix for observed SE
    B <- t(D1)%*%D1
    
    
    # Covariance
    cov_magis <- solve(A) %*% V %*% solve(A)
    cov_sand  <- solve(A) %*% B %*% solve(A)
    
    return(list(observed_se = sqrt(diag(cov_sand)),
                expected_se = sqrt(diag(cov_magis))))
    
  }else if(model=="2PL"){
    
    # Item Response Probability
    P<- 1/(1+exp((-1.7)*(a*(thetas-b))))
    
    # Calculate the residual for each item
    if(residual=="information"){
      # Information residual
      r<-a*(thetas-b)  
    }else if(residual=="standardized"){
      # Standardized residual
      r<-(dat-P)/sqrt(P*(1-P))  
    }else if(residual=="MSR"){
      # Modified standardized residual
      r<-(dat-P)/P
    }else if(residual=="DWMSR"){
      # Information residual and MSR
      r<-cbind(a*(thetas-b), (dat-P)/P)
    }else if(residual=="DWSTZ"){
      # Information residual and STZ
      r<-cbind(a*(thetas-b), (dat-P)/sqrt(P*(1-P)))
    }else{
      # If the residual is not identifiable
      stop("Cannot determine the residual type.")
    }
    
    # Calculate the weight based on the residual
    if(residual=="DWMSR"|residual=="DWSTZ"){
      # need to add the two weights for the dual methods
      if(weight.function =="Huber"){
        w<-huber(r[,1], tuning.par)+huber(r[,2], tuning.par)
      }else if(weight.function =="bisquare"){
        w<-bisquare(r[,1], tuning.par)+bisquare(r[,2], tuning.par)
      }else if(weight.function =="equal"){
        w<-1
      }
    }else{
      if(weight.function =="Huber"){
        w<-huber(r, tuning.par)
      }else if(weight.function =="bisquare"){
        w<-bisquare(r, tuning.par)
      }else if(weight.function =="equal"){
        w<-1
      }
    }
    
    # A value for expected and observed information (2nd derivative)
    A <- D^2 * a*w*P*(1 - P) 
    # V value for expected information (w* 2nd derivative)
    V <- D^2 * a*w^2*P*(1 - P)   
    
    # First derivative of log likelihood
    D1 <- a * c(D*(w*(dat - P)))
    # B value for observed SE
    B <- sum(D1^2)
    
    
    # Covariance
    cov_magis <- V/A^2
    cov_sand  <- B/A^2
    
    return(list(observed_se = sqrt(diag(cov_sand)),
                expected_se = sqrt(diag(cov_magis))))
    
  }else if(model=="GRM"){
    
    # Item response probability
    probs <- item.prob(thetas, "GRM", cbind(a, b))
    # Probability of responding in each category
    P<-probs$P[,,1]
    pstars<-cbind(rep(1, length(a)), probs$pstar[,,1], rep(0, length(a)))
    
    # Number of category thresholds
    nthresh<-ncol(b)
    
    # Calculate the residual for each item
    if(residual=="standardized"){
      
      expected.value<-P%*%matrix(c(1:(nthresh+1))) 
      response.variance<-apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*P
      # Standardized residual
      r<-(dat-expected.value)/sqrt(rowSums(response.variance))
      
    }else if(residual=="MSR"){
      
      expected.value<-P%*%matrix(c(1:(nthresh+1))) 
      
      # Modified standardized residual (does this even apply to polytomous?)
      r<-(dat-expected.value)/P[cbind(1:length(dat), dat)]
    }else{
      # If the residual is not identifiable
      stop("Residual type cannot be used with the GRM or is otherwise unidentifiable.")
    }
    
    # Calculate residual-based weight
    if(weight.function =="Huber"){
      w<-huber(r, tuning.par)
    }else if(weight.function =="bisquare"){
      w<-bisquare(r, tuning.par)
    }else if(weight.function =="equal"){
      w<-1
    }
    
    # Pstar for k-1
    P1<-pstars[cbind(1:length(a), dat)]
    # Pstar for k
    P2<-pstars[cbind(1:length(a), dat+1)]
    # P for k
    Pk<-P1-P2
    
    # A value for expected and observed information (2nd Derivative)
    A <- sum(D^2 * a^2*w*( (P1*(1-P1)*(1-P1-P1)-P2*(1-P2)*((1-P2)-P2))/Pk - (P1*(1-P1)-P2*(1-P2))^2/Pk^2 ) )
    # V value for expected information - here I have copied A and multiplied it by w again
    #V <- -sum(D^2 * a^2*w^2*( (P1*(1-P1)*(1-P1-P1)-P2*(1-P2)*((1-P2)-P2))/Pk - (P1*(1-P1)-P2*(1-P2))^2/Pk^2 ) )
    
    # First derivative of log likelihood
    D1 <- a * c(D*(w*(P1*(1-P1)-P2*(1-P2))/Pk))
    # B value for observed SE
    B <- sum(D1^2)
    
    #cov_magis<-V/A^2 # need to derive the ASE based on Magis SE
    cov_sand <- B/A^2
    
    return(list(observed_se = sqrt(cov_sand), residual = r#,expected_se = sqrt(cov_magis)
    ))
  }else if(model=="MGRM"){
    
    # Item Response Probability
    probs<-item.prob(t(thetas), "MGRM", cbind(a, d))
    P<-probs$P[,,1]
    pstars<-cbind(rep(1, nrow(a)), probs$pstar[,,1], rep(0, nrow(a)))
    
    # Number of category thresholds
    nthresh<-ncol(d)
    
    # Calculate the residual for each item
    if(residual=="standardized"){
      
      expected.value<-P%*%matrix(c(1:(nthresh+1))) 
      response.variance<-apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*P
      # Standardized residual
      r<-(dat-expected.value)/sqrt(rowSums(response.variance))
      
    }else if(residual=="MSR"){
      
      expected.value<-P%*%matrix(c(1:(nthresh+1))) 
      
      # Modified standardized residual (does this even apply to polytomous?)
      r<-(dat-expected.value)/P[cbind(1:length(dat), dat)]
      
    }else{
      # If the residual is not identifiable
      stop("Residual type cannot be used with the GRM or is otherwise unidentifiable.")
    }
    
    # Calculate residual-based weight
    if(weight.function =="Huber"){
      w<-huber(r, tuning.par)
    }else if(weight.function =="bisquare"){
      w<-bisquare(r, tuning.par)
    }else if(weight.function =="equal"){
      w<-1
    }
    
    # Pstar for k-1
    P1<-pstars[cbind(1:nrow(a), dat)]
    # Pstar for k
    P2<-pstars[cbind(1:nrow(a), dat+1)]
    # P for k
    Pk<-P1-P2
    
    # A matrix (KxK) for expected and observed information (Hessian)
    A <- D^2 * t(a*c(w*( (P1*(1-P1)*(1-P1-P1)-P2*(1-P2)*((1-P2)-P2))/Pk - (P1*(1-P1)-P2*(1-P2))^2/Pk^2 )) )%*%a
    # V matrix (KxK) for expected information - here I have copied A and multiplied it by w again
    #V <- -D^2 * t(a*w^2*( (P1*(1-P1)*(1-P1-P1)-P2*(1-P2)*((1-P2)-P2))/Pk - (P1*(1-P1)-P2*(1-P2))^2/Pk^2 ) )%*%a
    
    # First derivative of log likelihood
    D1 <- a * c(D*(w*(P1*(1-P1)-P2*(1-P2))/Pk))
    # K x K B Matrix for observed SE
    B <- t(D1)%*%D1
    
    #cov_magis<-V/A^2 # need to derive the ASE based on Magis SE
    cov_sand <- B/A^2
    
    return(list(observed_se = sqrt(diag(cov_sand))#,expected_se = sqrt(cov_magis)
    ))
  }
  
  
}

###### Bayesian SE #######
bayes.standard.error <- function(thetas, post.var, dat, a, b, weight.function = "equal", tuning.par=NULL, D = 1.7, eps = 1e-12) {
  # thetas: vector of posterior means of length l
  # post.var: vector of posterior variances (e.g., (standard.error)^2 from your MAP/EAP)
  # returns vector of Bayesian robust SEs according to Li & Rice (2023)
  if (length(thetas) != length(post.var)) stop("thetas and post.var must match length")
  l <- length(thetas)    
  # Number of category thresholds
  nthresh<-ncol(b)
  brse <- rep(NA, l)
  
  for (i in 1:l){
    
    # Item response probability
    probs <- item.prob(thetas, "GRM", cbind(a, b))
    # Probability of responding in each category
    P<-probs$P[,,1]
    pstars<-cbind(rep(1, length(a)), probs$pstar[,,1], rep(0, length(a)))
    
    
    # Calculate the residual for each item
    expected.value<-P%*%matrix(c(1:(nthresh+1))) 
    response.variance<-apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*P
    # Standardized residual
    r<-(dat-expected.value)/sqrt(rowSums(response.variance))
    
    # Calculate residual-based weight
    if(weight.function =="Huber"){
      w<-huber(r, tuning.par)
    }else if(weight.function =="bisquare"){
      w<-bisquare(r, tuning.par)
    }else if(weight.function =="equal"){
      w<-1
    }
    
    # Pstar for k-1
    P1<-pstars[cbind(1:length(a), dat)]
    # Pstar for k
    P2<-pstars[cbind(1:length(a), dat+1)]
    # P for k
    Pk<-P1-P2
    
    # A value for expected and observed information (2nd Derivative)
    A <- -sum(D^2 * a^2*w*( (P1*(1-P1)*(1-P1-P1)-P2*(1-P2)*((1-P2)-P2))/Pk - (P1*(1-P1)-P2*(1-P2))^2/Pk^2 ) )
    
    # First derivative of log likelihood
    D1 <- a * c(D*(w*(P1*(1-P1)-P2*(1-P2))/Pk))
    # B value for observed SE
    B <- sum(D1^2)
    
    if (!is.finite(post.var[i])) { brse[i] <- NA; next }
    if (!is.finite(A) || A <= 0 || !is.finite(B)) {
      brse[i] <- NA
      next
    }
    # approximate multiplicative correction factor
    correction <- (B / A)
    
    if (!is.finite(correction) || correction <= 0) { brse[i] <- NA; next }
    var_adjusted <- post.var[i] * correction
    brse[i] <- sqrt(var_adjusted)
  }
  return(list(brse=brse, residual = r))
}

###### MLE Function #####
theta.est.grm <- function(dat, a, b, iter=30, cutoff=0.01, init.val=0, weight.type="equal", tuning.par=NULL, D=1.7) {
  
  # Check if the turning parameter is given when the weight.type is not "normal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The tuning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  
  # Get dimensions of the input data
  l <- nrow(dat)  # number of subjects
  J <- ncol(dat)  # test length (number of items)
  nthresh <- ncol(b)  # number of threshold parameters
  
  # Initialize arrays for storing results
  theta.est2 <- standard.error <- matrix(data=NA, nrow=l)
  convergence <- matrix(0, nrow=l)
  theta.progression <- matrix(NA, nrow = l, ncol = iter)
  residual <- matrix(data=NA, nrow = l, ncol = J)
  
  # Loop to estimate theta for each subject
  for(i in 1:l){
    
    # Initialize theta value
    if(length(init.val) > 1){ #if more than one initial value specified
      theta <- init.val[i]
    } else {
      theta <- init.val
    }
    
    P0 <- 0
    
    # Iterative loop for maximum likelihood estimation of theta
    for (k in 1:iter){ #k iterations at maximum
      
      # Item response probabilities
      probs <- item.prob(theta, "GRM", cbind(a, b))
      
      # subset probabilities
      P.i<-probs$P[,,1]
      probs.resp<-P.i[cbind(1:J, dat[i,])]
      
      expected.value<-P.i%*%matrix(c(1:(nthresh+1))) 
      # Calculate standardized residual
      residual[i,]<-(dat[i,]-expected.value)/sqrt(rowSums(apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*P.i))  
      
      # Compute weighting term based on specified weight function (bisquare, Huber, equal)
      weighting.term <- NULL
      if (weight.type == "bisquare") {
        weighting.term <- bisquare(residual[i,], tuning.par)
      } else if (weight.type == "Huber") {
        weighting.term <- huber(residual[i,], tuning.par)
      } else if (weight.type == "equal"){
        weighting.term <- 1
      }
      
      # Check if weighting term is determined
      if (is.null(weighting.term)) {
        stop("Cannot determine the weighting function.")
      }
      
      # subset probabilities for derivatives
      pstars<-cbind(rep(1, J), probs$pstar[,,1], rep(0, J))
      
      ps0<-pstars[cbind(1:J, dat[i,])]
      qs0<-1-ps0
      ps1<-pstars[cbind(1:J, dat[i,]+1)]
      qs1<-1-ps1
      
      # Compute first and second derivatives of the log-likelihood
      D1 <- sum(D * a * weighting.term * (ps0 * qs0 - ps1 * qs1) / probs.resp)
      D2 <- sum(D^2 * a^2 * weighting.term * ((ps0 * qs0 * (qs0 - ps0) - ps1 * qs1 * (qs1 - ps1)) / probs.resp - (ps0 * qs0 - ps1 * qs1)^2 / probs.resp^2 ))
      
      # Check for NAs in the computation: if so, record nonconvergence
      if (is.na(theta - D1/D2)) {
        theta.est2[i] <- theta <- NA
        convergence[i,1] <- 1
        break
      }
      
      # Update theta based on Newton-Raphson method, using the first and second derivatives
      theta <- theta.progression[i,k] <- theta - D1/D2
      
      # Compute difference in log-likelihood for convergence criterion
      log_like <- sum(log(probs.resp)) - sum(log(P0))
      
      # Check for convergence: stop Newton-Raphson method if 
      # log-likelihood difference is less than cutoff
      if (abs(log_like) < cutoff){
        break
      }
      
      # Update initial probability (P0) for the next iteration
      P0 <- probs.resp
    }
    
    # Store final estimated theta
    theta.est2[i] <- theta
    
    se.int<-standard.errors(a, theta, dat[i,], b=b, residual = "standardized", weight.function = weight.type, tuning.par=tuning.par, model="GRM",  D=D)
    standard.error[i] <- se.int$observed_se
    residual[i,]<-se.int$residual
    
    # Handle cases where theta did not converge within the desired number of iterations
    if (k == iter) {
      theta.est2[i] <- standard.error[i] <- NA
      convergence[i,1] <- 1
    } else if (!is.na(theta) & theta < -3) {
      # Then check: if theta converged outside [-3, 3], replace it with -3 or 3 respectively
      theta <- -3
      
      se.int<-standard.errors(a, theta, dat[i,], b=b, residual = "standardized", weight.function = weight.type, tuning.par=tuning.par, model="GRM",  D=D)
      standard.error[i] <- se.int$observed_se
      residual[i,]<-se.int$residual
      theta.est2[i] <- theta
    } else if (!is.na(theta) & theta > 3) {
      theta <- 3
      se.int<-standard.errors(a, theta, dat[i,], b=b, residual = "standardized", weight.function = weight.type, tuning.par=tuning.par, model="GRM",  D=D)
      standard.error[i] <- se.int$observed_se
      residual[i,]<-se.int$residual
      theta.est2[i] <- theta
    }
  }
  
  # Return a list containing the estimated theta, binary indicator of nonconvergence, standard error, estimated theta over each iteration, and standardized residual
  return(list(theta = theta.est2, convergence = convergence, standard.error = standard.error, theta.progression = theta.progression, residual = residual))
}

###### MAP function #####
theta.map.grm <- function(dat, a, b, prior = c(0, 1), iter=30, cutoff=0.01, init.val=0, weight.type="equal", tuning.par=NULL, D=1.7) {
  
  # Check if the turning parameter is given when the weight.type is not "normal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The tuning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  
  # Get dimensions of the input data
  l <- nrow(dat)  # number of subjects
  J <- ncol(dat)  # test length (number of items)
  nthresh <- ncol(b)  # number of threshold parameters
  
  # prior indicates our prior values for bayesian estimation
  mu<-prior[1] # prior mean
  sigma2<-prior[2] # prior variance
  
  # Initialize arrays for storing results
  theta.est2 <- standard.error <- matrix(data=NA, nrow=l)
  convergence <- matrix(0, nrow=l)
  theta.progression <- matrix(NA, nrow = l, ncol = iter)
  residual <- matrix(data=NA, nrow = l, ncol = J)
  
  # Loop to estimate theta for each subject
  for(i in 1:l){
    
    # Initialize theta value
    if(length(init.val) > 1){ #if more than one initial value specified
      theta <- init.val[i]
    } else {
      theta <- init.val
    }
    
    P0 <- 0
    
    # Iterative loop for MAP estimation of theta
    for (k in 1:iter){ #k iterations at maximum
      
      # Item response probabilities 
      probs <- item.prob(theta, "GRM", cbind(a, b))
      
      # subset probabilities
      P.i<-probs$P[,,1]
      probs.resp<-P.i[cbind(1:J, dat[i,])]
      
      expected.value<-P.i%*%matrix(c(1:(nthresh+1))) 
      # Calculate standardized residual
      residual[i,]<-(dat[i,]-expected.value)/sqrt(rowSums(apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*P.i))  
      
      # Compute weighting term based on specified weight function (bisquare, Huber, equal)
      weighting.term <- NULL
      if (weight.type == "bisquare") {
        weighting.term <- bisquare(residual[i,], tuning.par)
      } else if (weight.type == "Huber") {
        weighting.term <- huber(residual[i,], tuning.par)
      } else if (weight.type == "equal"){
        weighting.term <- 1
      }
      
      # Check if weighting term is determined
      if (is.null(weighting.term)) {
        stop("Cannot determine the weighting function.")
      }
      
      # subset probabilities for derivatives
      pstars<-cbind(rep(1, J), probs$pstar[,,1], rep(0, J))
      
      ps0<-pstars[cbind(1:J, dat[i,])]
      qs0<-1-ps0
      ps1<-pstars[cbind(1:J, dat[i,]+1)]
      qs1<-1-ps1
      
      # Compute first and second derivatives of the log-likelihood
      D1 <- sum(D * a * weighting.term * (ps0 * qs0 - ps1 * qs1) / probs.resp)- (theta - mu)/sigma2
      D2 <- sum(D^2 * a^2 * weighting.term * ((ps0 * qs0 * (qs0 - ps0) - ps1 * qs1 * (qs1 - ps1)) / probs.resp - (ps0 * qs0 - ps1 * qs1)^2 / probs.resp^2 )) - 1/sigma2
      
      # Check for NAs in the computation: if so, record nonconvergence
      if (is.na(theta - D1/D2)) {
        theta.est2[i] <- theta <- NA
        convergence[i,1] <- 1
        break
      }
      
      # Update theta based on Newton-Raphson method, using the first and second derivatives
      theta <- theta.progression[i,k] <- theta - D1/D2
      
      # Compute difference in log-likelihood for convergence criterion
      log_like <- sum(log(probs.resp)) - sum(log(P0))
      
      # Check for convergence: stop Newton-Raphson method if 
      # log-likelihood difference is less than cutoff
      if (abs(log_like) < cutoff){
        break
      }
      
      # Update initial probability (P0) for the next iteration
      P0 <- probs.resp
    }
    
    # Store final estimated theta
    theta.est2[i] <- theta
    
    # Compute posterior sd
    post.sd <- sqrt(-1 / D2)
    #Compute sandwich standard error
    bse.int<-bayes.standard.error(theta, post.sd^2, dat[i,], a, b, weight.type, tuning.par, D = D, eps = 1e-12) 
    standard.error[i]<-bse.int$brse
    residual[i,]<-bse.int$residual
    
    # Handle cases where theta did not converge within the desired number of iterations
    if (k == iter) {
      theta.est2[i] <- standard.error[i] <- NA
      convergence[i,1] <- 1
    } else if (!is.na(theta) & theta < -3) {
      # Then check: if theta converged outside [-3, 3], replace it with -3 or 3 respectively
      theta <- -3
      # Compute posterior sd
      post.sd <- sqrt(-1 / D2)
      #Compute sandwich standard error 
      bse.int<-bayes.standard.error(theta, post.sd^2, dat[i,], a, b, weight.type, tuning.par, D = D, eps = 1e-12) 
      standard.error[i]<-bse.int$brse
      residual[i,]<-bse.int$residual
      
      theta.est2[i] <- theta
    } else if (!is.na(theta) & theta > 3) {
      theta <- 3
      # Compute posterior sd
      post.sd <- sqrt(-1 / D2)
      #Compute sandwich standard error 
      
      bse.int<-bayes.standard.error(theta, post.sd^2, dat[i,], a, b, weight.type, tuning.par, D = D, eps = 1e-12) 
      standard.error[i]<-bse.int$brse
      residual[i,]<-bse.int$residual
    }
  }
  
  # Return a list containing the estimated theta, binary indicator of nonconvergence, standard error, estimated theta over each iteration, and standardized residual
  return(list(theta = theta.est2, convergence = convergence, standard.error = standard.error, theta.progression = theta.progression, residual = residual))
}

###### EAP function #####
theta.eap.grm <- function(dat, a, b, prior=c(0,1), eap.quad.pts =seq(-4, 4, length.out = 41), iter=30, cutoff=0.01, init.val=0, weight.type="equal", tuning.par=NULL, D=1.7) {
  
  # Check if the turning parameter is given when the weight.type is not "normal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The tuning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  
  # Get dimensions of the input data
  l <- nrow(dat)  # number of subjects
  J <- ncol(dat)  # test length (number of items)
  nthresh <- ncol(b)  # number of threshold parameters
  
  # prior indicates our prior values for bayesian estimation
  mu<-prior[1] # prior mean
  sigma2<-prior[2] # prior variance
  
  # Initialize arrays for storing intermediate results
  theta.est2 <- standard.error <- matrix(data=NA, nrow=l)
  residual <- matrix(data=NA, nrow = l, ncol = J)
  theta.progression <- matrix(NA, nrow = l, ncol = iter)
  
  # Item response probabilities at each quadrature point
  probs.q <- item.prob(eap.quad.pts, "GRM", cbind(a, b))
  
  # Loop to estimate theta for each subject
  for(i in 1:l){
    P0<-0
    # Initialize theta value for residual
    if(length(init.val) > 1){ #if more than one initial value specified
      theta <- init.val[i]
    } else {
      theta <- init.val
    }
    
    for(k in 1:iter){
      # quadrature point probabilities
      probs.quad <- probs.q$P[cbind(seq_along(dat[i,]), dat[i,], rep(seq_len(dim(probs.q$P)[3]), each = length(dat[i,])))]
      dim(probs.quad) <- c(length(dat[i,]), dim(probs.q$P)[3])
      
      # Item response probabilities for person's residual
      probs<-item.prob(theta, "GRM", cbind(a, b))
      
      
      # subset probabilities for person residual
      P.i<-probs$P[,,1]
      probs.resp<-P.i[cbind(1:J, dat[i,])]
      
      expected.value<-P.i%*%matrix(c(1:(nthresh+1))) 
      # Calculate standardized residual
      residual[i,]<-(dat[i,]-expected.value)/sqrt(rowSums(apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*P.i))  
      
      # Compute weighting term based on specified weight function (bisquare, Huber, equal)
      if (weight.type == "bisquare") {
        weighting.term <- bisquare(residual[i,], tuning.par)
      } else if (weight.type == "Huber") {
        weighting.term <- huber(residual[i,], tuning.par)
      } else{
        weighting.term <- 1
      }
      
      # Prior density according to each density point
      f_x<-dnorm(eap.quad.pts, mu, sqrt(sigma2))
      
      # (Weighted) Likelihood according to person's data and each quad point
      likelihood<- apply(probs.quad, 2, function(x) prod((x^dat[i,])^weighting.term))
      
      # EAP theta estimate
      theta<-sum(eap.quad.pts*likelihood*f_x)/sum(likelihood*f_x)
      
      # Store final estimated theta
      theta.est2[i] <- theta.progression[i,k]<-theta
      
      # Compute difference in log-likelihood for convergence criterion
      log_like <- sum(log(probs.resp)) - sum(log(P0))
      
      # Check for convergence: stop Newton-Raphson method if 
      # log-likelihood difference is less than cutoff
      if (abs(log_like) < cutoff){
        break
      }
      
      # Update initial probability (P0) for the next iteration
      P0 <- probs.resp
    }  
    
    # Posterior standard deviation of EAP theta estimate according to Bock & Mislevy, 1989; Thissen et al., 1995
    post.sd<-sqrt(sum((eap.quad.pts - theta)^2*likelihood*f_x)/ sum(likelihood*f_x))
    
    #Compute sandwich standard error 
    bse.int<-bayes.standard.error(theta, post.sd^2, dat[i,], a, b, weight.type, tuning.par, D = D, eps = 1e-12) 
    standard.error[i]<-bse.int$brse
    residual[i,]<-bse.int$residual
  }
  
  # Return a list containing the estimated theta,  standard error, and standardized residual
  return(list(theta = theta.est2, standard.error = standard.error, theta.progression = theta.progression, residual = residual))
}

#' Response Probability Calculation (MIRT) -- deprecated (see item.prob)
#'
#' Calculate the item response probabilities given the item parameters and latent trait vector for person \emph{i}, which is
#' \eqn{\boldsymbol{\theta}_i = ({\theta}_{i1}, {\theta}_{i2}, ... {\theta}_{iL})'} for \emph{L} dimensions, according to the multidimensional IRT model (MIRT; McKinley & Reckase, 1983):
#' \deqn{P (X_{ij}=1 | \boldsymbol{\theta}_i) =  \frac{1}{1+ e^{-1.7(\boldsymbol{a_j}\boldsymbol\theta+d_j)}}} for \eqn{j=1,...,J} items.
#' An item \emph{j} has slope parameters \eqn{\boldsymbol{a}_j=(a_{1j}, a_{2j}, ..., a_{Lj})'} and intercept parameter \eqn{d_j}.
#' @param thetas A vector of length \emph{L} containing a subject's latent traits, where \emph{L} is the number of test dimensions
#' @param a A \eqn{J \times L} matrix containing fixed item slope parameters for \emph{L} dimensions and \emph{J} test items
#' @param d A vector of length \emph{J} containing fixed intercept parameters, for \emph{J} test items
#' @references McKinley, R. L., & Reckase, M. D. (1983, August). \emph{An Extension of the Two-Parameter Logistic Model to the Multidimensional Latent Space} (Research No. ONR83-2). Iowa City, IA: American College Testing Program.
#' @return P A vector of probabilities for endorsing each of \emph{J} test items
#' @return residual A vector of residuals detecting misfit between the subject and \emph{J} test items, according to the formula \eqn{r_{ij} = \textbf a_j\boldsymbol{\theta}_i + d_j} (see theta.est)
probs.calc <- function(thetas, a, d){
  thetas <- matrix(thetas, nrow = dim(a)[2], byrow = T)
  r <- a%*%thetas+d
  P<- 1/(1+exp((-1.7)*r))
  return(list(residual=r, P=matrix(P, ncol=1)))
}

#' GRM Response Probability \eqn{P^*} and Response Category Probabilities (P) Calculation (deprecated - see item.prob)
#'
#' Generate \eqn{P^*}, the threshold probability, and \eqn{P}, the probability of responding in each category for a vector of latent traits
#' using parameters needed for the GRM (Samejima, 1969).
#' The probability that a subject responds in or above a category \eqn{k} for item \eqn{j} is \eqn{P^*_{jk}(\theta) = \frac{1}{1+ e^{-a_j (\theta-b_{jk})}}},
#' for \eqn{K} categories and \eqn{K-1} threshold parameters  (\eqn{b_{j1}, ... b_{j,K-1}}), where \eqn{b_{jk}} separates response category \eqn{k} and \eqn{k+1} (\eqn{k=1,..K-1}) (Embretson & Reise, 2000).
#' \eqn{a_j} is the item discrimination parameter.  The probability of endorsing exactly category \eqn{k} is \eqn{P_{jk}(\theta) = P^*_{j,k}(\theta) - P^*_{j,k+1}(\theta),} where \eqn{P^*_{j1}(\theta) \equiv 1.0} and \eqn{P^*_{jK}(\theta) \equiv 0.0.}
#' @param thetas A vector of latent traits for \eqn{N} subjects
#' @param a A vector of discrimination parameters for \eqn{J} test items
#' @param b A matrix of category threshold parameters, with the number of rows corresponding to \eqn{J} test items and the number of columns corresponding to \eqn{K-1} thresholds between \eqn{K} categories
#' @references Embretson, S. E., & Reise, S. P. (2000). \emph{Item response theory for psychologists.} Mahwah, N.J: L. Erlbaum Associates.
#' @references Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. \emph{Psychometrika Monograph Supplement, 34} (4, Pt. 2), 100–100.
#' @return pstar A \eqn{J \times (K-1) \times N} array of category threshold probabilities for \eqn{J} items, \eqn{K} categories, and \eqn{N} subjects
#' @return P a \eqn{J \times K \times N} array of probabilities of responding in each of \eqn{K} categories for \eqn{J} items, across \eqn{N} subjects
#' @examples
#' thetas <- rnorm(15)  # Example latent traits for 15 subjects
#' a <- runif(10, 0.5, 1.5)  # Example discrimination parameters for 10 items
#' b <- t(apply(matrix(runif(10*4, -2.5,2.5), nrow = 10, ncol = 4), 1, sort))  # Example threshold parameters for 10 items and 4 thresholds (5 categories)
#' result <- probs.calc.grm(thetas, a, b)  # Calculate probabilities
probs.calc.grm<-function(thetas, a, b){
  # Number of subjects
  N<-length(thetas)
  # Number of items
  J<-nrow(b) 
  # Number of thresholds between categories
  nthresh<-ncol(b)
  a<-matrix(a, nrow=J, ncol=nthresh)
  
  # compute exponent
  exponent <- array(0, dim = c(J, nthresh, N))
  for(i in 1:N){
    exponent[,,i] <- a * as.matrix((thetas[i] - b))
  }
  # Calculate P* using the GRM formula
  pstar <- 1/(1+exp(-1.7*exponent))
  # Convert P* to category probabilities using the pstar_to_p function
  P<-pstar_to_p(pstar)
  
  # Return both P* and P
  return(list(pstar = pstar, P = P))
}

#' Item category response probability for MGRM (deprecated - see item.prob)

probs.mgrm<-function(thetas, a, b, D=1.7){
  
  thetas<-as.matrix(thetas)
  
  # Number of subjects
  N<-nrow(thetas)
  # Number of items
  J<-nrow(a) 
  # Number of thresholds between categories
  nthresh<-ncol(b)
  #Number of dimensions
  K<-ncol(a)
  
  # Generate threshold probabilities based on theta
  exponent <- array(dim = c(J,nthresh,N))
  for(i in 1:N){
    for(h in 1:nthresh){
      temp<-matrix(ncol=K, nrow=J)
      for(k in 1:K){
        temp[,k]<-thetas[i,k]-b[,h]
      }
      exponent[,h,i]<-rowSums(a*temp)
    }
  }  
  
  # Calculate P* using the GRM formula
  pstar <- 1/(1+exp(-D*exponent))
  # Convert P* to category probabilities using the pstar_to_p function
  P<-pstar_to_p(pstar)
  
  # Return both P* and P
  return(list(pstar = pstar, P = P, exponent = exponent))
}


#' Deprecated: see dat.gen() Likert-type data generation function
#'
#' Generate Likert-type data from probabilities of responding in each category
#' @param P A \eqn{J \times K \times N} array of probabilities of responding in each of \emph{K} categories over \emph{J} items for \emph{N} subjects
#' @return dat A \eqn{N \times J} matrix of randomly generated Likert-type data using the sample() function for \emph{N} subjects over \emph{J} items
#' @examples
#' thetas <- rnorm(15)  # Example latent traits for 15 subjects
#' a <- runif(10, 0.5, 1.5)  # Example discrimination parameters for 10 items
#' b <- t(apply(matrix(runif(10*4, -2.5,2.5), nrow = 10, ncol = 4), 1, sort))  # Example threshold parameters for 10 items and 4 thresholds (5 categories)
#' probs <- item.prob(thetas, "GRM", cbind(a, b))  # Calculate probabilities
#' data <- data.gen(probs$P)  # Generate Likert-type data

data.gen<-function(P){
  if(is.array(P)){ 
    # If more than one subject
    # Initialize matrix for generated data
    dat<-matrix(nrow = nrow(P), ncol = dim(P)[3])
    for(i in 1:dim(P)[3]){
      # Loop over each subject
      # For each item, sample a response category based on the probabilities
      dat[,i]<-apply(P[,,i], 1, function(x) sample(c(1:dim(P)[2]), 1, replace = T, prob=x))
    }
  }else{
    # If only one subject
    # For each item, sample a response category based on the probabilities
    dat<-apply(P, 1, function(x) sample(c(1:dim(P)[2]), 1, replace = T, prob=x))
  } 
  # Return the generated Likert-type data
  return(t(dat))
}



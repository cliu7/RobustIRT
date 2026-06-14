#' Convert category threshold probabilities to the probabilities of responding in each category
#'
#' Calculate \eqn{P}, the probabilities of responding in each category, from the \eqn{P^*} threshold values using the graded response model (GRM; Samejima, 1969).
#' The probability that a subject responds in or above a category \eqn{k} for item \eqn{j} is \eqn{P^*_{jk}(\theta) = \frac{1}{1+ e^{-a_j (\theta-b_{jk})}}},
#' for \eqn{K} categories and \eqn{K-1} threshold parameters  (\eqn{b_{j,1}, ..., b_{j,K-1}}), where \eqn{b_{j,k}} separates response category \eqn{k} and \eqn{k+1} (\eqn{k=1,...K-1}) (Embretson & Reise, 2000).
#' \eqn{a_j} is the item discrimination parameter.  The probability of endorsing exactly category \eqn{k} is \eqn{P_{jk}(\theta) = P^*_{j,k}(\theta) - P^*_{j,k+1}(\theta),} where \eqn{P^*_{j1}(\theta) \equiv 1.0} and \eqn{P^*_{jK}(\theta) \equiv 0.0.}
#' @references Embretson, S. E., & Reise, S. P. (2000). \emph{Item response theory for psychologists.} Mahwah, N.J: L. Erlbaum Associates.
#' @references Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. \emph{Psychometrika Monograph Supplement, 34} (4, Pt. 2), 100–100.
#' @param Pstar A \eqn{J \times K-1 \times N} array of \eqn{P^*} threshold probability values, for \eqn{K} categories, \eqn{J} items, and \eqn{N} subjects
#' @return The probabilities \eqn{P} of responding in each category.
#' @examples
#' # One-subject case
#' Pstar <- matrix(c(0.85, 0.50, 0.20,0.70, 0.40, 0.10,0.90, 0.60, 0.30), nrow = 3, byrow = TRUE)
#' pstar_to_p(Pstar)
#' 
#' # Multi-subject case
#' J <- 2   # items
#' K <- 4   # categories (0–3)
#' N <- 3   # persons
#' # Simulate P* values
#' Pstar <- array(runif(J * (K - 1) * N), dim = c(J, K - 1, N))
#' Pstar <- aperm(apply(Pstar, c(1, 3), sort, decreasing = TRUE), c(2, 1, 3))
#' pstar_to_p(Pstar)
#' @export

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
#'
#' Computes item response probabilities for select IRT models (1PL, Rasch, 2PL, MIRT, GRM, and MGRM), given ability and item parameters.
#' by constructing the appropriate linear predictors and applying the logistic function. Returns item response probabilities for dichotomous data or item category response probabilities for polytomous data.
#' @references Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee’s ability. In F. M. Lord & M. R. Novick (Eds.), Statistical Theories of Mental Test Scores (pp. 397–479). \emph{Addison‑Wesley}.
#' @references Lord, F. M., & Novick, M. R. (1968). Statistical theories of mental test scores. \emph{Addison-Wesley}.
#' @references Lord, F. M. (1980). Applications of item response theory to practical testing problems. \emph{Erlbaum}. https://doi.org/10.4324/9780203056615.
#' @references McKinley, R. L., & Reckase, M. D. (1983). An extension of the two‑parameter logistic model to the multidimensional latent space. \emph{Psychometrika}, 48(3), 369–382.
#' @references Mosteller, F., & Tukey, J. W. (1977). Data Analysis and Regression: A Second Course in Statistics. \emph{Addison‑Wesley}.
#' @references Muraki, E., & Engelhard, G. (1985). Full-information item factor analysis: Applications of EAP scores. \emph{Applied Psychological Measurement}, 9(4), 417–430
#' @references Rasch, G. (1960). Probabilistic models for some intelligence and attainment tests. \emph{Danish Institute for Educational Research}, 184.
#' @references Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. \emph{Psychometrika Monograph Supplement, 34} (4, Pt. 2), 100–100.
#' @param theta A numeric vector or matrix of latent trait values. 
#' @param ipars A matrix of item parameters. See examples for how to structure the columns of the matrix based on the model utilized.
#' @param model A character string specifying which IRT model to use. See details for formulae.
#' \itemize{
#'   \item \code{"Rasch"}: Allows item difficulty parameters to vary across items (Rasch, 1960).
#'   \item \code{"1PL"}: 1‑parameter logistic model with a common discrimination parameter and item‑specific difficulties. The Rasch model is the special case where all discriminations equal 1.
#'   \item \code{"2PL"}: 2‑parameter logistic model allowing both discrimination and difficulty to vary across items (Birnbaum, 1968).
#'   \item \code{"MIRT"}: Multidimensional extension of the 2PL model with item slope vectors across latent dimensions and item‑specific intercepts (McKinley & Reckase, 1983; Muraki & Engelhard, 1985).
#'   \item \code{"GRM"}: Graded response model for ordered polytomous items with item‑specific discrimination and ordered category thresholds (Samejima, 1969).
#'   \item \code{"MGRM"}: Multidimensional graded response model extending Samejima’s GRM to multiple latent dimensions, with slope vectors and ordered category thresholds (Muraki & Carlson, 1995).
#' }
#' @param D A positive scaling constant used for scaling the normal ogive model. Defaults to 1.7; alternatively is often set to 1.0.
#' @return For model accommodating dichotomous data ("Rasch", "1PL", "2PL", "MIRT"), returns an \eqn{N \times J} matrix of response probabilities \eqn{P(X = 1)}.
#' @return For models accommodating polytomous data ("GRM", "MGRM"), returns a list with:
#' \itemize{
#'   \item{pstar}: an array of cumulative probabilities \eqn{P^*(X \geq k)}.
#'   \item{P}: an array of category probabilities \eqn{P(X = k)}.
#' } 
#' @section IRT Models:
#'
#' \strong{Rasch}
#'
#' \deqn{
#'   P(X_{ij} = 1 \mid \theta_i) =
#'   \frac{1}{1 + e^{-(\theta_i - b_j)}}
#' }
#'
#' \strong{1PL}
#'
#' \deqn{
#'   P(X_{ij} = 1 \mid \theta_i) =
#'   \frac{1}{1 + e^{-Da(\theta_i - b_j)}}
#' }
#'
#' \strong{2PL}
#'
#' \deqn{
#'   P(X_{ij} = 1 \mid \theta_i) =
#'   \frac{1}{1 + e^{-Da_j(\theta_i - b_j)}}
#' }
#'
#' \strong{MIRT}
#'
#' \deqn{
#'   P(X_{ij} = 1 \mid \boldsymbol{\theta}_i) =
#'   \frac{1}{1 + e^{-D(\boldsymbol{a}_j\boldsymbol{\theta}_i + d_j)}}
#' }
#'
#' \strong{GRM}
#'
#' \deqn{
#'   P(X_{ij} = k \mid \theta_i) =
#'   \frac{1}{1 + e^{-Da_j(\theta_i - b_{jk})}}
#'   -
#'   \frac{1}{1 + e^{-Da_j(\theta_i - b_{j(k+1)})}}
#' }
#'
#' \strong{MGRM}
#'
#' \deqn{
#'   P(X_{ij} = k \mid \boldsymbol{\theta}_i) =
#'   \frac{1}{1 + e^{-D \sum_{l=1}^L a_{jl}(\boldsymbol{\theta}_{l} - d_{jk})}}
#'   -
#'   \frac{1}{1 + e^{-D \sum_{l=1}^L a_{jl}(\boldsymbol{\theta}_{l} - d_{j(k+1)})}}
#' }
#' 
#' @examples
#' # Rasch case
#' N <- 3
#' J <- 5
#' theta <- rnorm(N) # generate ability values
#' ipars <- rnorm(J)  # generate item difficulties
#' item.prob(theta, "Rasch", ipars)
                        
#' # 1PL case
#' N <- 4 # subjects
#' J <- 6 # items
#' theta <- rnorm(N) # generate ability values
#' ipars <- cbind(a = rep(1.2, J), # set item discrimination
#'                b = rnorm(J)) # generate item difficulties
#' item.prob(theta, "1PL", ipars)

#' # 2PL case
#' N <- 3 # subjects
#' J <- 5 # items
#' theta <- rnorm(N) # generate ability values
#' ipars <- cbind(a = runif(J, 0.5, 2), # generate item discrimination
#'                b = rnorm(J)) # generate item difficulty
#' item.prob(theta, "2PL", ipars) 
 
#' # MIRT case
#' N <- 2 # subjects
#' J <- 7 # items
#' L <- 3 # dimensions
#' theta <- matrix(rnorm(N * L), ncol = L) # N x L ability matrix
#' ipars <- cbind(matrix(runif(J * L, 0, 1), ncol = L), d = rnorm(J)) # generate slopes + intercept
#' item.prob(theta, "MIRT", ipars)

#' # GRM case
#' N <- 3 # subjects
#' J <- 5 # items
#' K <- 4 # categories
#' theta <- rnorm(N)
#' a <- runif(J, 0.5, 2) # J discriminations
#' b_raw <- matrix(rnorm(J * (K - 1)), nrow = J)  
#' b <- t(apply(b_raw, 1, sort)) # Sort thresholds within each item (row-wise)
#' ipars <- cbind(a, b) # Combine into a J x K matrix
#' item.prob(theta, "GRM", ipars)
 
#' # MGRM case
#' N <- 3 # subjects
#' J <- 5 # items
#' L <- 2 # dimensions
#' K <- 4 # categories
#' theta <- matrix(rnorm(N * L), ncol = L) 
#' a <- matrix(runif(J * L, 0, 1), ncol = L) # slopes 
#' b_raw <- matrix(rnorm(J * (K - 1)), nrow = J) # thresholds
#' b <- t(apply(b_raw, 1, sort))
#' ipars <- cbind(a, b)
#' item.prob(theta, "MGRM", ipars)
#' @export

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
  
  # Rasch predictor: theta - b
  # sapply loops over each person's theta value (Theta[,1]) 
  # For each x = Theta[n,1], compute x - ipars (vector of item difficulties) 
  # sapply returns J x N, then t() makes it N x J
   if(model=="RASCH"){
    ex<-t(sapply(Theta[,1], function(x) x-ipars))
  }
  
  # 1PL (a is constrained equal over items) or 2PL predictor: a * (theta - b)
  if(model=="1PL" | model=="2PL"){
    ex<-t(sapply(Theta[,1], function(x) ipars[,1]*(x-ipars[,2]))) # ipars[,1] = a_j (discrimination), ipars[,2] = b_j (difficulty)
  }
  
  # MIRT predictor: a·theta + d
  if(model=="MIRT"){
    a<-ipars[,1:L] # a: J x L matrix of discriminations
    d<-ipars[,L+1] # d: J-length vector of intercepts
    ex<-Theta %*% t(a)+ matrix(d, nrow=N, ncol=J, byrow=TRUE)
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
  if(model %in% c("RASCH", "1PL", "2PL", "MIRT")){
    return(P=invlogit(ex)) #Returns N x J matrix of P(X = 1)
  }
  
  # Apply logistic function and return probabilities for polytomous models
  if(model %in% c("GRM", "MGRM")){
    pstar<-invlogit(ex) # cumulative probabilities P*(X >= k)
    return(list(pstar=pstar, P=pstar_to_p(pstar))) #converts cumulative probabilities to category probabilities P(X = k)
  }
  
}
 
#' Residual Calculation
#' 
#' Computes the standardized and modified standardized residuals (MSRs; Yu & Cheng, 2019) for dichotomous and polytomous IRT models. For each observed response 
#' \eqn{y_j} to item \eqn{j}, the standardized residual and the MSR reflect the difference between the observed score and the model-implied expected value given the person’s estimated ability \eqn{\hat{\theta}}.
#' The basic standardized residual is given by \deqn{r_j = \frac{y_j - E(Y_j | \hat{\theta})}{\sqrt{\mathrm{Var}(Y_j | \hat{\theta})}}} where the expectation and variance are computed under the specified IRT model.
#' The MSR replaces the variance term with the conditional probability of the observed response category,
#' \eqn{P(y_j \mid \hat{\theta})}, yielding \deqn{r_j = \frac{y_j - E(Y_j | \hat{\theta})}{P(y_j|\hat{\theta})}}.
#' The information residual is available for dichotomous models only, capturing the difference between the person ability parameter and the difficulty parameter of an item, weighted by the item discrimination:
#' \eqn{a_j(\theta -b_j)} for unidimensional models and \eqn{\boldsymbol{a}_j \boldsymbol{\theta} +d_j)} for the MIRT model.
#' See \code{item.prob} for description of models and structure of \code{ipars}.
#' @references Yu, X & Cheng, Y. A change-point analysis procedure based on weighted residuals to detect back random responding. \emph{Psychological Methods} (Oct. 2019), pp. 658–674. DOI: 10.1037/met0000212.
#' @param theta An \eqn{N \times L} matrix of latent trait values, where \eqn{L} is the number of dimensions.
#' @param model Character string specifying the IRT model. See \code{item.prob} for supported models.
#' @param ipars A matrix or list of item parameters passed to \code{item.prob}. For dichotomous models, rows contain discrimination and difficulty parameters. For polytomous models, rows contain item discriminations and category thresholds.
#' @param dat An \eqn{N \times J} response matrix, with \eqn{N} respondents and \eqn{J} items.
#' @param residual Character string indicating which residual type to return: "standardized" or "msr". Defaults to both.
#' @param D Positive scaling constant for the normal-ogive approximation. Defaults to 1.7, otherwise often set to 1.0.

#' @return A list containing \eqn{N \times J} matrix for each residual specified:
#' \itemize{
#'   \item \code{standardized} {standardized residuals}
#'   \item \code{msr} {modified standardized residuals (Yu & Cheng, 2019)}
#'   \item \code{information} {information residual}
#' }
#' 
#' @examples
#' # 2PL model
#' dat <- matrix(c(1, 0, 1, 1), ncol = 1)
#' theta <- c(-1.0, 0.0, 0.5, 1.0)
#' ipars <- cbind(a = 1.0, b = -0.5)
#'
#' residual(theta, model = "2PL", ipars=ipars, dat=dat)
#'
#'
#' # MIRT model
#' dat <- matrix(c(1, 0, 1, 1, 0, 1, 0, 1), nrow = 4, byrow = TRUE)
#'
#' theta <- matrix(c(
#'   -1.0,  0,
#'    0,  0.5,
#'    1, -0.5,
#'    0.5, 1.0
#' ), ncol = 2, byrow = TRUE)
#'
#' ipars <- cbind(
#'   a1 = c(1.0, 2.0),
#'   a2 = c(0.5, 1.0),
#'   b  = c(-0.5, 0.5)
#' )
#'
#' residual(theta, model = "MIRT", ipars=ipars, dat=dat)
#'
#'
#' # GRM model
#' dat <- matrix(c(0,1,2, 1,2,3), nrow = 2, byrow = TRUE)
#' theta <- c(-0.5, 1.0)
#' 
#' ipars <- rbind(
#'   c(a = 1.0, b1 = -2.0, b2 = -1.0, b3 = 0.0),
#'   c(a = 0.5, b1 = -1.0, b2 =  0.0, b3 = 1.0)
#'
#' residual(theta, model = "GRM", ipars=ipars, dat=dat)
#' @export


residual<-function(theta, model, ipars, dat, resid = c("standardized", "msr", "information"), D=1.7){
    
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
  
  # Item (category) response probability
  probs<-item.prob(theta, model, ipars, D)
  
  out<-list() # initialize vector for storing output
  
  if(model %in% c("1PL", "2PL", "MIRT", "RASCH")){ # dichotomous data
    
    if("standardized" %in% resid){
      # standardized residual
      stz<-(dat-probs)/sqrt(probs*(1-probs))
      out$standardized<-stz
    }
    
    if("msr" %in% resid){
      # probability of observed response
      P.response <- ifelse(dat==1, probs, 1-probs)
      
      # modified standardized residual
      out$msr<-(dat-probs)/P.response
    }
    
    # information residual
    if("information" %in% resid){
      if(model=="MIRT"){
        a<-ipars[,1:L]
        d<-ipars[,L+1]
        info<-apply(a%*%t(theta), 2, function(x) x+d)
      }else{
        if(model=="RASCH"){
          ipars<-cbind(1, ipars)
        }
        a<-ipars[,1]
        b<-ipars[,2]
        info<-sapply(theta, function(x) a*(x-b))
      }
      out$information<-t(info)
    }
    
    
  }
  
  if(model %in% c("GRM", "MGRM")){ # polytomous data
    
    if("standardized" %in% resid){
      P<-probs$P
      
      # expected value of the response
      expected.val<- apply(P, 1, function(x) t(x) %*% matrix(1:K, K))
      
      # expected value of the squared response
      expected.val2<- apply(P, 1, function(x) t(x) %*% matrix(1:K, K)^2)
      
      # variance
      var.x<-expected.val2-expected.val^2
      
      # standardized residual
      stz<-(dat-expected.val)/sqrt(var.x)
      out$standardized<-stz
    }
    
    if("msr" %in% resid){
      # probability of observed response
      if(N==1){
        P.response <- P[cbind(1:J, c(dat))]
      }else{
        P.response <- t(sapply(seq_len(N), function(x) {P[cbind(seq_len(J), dat[x,], x)]}))
      }
      # modified standardized residual
      out$msr<-(dat-expected.val)/P.response
    }
  }
    
  return(out)
}   


#' Bisquare Weighting Function
#'
#' Calculate Tukey's bisquare weight (Mosteller & Tukey, 1977) given a residual and bisquare tuning parameter.
#' @param r A residual that measures the inconsistency of a response from the subject's assumed response model, on one item. Residuals of value NA are given a weight of 0.
#' @param B Bisquare tuning parameter. Larger values lead to less downweighting. For robust estimation, \code{B} is often set to 4.0 (Filonczuk & Cheng, 2025; Schuster & Yuan, 2011).
#' @references Filonczuk, A., & Cheng, Y. (2025). Robust estimation of the latent trait in graded response models. Behavior Research Methods, 57(1), 55. https://doi.org/10.3758/s13428-024-02574-2
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' @references Schuster, C., & Yuan, K.-H. (2011). Robust Estimation of Latent Ability in Item Response Models. \emph{Journal of Educational and Behavioral Statistics}, 36(6), 720–735. https://doi.org/10.3102/1076998610396890
#' @return Bisquare weight value. 
#' @examples 
#' 
#' # 1-person case
#' r <- c(-2, -1, 0, 1, 3)
#' B <- 4
#' bisquare(r, B)
#' 
#' # multi-person, multi-item case
#' r <- matrix(c( -2, -1, 0, 1, 
#'                1, 0.5, -0.2, 3, 
#'                0, 2.5, -3, 1 ), 
#'                nrow = 3, byrow = TRUE) 
#' B <- 4 
#' bisquare(r, B)
#' @export
bisquare<-function(r, B){
  w<-ifelse(is.nan(r), 0, 
            ifelse(abs(r) <= B, (1-(r/B)^2)^2, 0))
  
  return(w)
}

#' Huber Weighting Function
#'
#' Calculate the Huber weight (Huber, 1981) given a residual and Huber tuning parameter.
#' @param r A residual that measures the inconsistency of a response from the subject's assumed response model, on one item. Residuals of value NA are given a weight of 0.
#' @param H Huber tuning parameter. Larger values lead to less downweighting. For robust estimation, \code{H} is often set to 1.0 (Filonczuk & Cheng, 2025; Schuster & Yuan, 2011).
#' @references Filonczuk, A., & Cheng, Y. (2025). Robust estimation of the latent trait in graded response models. Behavior Research Methods, 57(1), 55. https://doi.org/10.3758/s13428-024-02574-2
#' @references Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250.
#' @references Schuster, C., & Yuan, K.-H. (2011). Robust Estimation of Latent Ability in Item Response Models. \emph{Journal of Educational and Behavioral Statistics}, 36(6), 720–735. https://doi.org/10.3102/1076998610396890
#' @return Huber weight value.
#' @examples 
#'
#' # 1-person case
#' r <- c(-3, -1, 0, 1.5, 4) 
#' H <- 1.5 
#' huber(r, H)
#' 
#' # multi-person, multi-item case
#' r <- matrix(c( -2, -0.5, 0, 2, 
#'                 3, 1.2, -1, 0.3 ), 
#'                 nrow = 2, byrow = TRUE) 
#' H <- 1 
#' huber(r, H)
#' @export

huber<-function(r, H){
  w<-ifelse(is.nan(r), 0, 
            ifelse(abs(r) <= H, 1, H/abs(r)))
  return(w)
}

#' Generates simulated item responses from model-implied response probabilities.
#'
#' Dichotomous items are sampled using Bernoulli trials, while polytomous (Likert-type) items are sampled using multinomial draws over ordered response categories. The function accepts three probability formats,
#' each corresponding to a different data‑generation scenario. 
#'
#' @param P A matrix or array of response probabilities. For polytomous data, each probability vector must sum to 1.
#' @param anchor Integer specifying the lowest category value. Typical values are 0 or 1; default is 0.
#' @param polytomous Is the data Likert-type? Must be specified as `TRUE` when P is a matrix of category response probabilities for polytomous data, containing data for one person. 
#' @details Let \eqn{N} denotes the number of persons, \eqn{J} the number of items, and \eqn{K} the number of response categories. The response probability input should be specified according to the corresponding scenario:
#' \itemize{
#'   \item Dichotomous data (\eqn{N \times J} matrix): Each entry \eqn{P[n, j]} gives the probability that person \eqn{n} answers item \eqn{j} correctly. Output is an \eqn{N \times J} matrix of 0/1 responses.
#'   \item Polytomous data for one person (\eqn{J \times K} matrix): Each row contains the category probabilities for one item, where \eqn{P[j, k]} is the probability of responding in category \eqn{k}. Output is a \eqn{1 \times J} matrix of integer category scores.
#'   \item Polytomous data for multiple persons (\eqn{J \times K \times N} array): Each slice \eqn{P[ , , n]} is a \eqn{J \times K} matrix of category response probabilities for person \eqn{n}. Output is an \eqn{N \times J} matrix of simulated category scores.
#' }
#' @return A matrix of simulated item responses, dependent on item response type.
#' \itemize{
#'   \item \code{Dichotomous}: \eqn{N \times J} matrix of 0/1 responses.
#'   \item \code{Polytomous}: \eqn{N \times J} matrix of integer category scores beginning at \code{anchor}.
#' }
#' @examples 
#' 
#' # Dichotomous Case with Bernoulli Sampling
#' # 5x4 matrix (5 persons, 4 items)
#' P_matrix <- matrix(c(0.2, 0.5, 0.8, 0.9,
#'                   0.1, 0.4, 0.7, 0.6,
#'                   0.3, 0.6, 0.9, 0.2,
#'                   0.5, 0.5, 0.5, 0.5,
#'                   0.9, 0.8, 0.4, 0.3),
#'                 nrow = 5, byrow = TRUE)
#'
#' dat.gen(P_matrix)
#' 
#' # Polytomous Case for One Person
#' #3x5 matrix of category probabilities (3 items, 5 response categories per item)
#' P_matrix <- matrix(c(0.05, 0.10, 0.20, 0.30, 0.35,
#'                             0.40, 0.30, 0.20, 0.05, 0.05,
#'                             0.10, 0.20, 0.40, 0.20, 0.10),
#'                           nrow = 3, byrow = TRUE)
#'
#' dat.gen(P_matrix, polytomous = TRUE, anchor = 0)
#' 
#' # Polytomous Case For Multiple Persons
#' # 4x4x10 array (4 items, 4 response categories, 10 persons)
#' P_array <- array(runif(4 * 4 * 10), dim = c(4, 4, 10))
#'
#' # Normalize each item–person probability vector
#' P_array <- apply(P_array, c(1, 3), function(x) x / sum(x))
#' P_array <- array(P_array, dim = c(4, 4, 10)) # restore array shape
#'
#' dat.gen(P_array, anchor = 0)
#'
#' @export
                          
dat.gen<-function(P, anchor = 0, polytomous = FALSE, seed=NULL){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(length(dim(P))==3){
    # If dealing with array of polytomous category probabilities
    out <- t(apply(P, c(1,3), function(p) sample(1:length(p), size = 1, prob = p))) - (1-anchor)
    
  }else if(polytomous==T){
    # If dealing with matrix of polytomous category probabilities for one subject
    out <- t(apply(P, 1, function(p) sample(1:length(p), size = 1, prob = p))) - (1-anchor)
    
  }else if(polytomous==F){
    # If dealing with matrix or vector of item success probabilities on dichotomous items
    U<-matrix(runif(length(P)), ncol = ncol(P), nrow = nrow(P))
    out <- ifelse(P>U, 1, 0)
  }
  return(out)
}

#' Standard error function
#' 
#' Computes standard errors accommodating robust procedures in multiple estimation frameworks.
#' Supported standard error types include:
#' \itemize{
#'   \item \code{Asymptotic SE} Information-based standard error incorporating weights from the robust estimation (Magis, 2014). Reduces to the expected Fisher information standard error when item weights are 1 or uniform across items.
#'   \item \code{Sandwich SE} The Fisher information-based standard error is weighted by a correction term to accounts for model misspecification (Huber, 1967; White, 1980).
#'   \item \code{Bayesian Posterior SD} Posterior standard deviation of the Bayesian estimate. The function currently only supports Bayesian standard deviations for the 2PL and GRM.
#'   \item \code{Bayesian Sandwich SD} The posterior standard deviation of the Bayesian estimate is weighted by a correction term to account for model misspecification (Li & Rice, 2023). The function currently only supports Bayesian standard errors for the 2PL and GRM.
#' }
#' The function accommodates robust weighting schemes (equal, Huber, bisquare) 
#' and supports MLE, MAP, and EAP estimation.
#' @param theta A numeric vector or matrix of latent trait values. For unidimensional models, a numeric vector of length \eqn{N}. For multidimensional models (MIRT, MGRM), an \eqn{N \times L} matrix.
#' @param ipars A matrix of item parameters, whose structure depends on the model. 
#' @param dat A \eqn{N \times J} matrix of polytomously-scored data (e.g., Likert-type) for \emph{N} subjects and \emph{J} items. Indexing begins at 0.
#' @param model Character string specifying the IRT model. See `item.prob` for supported models.
#' @param D A scaling constant. Defaults to 1.7; alternatively is often set to 1.0.
#' @param weight.type Type of weighting function to be used: "equal", "Huber", or "bisquare".
#' @param tuning.par Tuning parameter to be used with Huber or bisquare weights.
#' @param est.type Type of estimation to be used: "MLE", "MAP", or "EAP".
#' @param prior Numeric vector giving the mean and variance of the normal prior for MAP and EAP estimation types.
#' @param eap.quad.pts A numeric vector of quadrature points \eqn{\theta_q} used when computing EAP standard errors. Default is 41 equally spaced values across latent traits ranging from -4.0 to 4.0. Only used when "EAP" is included in \code{est.type}.
#' @details 
#' The function computes person-level standard errors for several IRT models by combining model-specific score functions with robust weighting and multiple
#' estimation frameworks (MLE, MAP, EAP). All SEs are derived from first- and second-order derivatives of the log-likelihood (or posterior) evaluated at each person's latent trait estimate.
#' 
#' For every person and item, the function first computes standardized residuals
#' and converts them into weights. These weights can be:
#' \itemize{
#'   \item \code{"equal"} — all responses weighted equally,
#'   \item \code{"Huber"} — down-weights large residuals,
#'   \item \code{"bisquare"} — strongly down-weights outliers.
#' }
#' The weights influence both the information-based SEs and the sandwich SEs.
#'
#' \strong{MLE standard errors.}  
#' For MLE, the function uses each model’s first and second derivatives to form:
#' \itemize{
#'   \item an expected information term (how much information the items provide),
#'   \item an empirical “sandwich” term (how variable the score function is).
#' }
#' The asymptotic SE uses only the information term, while the sandwich SE uses
#' both and is more robust to model–data misfit.
#'
#' \strong{MAP standard errors.}  
#' For MAP, the normal prior adds extra curvature to the information. This makes
#' the posterior SD smaller when the prior is strong. A Bayesian sandwich SE is
#' also returned by combining the posterior SD with the empirical variability.
#'
#' \strong{EAP standard errors.}  
#' For EAP, the function evaluates the likelihood at specific quadrature points
#' and computes the posterior variance directly from the weighted distribution.
#' A Bayesian sandwich SE is again produced by combining this posterior SD with
#' the empirical term.
#'
#' @references Huber, P. J. (1967). The Behavior of Maximum Likelihood Estimates Under Nonstandard Conditions. Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, 1, 221–233.
#' @references Li, K. Q., & Rice, K. M. (2023). A Bayesian “sandwich” for variance estimation (arXiv:2207.00100). arXiv. https://doi.org/10.48550/arXiv.2207.00100
#' @references Magis, D. (2014). On the asymptotic standard error of a class of robust estimators of ability in dichotomous item response models. British Journal of Mathematical and Statistical Psychology, 67(3), 430–450. https://doi.org/10.1111/bmsp.12027
#' @references White, H. (1980). A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity. Econometrica, 48(4), 817. https://doi.org/10.2307/1912934
#' @return A list whose elements depend on the specified model and `est.type`. Possible components include:
#' \itemize{
#'    \item \code{"asymptotic_MLE"} (information‑based SEs)
#'    \item \code{"sandwich_MLE"} (sandwich SEs)
#'    \item \code{"post_sd_MAP"} (posterior SD for MAP)
#'    \item \code{"sandwich_MAP"} (Bayesian sandwich SE for MAP)
#'    \item \code{"post_sd_EAP"} (posterior SD for EAP)
#'    \item \code{"sandwich_EAP"} (Bayesian sandwich SE for EAP)
#'    \item \code{"singular.matrix"} (indicator for singular information matrices in multidimensional data)
#' }    
#' @export
#' @examples 
#' library(mirt)
#' library(dplyr)
#' library(readr)
#' data(SAT12)
#' 
#' itemstats(SAT12, use_ts = FALSE)
#' # score the data (missing scored as 0)
#' head(SAT12)
#' dat <- key2binary(SAT12,
#'                  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' head(dat)
#' itemstats(dat)
#'
#' # score the data, missing (value of 8) treated as NA
#' SAT12missing <- SAT12
#' SAT12missing[SAT12missing == 8] <- NA
#' dat <- key2binary(SAT12missing,
#'                  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' head(dat)
#'
#' # potentially better scoring for item 32 (based on nominal model finding)
#' dat <- key2binary(SAT12,
#'                  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,3))
#'
#' # Rasch Model Example
#' Rasch_fit <- mirt(dat, 1, itemtype = "Rasch")
#' Rasch_theta <- fscores(Rasch_fit)
#' Rasch_ipars <- coef(Rasch_fit, IRTpars = TRUE, simplify = TRUE)$items[, "b", drop = FALSE]
#' Rasch_SE <- standard.errors(Rasch_theta, Rasch_ipars, dat, model = "Rasch")
#' 
#' # 1PL Model Example
#' fit_1PL <- mirt(dat, 1, itemtype = "1PL")
#' theta_1PL <- fscores(fit_1PL)
#' ipars_1PL <- cbind(a = rep(1, ncol(dat)), b = coef(fit_1PL, IRTpars = TRUE, simplify = TRUE)$items[, "b"])
#' SE_1PL <- standard.errors(theta_1PL, ipars_1PL, dat, model="1PL")
#' 
#' # 2PL Model Example
#' fit_2PL <- mirt(dat, 1, itemtype = "2PL")
#' theta_2PL <- fscores(fit_2PL)
#' ipars_2PL <- coef(fit_2PL, IRTpars = TRUE, simplify = TRUE)$items[, c("a1", "b")]
#' SE_2PL <- standard.errors(theta_2PL, ipars_2PL, dat, model = "2PL", est.type = c("MLE", "MAP"))
#' 
#' # MIRT Model Example
#' MIRT_fit <- mirt(dat, 2)
#' MIRT_theta <- fscores(MIRT_fit)
#' MIRT_ipars <- coef(MIRT_fit, IRTpars = TRUE, simplify = TRUE)$items
#' MIRT_SE <- standard.errors(MIRT_theta, MIRT_ipars, dat, model = "MIRT")
#' 
#' # GRM Model Example for Agreeableness
#' bfi <- read_csv("p234_bfi_demog_2020-04-24.csv")
#' A_items <- bfi %>% dplyr::select(starts_with("A"))
#' GRM_fit <- mirt(A_items, 1, itemtype = "graded")
#' GRM_theta <- fscores(GRM_fit)
#' GRM_ipars <- coef(GRM_fit, IRTpars = TRUE, simplify = TRUE)$items
#' GRM_SE <- standard.errors(theta = GRM_theta, ipars = GRM_ipars, dat = as.matrix(A_items), model = "GRM", est.type = c("MLE", "MAP", "EAP"), weight.type = "Huber", tuning.par = 1)
#' 
#' # MGRM Model Example
#' library(readr)
#' library(dplyr)
#' library(mirt)
#' bfi <- read_csv("p234_bfi_demog_2020-04-24.csv")
#' MGRM_items <- bfi %>% select(starts_with("t1_bfi"))
#' names(MGRM_items) <- gsub("t1_bfi_", "", names(MGRM_items))
#' names(MGRM_items)<- gsub("R", "", names(MGRM_items))
#'
#' MGRM_model <- "A = A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12
#'                C = C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12
#'                E = E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12
#'                N = N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12
#'                O = O1, O2, O3, O4, O5, O6, O7, O8, O9, O10, O11, O12"
#' MGRM_fit <- mirt(MGRM_items, model = MGRM_model, itemtype = "graded", method = "QMCEM", technical = list(NCYCLES = 500))
#' MGRM_theta <- fscores(MGRM_fit)
#' MGRM_ipars <- coef(MGRM_fit, IRTpars = TRUE, simplify = TRUE)$items
#' MGRM_dat <- as.matrix(MGRM_items)
#' MGRM_SE <- standard.errors(theta = MGRM_theta, ipars = MGRM_ipars, dat = MGRM_dat, model = "MGRM", est.type = "MLE", weight.type = "Huber", tuning.par  = 1)

                   
standard.errors<-function(theta, ipars, dat, model, D=1.7, weight.type = "equal", tuning.par = NULL, est.type = "MLE", prior=c(0,1), eap.quad.pts =seq(-4, 4, length.out = 41)){
  
  # Function to determine weights:
  weight.func <- function(r, weight.type2="equal", tuning.par2=NULL){
    if(weight.type2 == "equal") return(r*0+1)
    if(weight.type2 == "Huber") return(huber(r, tuning.par2))
    if(weight.type2 == "bisquare") return(bisquare(r, tuning.par2))
    if(length(weight.type2) == length(r)) return(weight.type2)
    stop("Invalid weight.type")
  }
  
  model<-toupper(model)
  
  N<-nrow(dat) # number of subjects
  J<-ncol(dat) # test length
  
  P<-item.prob(theta, model, ipars, D=D)
  r<-residual(theta, model, ipars, dat, resid = resid, D=D)
  w<- weight.func(r, weight.type, tuning.par)
  
  out<-list() # initialize vector for storing output
  
  if(model=="RASCH"){
    
    # first derivative
    D1<-t(sapply(1:N, function(x) D*(dat[x,]-P[x,])))
    
    # second derivative
    D2<-t(sapply(1:N, function(x) D^2*P[x,]*(1-P[x,])))
    
    A<- rowSums(w*D2) # expected info
    V<-rowSums(w^2*D2) # ASE numerator
    B <- rowSums((w*D1)^2) # sandwich B
    
    if("MLE"%in%est.type){
      out$asymptotic_MLE <- sqrt(V)/A
      out$sandwich_MLE <- sqrt(B)/A
    }
    
    if("MAP"%in%est.type){
      
      # Weighted Second Derivative
      D2.MAP<-rowSums(w*D2) + 1/prior[2]
      
      out$post_sd_MAP <- sqrt(1 / D2.MAP)
      out$sandwich_MAP <-out$post_sd_MAP*sqrt(B / A)
    }
    
    if("EAP"%in%est.type){
      # Prior density according to each density point
      f_x<-dnorm(eap.quad.pts, prior[1], sqrt(prior[2]))
      # Item response probabilities at each quadrature point
      probs.q <- item.prob(eap.quad.pts, "Rasch", ipars) 
      Q<-length(eap.quad.pts)
      
      # (Weighted) Likelihood according to person's data and each quad point
      likelihood<-matrix(NA, N, Q)
      for(i in 1:N){
        likelihood[i,]<- apply(probs.q, 1, function(x) prod((x^dat[i,]*(1-x)^(1-dat[i,]))^w[i,]))
      }
      
      # Posterior standard deviation of EAP theta estimate according to Bock & Mislevy, 1989; Thissen et al., 1995
      out$post_sd_EAP <-sapply(1:N, function(x) sqrt(sum((eap.quad.pts - theta[x])^2*likelihood[x,]*f_x)/ sum(likelihood[x,]*f_x)))
      out$sandwich_EAP <-out$post_sd_EAP*sqrt((B / A))
    } 

    return(out)
  }
  
  if(model=="1PL"){
    
    if(ncol(ipars)==2){
      a<-ipars[,1] # discrimination parameter
    }else{
      a<-rep(1, J)
    }
    
    # first derivative
    D1<-t(sapply(1:N, function(x) D*a*(dat[x,]-P[x,])))
    
    # second derivative
    D2<-t(sapply(1:N, function(x) D^2*a^2*P[x,]*(1-P[x,])))
    
    A<- rowSums(w*D2) # expected info
    V<-rowSums(w^2*D2) # ASE numerator
    B <- rowSums((w*D1)^2) # sandwich B
    
    
    if("MLE"%in%est.type){
      out$asymptotic_MLE <- sqrt(V)/A
      out$sandwich_MLE <- sqrt(B)/A
    }
    
    if("MAP"%in%est.type){
      
      # Weighted Second Derivative
      D2.MAP<-rowSums(w*D2) + 1/prior[2]
      
      out$post_sd_MAP <- sqrt(1 / D2.MAP)
      out$sandwich_MAP <-out$post_sd_MAP*sqrt(B / A)
    }
    
    if("EAP"%in%est.type){
      # Prior density according to each density point
      f_x<-dnorm(eap.quad.pts, prior[1], sqrt(prior[2]))
      # Item response probabilities at each quadrature point
      probs.q <- item.prob(eap.quad.pts, "1PL", ipars) 
      Q<-length(eap.quad.pts)
      
      # (Weighted) Likelihood according to person's data and each quad point
      likelihood<-matrix(NA, N, Q)
      for(i in 1:N){
        likelihood[i,]<- apply(probs.q, 1, function(x) prod((x^dat[i,]*(1-x)^(1-dat[i,]))^w[i,]))
      }
      
      # Posterior standard deviation of EAP theta estimate according to Bock & Mislevy, 1989; Thissen et al., 1995
      out$post_sd_EAP <-sapply(1:N, function(x) sqrt(sum((eap.quad.pts - theta[x])^2*likelihood[x,]*f_x)/ sum(likelihood[x,]*f_x)))
      out$sandwich_EAP <-out$post_sd_EAP*sqrt((B / A))
    } 
    
  }
  
  if(model=="2PL"){
    
    a<-ipars[,1] # discrimination parameter
    
    # first derivative
    D1<-t(sapply(1:N, function(x) D*a*(dat[x,]-P[x,])))
    
    # second derivative
    D2<-t(sapply(1:N, function(x) D^2*a^2*P[x,]*(1-P[x,])))
    
    A<- rowSums(w*D2) # expected info
    V<-rowSums(w^2*D2) # ASE numerator
    B <- rowSums((w*D1)^2) # sandwich B
    
    
    if("MLE"%in%est.type){
      out$asymptotic_MLE <- sqrt(V)/A
      out$sandwich_MLE <- sqrt(B)/A
    }
    
    if("MAP"%in%est.type){
      
      # Weighted Second Derivative
      D2.MAP<-rowSums(w*D2) + 1/prior[2]
      
      # Posterior standard deviation of MAP theta estimate
      out$post_sd_MAP <- sqrt(1 / D2.MAP)
      out$sandwich_MAP <-out$post_sd_MAP*sqrt(B / A)
    }
    
    if("EAP"%in%est.type){
      # Prior density according to each density point
      f_x<-dnorm(eap.quad.pts, prior[1], sqrt(prior[2]))
      # Item response probabilities at each quadrature point
      probs.q <- item.prob(eap.quad.pts, "2PL", ipars) 
      Q<-length(eap.quad.pts)
      
      # (Weighted) Likelihood according to person's data and each quad point
      likelihood<-matrix(NA, N, Q)
      for(i in 1:N){
        likelihood[i,]<- apply(probs.q, 1, function(x) prod((x^dat[i,]*(1-x)^(1-dat[i,]))^w[i,]))
      }
      
      # Posterior standard deviation of EAP theta estimate according to Bock & Mislevy, 1989; Thissen et al., 1995
      out$post_sd_EAP <-sapply(1:N, function(x) sqrt(sum((eap.quad.pts - theta[x])^2*likelihood[x,]*f_x)/ sum(likelihood[x,]*f_x)))
      out$sandwich_EAP <-out$post_sd_EAP*sqrt((B / A))
    } 
    
  }
  
  if(model=="MIRT"){
    if(!is.matrix(theta)) theta <- matrix(theta, nrow = 1)
    L<-ncol(theta)
    N<-nrow(theta)
    a<-ipars[,1:L]
    d<-ipars[,L+1]
    
    out_ase <- matrix(NA,N,L)
    out_sand  <- matrix(NA,N,L)
    out_singular <- matrix(0,N,L)
    
    for(i in 1:N){
      Pi<-P[i,]
      xi<-dat[i,]
      
      # Initialize LxL matrices
      A <- V <-B <- matrix(0,L,L)
      
      for(j in 1:nrow(a)){
        aj <- matrix(a[j, ], nrow = 1)
        wij<-w[i,j]
        gj <- matrix(D*wij*aj* (xi[j] - Pi[j]) , nrow=1) # jth contribution to the first derivative of the log likelihood
        B <- B + t(gj) %*% gj # sandwich B
        
        Ij<-D^2*(t(aj)%*% aj)*Pi[j]*(1-Pi[j]) # jth item information (2nd derivative)
        
        A <- A +wij*Ij # jth item expected info contribution
        V <- V + wij^2*Ij # ASE numerator contribution
        
      }
      
      if(length(A)==1){
        Ainv<-1/A
      } else if(det(A)<1e-12 || any(!is.finite(A))) {
        out_singular[i,]<-1
        out_ase[i,] <- NA
        out_sand[i,]  <- NA
      }else {
        Ainv <- solve(A)
        out_ase[i,] <- sqrt(diag(Ainv %*% V %*% Ainv))
        out_sand[i,]  <- sqrt(diag(Ainv %*% B %*% Ainv))
      }
      
    }
    
    out<-list(asymptotic_MLE = out_ase,
              sandwich_MLE = out_sand,
              singular.matrix = out_singular)
    
  }
  
  if(model=="GRM"){
    a<-ipars[,1] # discrimination parameters
    b<-ipars[,-1] # threshold parameters
    K<-ncol(b) # number of thresholds
    
    Pcat  <- P$P             
    Pstar <- array(NA, dim=c(J, K+2, N))
    Pstar[,1,] <- 1
    Pstar[,2:(K+1),] <- P$pstar
    Pstar[,K+2,] <- 0
    
    scores <- 1:K
    
    k_index <- t(dat)
    
    idx1 <- cbind(rep(1:J, N), as.vector(k_index), rep(1:N, each=J))
    idx2 <- cbind(rep(1:J, N), as.vector(k_index+1), rep(1:N, each=J))
    
    # Add boundary threshold probabilities
    P1 <- array(Pstar[idx1], dim=c(J,N))
    P2 <- array(Pstar[idx2], dim=c(J,N))
    Pk <- P1 - P2
    
    # First derivative
    D1 <- D*a*(P1*(1-P1) - P2*(1-P2)) / Pk
    
    # Second Derivative 
    D2 <- -D^2 * a^2 * ((P1*(1-P1)*(1-2*P1) - P2*(1-P2)*(1-2*P2))/Pk -
                          (P1*(1-P1) - P2*(1-P2))^2 / Pk^2)
    
    A <- colSums(t(w) * D2) # Expected information
    V <- colSums(t(w)^2 * D2) # ASE numerator
    B <- colSums((t(w) * D1)^2) # Sandwich B
    
    A[A <= 0 | !is.finite(A)] <- NA
    
    if("MLE"%in%est.type){
      out$asymptotic_MLE <- sqrt(V)/A
      out$sandwich_MLE <- sqrt(B)/A
    }
    
    if("MAP"%in%est.type){
      
      # Weighted Second Derivative
      D2.MAP<-A + 1/prior[2]
      
      out$post_sd_MAP <- sqrt(1 / D2.MAP)
      out$sandwich_MAP <-out$post_sd_MAP*sqrt(B / A)
    }
    
    if("EAP"%in%est.type){
      # Prior density according to each density point
      f_x<-dnorm(eap.quad.pts, prior[1], sqrt(prior[2]))
      # Item response probabilities at each quadrature point
      probs.q <- item.prob(eap.quad.pts, "GRM", ipars) 
      Q<-length(eap.quad.pts)
      
      # (Weighted) Likelihood according to person's data and each quad point
      likelihood<-matrix(NA, N, Q)
      for(i in 1:N){
        likelihood[i,]<- apply(probs.q, 1, function(x) prod((x^dat[i,]*(1-x)^(1-dat[i,]))^w[i,]))
      }
      
      # Posterior standard deviation of EAP theta estimate according to Bock & Mislevy, 1989; Thissen et al., 1995
      out$post_sd_EAP <-sapply(1:N, function(x) sqrt(sum((eap.quad.pts - theta[x])^2*likelihood[x,]*f_x)/ sum(likelihood[x,]*f_x)))
      out$sandwich_EAP <-out$post_sd_EAP*sqrt((B / A))
    } 
  }
  
  if(model=="MGRM"){
    if(!is.matrix(theta)) theta <- matrix(theta, nrow = 1)
    L <- ncol(theta)
    a<- ipars[, 1:L] # a: J x L matrix (first L columns: discrimination parameters)
    b<- ipars[, (L+1):ncol(ipars)] # b: J x K matrix (category threshold parameters)
    K <-ncol(b) # number of thresholds K
    
    # Calculate probabilities
    probs<-item.prob(theta, "MGRM", ipars)
    
    out_ase      <- matrix(NA, N, L)
    out_sand     <- matrix(NA, N, L)
    out_singular <- matrix(0,  N, L)
    
    for(i in 1:N){
      
      # Pstar: J x (K+2) boundary probabilities for person i
      Pstar_i <- matrix(NA, J, K + 2)
      Pstar_i[, 1] <- 1
      Pstar_i[, 2:(K+1)] <-probs$pstar[,, i]   # J x K slice
      Pstar_i[, K+2]<- 0
      
      Pcat_i <- probs$P[,,i]   # J x (K+1) category probabilities for person i
      
      # Accumulate A, V, B across items (LxL matrices)
      A_mat <- V_mat <- B_mat <- matrix(0, L, L)
      
      for(j in 1:J){
        aj  <- matrix(a[j, ], nrow = 1)   # 1 x L
        wij <- w[i, j]
        
        # Extract P* values
        ps0 <- Pstar_i[j, 1:(K+1)] 
        ps1 <- Pstar_i[j, 2:(K+2)]
        
        # category probs of length K+1
        Pk  <- Pcat_i[j, ]               
        
        qs0 <- 1 - ps0
        qs1 <- 1 - ps1
        
        # Numerators reused across derivatives
        num1 <- ps0 * qs0 - ps1 * qs1     # length K+1
        
        # First derivative: scalar (sum over categories), then outer to LxL 
        # d/dtheta log L_j = D * a_j . sum_k [ num1_k / Pk_k ]
        score_j <- D * sum(num1 / Pk)     # scalar score contribution
        
        g_j <- wij * t(aj) * score_j      # L x 1 gradient
        B_mat <- B_mat + g_j %*% t(g_j)
        
        # Second derivative: Expected Information (Fisher)
        # Note: the negative sign is absorbed so this is positive information
        num2 <- ps0*qs0*(qs0 - ps0) - ps1*qs1*(qs1 - ps1)  # = ps0*qs0*(1-2ps0) - ps1*qs1*(1-2ps1)
        
        # Expected info scalar (sum over k)
        Ij_scalar <- D^2 * sum(num2 / Pk - (num1^2) / Pk^2)
        # Note: this equals negative expected info for Fisher per item; we negate for positive info
        Ij_scalar <- -Ij_scalar
        
        Ij <- Ij_scalar * (t(aj) %*% aj)  # L x L item information
        
        A_mat <- A_mat + wij  * Ij
        V_mat <- V_mat + wij^2 * Ij
      }
      
      # Invert and compute SEs
      if(any(!is.finite(A_mat)) || any(!is.finite(V_mat))){
        out_singular[i, ] <- 1
        next 
      }
      
      det_A <- if(L > 1) det(A_mat) else A_mat[1,1]
      
      if(abs(det_A) < 1e-12){
        out_singular[i, ] <- 1
      } else {
        Ainv         <- if(L == 1) 1/A_mat else solve(A_mat)
        out_ase[i, ] <- sqrt(diag(Ainv %*% V_mat %*% Ainv))
        out_sand[i,] <- sqrt(diag(Ainv %*% B_mat %*% Ainv))
      }
    }
    
    out <- list(
      asymptotic_MLE   = out_ase,
      sandwich_MLE     = out_sand,
      singular.matrix = out_singular
    )
  }
  
  return(out)
}

#' Ability Estimation Function Using Robust Estimation
#' @param theta A numeric vector or matrix of latent trait values. For unidimensional models, a numeric vector of length \eqn{N}. For multidimensional models (MIRT, MGRM), an \eqn{N \times L} matrix.
#' @param ipars A matrix of item parameters, whose structure depends on the model. 
#' @param dat A \eqn{N \times J} matrix of polytomously-scored data (e.g., Likert-type) for \emph{N} subjects and \emph{J} items. Indexing begins at 0.
#' @param model Character string specifying the IRT model. See `item.prob` for supported models.
#' @param D A positive scaling constant used for scaling the normal ogive model. Defaults to 1.7; alternatively is often set to 1.0.
#' @param weight.type Type of weighting function to be used: "equal", "Huber", or "bisquare".
#' @param tuning.par Optional tuning parameter for Huber or bisquare weights.
#' @param est.type Type of estimation to be used: "MLE", "MAP", or "EAP".
#' @param init.val Initial value of the ability parameter for the Newton-Raphson method. Default is 0. The MLE is suggested as the initial value for robust estimation.
#' @param prior Numeric vector giving the mean and variance of the normal prior for MAP and EAP estimation types.
#' @param eap.quad.pts A numeric vector of quadrature points \eqn{\theta_q} used when computing EAP standard errors. These are typically q equally spaced values across the latent‑trait range. Only used when "EAP" is included in \code{est.type}.
#' @details The goal of robust estimation is to downweigh potentially aberrant responses to lessen their impact on the estimation of \eqn{\theta_i}. Robust estimates resist the harmful effects of response disturbances and tend to be less biased estimates of true ability than maximum likelihood estimates.
#'               Under a given item response model, let the probability of examinee \eqn{i} endorsing exactly category \eqn{k} on item \eqn{j} be denoted \eqn{P_{ijk}}.
#'               The contribution of item \emph{j} to the overall log-likelihood for one subject is weighted with a weight \eqn{\omega(r_{ij})} as a function of a residual \eqn{r_{ij}} for the item:
#'               \deqn{\sum^J_{j=1} \omega(r_{ij}) \sum^K_{k=1} u_{ijk}\text{log}P_{ijk} = 0 }
#'               \eqn{u_{jk}} is an indicator function: \deqn{u_{ijk} = \begin{cases}
#'                                                            1 & \text{if } X_{ij} = k; \\
#'                                                            0 & \text{otherwise}.
#'                                                            \end{cases} }
#'               The residual, which measures the inconsistency of a response from the subject's assumed response model, is \deqn{r_{ij} = \frac{1}{\sigma_{X_{ij}}}\left[X_{ij} - E(X_{ij}|\hat{\boldsymbol{\theta}}_i)\right]}.
#'               The difference in fit is determined between the observed response \eqn{X_{ij}} and expected score \eqn{E(X_{ij}|\hat{\boldsymbol{\theta}}_i) = \sum_{k=1}^KkP_{jk}(\hat{\boldsymbol{\theta}}_i)}, and scaled by the variance \eqn{\sigma_{X_{ij}}^2 = \sum_{k=1}^K (X_{ijk}-E[X_{ij}|\hat{\theta}_i])^2P_{jk}(\hat{\theta}_i).}
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
#' @return A list containing the following outputs:
#' \itemize{
#'   \item \code{theta} Ability estimates for \emph{N} subjects. NAs replace values that did not converge to any value. Estimates that converged to values less than -3.0 were replaced with -3.0, while estimates that converged to values greater than 3.0 were replaced with 3.0.
#'   \item \code{convergence} Indicators of convergence for \emph{N} subjects: a “0” indicates the value converged, while a “1” indicates the maximum likelihood estimation did not converge to any value.
#'   \item \code{standard.error} Standard errors of the theta estimates for \emph{N} subjects, given by the square root of the reciprocal of the Fisher information. NAs replace nonconverging values. 
#'   \item \code{theta.progression} A matrix with rows corresponding to each subject and columns corresponding to the number of iterations supplied to the input. Each column provides the updated theta estimate at each iteration of the Newton-Raphson algorithm until the change in log-likelihood for that subject reaches the cutoff value or the value is nonconverged (reaches infinite values).
#'   \item \code{residual} A \eqn{J \times N \times p} array containing residuals corresponding to the ability estimate for \emph{N} subjects respective to the \emph{J} test items at each iteration until convergence within maximum \emph{p} iterations, nonconvergence, or singular matrix is reached.
#' }
#' @export

robust.theta<-function(dat, ipars, model="2PL", D=1.7, dimen=NULL, weight.type = "equal", tuning.par = NULL, resid = "standardized", est.type = "MLE", init.val=0, prior=c(0,1), eap.quad.pts =seq(-4, 4, length.out = 41), iter=30, tol=0.01){
  
  model<-toupper(model)
  
  # Set checks to ensure proper inputs
  if(!(model%in%c("RASCH", "1PL", "2PL", "MIRT", "GRM", "MGRM"))){
    # accepts only one model
    stop(paste(model, " is not an accepted model."))
  }
  if(!(weight.type%in%c("equal", "Huber", "bisquare"))){
    # accepts only one weight type
    stop(paste(weight.type, " is not an accepted weight type."))
  }
  if(any(!(est.type%in%c("MLE", "MAP", "EAP")))){
    # accepts more than one estimation procedure
    stop(paste(est.type, " is not an accepted estimation procedure."))
  }
  if(weight.type != "equal"){
    if(is.null(tuning.par)){
      stop(paste("The turning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  
  # Function to determine weights:
  weight.func <- function(r, weight.type2="equal", tuning.par2=NULL){
    if(weight.type2 == "equal") return(rep(1, length(r)))
    if(weight.type2 == "Huber") return(huber(r, tuning.par2))
    if(weight.type2 == "bisquare") return(bisquare(r, tuning.par2))
    if(length(weight.type2) == length(r)) return(weight.type2)
    stop("Invalid weight.type")
  }
  
  # Function to get initial theta for subject i 
  get_init <- function(i, L){
    if(is.matrix(init.val)){ # if user specifies different initial values for each subject
      return(as.matrix(init.val[i, ]))      
    } else if(length(init.val) == L){ #if user specifies L values, the same initial values are used for each subject
      return(as.matrix(init.val)) 
    } else { #if user specifies one initial value, it will be used for all dimensions and each subject
      return(matrix(rep(init.val[1], L), nrow = L))
    }
  }
  
  # Extract data dimensions
  N<-nrow(dat) # number of subjects
  J<-ncol(dat) # number of items
  
  # for Bayesian estimates: extract priors
  if(any(est.type%in%c("MAP", "EAP"))){
    mu<-prior[1]
    sigma2<-prior[2]
  }
  
  # Initialize Output 
  residual.out <- matrix(NA, nrow = N, ncol = J)
  
  if(model%in%c("RASCH", "1PL", "2PL")){
    L<-1
    # Initialize Output 
    theta.prog <- array(NA, dim = c(N, iter, L))
    theta.out <- matrix(NA, nrow = N, ncol = L)
    se.out <- matrix(NA, nrow = N, ncol = L)
    convergence.out <- matrix(0, nrow = N, ncol = L)
    
    if(model=="RASCH"){
      ipars<- cbind(rep(1, J), ipars)
    }
    
    if("EAP"%in%est.type){
      f_x <- dnorm(eap.quad.pts, mu, sqrt(sigma2))
      probs.q <- item.prob(eap.quad.pts, model, ipars, D)
      Q <- length(eap.quad.pts)
      
      # Loop to estimate theta for each subject
      for(i in 1:N){
        P0<-0
        # Initialize theta value for residual
        theta<-get_init(i, L=1)
        
        for(k in 1:iter){
          # Residual Calculation
          r<-residual(theta, model, ipars, dat[i,], resid=resid, D=1.7)
          # Weight Calculation
          w.i<-weight.func(r=r, weight.type2=weight.type, tuning.par2=tuning.par)
          
          # (Weighted) Likelihood according to person's data and each quad point
          likelihood <- apply(probs.q, 1, function(x) prod((x^dat[i, ]*(1 - x)^(1 - dat[i, ]))^w.i))
          
          # EAP theta estimate
          theta<-sum(eap.quad.pts*likelihood*f_x)/sum(likelihood*f_x)
          
          # Store final estimated theta
          theta.new <- theta.prog[i,k,1]<-theta
          
          # Compute difference in log-likelihood for convergence criterion
          probs.i <- item.prob(theta.new, model, ipars, D)
          probs.resp <- ifelse(dat[i, ]==1, probs.i, 1-probs.i)
          log_like <- sum(log(probs.resp)) - sum(log(P0))
          
          P0 <- probs.resp
          theta <- theta.new
          
          if(k > 1 && abs(log_like) < tol){break}
        }
        
      }
    }else{
      
      
      # Newton-Raphson loop for MLE and MAP
      for(i in 1:N){
        
        theta <- get_init(i, 1)
        P0 <- 0
        
        for(k in 1:iter){
          
          #  Compute item (category) response probabilities 
          probs.i <- item.prob(theta, model, ipars, D)
          
        }
      }
    }
  }
  
  if(model%in%c("GRM")){
    
  }
  
  if(model%in%c("MIRT")){
    
  }
  
  if(model%in%c("MGRM")){
    
  }
  
}

#' Ability Estimation Function Using Robust Estimation (GRM)
#'
#' Calculate robust ability estimates using the GRM item response function with the given weight function, fixed item parameters, and item responses.
#' @param dat A \eqn{J \times N} matrix of polytomously-scored data (e.g., Likert-type) for \emph{J} items and \emph{N} subjects.
#' @param a Vector of slope parameters for \emph{J} items.
#' @param b A \eqn{J \times (K-1)} matrix of category threshold parameters for \emph{K} categories.
#' @param iter Max number of iterations. Default is 100.
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
#' @return A list containing the following outputs:
#' \itemize{
#'   \item \code{theta} Ability estimates for \emph{N} subjects. NAs replace values that did not converge to any value. Estimates that converged to values less than -3.0 were replaced with -3.0, while estimates that converged to values greater than 3.0 were replaced with 3.0.
#'   \item \code{convergence} Indicators of convergence for \emph{N} subjects: a “0” indicates the value converged, while a “1” indicates the maximum likelihood estimation did not converge to any value.
#'   \item \code{standard.error} Standard errors of the theta estimates for \emph{N} subjects, given by the square root of the reciprocal of the Fisher information. NAs replace nonconverging values. 
#'   \item \code{theta.progression} A matrix with rows corresponding to each subject and columns corresponding to the number of iterations supplied to the input. Each column provides the updated theta estimate at each iteration of the Newton-Raphson algorithm until the change in log-likelihood for that subject reaches the cutoff value or the value is nonconverged (reaches infinite values).
#'   \item \code{residual} A \eqn{J \times N \times p} array containing residuals corresponding to the ability estimate for \emph{N} subjects respective to the \emph{J} test items at each iteration until convergence within maximum \emph{p} iterations, nonconvergence, or singular matrix is reached.
#' }
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
#' dat<-dat.gen(probs$P, anchor=1)
#' 
#' # Make the data aberrant by reverse coding 20% items
#' ab.prop<-0.2
#' index<-sample(c(1:n), ab.prop*n)
#' ab.dat<-dat
#' ab.dat[, index]<-apply(matrix(dat[,index]), c(1,2), function(x) return(nthresh+2-x))
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
  resid <- matrix(data=NA, nrow = l, ncol = J)
  
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
      P.i<-probs$P
      probs.resp<-P.i[cbind(1:J, dat[i,])]
      
      expected.value<-P.i%*%matrix(c(1:(nthresh+1))) 
      # Calculate standardized residual
      resid[i,]<-(dat[i,]-expected.value)/sqrt(rowSums(apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*P.i))  
      
      # Compute weighting term based on specified weight function (bisquare, Huber, equal)
      weighting.term <- NULL
      if (weight.type == "bisquare") {
        weighting.term <- bisquare(resid[i,], tuning.par)
      } else if (weight.type == "Huber") {
        weighting.term <- huber(resid[i,], tuning.par)
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
    
    #se.int<-standard.errors(theta, cbind(a,b), dat[i,], "GRM",  weight.type = weight.type, tuning.par=tuning.par)
    #standard.error[i] <- se.int$observed_se
    #resid[i,]<-se.int$residual
    
    # Handle cases where theta did not converge within the desired number of iterations
    if (k == iter) {
      theta.est2[i] <-  NA
      convergence[i,1] <- 1
    } else if (!is.na(theta) & theta < -3) {
      # Then check: if theta converged outside [-3, 3], replace it with -3 or 3 respectively
      theta <- -3
      
      #se.int<-standard.errors(a, theta, dat[i,], b=b, residual = "standardized", weight.function = weight.type, tuning.par=tuning.par, model="GRM",  D=D)
      #standard.error[i] <- se.int$observed_se
      #residual[i,]<-se.int$residual
      theta.est2[i] <- theta
    } else if (!is.na(theta) & theta > 3) {
      theta <- 3
      #se.int<-standard.errors(a, theta, dat[i,], b=b, residual = "standardized", weight.function = weight.type, tuning.par=tuning.par, model="GRM",  D=D)
      #standard.error[i] <- se.int$observed_se
      #residual[i,]<-se.int$residual
      theta.est2[i] <- theta
    }
  }
  
  # Return a list containing the estimated theta, binary indicator of nonconvergence, standard error, estimated theta over each iteration, and standardized residual
  return(list(theta = theta.est2, convergence = convergence, standard.error = standard.error, theta.progression = theta.progression, residual = residual))
}
#' Robust Estimation of Item Parameters
#' 
#' @param dat A \eqn{N\times J} matrix of dichotomous response data
#' 
#' @references Hong, M., & Cheng, Y. (2019). Robust maximum marginal likelihood (RMML) estimation for item response theory models. Behavior Research Methods, 51(2), 573–588. https://doi.org/10.3758/s13428-018-1150-4
#' @export
#' @examples
#' # Load package and example data set 
#' library(mirt) 
#' data(Science) 
#' # Robust estimation of item parameters
#' robust.item(Science)
robust.item<-function(dat){
  
  # Initial Model Estimation 
  mod <- mirt(dat, 1, optimizer = 'NR') 
  # Person fit residual calculation
  per.fit <- personfit(mod, method = 'ML')$Zh 
  # Weight  
  weight <- pnorm(per.fit)*nrow(dat)/ sum(pnorm(per.fit)) 
  # Robust model estimation 
  robust.mod <- mirt(dat, 1, survey.weights = weight, optimizer = 'NR')
  
  return(robust.mod)
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
#' @return Histogram plot of residuals beneath a graph of the weight functions vs. the residuals.
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
#' @export

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
#' @export
                                       
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


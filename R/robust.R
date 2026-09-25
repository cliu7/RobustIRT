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
#' \eqn{a_j(\theta -b_j)} for unidimensional models and \eqn{\boldsymbol{a}_j' \boldsymbol{\theta} +d_j)} for the MIRT model.
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
#' dat <- matrix(c(4,1, 2,3), nrow = 2, byrow = TRUE)
#' theta <- c(-0.5, 1.0)
#' 
#' ipars <- rbind(
#'   c(a = 1.0, b1 = -2.0, b2 = -1.0, b3 = 0.0),
#'   c(a = 0.5, b1 = -1.0, b2 =  0.0, b3 = 1.0))
#'
#' grm.resid<-residual(theta, model = "GRM", ipars=ipars, dat=dat, resid="standardized")
#' @export


residual<-function(theta, model, ipars, dat=NULL, resid = c("standardized", "msr", "information"), D=1.7){
    
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
    K<-ncol(probs$P)
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
  
  if(length(out)==1){
    return(out[[1]])
  }else{
    return(out)
    
  }
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
#' Rasch_SE <- standard.errors(Rasch_theta, Rasch_ipars, dat, model = "RASCH")
#' 
#' # 2PL Model Example
#' fit_2PL <- mirt(dat, 1, itemtype = "2PL")
#' theta_2PL <- fscores(fit_2PL)
#' ipars_2PL <- coef(fit_2PL, IRTpars = TRUE, simplify = TRUE)$items[, c("a", "b")]
#' SE_2PL <- standard.errors(theta_2PL, ipars_2PL, dat, model = "2PL", est.type = c("MLE", "MAP"))
#' 
#' # MIRT Model Example
#' MIRT_fit <- mirt(dat, 2)
#' MIRT_theta <- fscores(MIRT_fit)
#' MIRT_ipars <- coef(MIRT_fit, IRTpars = TRUE, simplify = TRUE)$items
#' MIRT_ipars[32,"d"] <- -MIRT_ipars[32,"b"] 
#' MIRT_SE <- standard.errors(MIRT_theta, MIRT_ipars[,c("a1","a2","d")], dat, model = "MIRT")
#' 
#' # GRM Model Example for Agreeableness
#' bfi <- read_csv("p234_bfi_demog_2020-04-24.csv")
#' all_items <- bfi %>% dplyr::select(starts_with("t1_bfi_"))
#' colnames(all_items) <- sub("^t1_bfi_", "", colnames(all_items))
#' A_items <- all_items %>% dplyr::select(starts_with("A"))
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
#' MGRM_SE <- standard.errors(theta = MGRM_theta, ipars = MGRM_ipars, dat = MGRM_dat, model = "MGRM", est.type = "MLE", weight.type = "Huber", tuning.par = 1)


standard.errors<-function(theta, ipars, dat, model, D=1.7, weight.type = "equal", 
                          tuning.par = NULL, custom.weights = NULL, resid = "standardized", 
                          est.type = "MLE", prior=c(0,1), eap.quad.pts =seq(-4, 4, length.out = 41)){
  
  model<-toupper(model)
  est.type <- toupper(est.type)
  resid<- tolower(resid)
  
  # Item-level weight given residual vector
  compute.weights<- function(r_mat){
    if(weight.type == "equal") return(matrix(1, nrow(r_mat), ncol(r_mat)))
    if(weight.type == "Huber") return(huber(r_mat, tuning.par))
    if(weight.type == "bisquare") return(bisquare(r_mat, tuning.par))
    if(weight.type == "custom"){
      if(is.null(custom.weights)) stop("custom.weights must be supplied when weight.type = 'custom'.")
      return(custom.weights)
    } 
  }
  
  # Compute residual & weights 
  person.weights <- function(theta_mat, dat_mat, model_up, ipars_use){
    
    if(resid == "dual"){
      # dual: w(info) + w(msr)
      r.all <- residual(theta_mat, model_up, ipars_use, dat_mat, resid= c("information", "msr"), D = D)
      w1<-compute.weights(as.matrix(r.all$information))
      w2<-compute.weights(as.matrix(r.all$msr))
      return(w1+w2)
    }else{
      r.specific<- residual(theta_mat, model_up, ipars_use, dat_mat, resid= resid, D = D)
      return(compute.weights(as.matrix(r.specific)))
    }
  }
  
  dat<-as.matrix(dat)
  N<-nrow(dat) # number of subjects
  J<-ncol(dat) # test length
  
  mu<-prior[1]
  sigma2<-prior[2]
  P<- item.prob(theta, model, ipars, D=D)
  w<-person.weights(theta, dat, model, ipars)
  
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
      
      out$post_sd_MAP <- sqrt(1/D2.MAP)
      out$sandwich_MAP <-out$post_sd_MAP*sqrt(B/A)
    }
    
    if("EAP"%in%est.type){
      # Prior density according to each density point
      f_x<-dnorm(eap.quad.pts, prior[1], sqrt(prior[2]))
      # Item response probabilities at each quadrature point
      probs.q <- item.prob(eap.quad.pts, "RASCH", ipars, D) 
      Q<-length(eap.quad.pts)
      
      # (Weighted) Likelihood according to person's data and each quad point
      likelihood<-matrix(NA, N, Q)
      for(i in 1:N){
        likelihood[i,]<- apply(probs.q, 1, function(x) prod((x^dat[i,]*(1-x)^(1-dat[i,]))^w[i,]))
      }
      
      # Posterior standard deviation of EAP theta estimate according to Bock & Mislevy, 1989; Thissen et al., 1995
      out$post_sd_EAP <-sapply(1:N, function(x) sqrt(sum((eap.quad.pts - theta[x])^2*likelihood[x,]*f_x)/ sum(likelihood[x,]*f_x)))
      out$sandwich_EAP <-out$post_sd_EAP*sqrt((B/A))
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
      
      out$post_sd_MAP <- sqrt(1/D2.MAP)
      out$sandwich_MAP <-out$post_sd_MAP*sqrt(B/A)
    }
    
    if("EAP"%in%est.type){
      # Prior density according to each density point
      f_x<-dnorm(eap.quad.pts, prior[1], sqrt(prior[2]))
      # Item response probabilities at each quadrature point
      probs.q <- item.prob(eap.quad.pts, "1PL", ipars, D) 
      Q<-length(eap.quad.pts)
      
      # (Weighted) Likelihood according to person's data and each quad point
      likelihood<-matrix(NA, N, Q)
      for(i in 1:N){
        likelihood[i,]<- apply(probs.q, 1, function(x) prod((x^dat[i,]*(1-x)^(1-dat[i,]))^w[i,]))
      }
      
      # Posterior standard deviation of EAP theta estimate according to Bock & Mislevy, 1989; Thissen et al., 1995
      out$post_sd_EAP <-sapply(1:N, function(x) sqrt(sum((eap.quad.pts - theta[x])^2*likelihood[x,]*f_x)/ sum(likelihood[x,]*f_x)))
      out$sandwich_EAP <-out$post_sd_EAP*sqrt((B/A))
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
      out$post_sd_MAP <- sqrt(1/D2.MAP)
      out$sandwich_MAP <-out$post_sd_MAP*sqrt(B/A)
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
      out$sandwich_EAP <-out$post_sd_EAP*sqrt((B/A))
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
    pstar_arr <- P$pstar
    if(length(dim(pstar_arr)) == 2) pstar_arr <- array(pstar_arr, dim = c(dim(pstar_arr), 1))
    
    Pstar <- array(NA, dim = c(J, K+2, N))
    Pstar[,1,] <- 1
    Pstar[,2:(K+1),] <- pstar_arr
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
    D1 <- D*a*(P1*(1-P1) - P2*(1-P2))/Pk
    
    # Second Derivative 
    D2 <- D^2*a^2*((P1*(1-P1)*(1-2*P1) - P2*(1-P2)*(1-2*P2))/Pk -
                          (P1*(1-P1) - P2*(1-P2))^2/Pk^2)
    
    A <- colSums(t(w)*(-D2)) # Expected information
    V <- colSums(t(w)^2*(-D2)) # ASE numerator
    B <- colSums((t(w)*D1)^2) # Sandwich B
    
    A[A <= 0 | !is.finite(A)] <- NA
    
    if("MLE"%in%est.type){
      out$asymptotic_MLE <- sqrt(V)/A
      out$sandwich_MLE <- sqrt(B)/A
    }
    
    if("MAP"%in%est.type){
      
      # Weighted Second Derivative
      D2.MAP<-A + 1/prior[2]
      
      out$post_sd_MAP <- sqrt(1/D2.MAP)
      out$sandwich_MAP <-out$post_sd_MAP*sqrt(B/A)
    }
    
    if("EAP"%in%est.type){
      # Prior density according to each density point
      f_x<-dnorm(eap.quad.pts, prior[1], sqrt(prior[2]))
      # Item response probabilities at each quadrature point
      probs.q <- item.prob(eap.quad.pts, "GRM", ipars) 
      Pcat.q <- probs.q$P
      Q<-length(eap.quad.pts)
      
      # (Weighted) Likelihood according to person's data and each quad point
      likelihood<-matrix(NA, N, Q)
      for(i in 1:N){
        for(q in 1:Q){
          # JxK matrix of category probs for each quad point
          if(length(dim(Pcat.q))==3){
            Pq<- Pcat.q[,, q]   # J x K matrix for quad point q
          }else{
            Pq<- Pcat.q
          }
          # Probability of each person's observed response at quad point q
          P.obs <- Pq[cbind(1:J, dat[i,])]
          # Weighted likelihood contribution
          likelihood[i, q] <- prod(P.obs^w[i,])
        }
      }
      
      # Posterior standard deviation of EAP theta estimate according to Bock & Mislevy, 1989; Thissen et al., 1995
      out$post_sd_EAP <-sapply(1:N, function(x) sqrt(sum((eap.quad.pts - theta[x])^2*likelihood[x,]*f_x)/ sum(likelihood[x,]*f_x)))
      out$sandwich_EAP <-out$post_sd_EAP*sqrt((B/A))
    } 
  }
  
  if(model=="MGRM"){
    if(!is.matrix(theta)){
      theta <- matrix(theta, nrow = 1)
    }
    L <- ncol(theta)
    a<- ipars[, 1:L] # JxL matrix of discrimination parameters
    b<- ipars[, (L+1):ncol(ipars)] # JxK matrix of category thresholds
    K <-ncol(b) # number of thresholds
    
    # Calculate probabilities
    probs<-item.prob(theta, "MGRM", ipars)
    
    out_ase <- matrix(NA, N, L)
    out_sand <- matrix(NA, N, L)
    out_singular <- matrix(0,  N, L)
    
    for(i in 1:N){
      
      # Pstar: J x (K+2) boundary probabilities for person i
      Pstar_i <- matrix(NA, J, K + 2)
      Pstar_i[, 1] <- 1
      Pstar_i[, 2:(K+1)] <-probs$pstar[,, i]
      Pstar_i[, K+2]<- 0
      
      # Extract boundary probabilities
      ps0<- Pstar_i[cbind(1:J, dat[i, ])]
      ps1<- Pstar_i[cbind(1:J, dat[i, ] + 1)]
      qs0 <-1-ps0
      qs1 <-1-ps1
      Pk <- pmax(ps0-ps1, 1e-12)
      
      # Add A, V, B across items
      A_mat <- V_mat <- B_mat <- matrix(0, L, L)
      for(j in 1:J){
        aj <- matrix(a[j, ], ncol = 1) 
        wij <- w[i, j]
        
        # Numerators
        num1<- ps0[j]*qs0[j] - ps1[j]*qs1[j]
        num2<- ps0[j]*qs0[j]*(qs0[j]-ps0[j]) - ps1[j]*qs1[j]*(qs1[j]-ps1[j])
        
        # First derivative contribution for this item
        score_j <- D*num1/Pk[j] # scalar score contribution
        g_j <- wij*aj*score_j # L x 1 gradient
        
        # Item information 
        info_sc <- -(D^2*(num2/Pk[j]-(num1/Pk[j])^2))
        Ij<- info_sc*(aj %*%t(aj)) 
        
        A_mat<- A_mat+wij*Ij
        V_mat<- V_mat+wij^2*Ij
        B_mat<- B_mat+g_j%*%t(g_j)
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
                               
#' Robust Latent Trait Estimation
#'
#' Estimates latent trait(s) under robust weighting for the Rasch, 1PL, 2PL, MIRT,
#' GRM, and MGRM item response models. The function supports robust ML estimation 
#' for all models, MAP and EAP estimation for unidimensional models, and WLE 
#' estimation for the Rasch-2PL models. Standard errors are returned
#' alongside the ability estimates.
#' 
#' @param dat An \eqn{N \times J} matrix of item responses (\eqn{N} subjects,
#'   \eqn{J} items). Responses for dichotomous models must be 0/1. Responses for
#'   polytomous models must be integer-valued starting at 1. Missing data is not supported.
#' @param ipars A matrix of item parameters structured identically to the
#'   \code{ipars} argument of \code{\link{item.prob}}.
#' @param model The IRT model to be used: \code{"Rasch"}, \code{"1PL"}, \code{"2PL"},
#'   \code{"MIRT"}, \code{"GRM"}, or \code{"MGRM"}. See \code{\link{item.prob}}
#'   for details on each model.
#' @param D Scaling constant for the normal ogive model. Default 1.7; often set 
#'   to 1.0 for the logistic scale.
#' @param residual The residual used to compute weights:
#'   \code{"standardized"} (default), \code{"msr"}, \code{"information"}, or \code{"dual"}.
#' @param weight.type Weighting scheme: \code{"equal"} (default), \code{"Huber"},
#'   \code{"bisquare"}, or \code{"custom"}.
#' @param tuning.par Tuning parameter for Huber or bisquare weights.  Required when
#'   \code{weight.type} is not \code{"equal"} or \code{"custom"}.
#' @param custom.weights An \eqn{N \times J} numeric matrix of user-specified weights
#'   in \eqn{[0, 1]}.  Required when \code{weight.type = "custom"}.
#' @param est.type Estimation methods: \code{"MLE"}, \code{"MAP"}, \code{"EAP"}, \code{"WLE"}.  
#'   Can specify more than one in a vector (e.g., \code{c("MLE","MAP")}) for supported models.
#' @param dimen Numeric integer specifying the number of dimensions for the IRT model.
#'   Required when using the MGRM.
#' @param prior Length-2 numeric vector \code{c(mu, sigma2)} for the normal prior
#'   used in MAP and EAP.  Default is the standard normal \code{c(0, 1)}.
#' @param eap.quad.pts Numeric vector of quadrature points for EAP. Default is 41
#'   equally spaced points in the interval \eqn{[-4, 4]}.
#' @param up.bound Upper bound for truncating \eqn{\hat{\theta}}. Any converged 
#'   estimate greater than \code{up.bound} is replaced by the specified value. Default is 3.0.
#' @param low.bound Lower bound for truncating \eqn{\hat{\theta}}. Any converged 
#'   estimate less than \code{low.bound} is replaced by the specified value. Default is -3.0.
#' @param init.val Starting value(s) for the Newton–Raphson algorithm. May be a
#'   scalar (applied to all subjects and dimensions), a length-\eqn{L} vector for \eqn{L} dimensions
#'   (applied to all subjects), or an \eqn{N \times L} matrix.  Default is 0.
#' @param iter Maximum number of Newton–Raphson iterations.  Default is 30.
#' @param tol Convergence tolerance on the absolute change in log-likelihood.
#'   Default is 0.01.
#' @param score.tol Convergence tolerance on the absolute change in the score function.
#'   Default is \code{tol}.
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
#'               Convergence occurs when the absolute change in the log-likelihood is below \code{tol} AND the absolute change in the score function is below \code{score.tol} between two iterations, not exceeding the maximum number of iterations specified in \code{iter}.
#' 
#' @section Estimation types:
#' \itemize{
#'   \item \code{"MLE"} — Maximum likelihood estimation with the Newton–Raphson algorithm.
#'   \item \code{"MAP"} — Maximum a posteriori (Bayesian mode) with the Newton–Raphson algorithm
#'         under a normal prior. Available for Rasch, 1PL, 2PL, and GRM only.
#'   \item \code{"EAP"} — Expected a posteriori under a normal prior using Gauss–Hermite-style quadrature.
#'         Available for Rasch, 1PL, and 2PL only.
#'   \item \code{"WLE"} - Weighted likelihood estimation (Warm, 1989) with the Newton-Raphson algorithm,
#'         equivalent to the MAP estimate under Jeffreys' prior. Available for Rasch, 1PL, and 2PL only.
#' }
#'
#' @section Residuals:
#' \itemize{
#'   \item \code{"information"} — Information residual \eqn{r_{ij} = a_j(\theta_i - b_j)}
#'         or \eqn{r_{ij} = \boldsymbol{a}_j'\boldsymbol{\theta}_i + d_j} in the 
#'         multidimensional case. Available only for the Rasch, 1PL, 2PL, and 
#'         MIRT models. The information residual detects anomalies in the predictor 
#'         space (i.e., \eqn{\boldsymbol{\theta}}) to produce Mallows-class weights) 
#'         and supports the asymptotic standard error (Magis, 2014).
#'   \item \code{"standardized"} — Pearson residual
#'         \eqn{z_j = (x_{ij} - E[X_{ij}|\theta_i]) / \sqrt{Var[X_{ij}|\theta_i]}}.
#'         Available for all models. The standardized residual detects anomalies 
#'         in the observed response space. 
#'   \item \code{"msr"} — Modified standardized residual
#'         \eqn{M_{ij} = (x_{ij} - E[X_{ij}|\theta_i]) / P(x_{ij}|\theta_i)}.
#'         Available for all models. The modified standardized residual detects 
#'         anomalies in the observed response space. 
#'   \item \code{"dual"} — Dual-weight estimator (Chen et al., 2025). Computes
#'         \eqn{w(r_{1ij}) + w(r_{2ij})} where \eqn{r_{1ij}} is the information residual
#'         and \eqn{r_{2ij}} is the MSR. Available for the Rasch,
#'         1PL, and 2PL models only.
#' }
#'
#' @section Weight functions:
#' \itemize{
#'   \item \code{"equal"} — All weights are 1 (e.g., non-robust MLE/MAP/EAP).
#'   \item \code{"Huber"} — Huber (1981) weight. Requires \code{tuning.par}, recommended at 1.
#'   \item \code{"bisquare"} — Tukey bisquare weight (Mosteller & Tukey, 1977). 
#'         Requires \code{tuning.par}, recommended at 4.
#'   \item \code{"custom"} — User-supplied \eqn{N \times J} matrix of weights 
#'         specified in \code{custom.weights}.
#' }
#'
#' @section Standard errors:
#' Every returned estimation type includes a sandwich standard error or deviation
#' of the ability estimate, along with the traditionally reported value:
#' \describe{
#'   \item{\code{ase}}{Asymptotic standard error that accounts for robust weights (Magis, 2014). For Rasch/
#'         1PL/2PL with \code{residual = "information"} this is the exact asymptotic standard error;
#'         otherwise it is the general weighted-information standard error
#'         \eqn{\sqrt{\sum w^2 I_j} / \sum w I_j}.}
#'   \item{\code{post_sd}}{Posterior SD for Bayesian estimators, MAP and EAP.}
#'   \item{\code{sandwich}}{Huber-White sandwich standard error (White, 1980). The asymptotic 
#'         standard error is weighted by a factor to account for model misspecification. For 
#'        Bayesian estimators MAP and EAP, the sandwich version of the posterior 
#'        standard deviation is provided (Li & Rice, 2023).}
#' }
#'
#' @references Chen, J., Cheng, Y., & Li, Z. (2025). Dual-weight robust estimation of
#'   latent ability. \emph{Journal of Educational and Behavioral Statistics}.
#' @references Embretson, S. E., & Reise, S. P. (2000). \emph{Item response theory for psychologists.} Mahwah, N.J: L. Erlbaum Associates.
#' @references Filonczuk, A., & Cheng, Y. (2025). Robust estimation of the latent trait
#'   in graded response models. \emph{Behavior Research Methods, 57}(1), 55.
#' @references Huber, P. (1981). \emph{Robust Statistics}. John Wiley & Sons, Inc.
#' @references Li, K. Q., & Rice, K. M. (2023). A Bayesian “sandwich” for variance estimation. arXiv. 
#' @references Magis, D. (2014). On the asymptotic standard error of a class of robust
#'   estimators of ability in dichotomous item response models.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 67}(3), 430–450.
#' @references Maeda, H., & Zhang, B. (2020). Bayesian Extension of Biweight and 
#'   Huber Weight for Robust Ability Estimation. \emph{Journal of Educational Measurement}, 57(1), 51–70. 
#' @references Mislevy, R. J., & Bock, R. D. (1982). Biweight estimates of latent
#'   ability. \emph{Educational and Psychological Measurement, 42}(3), 725–737.
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}.
#'   Addison-Wesley Publishing Company.
#' @references Schuster, C., & Yuan, K.-H. (2011). Robust estimation of latent ability in
#'   item response models. \emph{Journal of Educational and Behavioral Statistics,
#'   36}(6), 720–735.
#' @references White, H. (1980). A Heteroskedasticity-Consistent Covariance Matrix 
#'   Estimator and a Direct Test for Heteroskedasticity. \emph{Econometrica, 48}(4), 817–838.
#' @references Yu, X., & Cheng, Y. (2019). A change-point analysis procedure based on
#'   weighted residuals to detect back random responding.
#'   \emph{Psychological Methods, 24}(5), 658–674.
#' @references Warm, T.A. Weighted likelihood estimation of ability in item 
#'   response theory. Psychometrika 54, 427–450 (1989).
#'
#'
#' @return A named list. Elements present depend on \code{est.type}:
#' \describe{
#'   \item{\code{theta_MLE}}{(\eqn{N \times L}) Robust MLE estimates.}
#'   \item{\code{ase_MLE}}{(\eqn{N \times L}) Asymptotic standard error or Fisher-information based 
#'     standard error depending on if the conditions presented in Magis (2014) are met
#'     (see Details).}
#'   \item{\code{sandwich_MLE}}{(\eqn{N \times L}) Huber-White sandwich standard error.}
#'   \item{\code{convergence_MLE}}{(\eqn{N \times L}) Nonconvergence indicators: 0 = converged, 1 = did not converge, 2 = converged outside of bounds \code{low.bound} and \code{up.bound}. }
#'   \item{\code{theta_MAP}}{(\eqn{N \times 1}) Robust MAP estimates.
#'     Rasch/1PL/2PL/GRM only.}
#'   \item{\code{post_sd_MAP}}{(\eqn{N \times 1}) Posterior standard deviation of the MAP.}
#'   \item{\code{sandwich_MAP}}{(\eqn{N \times 1}) Sandwich equivalent of the posterior standard deviation.}
#'   \item{\code{convergence_MAP}}{(\eqn{N \times 1}) Nonconvergence indicators: 0 = converged, 1 = did not converge, 2 = converged outside of bounds \code{low.bound} and \code{up.bound}}
#'   \item{\code{theta_EAP}}{(\eqn{N \times 1}) Robust EAP estimates.
#'     Rasch/1PL/2PL only.}
#'   \item{\code{post_sd_EAP}}{(\eqn{N \times 1}) Posterior standard deviation of the EAP.}
#'   \item{\code{sandwich_EAP}}{(\eqn{N \times 1}) Sandwich equivalent of the posterior standard deviation.}
#'   \item{\code{weights}}{(\eqn{N \times J}) Item weights at the iteration of convergence.}
#'   \item{\code{residuals}}{(\eqn{N \times J}) Residuals used for weighting 
#'     at the iteration of convergence.}
#' }
#' @examples
#' ## 2PL MLE + MAP with Huber weights (information residual)
#' set.seed(25)
#' N<-50 
#' J<-20
#' ipars<-cbind(a = runif(J, 0.8, 2), b = rnorm(J))
#' theta_true<-rnorm(N)
#' P<-item.prob(theta_true, "2PL", ipars)
#' dat<-dat.gen(P)
#' out<-robust.theta(dat, ipars, model = "2PL", resid = "information", weight.type = "Huber", 
#'                     tuning.par = 1, est.type = c("MLE", "MAP"))
#' head(out$theta_MLE)
#' head(out$ase_MLE)
#' head(out$post_sd_MAP)
#'
#' ## GRM MLE with standardized residual + bisquare weights 
#' set.seed(24)
#' J<-15
#' K<-5 
#' N<-100
#' a<-runif(J, 0.8, 2)
#' b<-t(apply(matrix(rnorm(J*(K-1)), J), 1, sort))
#' ip<-cbind(a, b)
#' th<-rnorm(N)
#' P<-item.prob(th, "GRM", ip)
#' dat<-dat.gen(P$P, anchor = 1)
#' out_grm<-robust.theta(dat, ip, model = "GRM", resid = "standardized", weight.type = "bisquare", tuning.par = 4, est.type = "MLE")
#' head(out_grm$theta_MLE)
#'
#' ## MIRT MLE with information residual Error 
#' set.seed(26)
#' L<-2
#' J<-50
#' N<-80
#' a<-matrix(runif(J*L, 0.5, 1.5), J, L)
#' d<-rnorm(J)
#' ip<-cbind(a, d)
#' th<-matrix(rnorm(N*L), N, L)
#' P<-item.prob(th, "MIRT", ip)
#' dat<-dat.gen(P)
#' out_mirt<-robust.theta(dat, ip, model = "MIRT", resid = "information", 
#'                           weight.type = "Huber", tuning.par = 1, est.type = "MLE", score.tol = 0.1)
#' head(out_mirt$theta_MLE)
#' 
#' ## MGRM MLE with standardized residual 
#' set.seed(27)
#' L<-2
#' K<-5
#' J<-50
#' N<-20
#' a<-matrix(runif(J*L, 0.5, 1.5), J, L)
#' d<-t(apply(matrix(rnorm(J*(K-1)), J), 1, sort))
#' ip<-cbind(a, d)
#' th<-matrix(rnorm(N*L), N, L)
#' P<-item.prob(th, "MGRM", ip)$P
#' dat<-dat.gen(P, anchor=1)
#' out_mgrm<-robust.theta(dat, ip, model = "MGRM", resid = "standardized", dimen=2,
#'                           weight.type = "Huber", tuning.par = 1, est.type = "MLE")
#' head(out_mgrm$theta_MLE)
#' 
#' ## WLE Example: 2PL 
#' set.seed(20)
#' thetas<-c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
#' a<-runif(15, 0.8, 2)
#' b<-rnorm(15)
#' P<-item.prob(thetas, "2PL", cbind(a, b))
#' dat<-dat.gen(P, anchor = 0)
#' out_wle<-robust.theta(dat, cbind(a, b), model = "2PL", resid = "information", 
#'                           weight.type = "Huber", tuning.par = 1, est.type = "WLE")
#' head(out_wle$theta_WLE)
#' @export

robust.theta<-function(dat, ipars, model= "2PL", D = 1.7, resid = "standardized", 
                         weight.type = "equal", tuning.par = NULL, custom.weights = NULL,
                         est.type = "MLE", dimen = NULL, prior = c(0, 1), eap.quad.pts = seq(-4, 4, length.out = 41),
                         init.val = 0, up.bound = 3, low.bound = -3, iter = 30, tol = 0.01, score.tol = tol){
  
  ##### Ensure proper input #####
  
  model<-toupper(model)
  est.type<-toupper(est.type)
  resid<-tolower(resid)
  
  dich.mods<-c("RASCH", "1PL", "2PL")
  
  # Set checks to ensure proper inputs
  if(!(model %in% c(dich.mods, "MIRT", "GRM", "MGRM")))
    stop(paste(model, "is not a supported model." ))
  
  if(!(weight.type %in% c("equal", "Huber", "bisquare", "custom")))
    stop(paste(weight.type, "is not a supported weight.type."))
  
  if(weight.type %in% c("Huber", "bisquare") && is.null(tuning.par))
    stop(paste("tuning.par must be supplied when weight.type = ", weight.type))
  
  if(weight.type == "custom" && is.null(custom.weights))
    stop("custom.weights must be supplied when weight.type = 'custom'.")
  
  if(!(resid %in% c("standardized", "msr", "information", "dual")))
    stop(paste(resid, "is not a supported residual type."))
  
  if(resid %in% c("information", "dual") && !(model %in% c(dich.mods, "MIRT")))
    stop("The 'information' and 'dual' residuals are only available for Rasch, 1PL, and 2PL models.")
  
  if(any(!(est.type %in% c("MLE", "MAP", "EAP", "WLE"))))
    stop("est.type must be one or more of 'MLE', 'MAP', 'EAP', 'WLE'.")
  
  if("EAP" %in% est.type && !(model %in% dich.mods))
    stop("EAP estimation is only available for Rasch, 1PL, and 2PL models.")
  
  if("MAP" %in% est.type && !(model %in% dich.mods))
    stop("MAP estimation is only available for Rasch, 1PL, and 2PL models.")
  
  if("WLE" %in% est.type && !(model %in% dich.mods))
    stop("WLE estimation is only available for Rasch, 1PL, and 2PL models.")
  
  # Initialize variables
  dat<-as.matrix(dat)
  N<-nrow(dat)
  J<-ncol(dat)
  
  if(weight.type == "custom" && (dim(custom.weights)!=dim(dat)))
    stop("Custom weights must be a matrix the same dimension as the data")
  
  mu<-prior[1]
  sigma2<-prior[2]
  
  ##### Internal functions #####
  
  # Item-level weight given residual vector
  compute.weights<-function(r_mat){
    if(weight.type == "equal") return(matrix(1, nrow(r_mat), ncol(r_mat)))
    if(weight.type == "Huber") return(huber(r_mat, tuning.par))
    if(weight.type == "bisquare") return(bisquare(r_mat, tuning.par))
    if(weight.type == "custom"){
      if(is.null(custom.weights)) stop("custom.weights must be supplied when weight.type = 'custom'.")
      return(custom.weights)
    } 
  }
  
  # Compute residual & weights for person i at current theta_i
  person.weights<-function(theta_i, dat_i, model_up, ipars_use){
    
    if(resid == "dual"){
      # dual: w(info) + w(msr)
      r.all<-residual(theta_i, model_up, ipars_use, matrix(dat_i, nrow = 1), resid= c("information", "msr"), D = D)
      w1<-compute.weights(matrix(r.all$information, nrow=1))
      w2<-compute.weights(matrix(r.all$msr, nrow=1))
      return(w1+w2)
    }else{
      r.specific<- residual(theta_i, model_up, ipars_use, matrix(dat_i, nrow = 1), resid= resid, D = D)
      return(compute.weights(as.matrix(r.specific, nrow=1)))
    }
  }
 
  # Get initial theta for subject i 
  get_init<-function(i, L){
    if(is.matrix(init.val)){ # if user specifies different initial values for each subject
      return(as.matrix(init.val[i,]))      
    }else if(length(init.val) == L){ #if user specifies L values, the same initial values are used for each subject
      return(as.matrix(init.val)) 
    }else{ #if user specifies one initial value, it will be used for all dimensions and each subject
      return(matrix(rep(init.val[1], L), nrow = L))
    }
  }
  
  # Unidimensional NR step for Rasch/1PL/2PL [returns list(D1, D2)]
  # If bayes=T, DI / and D2 are of the log posterior
  dich_derivs<-function(theta_i, dat_i, ipars_use, w_i, bayes = FALSE, wle = FALSE){
    P_i<-as.numeric(item.prob(theta_i, model, ipars_use, D))
    a_j<-ipars_use[, 1]
    D1<-sum(w_i*D*a_j*(dat_i - P_i))
    D2<-sum(w_i*D^2*a_j^2*P_i*(1 - P_i)) # this is the same as information, or I in Warm's WLE
    if(bayes){
      D1<-D1 - (theta_i - mu)/sigma2
      D2<-D2 + 1/sigma2
    }
    
    if(wle){
      Q_i = 1-P_i
      P.p<-D*a_j*P_i*Q_i # P prime
      P.pp<-D^2*a_j^2*P_i*Q_i*(1-2*P_i) #P prime prime
      
      # Bias correction numerator and its derivative
      J.i <- sum(w_i*P.p*P.pp / (P_i*Q_i))
      J.p <- sum(w_i*D^4*a_j^4*P_i*Q_i*(1 - 6*P_i*Q_i))
      
      D1 <- D1 + J.i/(2*D2)
      D2 <- D2-(D2*J.p-J.i^2) / (2*D2^2) # make sure this is consistent with later negation
    }
    list(D1 = D1, D2 = -D2, P_i = P_i, info_j = D^2*a_j^2*P_i*(1 - P_i))
  }
  
  # GRM NR step
  # If bayes=T, DI / and D2 are of the log posterior
  grm_derivs<-function(theta_i, dat_i, ipars_use, w_i, bayes = FALSE) {
    a<-ipars_use[, 1]
    b<-ipars_use[, -1, drop = FALSE]
    prb<-item.prob(theta_i, "GRM", ipars_use, D)
    K<-ncol(b)
    Pcat<-prb$P   
    pstar_ext<-cbind(1, prb$pstar[,,1], 0)  
    if(length(dim(pstar_ext)) == 3) pstar_ext<-pstar_ext[,, 1]
    
    ps0<-pstar_ext[cbind(seq_len(J), dat_i)]
    ps1<-pstar_ext[cbind(seq_len(J), dat_i + 1)]
    qs0<-1-ps0
    qs1<-1-ps1
    Pk<-ps0-ps1
    
    D1<-sum(w_i * D * a * (ps0 * qs0 - ps1 * qs1) / Pk)
    D2<-sum(w_i * D^2 * a^2 * ((ps0*qs0*(qs0-ps0) - ps1*qs1*(qs1-ps1)) / Pk -
                                   (ps0*qs0 - ps1*qs1)^2 / Pk^2))
    # For SE: scored expected info per item (positive)
    info_j<--D^2 * a^2 * ((ps0*qs0*(qs0-ps0) - ps1*qs1*(qs1-ps1)) / Pk -
                              (ps0*qs0 - ps1*qs1)^2 / Pk^2)
    
    if (bayes) {
      D1<-D1-(theta_i-mu)/sigma2
      D2<-D2 + 1/sigma2
    }
    list(D1 = D1, D2 = D2, Pk = Pk, ps0 = ps0, ps1 = ps1, info_j = info_j)
  }
  
  ##### More preparation #####
  
  ipars_use<-ipars
  if(model == "RASCH"){
    if(is.vector(ipars_use) || ncol(ipars_use) == 1)
      ipars_use<-cbind(rep(1, J), as.numeric(ipars_use))
  }
 
  # Extract dimension
  if(model == "MIRT"){
    L<-ncol(ipars_use) - 1   # last column is d
  }else if(model == "MGRM"){
    L<-dimen
  }else{
    L<-1
  }
  
  # For storing output
  out<-list()
  
  ##### MLE (Rasch, 1PL, 2PL, GRM, MIRT, MGRM) ####
  if("MLE" %in% est.type){
    
    # Unidimensional dichotomous models
    if(model %in% dich.mods){
      
      # Initialize output for storing
      theta_mle<-rep(NA, N)
      conv_mle<-rep(0, N)
      w_final<-matrix(NA, N, J)
      r_final<-matrix(NA, N, J)
      
      for(i in 1:N){
        th_i<-get_init(i, 1)
        P0<-0
        
        for(k in 1:iter){
          w_i<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
          dv<-dich_derivs(th_i, dat[i,], ipars_use, w_i)
          
          if(!is.finite(dv$D1) || !is.finite(dv$D2) || abs(dv$D2) < 1e-14){
            conv_mle[i]<-1
            break
          }
          
          th_new<-th_i-dv$D1/dv$D2
          
          if(!is.finite(th_new)){ 
            conv_mle[i]<-1
            break 
          }
          
          # Change in log likelihood for convergence criteria 1
          log_like<-sum(log(pmax(ifelse(dat[i,] == 1, dv$P_i, 1 - dv$P_i), 1e-12))) - sum(log(pmax(P0, 1e-12)))
          P0<- ifelse(dat[i,] == 1, dv$P_i, 1 - dv$P_i)
          th_i<-th_new
          
          # End NR algorithm if all 3 convergence criteria are met
          if(k > 1 && abs(log_like) < tol && abs(dv$D1) < score.tol){break}
          
          # Check for reaching max # of iterations
          if(k == iter){
            conv_mle[i]<-1
            th_i<-NA
          }
            
        }
        #Check for convergence outside bounds (convergence code = 2)
        if(!is.na(th_i) & th_i<low.bound){
          th_i<-low.bound
          conv_mle[i]<-2
        }else if(!is.na(th_i) & th_i > up.bound){
          th_i<-up.bound
          conv_mle[i]<-2
        }
        
        # Final weights / residuals at converged theta
        theta_mle[i]<-th_i
        if(is.na(th_i)){
          r_final[i,]<-w_final[i,]<-rep(NA, J)
        }else{
          w_final[i,]<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
          r_final[i,]<-residual(th_i, model, ipars_use, matrix(dat[i,], 1), resid=resid, D=D)
        }
      } # end person loop
      
      out$theta_MLE<- matrix(theta_mle, ncol = 1)
      
      # Compute SEs for converged subjects
      conv_idx<-which(conv_mle!=1 & !is.na(theta_mle))
      
      ase_MLE<-rep(NA, N)
      sand_MLE<-rep(NA, N)
      
      if(length(conv_idx) > 0){
        se.all<-standard.errors(theta_mle[conv_idx], ipars_use, dat[conv_idx, , drop = FALSE], model, D, 
                                  weight.type, tuning.par, custom.weights, resid, "MLE", prior, eap.quad.pts)
        ase_MLE[conv_idx]<-se.all$asymptotic_MLE
        sand_MLE[conv_idx]<-se.all$sandwich_MLE
      }
      
      out$ase_MLE<-matrix(ase_MLE, ncol=1)
      out$sandwich_MLE<-matrix(sand_MLE, ncol=1)
      out$convergence_MLE<-matrix(conv_mle, ncol = 1)
      out$weights_MLE<-w_final
      out$residuals_MLE<-r_final
    }
    
    ### GRM (unidimensional)
    if(model == "GRM"){
      a_grm<-ipars_use[, 1]
      b_grm<-ipars_use[, -1, drop = FALSE]
      
      theta_mle <-rep(NA, N)
      conv_mle <-rep(0, N)
      w_final<-matrix(NA, N, J)
      r_final<-matrix(NA, N, J)
      
      for(i in 1:N){
        th_i<-get_init(i, 1)
        P0<-0
        
        for(k in 1:iter){
          w_i<-as.numeric(person.weights(th_i, dat[i,], "GRM", ipars_use))
          gd<-grm_derivs(th_i, dat[i,], ipars_use, w_i)
          
          if(!is.finite(gd$D1) || !is.finite(gd$D2) || abs(gd$D2) < 1e-14){
            conv_mle[i]<-1
            break
          }
          
          th_new<-th_i- gd$D1/gd$D2
          if(!is.finite(th_new)){ 
            conv_mle[i]<-1
            break 
          }
          
          Pk_log<-log(pmax(gd$Pk, 1e-12))
          log_like<-sum(Pk_log) - sum(log(pmax(P0, 1e-12)))
          P0<-gd$Pk
          th_i<-th_new
          
          # End NR algorithm if all 3 convergence criteria are met
          if(k > 1 && abs(log_like) < tol && abs(gd$D1) < score.tol){break}
          
          # Check for reaching max # of iterations 
          if(k == iter){
            conv_mle[i]<-1
            th_i<-NA
          }
        }
        
        #Check for convergence outside bounds (convergence code = 2)
        if(!is.na(th_i) & th_i<low.bound){
          th_i<-low.bound
          conv_mle[i]<-2
        }else if(!is.na(th_i) &  th_i > up.bound){
          th_i<-up.bound
          conv_mle[i]<-2
        }
        
        theta_mle[i]<-th_i
        # Final weights / residuals at converged theta
        if(is.na(th_i)){
          r_final[i,]<-w_final[i,]<-rep(NA, J)
        }else{
          w_final[i,]<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
          r_final[i,]<-residual(th_i, model, ipars_use, matrix(dat[i,], 1), resid=resid, D=D)
          
        }
      }
      
      out$theta_MLE<-matrix(theta_mle, ncol = 1)
      
      # Compute SEs for converged subjects
      conv_idx<-which(conv_mle!=1& !is.na(theta_mle))
      
      ase_MLE<-rep(NA, N)
      sand_MLE<-rep(NA, N)
      
      if(length(conv_idx) > 0){
        se.all<-standard.errors(theta_mle[conv_idx], ipars_use, dat[conv_idx, , drop = FALSE], model, D, 
                                  weight.type, tuning.par, custom.weights, resid, "MLE", prior, eap.quad.pts)
        ase_MLE[conv_idx]<-se.all$asymptotic_MLE
        sand_MLE[conv_idx]<-se.all$sandwich_MLE
      }
      
      out$ase_MLE = matrix(ase_MLE, ncol = 1)
      out$sandwich_MLE = matrix(sand_MLE, ncol = 1)
      out$convergence_MLE<-matrix(conv_mle, ncol = 1)
      out$weights_MLE<-w_final
      out$residuals_MLE<-r_final
    }
    
    ### MIRT
    if(model == "MIRT"){
      a_m<-ipars_use[,1:L, drop = FALSE]
      d_m<-ipars_use[,L+1]
      
      theta_mle<-matrix(NA, N, L)
      conv_mle<-matrix(0, N, L)
      sing_flag<-matrix(0, N, L)
      w_final<-matrix(NA, N, J)
      r_final<-matrix(NA, N, J)
      
      for(i in 1:N){
        th_i<-get_init(i, L)
        P0<-0
        
        for(k in 1:iter){
          P_i<- as.numeric(item.prob(matrix(th_i, 1), "MIRT", ipars_use, D))
          w_i<- as.numeric(person.weights(matrix(th_i, 1), dat[i,], "MIRT", ipars_use))
          
          # Gradient (L x 1) and Hessian (L x L)
          grad<-matrix(0, L, 1)
          Hess<-matrix(0, L, L)
          for(j in 1:J){
            aj<-matrix(a_m[j,], L, 1)
            pj<-P_i[j]
            wij<-w_i[j]
            grad<-grad + wij*D*aj*(dat[i, j] - pj)
            Hess<-Hess + wij*D^2*(aj %*% t(aj))*pj*(1 - pj)
          }
          
          if(any(!is.finite(grad)) || any(!is.finite(Hess))){
            conv_mle[i,]<-1
            break
          }
          
          det_H<-ifelse(L == 1, Hess[1,1], det(Hess))
          
          if(abs(det_H) < 1e-12){ 
            conv_mle[i,]<-1
            break 
          }
          
          if(L==1){
            H_inv<-matrix(1 / Hess[1,1])
          }else{
            H_inv<-solve(Hess)
          }
          
          th_new<-th_i + as.numeric(H_inv %*% grad)  # grad is +, Hessian is neg expected info
          
          if(any(!is.finite(th_new))){ 
            conv_mle[i,]<-1
            break 
          }
          
          log_like<-sum(log(pmax(ifelse(dat[i,] == 1, P_i, 1 - P_i), 1e-12))) -sum(log(pmax(P0, 1e-12)))
          P0<-ifelse(dat[i,] == 1, P_i, 1 - P_i)
          th_i<-th_new
          
          #  End NR algorithm if all 3 convergence criteria are met
          # L-infinity norm employed for multidimensional theta vector
          if(k > 1 && abs(log_like) < tol && max(abs(grad))<score.tol){break}
          
          # Check if convergence within # of iterations
          if(k == iter){conv_mle[i,]<-1}
        }
        
        if(any(conv_mle[i,] == 1)){
            theta_mle[i,]<-rep(NA, L)
          }else if(any(th_i>up.bound) || any(th_i<low.bound)){
            # Replace any converged estimates outside the bounds
            theta_mle[i,]<-pmax(low.bound, pmin(up.bound, th_i))
            conv_mle[i,]<-2
          }else{
            theta_mle[i,]<-th_i
          }
        
        # Compute item-level weights and residuals
        if(any(is.na(th_i))){
          r_final[i,]<-w_final[i,]<-rep(NA, J)
        }else{
          w_final[i,]<-as.numeric(person.weights(matrix(th_i, 1), dat[i,], "MIRT", ipars_use))
          r_final[i,]<-as.numeric(residual(matrix(th_i, 1), "MIRT", ipars_use,
                                           matrix(dat[i,], 1), resid, D=D))
          
        }
      }
      
      out$theta_MLE<-theta_mle
      
      # Compute SEs for converged subjects
      conv_idx<-which(apply(conv_mle, 1, function(x) all(x != 1)) & 
                          apply(theta_mle, 1, function(x) all(!is.na(x))))
      
      ase_MLE<-matrix(NA, N, L)
      sand_MLE<-matrix(NA, N, L)
      
      if(length(conv_idx) > 0){
        se.all<-standard.errors(theta_mle[conv_idx,], ipars_use, dat[conv_idx, , drop = FALSE], model, D, 
                                  weight.type, tuning.par, custom.weights, resid, "MLE", prior, eap.quad.pts)
        ase_MLE[conv_idx,]<-se.all$asymptotic_MLE
        sand_MLE[conv_idx,]<-se.all$sandwich_MLE
      }
      out$ase_MLE = ase_MLE
      out$sandwich_MLE = sand_MLE
      out$singular_flag = sing_flag
      out$convergence_MLE<-conv_mle
      out$weights_MLE<-w_final
      out$residuals_MLE<-r_final
    }
    
    ### MGRM
    if(model == "MGRM"){
      
      a_mg<-ipars_use[, 1:L, drop = FALSE]
      b_mg<-ipars_use[, (L+1):ncol(ipars_use), drop = FALSE]
      K_mg<-ncol(b_mg)
      
      theta_mle<-matrix(NA, N, L)
      conv_mle<-matrix(0, N, L)
      sing_flag<-matrix(0, N, L)
      w_final<-matrix(NA, N, J)
      r_final<-matrix(NA, N, J)
      
      for(i in 1:N){
        th_i<- get_init(i, L)
        P0<- 0
        
        for(k in 1:iter){
          prb<-item.prob(matrix(th_i, 1), "MGRM", ipars_use, D)
          
          # Pstar extended J x (K_mg+2)
          Pstar_i<- matrix(NA, J, K_mg + 2)
          Pstar_i[,1]<-1
          Pstar_i[,2:(K_mg+1)] <-prb$pstar[,, 1]
          Pstar_i[,K_mg+2]<-0
          
          Pcat_i<-prb$P   # J x (K_mg+1)
          
          ps0<-Pstar_i[cbind(1:J, dat[i,])]
          ps1<-Pstar_i[cbind(1:J, dat[i,] + 1)]
          qs0<-1 - ps0
          qs1<-1 - ps1
          Pk<-pmax(ps0 - ps1, 1e-12)
          
          w_i<-as.numeric(person.weights(matrix(th_i, 1), dat[i,], "MGRM", ipars_use))
          
          # Gradient and Hessian
          grad<-matrix(0, L, 1)
          Hess<-matrix(0, L, L)
          for(j in 1:J){
            aj<-matrix(a_mg[j,], L, 1)
            wij<-w_i[j]
            num1<-ps0[j]*qs0[j] - ps1[j]*qs1[j]
            num2<-ps0[j]*qs0[j]*(qs0[j]-ps0[j]) - ps1[j]*qs1[j]*(qs1[j]-ps1[j])
            score_j<-D*num1 / Pk[j]
            grad<-grad + wij*aj*score_j
            info_sc<-D^2*(num2 / Pk[j] - (num1 / Pk[j])^2)
            Ij_mat<-(-info_sc)*(aj %*% t(aj))
            Hess<-Hess + wij*Ij_mat
          }
          
          if(any(!is.finite(grad)) || any(!is.finite(Hess))){ 
            conv_mle[i,]<-1
            break 
          }
          
          # Checks for singular matrix
          det_H<-ifelse(L == 1, Hess[1,1], det(Hess))
          if(abs(det_H) < 1e-12){ 
            sing_flag[i,]<-1
            break 
          }
          
          # More checks for singular matrix and update theta estimate
          if(L == 1){
            H_inv<-matrix(1 / Hess[1,1])
            th_new<-th_i + as.numeric(H_inv %*% grad)
            
          }else{
            chck<-try(solve(Hess), silent = TRUE)
            if(!inherits(chck, "try-error")){
              H_inv<-solve(Hess)
              th_new<-th_i + as.numeric(H_inv %*% grad)
            }else{
              th_new<-NA
              sing_flag[i,]<-1
              break
            }
              
          }
          
          if(any(!is.finite(th_new))){ 
            conv_mle[i,]<-1
            break 
          }
          
          log_like<-sum(log(Pk)) - sum(log(pmax(P0, 1e-12)))
          P0<-Pk
          th_i<-th_new
          
          #  End NR algorithm if all 3 convergence criteria are met
          # L-infinity norm employed for multidimensional theta vector
          if(k > 1 && abs(log_like) < tol && max(abs(grad))<score.tol){break}
          
          # Check if convergence within # of iterations
          if(k == iter){conv_mle[i,]<-1}
        }
        
        if(any(conv_mle[i,] == 1)){
          theta_mle[i,]<-rep(NA, L)
        }else if(any(th_i>up.bound) || any(th_i<low.bound)){
          # Replace any converged estimates outside the bounds
          theta_mle[i,]<-pmax(low.bound, pmin(up.bound, th_i))
          conv_mle[i,]<-2
        }else{
          theta_mle[i,]<-th_i
        }
        
        # Calculate item-level weights and residuals
        if(any(is.na(th_i))){
          r_final[i,]<-w_final[i,]<-rep(NA, J)
        }else{
          w_final[i,]<-as.numeric(person.weights(matrix(th_i, 1), dat[i,], "MGRM", ipars_use))
          r_final[i,]<-as.numeric(residual(matrix(th_i, 1), "MGRM", ipars_use,
                                           matrix(dat[i,], 1), resid, D=D))
        }
      }
      
      out$theta_MLE<-theta_mle
      
      # Compute SEs for converged subjects
      conv_idx<-which(apply(conv_mle, 1, function(x) all(x !=1)) & 
                          apply(theta_mle, 1, function(x) all(!is.na(x))))
      
      ase_MLE<-matrix(NA, N, L)
      sand_MLE<-matrix(NA, N, L)
      
      if(length(conv_idx) > 0){
        # make sure formatting of thetas is proper
        if(length(conv_idx)==1){
          th_converged<-matrix(theta_mle[conv_idx,], nrow = 1)
        }else{
          th_converged<-theta_mle[conv_idx,]
        }
      
        se.all<-standard.errors(th_converged, ipars_use, dat[conv_idx, , drop = FALSE], model, D, 
                                  weight.type, tuning.par, custom.weights, resid, "MLE", prior, eap.quad.pts)
        ase_MLE[conv_idx,]<-se.all$asymptotic_MLE
        sand_MLE[conv_idx,]<-se.all$sandwich_MLE
      }
      out$ase_MLE = ase_MLE
      out$sandwich_MLE = sand_MLE
      out$singular_flag = sing_flag
      out$convergence_MLE<-conv_mle
      out$weights_MLE<-w_final
      out$residuals_MLE<-r_final
    }
    
  } # end MLE block
  
  ##### MAP (Rasch, 1PL, 2PL, GRM) #####
  if("MAP" %in% est.type){
    
    deriv_fn<-if(model=="GRM"){grm_derivs}else{dich_derivs}
    
    theta_map<-rep(NA, N)
    conv_map<-rep(0, N)
    w_final<-matrix(NA, N, J)
    r_final<-matrix(NA, N, J)
    
    for(i in 1:N){
      th_i<-get_init(i, 1)
      P0<-0
      
      for(k in 1:iter){
        w_i<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
        dv<-deriv_fn(th_i, dat[i,], ipars_use, w_i, bayes = TRUE)
        
        if(!is.finite(dv$D1) || !is.finite(dv$D2) || abs(dv$D2) < 1e-14){
          conv_map[i]<-1
          break
        }
        
        th_new<-th_i-dv$D1 / dv$D2
        if(!is.finite(th_new)){ 
          conv_map[i]<-1
          break 
        }
        
        # Convergence check with change in log-posterior
        P_cur<-if(model == "GRM"){
          pmax(dv$Pk, 1e-12)
        }else{
          pmax(ifelse(dat[i,] == 1, dv$P_i, 1 - dv$P_i), 1e-12)
        }
        log_like<-sum(log(P_cur)) - sum(log(pmax(P0, 1e-12)))
        P0<-P_cur
        th_i<-th_new
        
        # End NR algorithm if all 3 convergence criteria are met
        if(k > 1 && abs(log_like) < tol && abs(dv$D1)<score.tol){break}
        
        # Check for reaching max # of iterations
        if(k == iter){
           conv_map[i]<-1
           th_i<-NA
        }
      }
    #Check for convergence outside bounds (convergence code = 2)
    if(!is.na(th_i) & th_i<low.bound){
      th_i<-low.bound
      conv_mle[i]<-2
    }else if(!is.na(th_i) & th_i > up.bound){
      th_i<-up.bound
      conv_mle[i]<-2
    }
      # Calculate final item-level weights and residuals
      if(is.na(th_i)){
        r_final[i,]<-w_final[i,]<-rep(NA, J)
      }else{
        w_final[i,]<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
        r_final[i,]<-as.numeric(residual(matrix(th_i, 1), model, ipars_use, matrix(dat[i,], 1), resid, D=D))
      }
      
      theta_map[i]<-th_i
    }
    
    out$theta_MAP<-matrix(theta_map, ncol = 1)
    
    # Compute SEs for converged subjects
    conv_idx<-which(conv_map!=1 & !is.na(theta_map))
    post_sd_MAP<-rep(NA, N)
    sand_MAP<-rep(NA, N)
    
    if(length(conv_idx) > 0){
      se.all<-standard.errors(theta_map[conv_idx], ipars_use, dat[conv_idx, , drop = FALSE], model, D, 
                                weight.type, tuning.par, custom.weights, resid, "MAP", prior, eap.quad.pts)
      post_sd_MAP[conv_idx]<-se.all$post_sd_MAP
      sand_MAP[conv_idx]<-se.all$sandwich_MAP
    }
    out$post_sd_MAP<-matrix(post_sd_MAP, ncol = 1)
    out$sandwich_MAP <-matrix(sand_MAP, ncol = 1)
    
    out$convergence_MAP<-matrix(conv_map, ncol = 1)
    out$weights_MAP<-w_final
    out$residuals_MAP<-r_final
  } # end MAP block
  
  ##### EAP (Rasch, 1PL, 2PL) #####
  if("EAP" %in% est.type){
    
    f_x<-dnorm(eap.quad.pts, mu, sqrt(sigma2))
    probs_q<-item.prob(eap.quad.pts, model, ipars_use, D) 
    Q<-length(eap.quad.pts)
    
    theta_eap<-rep(NA, N)
    w_final<-matrix(NA, N, J)
    r_final<-matrix(NA, N, J)
    
    for(i in 1:N){
      th_i<-get_init(i, 1)[1]
      P0<-0
      
      for(k in 1:iter){
        w_i<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
        
        # Weighted likelihood at each quadrature point
        lik_q<-apply(probs_q, 1, function(pj)
          prod((pj^dat[i,]*(1 - pj)^(1 - dat[i,]))^w_i))
        
        denom<-sum(lik_q*f_x)
        if(!is.finite(denom) || denom == 0) break
        
        th_new<-sum(eap.quad.pts*lik_q*f_x) / denom
        
        P_cur<-pmax(ifelse(dat[i,] == 1,
                                item.prob(th_i, model, ipars_use, D),
                                1 - item.prob(th_i, model, ipars_use, D)), 1e-12)
        log_like<-sum(log(P_cur)) - sum(log(pmax(P0, 1e-12)))
        P0<-P_cur
        th_i<-th_new
        
        if(k > 1 && abs(log_like) < tol) break
      }
      
      w_final[i,]<-w_i_f<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
      r_final[i,]<-as.numeric(residual(matrix(th_i, 1), model, ipars_use, matrix(dat[i,], 1), resid, D=D))
      lik_f<-apply(probs_q, 1, function(pj) prod((pj^dat[i,]*(1 - pj)^(1 - dat[i,]))^w_i_f))
      denom_f<-sum(lik_f*f_x)
      
      if(is.finite(denom_f) && denom_f > 0){
        th_eap_i<-sum(eap.quad.pts*lik_f*f_x) / denom_f
        
        theta_eap[i]<-pmax(low.bound, pmin(up.bound, th_i))
      }
    }
    
    out$theta_EAP<-matrix(theta_eap, ncol = 1)
    
    # Compute SEs for converged subjects
    conv_idx<-which(conv_map==0 & !is.na(theta_eap))
    ase_EAP<-rep(NA, N)
    sand_EAP<-rep(NA, N)
    
    if(length(conv_idx) > 0){
      se.all<-standard.errors(theta_eap[conv_idx], ipars_use, dat[conv_idx, , drop = FALSE], model, D, 
                                weight.type, tuning.par, custom.weights, resid, "EAP", prior, eap.quad.pts)
      ase_EAP[conv_idx]<-se.all$asymptotic_EAP
      sand_EAP[conv_idx]<-se.all$sandwich_EAP
    }
    out$post_sd_EAP<-matrix(ase_EAP, ncol = 1)
    out$sandwich_EAP <-matrix(sand_EAP, ncol = 1)
    
    out$weights_EAP<-w_final
    out$residuals_EAP<-r_final
  } # end EAP block
  
  ##### WLE (Rasch, 1PL, 2PL) #####
  if("WLE" %in% est.type){
    
    # Unidimensional dichotomous models
    if(model %in% dich.mods){
      
      # Initialize output for storing
      theta_wle<-rep(NA, N)
      conv_wle<-rep(0, N)
      w_final<-matrix(NA, N, J)
      r_final<-matrix(NA, N, J)
      
      for(i in 1:N){
        th_i<-get_init(i, 1)
        P0<- item.prob(th_i, model, ipars_use, D)
        
        for(k in 1:iter){
          w_i<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
          dv<-dich_derivs(th_i, dat[i,], ipars_use, w_i, wle=TRUE)
          
          if(!is.finite(dv$D1) || !is.finite(dv$D2) || abs(dv$D2) < 1e-14){
            conv_wle[i]<-1
            break
          }
          
          th_i<-th_i-dv$D1/dv$D2
          
          if(!is.finite(th_i)){ 
            conv_wle[i]<-1
            break 
          }
          
          P1<-item.prob(th_i, model, ipars_use, D)
          # Compute difference in log-likelihood for convergence criterion
          log_like <- sum(dat[i,]*log(P1)+(1-dat[i,])*log(1-P1)) - sum(dat[i,]*log(P0)+(1-dat[i,])*log(1-P0))
          
          # Check for NaNs in log_like: if so, record nonconvergence
          if (is.nan(log_like)) {
            th_i<- NA
            conv_wle[i] <- 1
            break
          }
          
          # End NR algorithm if all 3 convergence criteria are met
          if(k > 1 && abs(log_like) < tol && abs(dv$D1) < score.tol){break}
          
          # Check for reaching max # of iterations
          if(k == iter){
            conv_wle[i]<-1
            th_i<-NA
          }
          
          # Update initial probability (P0) for the next iteration
          P0<-P1
        
        }
        #Check for convergence outside bounds (convergence code = 2)
        if(!is.na(th_i) & th_i<low.bound){
          th_i<-low.bound
          conv_wle[i]<-2
        }else if(!is.na(th_i) & th_i > up.bound){
          th_i<-up.bound
          conv_wle[i]<-2
        }
        
        # Final weights / residuals at converged theta
        theta_wle[i]<-th_i
        if(is.na(th_i)){
          r_final[i,]<-w_final[i,]<-rep(NA, J)
        }else{
          w_final[i,]<-as.numeric(person.weights(th_i, dat[i,], model, ipars_use))
          r_final[i,]<-residual(th_i, model, ipars_use, matrix(dat[i,], 1), resid=resid, D=D)
        
        }
        
        
        
      } # end person loop
      
      out$theta_WLE<- matrix(theta_wle, ncol = 1)
      
      # Compute SEs for converged subjects
      conv_idx<-which(conv_wle!=1 & !is.na(theta_wle))
      
      ase_WLE<-rep(NA, N)
      sand_WLE<-rep(NA, N)
      
      if(length(conv_idx) > 0){
        
        # need to update standard errors to include WLE sandwich SEs
        se.all<-standard.errors(theta_wle[conv_idx], ipars_use, dat[conv_idx, , drop = FALSE], model, D, 
                                weight.type, tuning.par, custom.weights, resid, "MLE", prior, eap.quad.pts)
        ase_WLE[conv_idx]<-se.all$asymptotic_MLE
        sand_WLE[conv_idx]<-se.all$sandwich_MLE
      }
      
      out$ase_WLE<-matrix(ase_WLE, ncol=1)
      out$sandwich_WLE<-matrix(sand_WLE, ncol=1)
      out$convergence_WLE<-matrix(conv_wle, ncol = 1)
      out$weights_WLE<-w_final
      out$residuals_WLE<-r_final
    }
    
    
    
  } # end WLE block
  
  return(out)
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
#' The robust marginal maximum likelihood (RMML) estimator (Hong & Cheng, 2019) provides estimates of item parameters for various IRT models when both item and latent trait parameters are unknown. The robust algorithm is designed to reduce the bias in item parameter estimates due to aberrant responding.
#' A weight, \eqn{w_i}, is assigned to each examinee such that smaller weights reflect larger model misfit for the examinee. The weights are then applied to the likelihood according to 
#' \eqn{L_w(\boldsymbol{x}|\boldsymbol{\gamma}) = \prod_{i=1}^N L(\boldsymbol{x}_i|\boldsymbol{\gamma})^{w_i},}
#' where \eqn{L(\boldsymbol{x}_i|\boldsymbol{\gamma})} is the likelihood of the response vector \eqn{\boldsymbol{x}_i} for examinee \eqn{i} (\eqn{i=1,...,N}) given \eqn{\boldsymbol{\gamma},} the matrix of item parameters for the test. The EM algorithm proceeds, maximizing this weighted likelihood.
#' 
#' Hong and Cheng (2019) derive weights based on the person-fit statistic \eqn{l_z^{`*`}}, where lower values indicate greater misfit between the response pattern and the model (Snijders, 2001). Since \eqn{l_z^{`*`}} follows the asymptotic standard normal distribution, weights were derived as the normalized \emph{p}-value of the residual under a one-sided alternative hypothesis, e.g.,
#' \eqn{w_i = \frac{\Phi(l^{`*`}_{zi})}{\sum_{i=1}^N \Phi(l^{`*`}_{zi})},}
#' where \eqn{\Phi(\cdot)} is the CDF of the standard normal distribution. Note that the user may define different weights with the argument \code{survey.weights}.
#' 
#' @param dat A \eqn{N\times J} matrix of response data
#' @param survey.weights Optional argument to assign user-specified weights. Default weights are calculated based on \eqn{l^{`*`}_{zi}.}
#' @param ... Additional arguments to be passed to \code{mirt()}
#' 
#' @references Hong, M., & Cheng, Y. (2019). Robust maximum marginal likelihood (RMML) estimation for item response theory models. Behavior Research Methods, 51(2), 573–588. https://doi.org/10.3758/s13428-018-1150-4
#' @references Snijders, T. A. (2001). Asymptotic null distribution of person fit statistics with estimated person parameter. Psychometrika, 66(3), 331-342.
#' 
#' @return An object of class \code{mirt}
#' 
#' @examples
#' # Load package and example data set 
#' library(mirt) 
#' data(Science) 
#' # Robust estimation of item parameters
#' robust.item(Science)
#' 
#' Multidimensional GRM example: 5 Factors
#' data(BFI2)
#' robust.item(BFI2[,20:79], model=5, TOL=0.001, method="QMCEM")
#' @export
robust.item<-function(dat, survey.weights=NULL, ...){
  
  if(is.null(survey.weights)){
      
    # Initial Model Estimation 
    mod <- mirt(dat, ...) 
    # Person fit residual calculation
    per.fit <- personfit(mod, method = 'ML')$Zh 
    # Weight  
    weight <- pnorm(per.fit)*nrow(dat)/ sum(pnorm(per.fit)) 
    # Robust model estimation 
    robust.mod <- mirt(dat, survey.weights=weight, ...)
  }else{
    robust.mod <- mirt(dat, survey.weights=survey.weights, ...)
  }
  
  return(robust.mod)
}

#' Robust Estimation of Working Speed
#'
#' To reduce the influence of spuriously slow or spuriously fast response times on the estimate of an examinee's working speed (\eqn{\tau}), robust M-estimates can be obtained (Hong et al., 2020).
#' Under the log-normal model (van der Linden, 2006), the log response time (\eqn{t_j}) for an examinee on item \eqn{j} (\eqn{j=1,2,...J}), follows 
#' \eqn{\log(t_j) | \tau ~ N(\Beta_j - \tau, \frac{1}{\alpha_j^2})},
#' assuming independence of log item response times for the examinee. The probability density function for the model is then given by
#' \eqn{f(t_i|\tau) = \frac{\alpha_j}{t_j \sqrt{2\pi}} \exp{(-\frac{1}{2}(a_j(\log t_j - (\Beta_j - \tau)))^2)}}
#' The robust maximum likelihood estimate of \eqn{\tau} is 
#' \eqn{\hat{\tau}^{RML} = \frac{\sum_{j=1}^J w(r_j) \alpha_j^2 (\Beta_j - \log t_j)}{\sum_{j=1}^J w(r_j) \alpha_j^2 }}
#' where weights \eqn{w(\cdot)} are defined by either the Huber (Huber, 1981) or bisquare (Mosteller & Tukey, 1977) weight functions (see \code{huber()} and \code{bisquare()} functions for more detail).
#' A residual capturing the difference between the observed and model-expected response time is calculated with \eqn{r_j = \alpha_j(\log(t_{ij}) - \Beta_j - \tau_i)}.
#' Although \eqn{\hat{\tau}^{RML}} yields a closed-form solution, the residuals and corresponding weights are updated based on the previously estimated \eqn{\hat{\tau}^{RML}} in an iterative manner, with the initial value set at the maximum likelihood estimate (where all weights are fixed at 1).
#' Convergence occurs when the absolute change in \eqn{\hat{\tau}^{RML}} between two iterations is within the specified tolerance.
#' 
#' @param dat A \eqn{N\times J} matrix of numerical response time data
#' @param ipars Item-level parameters for the log-normal model, structured with the first column containing discrimination parameters \eqn{\alpha_j} and the second column containing time-intensity parameters \eqn{\Beta_j} for \eqn{j=1,...,J}.
#' @param weight.type Weighting scheme: \code{"equal"} (default), \code{"Huber"},
#'   \code{"bisquare"}, or \code{"custom"}.
#' @param tuning.par Tuning parameter for Huber or bisquare weights.  Required when
#'   \code{weight.type} is not \code{"equal"} or \code{"custom"}.
#' @param custom.weights An \eqn{N \times J} numeric matrix of user-specified weights
#'   in \eqn{[0, 1]}. Required when \code{weight.type = "custom"}.
#' @param iter Maximum number of iterations. Default is 30.
#' @param tol Convergence tolerance on the absolute change in working speed estimate between two iterations.
#'   Default is 0.0001.
#' @details The working speed parameter, \eqn{tau}, captures how "speedily" a subject responds. 
#'          In effect, a larger \eqn{tau} corresponds to quicker responses, while lower \eqn{tau}s correspond to slower responses.
#'          When an anomalous handful of response times are spuriously quicker than others, the MLE of \eqn{tau} may become inflated (suggesting more speediness), and the robust estimate counteracts this bias, yielding a smaller estimate.
#'          Conversely, spuriously slower responses may yield a smaller MLE than the true working speed, while the robust estimate works to increase the MLE closer to the truth.
#'          For identifiability purposes, \eqn{tau} is assumed to have a mean of 0.
#' 
#' @references Huber, P. (1981). \emph{Robust Statistics}. John Wiley & Sons, Inc.
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}.
#'   Addison-Wesley Publishing Company.
#' @return A list containing
#' \describe{
#'   \item{\code{tau}}{ (\eqn{N \times 1}) matrix of working speed estimates.}
#'   \item{\code{residuals}}{ (\eqn{N \times J}) matrix of item-level residuals (\eqn{r_j}) capturing the difference between the observed and model-expected response time for each individual.}
#'   \item{\code{nonconvergence}}{ (\eqn{N \times 1}) Nonconvergence indicator (0 = converged, 1 = not converged).}
#' }
#' @examples
#' 
#' # Example with Canadian PISA data
#' library(pisaRT)
#' dat<- pisaW[,26:37]
#' # Estimate item parameters for the log-normal model under normal responding
#' library(LNIRT)
#' mod.ln <- LNRT(RT = dat, data = dat)$Post.Means
#' ipars<-cbind(mod.ln$Time.Discrimination, mod.ln$Time.Intensity)
#' # Reduce response time by 1/4 for people as if they had item preknowledge on items 2, 4, 6, 8
#' dat[,c(2, 4, 6, 8)]<-dat[,c(2, 4, 6, 8)] - log(4)
#' 
#' # Robust estimation of tau with aberrant data
#' out.rt.ab<-robust.rt(dat = exp(dat), ipars, weight.type = "bisquare", tuning.par = 4)
#' 
#' # Robust estimation of tau with normal response data
#' out.rt<-robust.rt(dat = exp(dat), ipars, weight.type = "equal")
#' head(cbind(out.rt.ab$tau, out.rt$tau))
#' 
#' # Robust estimates are lower, while MLEs suggest speedier responses
#' 
#' @export

robust.rt<-function(dat, ipars, weight.type = "equal", tuning.par = NULL, custom.weights = NULL, iter=30, tol=0.0001){
  
  ##### Ensure proper input #####
  if(!(weight.type %in% c("equal", "Huber", "bisquare", "custom")))
    stop(paste(weight.type, "is not a supported weight.type."))
  
  if(weight.type %in% c("Huber", "bisquare") && is.null(tuning.par))
    stop(paste("tuning.par must be supplied when weight.type = ", weight.type))
  
  if(weight.type == "custom" && is.null(custom.weights))
    stop("custom.weights must be supplied when weight.type = 'custom'.")
  
  ##### Internal notation and functions #####
  J<-ncol(dat)
  N<-nrow(dat)
  
  alphas<-ipars[,1]
  betas<-ipars[,2]
  
  # Item-level weight given residual vector
  compute.weights<-function(r_mat){
    if(weight.type == "equal") return(matrix(1, nrow(r_mat), ncol(r_mat)))
    if(weight.type == "Huber") return(huber(r_mat, tuning.par))
    if(weight.type == "bisquare") return(bisquare(r_mat, tuning.par))
    if(weight.type == "custom"){
      if(is.null(custom.weights)) stop("custom.weights must be supplied when weight.type = 'custom'.")
      return(custom.weights)
    } 
  }
  
  # Matrices for output
  tau<-matrix(NA, N)
  convergence<-matrix(1, N)
  residual_1<-matrix(NA, N, J)
  
  ##### Begin Computation #####
  for(i in 1:N){
    dat_i<-dat[i,]
    tau0<- sum(alphas^2*(betas-log(dat_i)))/sum(alphas^2)
    
    for(m in 1:iter){
      res<-alphas*(log(dat_i)-(betas-tau0))
      wgts<-compute.weights(as.matrix(res))
      tau1<- sum(wgts*alphas^2*(betas-log(dat_i)))/sum(wgts*alphas^2)
      
      if(abs(tau1-tau0)<tol){
        convergence[i,]<-0
        break
      }
      tau0<-tau1
    }
    tau[i]<-tau1
    residual_1[i,]<-as.matrix(alphas*(log(dat_i)-(betas-tau0)))
    
  } # end person loop of computations
  
  return(list(tau = tau, residuals = residual_1, nonconvergence = convergence))
}

#' Plot histogram of residuals along plot of the weight function
#'
#' Plot a histogram of residuals along the graph of the weighting function (dependent on the tuning parameter) as a function of the residual
#' @param r A vector of residuals
#' @param H Huber tuning parameter
#' @param B Bisquare tuning parameter
#' @param x.axis A vector containing two values: the lower and upper limits for the x-axis. Default is the minimum and maximum values of the residuals supplied in \code{r}, excluding infinite values.
#' @details This function allows the user to visualize the amount of downweighting applied to their data in robust latent trait estimation based on the tuning parameter and item-level residuals, in order to choose a tuning parameter that suits their data well.
#'               The values \eqn{H=1} and \eqn{B=4} are recommended for Huber (Huber, 1981) and bisquare (Mosteller & Tukey, 1977) weighting, respectively, but the amount of downweighting can be made more conservative or more liberal depending on the research goals.
#'               To increase the downweighting applied in estimation, use a smaller tuning parameter. To decrease the amount of downweighting, apply a greater tuning parameter.
#'               The histogram of residuals is plotted below a plot of the weight curve(s) at the tuning parameter. If \emph{H} is supplied to the function, the Huber curve is plotted, and, if \emph{B} is supplied, Tukey's bisquare weight is plotted. If both \emph{H} and \emph{B} are supplied, the weights are plotted simultaneously.
#'               If \emph{H} is supplied, dotted vertical lines at \emph{H} and \emph{-H} reveal the cut between data that is downweighted (\eqn{|r_ij|>H}) and data that is not downweighted (\eqn{|r_ij|<H}), e.g., receives full weighting (\eqn{w_ij=1.0}).
#'               If \emph{B} is supplied, dotted vertical lines at \emph{B} and \emph{-B} reveal the cut between data that receives some weight (\eqn{|r_ij|<B}) and data that is removed (\eqn{|r_ij|>B}), e.g., receives a weight of 0.
#'               If no tuning parameter is supplied, the histogram of residuals is provided.
#' @references Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' @return Histogram plot of residuals beneath a graph of the weight functions vs. the residuals.
#' @examples
#' 
#' ### Visualize Huber weight function
#' ##  Use residual() function
#' thetas<-c(-2, -1, 0, 1, 2)
#' ipars<-data.frame(a = c(1, 1.2, 1.5, 1, 1.2, 1.5, 1, 1.2, 1.5), b = c(-2, -1.5, -1, -.5, 0, .5, 1, 1.5, 2))
#' info_r<-residual(thetas, "2PL", ipars, resid = "information")
#' choose.tuco(info_r, H=1)
#' 
#' 
#' ### Compare Huber and bisquare weights
#' ##  Use robust.theta() output
#' data(BFI2)
#' dat<-BFI2[,20:79]
#' specific<- rep(1:5, 12) #which factor each item loads on
#' mod<-mirt(dat, specific, TOL=0.001, method="QMCEM")
#' ipars<-coef(mod, simplify=T)$items
#' ipars<-cbind(ipars[,1:5], -ipars[,6:9]/rowSums(ipars[,1:5])) # reparameterize
#' out <- robust.theta(dat, ipars, model= "MGRM", resid = "standardized", dimen=5, D=1)
#' 
#' choose.tuco(r=out$residuals_MLE, H=1, B=4, c(-6, 6))
#' @export

choose.tuco<-function(r, H=NULL, B=NULL, x.axis=NULL){
  
  # r is a vector or matrix of residuals
  residuals<-data.frame(Residual =c(r[!is.infinite(r)& !is.nan(r)]))
  
  if(is.null(x.axis)){
    x.axis<- c(min(residuals$Residual, na.rm =T), max(residuals$Residual, na.rm =T))
  }
  # Histogram of residuals
  hist.out<-ggplot(residuals, aes(x=Residual))+
    geom_histogram(aes(y = after_stat(density)), bins=50)+ 
    xlim(x.axis[1], x.axis[2])
    ylab("Density")
  weight.out<-ggplot()+
    xlim(x.axis[1], x.axis[2] ) + ylab("Weight")+
    theme(legend.position = c(.9, .74))+
    ggtitle("Weights Applied in Estimation")
  
  if(!is.null(H)){
    hist.out<- hist.out + 
      geom_vline(xintercept = -H, linetype="dashed", color = "firebrick")+ 
      geom_vline(xintercept = H, linetype="dashed", color = "firebrick")
    weight.out<-weight.out+
      stat_function(fun=function(x) huber(x, H), aes(colour = "Huber"))+
      geom_vline(xintercept = -H, linetype="dashed", color = "firebrick")+ 
      geom_vline(xintercept = H, linetype="dashed", color = "firebrick")
  }
  
  if(!is.null(B)){
    hist.out<- hist.out +  
      geom_vline(xintercept = -B, linetype="dashed", color = "darkcyan")+ 
      geom_vline(xintercept = B, linetype="dashed", color = "darkcyan")
    weight.out<-weight.out+
      stat_function(fun=function(x) bisquare(x, B),  aes(colour = "Bisquare"))+ 
      geom_vline(xintercept = -B, linetype="dashed", color = "darkcyan")+ 
      geom_vline(xintercept = B, linetype="dashed", color = "darkcyan")
  }
  
  if(!is.null(H) & !is.null(B)){
    
    weight.out<-weight.out+
      scale_color_manual(name = "Function", breaks=c('Bisquare', 'Huber'), values=c('Bisquare'="darkcyan", 'Huber'='firebrick')) 
    return(do.call(ggarrange, c(list(weight.out+xlab(NULL), hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  
  }else if(is.null(H) & !is.null(B)){
    
    weight.out<-weight.out+
      scale_color_manual(name = "Function", breaks=c('Bisquare'), values=c('Bisquare'="darkcyan"))
    return(do.call(ggarrange, c(list(weight.out+xlab(NULL), hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  
  }else if(is.null(B) & !is.null(H)){
    
    weight.out<-weight.out+
      scale_color_manual(name = "Function", breaks=c('Huber'), values=c('Huber'='firebrick')) 
    return(do.call(ggarrange, c(list(weight.out+xlab(NULL), hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  
  }else{ # when no tuning parameter is specified, only return histogram of residuals
    return(hist.out+ggtitle("Histogram of Residuals"))
  }
}

#' Plot to compare robust latent trait estimates with non-robust estimates
#'
#' Robust latent trait estimates are plotted against the corresponding non-robust latent trait estimates for each subject and compared against the identity line for an aggregate sample-level analysis. Summary statistics, including the Euclidean distance between the two estimates for the individual, are provided for a subject-level analysis.
##' @param dat An \eqn{N \times J} matrix of item responses (\eqn{N} subjects,
#'   \eqn{J} items). Responses for dichotomous models must be 0/1. Responses for
#'   polytomous models must be integer-valued starting at 1. Missing data is not supported.
#' @param ipars A matrix of item parameters structured identically to the
#'   \code{ipars} argument of \code{\link{item.prob}}.
#' @param model The IRT model to be used: \code{"Rasch"}, \code{"1PL"}, \code{"2PL"},
#'   \code{"MIRT"}, \code{"GRM"}, or \code{"MGRM"}. See \code{\link{item.prob}}
#'   for details on each model.
#' @param weight.type Weighting scheme for the robust estimates: \code{"Huber"},
#'   \code{"bisquare"}, or \code{"custom"}.
#' @param tuning.par Tuning parameter for Huber or bisquare weights. Required when
#'   \code{weight.type} is \code{"Huber"} or \code{"bisquare"}.
#' @param same.plot.dim Should estimates across all dimensions (if the model is multidimensional) be aggregated and plotted jointly on one plot? Default is \code{FALSE}.
#' @param ids Optional vector of identification strings to link to summary statistics output.
#' @param ab.ids Optional vector of identification strings for subjects who have been flagged or suspected of anomalous response behavior. The latent trait estimates for these subjects are plotted in red. If no IDs are supplied in \code{ids}, the \code{ab.ids} should correspond to the row number of the flagged subjects in \code{dat}.
#' @param ... Additional arguments to be passed to \code{robust.theta()}
#' @details When the data does not contain anomalous responses, robust estimates should not differ greatly from the non-robust estimate for a given estimation method. As such, when the robust estimates are plotted against their non-robust counterpart, the data should roughly follow the identity line.
#'          When this pattern deviates from the identity line, aberrant responses may be present.
#'          For instance, if the robust estimate tends to be larger in magnitude than the non-robust estimate, it is likely that the algorithm is reducing the impact of noise such as random responses.
#'          This plot can reveal aggregate sample-level trends, as well as highlight extreme anomalies (e.g., a respondent whose robust estimate deviates extremely from its non-robust estimate, moreso than the remaining sample), which can be further checked in the \code{Summary Statistics} output.
#'          Under \code{Summary Statistics}, the Euclidean distance between the robust and non-robust latent trait vector estimates is reported, along with the corresponding estimates and the raw data. 
#'          The table is sorted by decreasing Euclidean distance to highlight subjects with the most extreme discrepancies, including those that may appear anomalous on the plot. Investigation into the response patterns of these individuals may provide insight into why these discrepancies occur.  
#'          Note that person fit statistics are preferred for the purpose of detecting aberrant respondents. A larger Euclidean distance may not necessarily mean that the response pattern contains anomalies. However, the discrepancy between the robust and non-robust latent trait estimates may highlight particular individuals for further investigation and/or confirmation that robust estimation is effective in the sample.
#'          
#' @return `Summary Statistics` A data frame containing the ID (if supplied in \code{ids}), Euclidean distance between the robust and non-robust latent trait estimate(s), robust latent trait estimate(s), non-robust latent trait estimate(s), and response vector for each respondent. The subjects are organized by greatest to least Euclidean distance.
#' @return `plot` If \code{same.plot.dim = TRUE}, a plot of robust estimates against non-robust latent trait estimates across subjects and dimensions. If \code{same.plot.dim=FALSE}, a different plot is returned for each dimension, labeled as Dimension 1, Dimension 2, etc. 
#' @examples
#' 
#' # 2Pl Example
#' 
#' # MGRM Example: Different Plots
#' library(mirt)
#' data(BFI2)
#' dat<-BFI2[,20:79]
#' specific<- rep(1:5, 12) #which factor each item loads on
#' mod<-mirt(dat, specific, TOL=0.001, method="QMCEM")
#' ipars<-coef(mod, simplify=T)$items
#' ipars<-cbind(ipars[,1:5], -ipars[,6:9]/rowSums(ipars[,1:5])) # reparameterize
#' out <- robust.theta.comparisons(dat, ipars, model= "MGRM", weight.type = "bisquare", tuning.par = 4, resid = "standardized", dimen=5, D=1)
#' head(out$`Summary Statistics`)
#' out$`Dimension 1`
#' out$`Dimension 2`
#' 
#' @export
   
robust.theta.comparisons<-function(dat, ipars, model, weight.type, tuning.par, est.type = "MLE", same.plot.dim = FALSE, ids = NULL, ab.ids = NULL, ...){
  
  th.est.type<-paste0("theta_", est.type)
  ##### Latent trait estimation
  
  # Equal weight
  th.eq<-robust.theta(dat, ipars, model, weight.type = "equal", tuning.par = NULL, ...)
  
  if(weight.type == "custom"){
    # Custom weights
    th.rob<-robust.theta(dat, ipars, model, weight.type = "custom", tuning.par = NULL, ...)
    
  }else{
    # Robust estimation
    th.rob<-robust.theta(dat, ipars, model, weight.type = weight.type, tuning.par = tuning.par, ...)
    
  }
  
  ##### Create summary statistic table
  
  #Euclidean distance between robust and nonrobust estimate
  euc.dist <- sqrt(rowSums(as.matrix((th.eq[[th.est.type]] - th.rob[[th.est.type]])^2)))
  
  sum.stats<-data.frame(Distance=euc.dist,
                        Robust=th.rob[[th.est.type]],
                        Equal=th.eq[[th.est.type]],
                        dat)
  
  # If IDs are supplied, create a column with them
  if(!is.null(ids)){
    sum.stats<-cbind(ID = ids, sum.stats)
  }else{
    ids<-1:nrow(dat)
  }
  
  # Order by Euclidean distance
  sum.stats<-sum.stats[order(sum.stats$Distance, decreasing = T),]
  
  ##### Create plots
  
  if(same.plot.dim){ # If estimates for all dimensions of the latent trait are to be plotted jointly on one plot
    df<-data.frame(x=c(th.eq[[th.est.type]]), y=c(th.rob[[th.est.type]]))
    
    if(!is.null(ab.ids)){
      df<-cbind(df, Status = ifelse(ids%in%ab.ids, "Flagged", "Not Flagged"))
    } else{
      df<-cbind(df, Status = "Not Flagged")
    }
    
    out<-vector("list", 2)
    names(out) <- c("plot", "Summary Statistics")
    
    p<-ggplot(df, aes(x=x, y=y, color = Status))+
      geom_point()+
      theme_classic()+
      scale_colour_manual(values = setNames(c('red','black'),c("Flagged", "Not Flagged")))+
      geom_abline(slope = 1, intercept = 0)+
      ggtitle("Robust vs. Non-Robust Estimates")+
      labs(x="Equal-Weight Estimates", y="Robust Estimates")
    
    if(is.null(ab.ids)){
      p<-p+guides(color = "none")
    }
    
    out[[1]]<-p
    
  }else{ # If estimates are to be plotted on separate plots for different latent trait dimensions
    
    dimen<-ncol(th.rob[[th.est.type]])
    out<-vector("list", dimen+1)
    names(out) <- c(paste0("Dimension ", 1:dimen), "Summary Statistics")
    
    for(i in 1:dimen){
      df<-data.frame(x=th.eq[[th.est.type]][,i], y=th.rob[[th.est.type]][,i])
      
      if(!is.null(ab.ids)){
        df<-cbind(df, Status = ifelse(ids%in%ab.ids, "Flagged", "Not Flagged"))
      } else{
        df<-cbind(df, Status = "Not Flagged")
      }
      
      p<-ggplot(df, aes(x=x, y=y, color = Status))+
        geom_point()+
        theme_classic()+
        scale_colour_manual(values = setNames(c('red','black'),c("Flagged", "Not Flagged")))+
        geom_abline(slope = 1, intercept = 0)+
        ggtitle(paste0("Robust vs. Non-Robust Estimates: Dimension ", i))+
        labs(x="Equal-Weight Estimates", y="Robust Estimates")
      
      if(is.null(ab.ids)){
        p<-p+guides(color = "none")
      }
      
      out[[i]] <- p
    }
  }
  
  out[[length(out)]]<-sum.stats
  
  return(out)
}

#' Comparing Robust with MML Item Parameter Estimates 
#' 
#' The following function returns a plot comparing the robust estimate (Hong & Cheng, 2019) of an item parameter against the corresponding marginal maximum likelihood (MML) estimate for all items.
#' The item parameters are labeled according to their label in \code{mirt()}.
#' 
#' @param dat A \eqn{N\times J} matrix of response data
#' @param ... Additional arguments to be passed to \code{mirt()}
#' 
#' @references Hong, M., & Cheng, Y. (2019). Robust maximum marginal likelihood (RMML) estimation for item response theory models. Behavior Research Methods, 51(2), 573–588. https://doi.org/10.3758/s13428-018-1150-4
#' @details When the model fits the data, item parameters estimated with RMML should be approximately equivalent to those estimated with MML, thereby roughly following the identity line. 
#'          When the two estimates differ greatly for a parameter on an item, aberrant responding may be present such that its effects are mitigated by robust weighting. 
#'          The plots returned by \code{robust.item.plots()} can therefore be used to spot items that are greatly affected by robust estimation, for further investigation.
#' 
#' @return A list containing one plot for each item parameter, structured as \code{`Parameter Estimates for a1`}, \code{`Parameter Estimates for d`}, etc., in addition to
#' \describe{
#'   \item{\code{`MML Item Parameter Estimates`}}{ Matrix of item parameter estimates obtained with MML as defined by \code{mirt()}, each row corresponding to an item.}
#'   \item{\code{`RMML Item Parameter Estimates`}}{ Matrix of item parameter estimates obtained with robust MML as defined by \code{mirt()}, each row corresponding to an item.}
#' }
#' 
#' 
#' @examples
#' 
#' # SAT12 Data: 3-Dimensional MIRT Model
#' SAT12[SAT12 == 8] <- NA #set 8 as a missing value
#' # Correct answer key
#' library(mirt)
#' key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5) 
#' scoredSAT12 <- key2binary(SAT12, key)
#' specific <- c(2, 3, 2, 3, 3, 2, 1, 2, 1, 1, 1, 3, 1, 3, 1, 2, 1, 1, 3, 3, 1, 1, 3, 1, 3, 3, 1, 3, 2, 3, 1,2) #which factor each item loads on
#' 
#' plots_out <- robust.item.plots(scoredSAT12, model=specific)
#' 
#' # View three discrimination parameters and intercept
#' plots_out$`Parameter Estimates for a1`
#' plots_out$`Parameter Estimates for a2`
#' plots_out$`Parameter Estimates for a3`
#' plots_out$`Parameter Estimates for d`
#' 
#' @export

robust.item.plots<-function(dat, ...){
  
  # Initial Model Estimation 
  mod <- mirt(dat, ...) 
  mml.ipars<-coef(mod, simplify = T)$items
  # Person fit residual calculation
  per.fit <- personfit(mod, method = 'ML')$Zh 
  # Weight  
  weight <- pnorm(per.fit)*nrow(dat)/ sum(pnorm(per.fit)) 
  
  # Robust model estimation 
  robust.mod <- mirt(dat, survey.weights=weight, ...)
  rmml.ipars<-coef(robust.mod, simplify = T)$items
  plot.titles<-paste("Parameter Estimates for", colnames(rmml.ipars))
  n.ipar<-ncol(rmml.ipars)
  out<-vector("list", n.ipar+2)
  names(out) <- c(plot.titles, "MML Item Parameter Estimates", "RMML Item Parameter Estimates")
  
  for(i in 1:n.ipar){
    df<-data.frame(x=mml.ipars[,i], y=rmml.ipars[,i])
    p<-ggplot(df, aes(x=x, y=y))+
      geom_point()+
      theme_classic()+
      geom_abline(slope = 1, intercept = 0)+
      ggtitle(plot.titles[i])+
      labs(x="MML Estimates", y="Robust Estimates")
    out[[i]] <- p
  }
  out[[n.ipar+1]]<-mml.ipars
  out[[n.ipar+2]]<-rmml.ipars
  return(out)
}

#' Histogram of RMMLE Weights 
#' 
#' The following function returns a histogram of weights computed for robust MML estimates of item parameters (Hong & Cheng, 2019). 
#' Weights are derived based on the person-fit statistic \eqn{l_z^{`*`}}, where lower values indicate greater misfit between the response pattern and the model (Snijders, 2001). 
#' Since \eqn{l_z^{`*`}} follows the asymptotic standard normal distribution, weights are the normalized \emph{p}-value of the residual under a one-sided alternative hypothesis, e.g.,
#' \eqn{w_i = \frac{\Phi(l^{`*`}_{zi})}{\sum_{i=1}^N \Phi(l^{`*`}_{zi})},}
#' where \eqn{\Phi(\cdot)} is the CDF of the standard normal distribution. 
#' 
#' @param dat A \eqn{N\times J} matrix of response data
#' @param freq Logical; if \code{TRUE}, the histogram graphic is a representation of frequencies; if \code{FALSE}, probability densities are plotted (so that the histogram has a total area of one). Defaults to \code{TRUE}.
#' @param ... Additional arguments to be passed to \code{mirt()}
#' 
#' @references Hong, M., & Cheng, Y. (2019). Robust maximum marginal likelihood (RMML) estimation for item response theory models. Behavior Research Methods, 51(2), 573–588. https://doi.org/10.3758/s13428-018-1150-4
#' @references Snijders, T. A. (2001). Asymptotic null distribution of person fit statistics with estimated person parameter. Psychometrika, 66(3), 331-342.
#' @references Sinharay, S. (2016). The choice of the ability estimate with asymptotically correct standardized person-fit statistics. British Journal of Mathematical and Statistical Psychology, 69(2), 175–193. https://doi.org/10.1111/bmsp.12067.
#' @details The IRT model is estimated with, first, MML and, second, RMML. Four histograms are returned. Two histograms display the distribution of person-fit residuals, \eqn{l^{`*`}_{zi}}, 
#'          across the sample. In the first of these histograms, \eqn{l^{`*`}_{zi}} is calculated from item parameters estimated with MML, while the second histogram displays the residuals calculated using robust MML item parameter estimates.
#'          As RMML reduces the bias in item parameter estimates due to aberrant responding, it is expected that person-fit residuals consequently are more accurate (e.g., see ).  
#'          Likewise, the latter two histograms display the distribution of person-specific weights, with the first graph utilizing MML and the second graph utilizing robust MML item parameter estimates.
#'          Under the null hypothesis, the p-values for a continuous test statistic are said to follow a standard normal distribution. 
#'          Given that the weights recommended for robust MML estimation are proportional to the p-values of \eqn{l^{`*`}_{zi}} according to the standard normal distribution, it is expected that, under model fit, the weights produce an approximately uniform distribution.
#'          Deviations from the uniform distribution, under appropriately large sample sizes, may provide insight into the behavior of the examinees, e.g., a large uptick near 0 may indicate prevalent aberrant responding.
#'          As with the residuals, these weights are expected to become more accurate after robust weighting and may reveal anomalous patterns more closely.
#'          Note that the distribution of weights may differ from \eqn{U(0,1)} due to the inverse weighting of the p-value by the sum of weights across the sample.
#' 
#' @return List with four histograms
#'  \describe{
#'    \item{\code{Person-Fit Indices (MMLE)}}{Histogram of person-fit residuals \eqn{l_z^{`*`}} across individuals based on the IRT model estimated with MML.}
#'    \item{\code{Person-Fit Indices (RMMLE)}}{Histogram of person-fit residuals \eqn{l_z^{`*`}} across individuals based on the IRT model estimated with robust MML.}
#'    \item{\code{Subject-Level Weights (MMLE)}}{(\eqn{N \times L}) Histogram of person-level weights based on the IRT model estimated with MML.}
#'    \item{\code{Subject-Level Weights (RMMLE)}}{(\eqn{N \times L}) Histogram of person-level weights based on the IRT model estimated with robust MML.}
#' }
#' @examples
#' 
#' # SAT12 Data: 3-Dimensional MIRT Model
#' SAT12[SAT12 == 8] <- NA #set 8 as a missing value
#' # Correct answer key
#' library(mirt)
#' key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5) 
#' scoredSAT12 <- key2binary(SAT12, key)
#' specific <- c(2, 3, 2, 3, 3, 2, 1, 2, 1, 1, 1, 3, 1, 3, 1, 2, 1, 1, 3, 3, 1, 1, 3, 1, 3, 3, 1, 3, 2, 3, 1,2) #which factor each item loads on
#' 
#' hist_out <- hist.weights(scoredSAT12, model=specific)
#' 
#' # Compare weights before and after robust estimation
#' hist_out$`Subject-Level Weights (MMLE)` 
#' hist_out$`Subject-Level Weights (RMMLE)`
#' 
#' @export

hist.weights<-function(dat, freq = TRUE, ...){
  
  # Initial Model Estimation 
  mod <- mirt(dat, ...) 
  # Person fit residual calculation
  per.fit <- personfit(mod, method = 'ML')$Zh 
  # Weight  
  weight <- pnorm(per.fit)*nrow(dat)/ sum(pnorm(per.fit)) 
  
  # Robust model estimation 
  robust.mod <- mirt(dat, survey.weights=weight, ...)
  # Person fit residual calculation
  per.fit2 <- personfit(robust.mod, method = 'ML')$Zh 
  # Weight  
  weight2 <- pnorm(per.fit2)*nrow(dat)/ sum(pnorm(per.fit2)) 
  
  return(list(`Person-Fit Indices (MMLE)` = hist(per.fit, main = "Histogram of Person-Fit Residuals (MMLE)", xlab = "Residuals", freq=freq)$hist, 
              `Person-Fit Indices (RMMLE)` = hist(per.fit2, main = "Histogram of Person-Fit Residuals (RMMLE)", xlab = "Residuals", freq=freq)$hist, 
              `Subject-Level Weights (MMLE)` = hist(weight, main = "Histogram of Subject-Level Weights (MMLE)", xlab = "Weights", freq=freq)$hist, 
              `Subject-Level Weights (RMMLE)` = hist(weight2, main = "Histogram of Subject-Level Weights (RMMLE)", xlab = "Weights", freq=freq)$hist))
  
}
                                    

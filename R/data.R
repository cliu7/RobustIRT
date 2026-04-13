#' SAT12 data
#' 
#' @description
#' Data obtained from the TESTFACT (Woods et al., 2003) manual, with 32 response pattern scored items for a grade 12 science assessment test (SAT) measuring topics of chemistry, biology, and physics. 
#' The scoring key for these data is [1, 4, 5, 2, 3, 1, 2, 1, 3, 1, 2, 4, 2, 1, 5, 3, 4, 4, 1, 4, 3, 3, 4, 1, 3, 5, 1, 3, 1, 5, 4, 5], respectively. 
#' However, careful analysis using the nominal response model suggests that the scoring key for item 32 may be incorrect, and should be changed from 5 to 3.
#'
#' @format A data frame with 600 examinees and 32 scored item responses. Each column corresponds to one item, coded as integers 1–5 representing the selected response option.
#' 
#' @author(s)
#' Phil Chalmers rphilip.chalmers@gmail.com
#'
#' @source \url{https://philchalmers.github.io/mirt/html/SAT12.html}
#' 
#' @references Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory Package for the R Environment. Journal of Statistical Software, 48(6), 1-29. \url{doi: 10.18637/jss.v048.i06}
#' @references McKinley, R. L., & Reckase, M. D. (1983). An extension of the two‑parameter logistic model to the multidimensional latent space. \emph{Psychometrika}, 48(3), 369–382.
#' @references Muraki, E., & Engelhard, G. (1985). Full-information item factor analysis: Applications of EAP scores. \emph{Applied Psychological Measurement}, 9(4), 417–430
#' @references Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., & Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics, and Full-information Item Factor Analysis [Computer software]. Lincolnwood, IL: Scientific Software International.
#' 
#' @usage data(SAT12)
#'
#' @examples
#' Run examples
#'
#' ## Not run: 
#'
#' itemstats(SAT12, use_ts = FALSE)
#'
#' #  score the data (missing scored as 0)
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
#' ## End(Not run)
"SAT12"
#'
#' BFI2
#'
#' Polytomous data from the Big Five Inventory-2
#'
#' Data obtained from the Big Five Inventory-2 (BFI-2; Soto & John, 2017), collected from adolescents ages 14 to 17 years old enrolled in a high school AP Statistics course (Ober et al., 2021). The inventory contains 60 items on a 5-point Likert scale, where exactly 12 items load onto each of the five personality factors in the BFI-2 (Extraversion, Agreeableness, Conscientiousness, Negative Emotionality, and Open-Mindedness).  Data for 838 subjects was collected, after removing four subjects with missing data. No missing data remains. Variables are labeled by the first letter of the personality factor (e.g., E = Extraversion, A = Agreeableness, etc.) and item number 1 through 12. Items marked with a "R" were negatively worded on the assessment and responses have been computationally reverse coded.
#'
#' @references Ober, T. M., Cheng, Y., Jacobucci, R., & Whitney, B. M. (2021, January). Examining the factor structure of the Big Five Inventory-2 personality domains with an adolescent sample. \emph{Psychological Assessment, 33}(1), 14–28. \url{doi:10.1037/pas0000962}
#' @references Soto, C. J., & John, O. P. (2017). The next Big Five Inventory (BFI-2): Developing and assessing a hierarchical model with 15 facets to enhance bandwidth, fidelity, and predictive power. \emph{Journal of Personality and Social Psychology, 113}, 117–143. \url{http://dx.doi.org/10.1037/pspp0000096}
#' @source \url{https://osf.io/awvnd}
"BFI2"

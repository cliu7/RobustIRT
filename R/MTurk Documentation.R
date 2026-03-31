# mydata <- readr::read_csv("data_amt.csv")
# usethis::use_data(mydata, overwrite = TRUE)

#' MTurk Comparison Data
#' 
#' @Description
#' A dataset containing 26 survey items collected from Amazon Mechanical Turk (AMT). Additional demographic variables are also included.
#' Respondents were recruited from MTurk and were paid $0.20 to participate. To qualify, workers were required to have completed at least 100 prior HITs  
#' with an acceptance rate of 97% or higher.
#'
#' @format A data frame with 1402 rows (participants) and 34 columns (survey items and follow-up information). 
#' \describe{
#'   \item{country}{The country the user's network connection was based.}
#'   \item{engnat}{Whether the respondent reported English as their native language (0 = no, 1 = yes).}
#'   \item{age}{Age in years.}
#'   \item{gender}{Gender indicator.)}
#'   \item{Q1}{I am tall.}
#'   \item{Q2}{I am short.}
#'   \item{Q3}{I have to stand on a stool to reach tall kitchen shelves.}
#'   \item{Q4}{I have to stand in the back in group photos to not cover up other people.}
#'   \item{Q5}{I hit my head on low ceilings.}
#'   \item{Q6}{I rarely meet people with more height than me.}
#'   \item{Q7}{I'm kind of a midget.}
#'   \item{Q8}{Airplane seats never have enough room for my long legs.}
#'   \item{Q9}{When I hug people, my head is underneath their chin.}
#'   \item{Q10}{I have gangly limbs.}
#'   \item{Q11}{I have been sent to the hospital by an electric shock.}
#'   \item{Q12}{I own a goat.}
#'   \item{Q13}{I know the 'happy birthday to you..' song.}
#'   \item{Q14}{I have been asked for money by beggars.}
#'   \item{Q15}{I prefer to play it safe and avoid danger.}
#'   \item{Q16}{I prefer variety to routine.}
#'   \item{Q17}{I rarely clean house.}
#'   \item{Q18}{I rarely complain.}
#'   \item{Q19}{I rarely overindulge.}
#'   \item{Q20}{I accept what others say.}
#'   \item{Q21}{I enjoy being part of a loud crowd.}
#'   \item{Q22}{I offend no one.}
#'   \item{Q23}{I see that nobody gets left out.}
#'   \item{Q24}{I try not to deceive others.}
#'   \item{Q25}{I try out new things.}
#'   \item{Q26}{I will push people around to get what I want.}
#'   \item{feet}{Height reported in feet (if the respondent used imperial units).}
#'   \item{inch}{Height reported in inches (if the respondent used imperial units).}
#'   \item{cm}{Height reported in centimeters (if the respondent used metric units).}
#'   \item{submittime}{The time (PST) the survey was submitted.}
#'
#' @details
#' This dataset is included to demonstrate robust estimation for various item response theory (IRT) models, particularly in contexts where aberrant responding may occur. 
#' The item data was collected using a 5-point Likert scale, with the following labels: 1 = Strongly disagree, 2 = Disagree, 3 = Neither agree nor disagree, 4 = Agree, and 5 = Strongly agree. A 0 indicates no response.
#'
#' @source Open Psychometrics (2019, December 29). A quality comparison of data collected on this website to data collected on Amazon Mechanical Turk. https://openpsychometrics.org/_rawdata/validity/
#' 
#' @usage data(MTurk)
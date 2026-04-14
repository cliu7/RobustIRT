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
"SAT12"
#' BFI2 data
#' 
#' Polytomous data from the Big Five Inventory-2
#' 
#' @description 
#' Data obtained from the Big Five Inventory-2 (BFI-2; Soto & John, 2017), collected from adolescents ages 14 to 17 years old enrolled in a high school AP Statistics course (Ober et al., 2021). 
#' The inventory contains 60 items on a 5-point Likert scale, where exactly 12 items load onto each of the five personality factors in the BFI-2 (Extraversion, Agreeableness, Conscientiousness, Negative Emotionality, and Open-Mindedness).  
#' Data for 838 subjects is contained, after removing four subjects with missing data. No missing data remains. 
#' Variables are labeled by the first letter of the personality factor (e.g., E = Extraversion, A = Agreeableness, etc.) and item number 1 through 12.
#' Items marked with an 'R' are negatively worded and were properly reverse coded for this dataset.
#' 
#' \itemize{
#'   \item{\code{id}}{Participant ID}.
#'   \item{\code{pilot}}{Data collection phase}.
#'   \item{\code{current_age_years_n}}{Current age in years}.
#'   \item{\code{biological_sex_c}}{Biological sex}.
#'   \item{\code{expected_education_c}}{Highest level of expected education}.
#'   \item{\code{parent_education_c}}{Highest level of education obtained by parents}.
#'   \item{\code{free_reduced_lunch_c}}{Indicator of eligibility for free or reduced lunch}.
#'   \item{\code{race_AmerIn1}}{American Indian or Alaska Native (1 = Yes).}
#'   \item{\code{race_Asian1}}{Asian (1 = Yes).}
#'   \item{\code{race_Black1}}{Black or African American (1 = Yes).}
#'   \item{\code{race_Hawa1}}{Native Hawaiian or Pacific Islander (1 = Yes).}
#'   \item{\code{race_White1}}{White (1 = Yes).}
#'   \item{\code{race_MexAmer1}}{Mexican American (1 = Yes).}
#'   \item{\code{race_Puerto1}}{Puerto Rican (1 = Yes).}
#'   \item{\code{race_HispOther1}}{Other Hispanic origin (1 = Yes).}
#'   \item{\code{race_HispAny1}}{Any Hispanic origin (1 = Yes).}
#'   \item{\code{race_Other1}}{Other race (1 = Yes).}
#'   \item{\code{race_PreferNot1}}{Prefer not to respond (1 = Yes).}
#'   \item{\code{race_Multi1}}{Multiracial (1 = Yes).}
#'   \item{\code{t1_bfi_E1}}{Is outgoing, sociable.}.
#'   \item{\code{t1_bfi_A1}}{Is compassionate.}.
#'   \item{\code{t1_bfi_C1R}}{Tends to be disorganized. [R]}.
#'   \item{\code{t1_bfi_N1R}}{Is relaxed, handles stress well. [R]}.
#'   \item{\code{t1_bfi_O1R}}{Has few artistic interests. [R]}.
#'   \item{\code{t1_bfi_E2}}{Has an assertive personality.}.
#'   \item{\code{t1_bfi_A2}}{Treats others with respect.}.
#'   \item{\code{t1_bfi_C2R}}{Tends to be lazy. [R]}.
#'   \item{\code{t1_bfi_N2R}}{Stays optimistic after experiencing a setback. [R]}.
#'   \item{\code{t1_bfi_O2}}{Is curious about many different things.}.
#'   \item{\code{t1_bfi_E3R}}{Rarely feels excited or eager. [R]}.
#'   \item{\code{t1_bfi_A3R}}{Tends to find fault with others. [R]}.
#'   \item{\code{t1_bfi_C3}}{Is dependable, steady.}.
#'   \item{\code{t1_bfi_N3}}{Is moody, has up and down mood swings.}.
#'   \item{\code{t1_bfi_O3}}{Is inventive, finds clever ways to do things.}.
#'   \item{\code{t1_bfi_E4R}}{Tends to be quiet. [R]}.
#'   \item{\code{t1_bfi_A4R}}{Feels little sympathy for others. [R]}.
#'   \item{\code{t1_bfi_C4}}{Is systematic, likes to keep things in order.}.
#'   \item{\code{t1_bfi_N4}}{Can be tense.}.
#'   \item{\code{t1_bfi_O4}}{Is fascinated by art, music, or literature.}.
#'   \item{\code{t1_bfi_E5}}{Is dominant, acts as a leader.}.
#'   \item{\code{t1_bfi_A5R}}{Starts arguments with others. [R]}.
#'   \item{\code{t1_bfi_C5R}}{Has difficulty getting started on tasks. [R]}.
#'   \item{\code{t1_bfi_N5R}}{Feels secure, comfortable with self. [R]}.
#'   \item{\code{t1_bfi_O5R}}{Avoids intellectual, philosophical discussions. [R]}.
#'   \item{\code{t1_bfi_E6R}}{Is less active than other people. [R]}.
#'   \item{\code{t1_bfi_A6}}{Has a forgiving nature.}.
#'   \item{\code{t1_bfi_C6R}}{Can be somewhat careless. [R]}.
#'   \item{\code{t1_bfi_N6R}}{Is emotionally stable, not easily upset. [R]}.
#'   \item{\code{t1_bfi_O6R}}{Has little creativity. [R]}.
#'   \item{\code{t1_bfi_E7R}}{Is sometimes shy, introverted. [R]}.
#'   \item{\code{t1_bfi_A7}}{Is helpful and unselfish with others.}.
#'   \item{\code{t1_bfi_C7}}{Keeps things neat and tidy.}.
#'   \item{\code{t1_bfi_N7}}{Worries a lot.}.
#'   \item{\code{t1_bfi_O7}}{Values art and beauty.}.
#'   \item{\code{t1_bfi_E8R}}{Finds it hard to influence people. [R]}.
#'   \item{\code{t1_bfi_A8R}}{Is sometimes rude to others. [R]}.
#'   \item{\code{t1_bfi_C8}}{Is efficient, gets things done.}.
#'   \item{\code{t1_bfi_N8}}{Often feels sad.}.
#'   \item{\code{t1_bfi_O8}}{Is complex, a deep thinker.}.
#'   \item{\code{t1_bfi_E9}}{Is full of energy.}.
#'   \item{\code{t1_bfi_A9R}}{Is suspicious of others' intentions. [R]}.
#'   \item{\code{t1_bfi_C9}}{Is reliable, can always be counted on.}.
#'   \item{\code{t1_bfi_N9R}}{Keeps their emotions under control. [R]}.
#'   \item{\code{t1_bfi_O9R}}{Has difficulty imagining things. [R]}.
#'   \item{\code{t1_bfi_E10}}{Is talkative.}.
#'   \item{\code{t1_bfi_A10R}}{Can be cold and uncaring. [R]}.
#'   \item{\code{t1_bfi_C10R}}{Leaves a mess, doesn't clean up. [R]}.
#'   \item{\code{t1_bfi_N10R}}{Rarely feels anxious or afraid. [R]}.
#'   \item{\code{t1_bfi_O10R}}{Thinks poetry and plays are boring. [R]}.
#'   \item{\code{t1_bfi_E11R}}{Prefers to have others take charge. [R]}.
#'   \item{\code{t1_bfi_A11}}{Is polite, courteous to others.}.
#'   \item{\code{t1_bfi_C11}}{Is persistent, works until the task is finished.}.
#'   \item{\code{t1_bfi_N11}}{Tends to feel depressed, blue.}.
#'   \item{\code{t1_bfi_O11R}}{Has little interest in abstract ideas. [R]}.
#'   \item{\code{t1_bfi_E12}}{Shows a lot of enthusiasm.}.
#'   \item{\code{t1_bfi_A12}}{Assumes the best about people.}.
#'   \item{\code{t1_bfi_C12R}}{Sometimes behaves irresponsibly. [R]}.
#'   \item{\code{t1_bfi_N12}}{Is temperamental, gets emotional easily.}.
#'   \item{\code{t1_bfi_O12}}{Is original, comes up with new ideas.}.
#' }
#'  
#' @details The BFI2 is a 60‑item personality assessment measuring the Big Five domains: Agreeableness, Conscientiousness, Extraversion, Negative Emotionality, and Open‑Mindedness. 
#' Each domain contains 12 items, organized into 15 facets (e.g., Compassion, Organization, Anxiety). Items are scored on a 1–5 Likert scale and include both positively and negatively keyed statements.
#'
#' @source \url{https://osf.io/awvnd}
#' 
#' @references Ober, T. M., Cheng, Y., Jacobucci, R., & Whitney, B. M. (2021, January). Examining the factor structure of the Big Five Inventory-2 personality domains with an adolescent sample. \emph{Psychological Assessment, 33}(1), 14–28. \url{doi:10.1037/pas0000962}
#' @references Soto, C. J., & John, O. P. (2017). The next Big Five Inventory (BFI-2): Developing and assessing a hierarchical model with 15 facets to enhance bandwidth, fidelity, and predictive power. \emph{Journal of Personality and Social Psychology, 113}, 117–143. \url{http://dx.doi.org/10.1037/pspp0000096}
#' 
#' @usage data(BFI2)
"BFI2"
#' MTurk Height Data
#' 
#' @description
#' A dataset containing 26 survey items collected from Amazon Mechanical Turk (AMT). Additional demographic variables are also included.
#' Respondents were recruited from MTurk and were paid $0.20 to participate. To qualify, workers were required to have completed at least 100 prior human intelligence tasks (HITs)  
#' with an acceptance rate of 97% or higher.
#'
#' @format 
#' A data frame with 1402 rows (participants) and 34 columns (survey items and follow-up information). 
#' \itemize{
#'   \item{\code{country}}{The country the user's network connection was based.}
#'   \item{\code{engnat}}{Whether the respondent reported English as their native language (0 = no, 1 = yes).}
#'   \item{\code{age}}{Age in years.}
#'   \item{\code{gender}}{Gender indicator.)}
#'   \item{\code{Q1}}{I am tall.}
#'   \item{\code{Q2}}{I am short.}
#'   \item{\code{Q3}}{I have to stand on a stool to reach tall kitchen shelves.}
#'   \item{\code{Q4}}{I have to stand in the back in group photos to not cover up other people.}
#'   \item{\code{Q5}}{I hit my head on low ceilings.}
#'   \item{\code{Q6}}{I rarely meet people with more height than me.}
#'   \item{\code{Q7}}{I'm kind of a midget.}
#'   \item{\code{Q8}}{Airplane seats never have enough room for my long legs.}
#'   \item{\code{Q9}}{When I hug people, my head is underneath their chin.}
#'   \item{\code{Q10}}{I have gangly limbs.}
#'   \item{\code{Q11}}{I have been sent to the hospital by an electric shock.}
#'   \item{\code{Q12}}{I own a goat.}
#'   \item{\code{Q13}}{I know the 'happy birthday to you..' song.}
#'   \item{\code{Q14}}{I have been asked for money by beggars.}
#'   \item{\code{Q15}}{I prefer to play it safe and avoid danger.}
#'   \item{\code{Q16}}{I prefer variety to routine.}
#'   \item{\code{Q17}}{I rarely clean house.}
#'   \item{\code{Q18}}{I rarely complain.}
#'   \item{\code{Q19}}{I rarely overindulge.}
#'   \item{\code{Q20}}{I accept what others say.}
#'   \item{\code{Q21}}{I enjoy being part of a loud crowd.}
#'   \item{\code{Q22}}{I offend no one.}
#'   \item{\code{Q23}}{I see that nobody gets left out.}
#'   \item{\code{Q24}}{I try not to deceive others.}
#'   \item{\code{Q25}}{I try out new things.}
#'   \item{\code{Q26}}{I will push people around to get what I want.}
#'   \item{\code{feet}}{Height reported in feet (if the respondent used imperial units).}
#'   \item{\code{inch}}{Height reported in inches (if the respondent used imperial units).}
#'   \item{\code{cm}}{Height reported in centimeters (if the respondent used metric units).}
#'   \item{\code{submittime}}{The time (PST) the survey was submitted.}
#'}
#' @details 
#' The item data followed a 5-point Likert scale, with the following labels: 1 = Strongly disagree, 2 = Disagree, 3 = Neither agree nor disagree, 4 = Agree, and 5 = Strongly agree. A 0 indicates no response.
#' The first four items include demographic information. Q1 through Q10 measure height, while Q13 through Q26 measure general personality characteristics. 
#' The last four items ask users to report their hight and contain a record of submission time. 
#' Q11 and Q12 are bogus items, designed to detect inattentive responders, or those who endorse the item although it is unlikely to be true.
#' Due to factors such as disengagement and low motivation that are common on low stakes assessment, MTurk data may be especially susceptible to aberrant response behavior. 
#' When fitting an item response theory model, robust estimation may be applied to reduce bias in the latent trait estimate for respondents exhibiting aberrant responses.
#' However, when a subject's response pattern contains aberrant responses for the majority of or all items (e.g., "bot" responding), robust estimation will not be effective in mitigating the aberrances and unable to reduce the bias in the latent trait.
#'
#' @source Open Psychometrics (2019, December 29). A quality comparison of data collected on this website to data collected on Amazon Mechanical Turk. https://openpsychometrics.org/_rawdata/validity/
#' 
#' @usage data(MTurkHeight)
"MTurkHeight"

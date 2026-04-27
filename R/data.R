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



#'  ################ Simulating Datasets For Aberrant Responding ###################

set.seed(123) # set the seed to ensure reproducibility

N <- 1000 # number of individuals
J <- 30 # number of items
theta <- rnorm(N, 0, 1) # generate theta values from normal distribution

a <- runif(J, 1.0, 2.0) # generate discrimination values between 1.0 and 2.0
b <- rnorm(J) # generate difficulty values


#' Suboptimal Responding with the 2PL Model
#' 
#' @description 
#' Simulates dichotomous item responses under a 2PL model in which a subset of
#' individuals respond suboptimally by having their latent trait values reduced.
#' @details
#' Extracts item parameters from uniform (discrimination) and normal
#' (difficulty) distributions. Forty percent of individuals are randomly selected
#' to respond suboptimally. Their latent ability values are shifted downward by
#' one standard unit before computing 2PL response probabilities. Responses are
#' generated using data from the GDINA package.
#' @param N Number of individuals.
#' @param J Number of items.
#' @param theta Vector of latent ability traits.
#' @param a Item discrimination parameters.
#' @param b Item difficulty parameters.
#' @param subopt_ids Individuals responding suboptimally.
#' @param prob_2pl Function computing 2PL probabilities.
#' @param dat Simulated response matrix.
#' @references C. Schuster and K.-H. Yuan. Robust Estimation of Latent Ability in Item Response Models. Journal of Educational and Behavioral Statistics, 36(6):720–735, Dec. 2011. ISSN 1076-9986, 1935-1054. doi: 10.3102/1076998610396890. URL http://journals.sagepub.com/doi/10.3102/1076998610396890.

# 1. Identify 40% of examinees to respond suboptimally
n_subopt <- floor(0.40 * N)
subopt_ids <- sample(1:N, n_subopt)

# 2. Create updated theta vector
theta_new <- theta
theta_new[subopt_ids] <- theta_new[subopt_ids] - 1

# 2PL probabilities function
prob_2pl <- function(theta, a, b) {
  exp(1.7 * a * (theta - b)) / (1 + exp(1.7 * a * (theta - b)))
}

# Create a matrix with 2PL generated probabilities
prob_matrix <- sapply(1:J, function(j) {
  prob_2pl(theta_new, a[j], b[j])
})

# Generate the data using the dat.gen function and the probabilities above
dat <- dat.gen(prob_matrix) 

save(simdat=dat, aberrant_ids=subopt_ids, ipars=cbind(a,b), file = "simSuboptimal2PL.RData")


#' Back Random Responding (BRR) with the 2PL Model
#' 
#' @description
#' Simulates dichotomous responses where a subset of items at the end of the test
#' are replaced with random responses for a subset of individuals.
#' @details
#' Responses are first generated under a standard 2PL model. The last 40% of
#' items are designated as BRR items. A prevalence rate determines which
#' individuals exhibit BRR. For these individuals, responses to BRR items are
#' overwritten with probability 0.20.
#' @param prev_rate Proportion of individuals showing BRR.
#' @param last_items Items overwritten with random responses.
#' @param brr_ids Individuals exhibiting BRR.
#' @param brr_prob The response probability applied to BRR items for affected individuals.
#' @references Clark, M. E., Gironda, R. J., & Young, R. W. (2003). Detection of back random responding: Effectiveness of MMPI-2 and Personality Assessment Inventory validity indices. Psychological Assessment, 15(2), 223–234. https://doi.org/10.1037/1040-3590.15.2.223
#' @references Yu, X., & Cheng, Y. (2019). A change-point analysis procedure based on weighted residuals to detect back random responding. Psychological Methods, 24(5), 658–674. https://doi.org/10.1037/met0000212
#' 
theta <- rnorm(N, 0, 1)
a <- runif(J, 1.0, 2.0) # generate discrimination values between 1.0 and 2.0
b <- rnorm(J) # generate difficulty values

# Create a matrix with 2PL generated probabilities
prob_matrix <- sapply(1:J, function(j) {
  prob_2pl(theta, a[j], b[j])
})

# Generate the data using the dat.gen function and the probabilities above
dat <- dat.gen(prob_matrix) 

# Identify the last 40% of items
last_items <- tail(1:J, floor(0.40 * J))

# Probability of BRR
brr_prob <- 0.20

# Overwrite those items with BRR responses
for (j in last_items) {
  dat[, j] <- rbinom(N, 1, brr_prob)
}

# subset of people, prevalence rate
prev_rate <- 0.20
brr_ids <- sample(1:N, size = floor(prev_rate * N))

simdat = dat
aberrant_ids = brr_ids
ipars = cbind(a, b)

save(simdat=dat, aberrant_ids=brr_ids, ipars=cbind(a, b), file = "simBRR2PL.RData")

#' Cheating with the 2PL Model
#' 
#' @description
#' Simulates dichotomous responses where low‑ability individuals are given perfect
#' scores on the most difficult items.
#' @details
#' After generating 2PL responses, the lowest 20% of individuals are identified. 
#' The most difficult 20% of items are also identified. For these individual-item combinations, 
#' responses are overwritten with correct answers to reflect cheating behavior.
#' @param low_theta_ids Individuals with lowest ability.
#' @param dif_items Most difficult items.
#' @param dat Response matrix with cheating applied.
#' @references Wang, C., Xu, G., Shang, Z., & Kuncel, N. (2018). Detecting Aberrant Behavior and Item Preknowledge: A Comparison of Mixture Modeling Method and Residual Method. Journal of Educational and Behavioral Statistics, 43(4), 469–501. https://doi.org/10.3102/1076998618767123
#' 
theta <- rnorm(N, 0, 1)
a <- runif(J, 1.0, 2.0) # generate discrimination values between 1.0 and 2.0
b <- rnorm(J) # generate difficulty values

# Create a matrix with 2PL generated probabilities
prob_matrix <- sapply(1:J, function(j) {
  prob_2pl(theta, a[j], b[j])
})

# Generate the data using the dat.gen function and the probabilities above
dat <- dat.gen(prob_matrix) 

# Identify the 20% lowest thetas (lowest ability)
n_low <- floor(0.20 * N)
low_theta_ids <- order(theta)[1:n_low]

# Identify the 20% most difficult items
n_dif <- floor(0.20 * J)
dif_items <- order(b, decreasing = TRUE)[1:n_dif]

# Overwrite those items with 
for (j in dif_items) {
  dat[low_theta_ids, j] <- rbinom(n_low, 1, 1)
}

save(simdat=dat, aberrant_ids=low_theta_ids, ipars=cbind(a,b), file = "simCheating2PL.RData")

#' Warm-up Responding with the 2PL Model
#' @description
#' Simulates individuals who perform poorly on early items due to warm-up effects.
#' @details
#' Responses are generated under a 2PL model. A random 20% of individuals are
#' selected to exhibit warm‑up behavior. For these individuals, responses to the
#' first 30% of items are overwritten with incorrect responses.
#' @param warmup_ids Individuals showing warm‑up behavior.
#' @param warmup_items Early items overwritten with zeros.
#' @references Meijer, R. R. (2002). Outlier Detection in High-Stakes Certification Testing. Journal of Educational Measurement, 39(3), 219–233. https://doi.org/10.1111/j.1745-3984.2002.tb01175.x
#' 
theta <- rnorm(N, 0, 1)
a <- runif(J, 1.0, 2.0) # generate discrimination values between 1.0 and 2.0
b <- rnorm(J) # generate difficulty values

# Create a matrix with 2PL generated probabilities
prob_matrix <- sapply(1:J, function(j) {
  prob_2pl(theta, a[j], b[j])
})

# Generate the data using the dat.gen function and the probabilities above
dat <- dat.gen(prob_matrix) 

# Identify 30% of participants to "warm-up"
prev_rate <- 0.20
warmup_ids <- sample(1:N, size = floor(prev_rate * N)) #MAKE RANDOM

# Identify the first 30% of items
n_items_warmup <- floor(0.30 * J)
warmup_items <- 1:n_items_warmup

# Overwrite those responses as incorrect (0)
for (j in warmup_items) {
  dat[warmup_ids, j] <- 0
}

save(simdat=dat, aberrant_ids=warmup_ids, ipars=cbind(a, b), file = "simWarmUp2PL.RData")

#' Improper Reverse Coding with the Graded Response Model (GRM)
#' 
#' @description
#' Generates polytomous responses under a graded response model (GRM) and applies
#' improper reverse coding to a subset of individuals and items.
#' @details
#' Item thresholds are generated and then sorted within each item to maintain the
#' required ordering. GRM category probabilities are obtained. A random 20% of individuals 
#' and 30% of items are selected, and responses are recoded using \eqn{X' = K + 1 - X}, 
#' yielding improperly reversed categories.
#' @param K Number of response categories.
#' @param irc_ids Individuals with IRC behavior.
#' @param irc_items Items subjected to reverse coding.
#' @references Hughes, G. D. (2009). The Impact of Incorrect Responses to Reverse-Coded Survey Items. Research in the Schools, 14
#' 
K <- 5  # maximum Likert category

theta <- rnorm(N, 0, 1)
a <- runif(J, 1.0, 2.0) # generate discrimination values between 1.0 and 2.0
b_raw <- matrix(rnorm(J + (K-1)), nrow = J)
b <- t(apply(b_raw, 1, sort)) # Sort thresholds within each item (row-wise)

# Create a matrix with GRM generated probabilities
prob_matrix <- item.prob(theta, model = "GRM", ipars = cbind(a, b), D=1.7)

# Generate the data using the dat.gen function and the probabilities above
dat <- dat.gen(prob_matrix$P, anchor=1, polytomous = TRUE) 

# Identify 30% of items
n_irc <- floor(0.30 * J)
irc_items <- sample(1:J, n_irc)

# Do this for percentage of people
prev_rate <- 0.20
irc_ids <- sample(1:N, size = floor(prev_rate * N))

# Apply improperly reverse coded responding: observed = K + X - 1
for (j in irc_items) {
  dat[irc_ids, j] <- K + 1 - dat[irc_ids, j]
}

save(simdat=dat, aberrant_ids=irc_ids, ipars=cbind(a, b), file = "simIRCGRM.RData")

#' Back Random Responding with the Multidimensional Graded Response Model (MGRM)
#' 
#' @description
#' Generates multidimensional polytomous responses under an MGRM and applies
#' random responding to a subset of individuals and items.
#' @details
#' A simple‑structure loading matrix assigns each item to one dimension. After
#' generating MGRM responses, 40% of items and 20% of individuals are selected.
#' For these individual–item combinations, responses are replaced with uniformly random
#' category selections.
#' @param L Number of latent dimensions.
#' @param brrmgrm_items Items with BRR.
#' @param brrmgrm_ids Individuals exhibiting BRR.
#' @param brrmgrm_prob The probability assigned to each response category under BRR.
#' @references Clark, M. E., Gironda, R. J., & Young, R. W. (2003). Detection of back random responding: Effectiveness of MMPI-2 and Personality Assessment Inventory validity indices. Psychological Assessment, 15(2), 223–234. https://doi.org/10.1037/1040-3590.15.2.223
#' @references Yu, X., & Cheng, Y. (2019). A change-point analysis procedure based on weighted residuals to detect back random responding. Psychological Methods, 24(5), 658–674. https://doi.org/10.1037/met0000212
#' 
N <- 2000 # number of individuals
J <- 50 # number of items
K <- 5  # maximum Likert category
L <- 3 # number of dimensions
theta <- matrix(rnorm(N * L), ncol = L)

# Each item loads on exactly one dimension (simple structure)
a <- matrix(0, nrow = J, ncol = L)

# Assign each dimension an equal number of items
items_per_dim <- ceiling(J / L)
dim_assign <- rep(1:L, each = items_per_dim)[1:J]

for (j in 1:J) {a[j, dim_assign[j]] <- runif(1, 1.0, 2.0)}
b_raw <- matrix(rnorm(J * (K - 1)), nrow = J)
b <- t(apply(b_raw, 1, sort))

ipars <- cbind(a, b)
prob_matrix <- item.prob(theta, model = "MGRM", ipars = ipars, D = 1.7)
dat <- dat.gen(prob_matrix$P, anchor = 1, polytomous = TRUE)

#Apply back random responding
n_brrmgrm_items <- floor(0.40 * J)
brrmgrm_items <- sample(1:J, n_brrmgrm_items)

# Choose 20% of people to respond randomly
prev_rate <- 0.20
brrmgrm_ids <- sample(1:N, size = floor(prev_rate * N))

# Random responding probability
brrmgrm_prob <- 1 / K   # uniform over categories

for (j in brrmgrm_items) {
  dat[brrmgrm_ids, j] <- sample(1:K, size = length(brrmgrm_ids), replace = TRUE)
}

save(simdat=dat, aberrant_ids=brrmgrm_ids, ipars, file = "simBRRMGRM.RData")

#' Cheating with the Multidimensional Item Response Theory Model (MIRT)
#' 
#' @description
#' Generates polytomous responses under a multidimensional graded response model
#' and applies cheating behavior to low‑ability individuals on the hardest items.
#' @details
#' Individuals falling in the lowest 20% on the first dimension are identified.
#' Item difficulty is computed as the mean of the category thresholds, and the
#' hardest 20% of items are selected. For these examinees and items, responses
#' are replaced with the highest category.
#' @param low_theta_mirt_ids Individuals who are cheating.
#' @param dif_items Hardest items.
#' @references Cui, Y., & Li, J. (2015). Evaluating Person Fit for Cognitive Diagnostic Assessment. Applied Psychological Measurement, 39(3), 223–238. https://doi.org/10.1177/0146621614557272
#' @references Meijer, R. R. (1996). Person-Fit Research: An Introduction. Applied Measurement in Education, 9(1), 3–8. https://doi.org/10.1207/s15324818ame0901

# J = 50 items, N = 2000, L = 2, not simple structure

N <- 2000 # number of individuals
J <- 50 # number of items
K <- 5  # maximum Likert category
L <- 2 # number of dimensions
theta <- matrix(rnorm(N * L), ncol = L)

a <- matrix(runif(J * L, 1.0, 2.0), nrow = J, ncol = L)
b_raw <- matrix(rnorm(J * (K - 1)), nrow = J)
b <- t(apply(b_raw, 1, sort))
ipars <- cbind(a, b)

prob_matrix <- item.prob(theta, model = "MGRM", ipars = ipars, D = 1.7)
dat <- dat.gen(prob_matrix$P, anchor = 1, polytomous = TRUE)

# Apply cheating

# Identify lowest 20% of examinees by ability on dimension 1
n_low <- floor(0.20 * N)
low_theta_mirt_ids <- order(theta[, 1])[1:n_low]

# Identify hardest 20% of items (highest average threshold)
item_difficulty <- rowMeans(b)
n_dif <- floor(0.20 * J)
dif_items <- order(item_difficulty, decreasing = TRUE)[1:n_dif]

# Overwrite those responses with perfect scores (cheating)
for (j in dif_items) {
  dat[low_theta_mirt_ids, j] <- K   # highest category
}

save(simdat=dat, aberrant_ids=low_theta_mirt_ids, ipars=cbind(a, b), file = "simCheatingMirt.RData")


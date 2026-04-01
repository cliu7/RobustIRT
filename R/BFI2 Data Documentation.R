#' BFI2 data
#' 
#' Polytomous data from the Big Five Inventory-2
#' 
#' @description Data obtained from the Big Five Inventory-2 (BFI-2; Soto & John, 2017), collected from adolescents ages 14 to 17 years old enrolled in a high school AP Statistics course (Ober et al., 2021). 
#' The inventory contains 60 items on a 5-point Likert scale, where exactly 12 items load onto each of the five personality factors in the BFI-2 (Extraversion, Agreeableness, Conscientiousness, Negative Emotionality, and Open-Mindedness).  
#' Data for 838 subjects was collected, after removing four subjects with missing data. No missing data remains. 
#' Variables are labeled by the first letter of the personality factor (e.g., E = Extraversion, A = Agreeableness, etc.) and item number 1 through 12. 
#' 
#' @format A data frame with 600 examinees and 32 scored item responses. Each column corresponds to one item, coded as integers 1–5 representing the selected response option.
#' \describe{
#'   \item{id}{Participant ID}.
#'   \item{pilot}{Data collection phase}.
#'   \item{current_age_years_n}{Current age in years}.
#'   \item{biological_sex_c}{Biological sex}.
#'   \item{expected_education_c}{Highest level of expected education}.
#'   \item{parent_education_c}{Highest level of education obtained by parents}.
#'   \item{free_reduced_lunch_c}{Indicator of eligibility for free or reduced lunch}.
#'   \item{race_AmerIn1}{American Indian or Alaska Native (1 = Yes).}
#'   \item{race_Asian1}{Asian (1 = Yes).}
#'   \item{race_Black1}{Black or African American (1 = Yes).}
#'   \item{race_Hawa1}{Native Hawaiian or Pacific Islander (1 = Yes).}
#'   \item{race_White1}{White (1 = Yes).}
#'   \item{race_MexAmer1}{Mexican American (1 = Yes).}
#'   \item{race_Puerto1}{Puerto Rican (1 = Yes).}
#'   \item{race_HispOther1}{Other Hispanic origin (1 = Yes).}
#'   \item{race_HispAny1}{Any Hispanic origin (1 = Yes).}
#'   \item{race_Other1}{Other race (1 = Yes).}
#'   \item{race_PreferNot1}{Prefer not to respond (1 = Yes).}
#'   \item{race_Multi1}{Multiracial (1 = Yes).}
#'   \item{t1_bfi_E1}{Is outgoing, sociable.}.
#'   \item{t1_bfi_A1}{Is compassionate.}.
#'   \item{t1_bfi_C1R}{Tends to be disorganized. [R]}.
#'   \item{t1_bfi_N1R}{Is relaxed, handles stress well. [R]}.
#'   \item{t1_bfi_O1R}{Has few artistic interests. [R]}.
#'   \item{t1_bfi_E2}{Has an assertive personality.}.
#'   \item{t1_bfi_A2}{Treats others with respect.}.
#'   \item{t1_bfi_C2R}{Tends to be lazy. [R]}.
#'   \item{t1_bfi_N2R}{Stays optimistic after experiencing a setback. [R]}.
#'   \item{t1_bfi_O2}{Is curious about many different things.}.
#'   \item{t1_bfi_E3R}{Rarely feels excited or eager. [R]}.
#'   \item{t1_bfi_A3R}{Tends to find fault with others. [R]}.
#'   \item{t1_bfi_C3}{Is dependable, steady.}.
#'   \item{t1_bfi_N3}{Is moody, has up and down mood swings.}.
#'   \item{t1_bfi_O3}{Is inventive, finds clever ways to do things.}.
#'   \item{t1_bfi_E4R}{Tends to be quiet. [R]}.
#'   \item{t1_bfi_A4R}{Feels little sympathy for others. [R]}.
#'   \item{t1_bfi_C4}{Is systematic, likes to keep things in order.}.
#'   \item{t1_bfi_N4}{Can be tense.}.
#'   \item{t1_bfi_O4}{Is fascinated by art, music, or literature.}.
#'   \item{t1_bfi_E5}{Is dominant, acts as a leader.}.
#'   \item{t1_bfi_A5R}{Starts arguments with others. [R]}.
#'   \item{t1_bfi_C5R}{Has difficulty getting started on tasks. [R]}.
#'   \item{t1_bfi_N5R}{Feels secure, comfortable with self. [R]}.
#'   \item{t1_bfi_O5R}{Avoids intellectual, philosophical discussions. [R]}.
#'   \item{t1_bfi_E6R}{Is less active than other people. [R]}.
#'   \item{t1_bfi_A6}{Has a forgiving nature.}.
#'   \item{t1_bfi_C6R}{Can be somewhat careless. [R]}.
#'   \item{t1_bfi_N6R}{Is emotionally stable, not easily upset. [R]}.
#'   \item{t1_bfi_O6R}{Has little creativity. [R]}.
#'   \item{t1_bfi_E7R}{Is sometimes shy, introverted. [R]}.
#'   \item{t1_bfi_A7}{Is helpful and unselfish with others.}.
#'   \item{t1_bfi_C7}{Keeps things neat and tidy.}.
#'   \item{t1_bfi_N7}{Worries a lot.}.
#'   \item{t1_bfi_O7}{Values art and beauty.}.
#'   \item{t1_bfi_E8R}{Finds it hard to influence people. [R]}.
#'   \item{t1_bfi_A8R}{Is sometimes rude to others. [R]}.
#'   \item{t1_bfi_C8}{Is efficient, gets things done.}.
#'   \item{t1_bfi_N8}{Often feels sad.}.
#'   \item{t1_bfi_O8}{Is complex, a deep thinker.}.
#'   \item{t1_bfi_E9}{Is full of energy.}.
#'   \item{t1_bfi_A9R}{Is suspicious of others' intentions. [R]}.
#'   \item{t1_bfi_C9}{Is reliable, can always be counted on.}.
#'   \item{t1_bfi_N9R}{Keeps their emotions under control. [R]}.
#'   \item{t1_bfi_O9R}{Has difficulty imagining things. [R]}.
#'   \item{t1_bfi_E10}{Is talkative.}.
#'   \item{t1_bfi_A10R}{Can be cold and uncaring. [R]}.
#'   \item{t1_bfi_C10R}{Leaves a mess, doesn't clean up. [R]}.
#'   \item{t1_bfi_N10R}{Rarely feels anxious or afraid. [R]}.
#'   \item{t1_bfi_O10R}{Thinks poetry and plays are boring. [R]}.
#'   \item{t1_bfi_E11R}{Prefers to have others take charge. [R]}.
#'   \item{t1_bfi_A11}{Is polite, courteous to others.}.
#'   \item{t1_bfi_C11}{Is persistent, works until the task is finished.}.
#'   \item{t1_bfi_N11}{Tends to feel depressed, blue.}.
#'   \item{t1_bfi_O11R}{Has little interest in abstract ideas. [R]}.
#'   \item{t1_bfi_E12}{Shows a lot of enthusiasm.}.
#'   \item{t1_bfi_A12}{Assumes the best about people.}.
#'   \item{t1_bfi_C12R}{Sometimes behaves irresponsibly. [R]}.
#'   \item{t1_bfi_N12}{Is temperamental, gets emotional easily.}.
#'   \item{t1_bfi_O12}{Is original, comes up with new ideas.}.
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
#'
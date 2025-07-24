### doc-data.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 13 2025 (11:30) 
## Version: 
## Last-Updated: jul 22 2025 (09:44) 
##           By: Brice Ozenne
##     Update #: 62
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * CHARM
#' @title RCT In Chronic Heart Failure Assessing an Inhibitor of the Renin-Angiotensin System. 
#' @name CHARM
#' @rdname data-CHARM
#'
#' @description Simulated data inspired from the CHARM-Preserved Trial comparing cardiovascular death or admission to hospital for chronic heart failure (CHF) among patients with CHF treated with candesartan or placebo.
#'
#' \itemize{
#' \item \code{id}: study participant.
#' \item \code{treatment}: treatment arm: candesartan (T) or placebo (C).
#' \item \code{Mortality}: time from inclusion (say diagnosis) to death or loss to follow-up.
#' \item \code{statusMortality}: event triggering the end of follow-up: death (1), drop-out (0).
#' \item \code{Hospitalization}: time from inclusion (say diagnosis) to hospitalization or end of follow-up.
#' \item \code{statusHospitalization}: hospitalization (1) or end of follow-up (0).
#' \item \code{Composite}: time to hospitalization, death, or loss to follow-up, whichever comes first.
#' \item \code{statusComposite}: hospitalization or death (1), drop-out (0).
#' }
#' 
#' @author Johan Verbeeck
#' @docType data
#' @usage data(CHARM, package = "BuyseTest")
#' @references Yusuf Salim, et al. "Effects of candesartan in patients with chronic heart failure and preserved left-ventricular ejection fraction: the CHARM-Preserved Trial". The Lancet (2003) 9386(362):777-781.
#' @keywords datasets
"CHARM"

## * EB
#' @title Rare disease trial
#' @name EB
#' @rdname data-EB
#'
#' @description 16 pediatric subjects suffering from a aare skin disease (epidermolysis bullosa simplex) treated with placebo or diacerin cream in a longitudinal cross-over trial (14 paired).
#'
#' \itemize{
#' \item \code{Id}: study participant.
#' \item \code{Time}: end of the 4 weeks treatment period (t4) or 12 weeks after treatment ended (t12)
#' \item \code{treatment}: treatment arm: diacerin cream (V) or placebo (P).
#' \item \code{StdDiffCount}: blister uncertainty
#' \item \code{Bin}: indicator of 40% reduction from baseline in the number of blisters within the treated areas.
#' \item \code{DiffQoL}: Change in quality of life since baseline.
#' \item \code{period}: same as time but coded as numeric: 1 for t4 and 2 for t12
#' }
#' 
#' @docType data
#' @usage data(CHARM, package = "BuyseTest")
#' @references Wally et al. "Diacerein orphan drug development for epidermolysis bullosa simplex: A phase 2/3 randomized, placebo-controlled, double-blind clinical trial". Journal of American Academy of Dermatology (2018) 78(5):892-901. https://doi.org/10.1016/j.jaad.2018.01.019.
#' @keywords datasets
"EB"

## * Profil
#' @title N-of-1 trials with On-demand Sildenafil as a Treatment for Raynaud Phenomenon
#' @name profil
#' @rdname data-profil
#'
#' @description Data from the Profil trial, a multiple crossover study where 38 patients had repeat blocks.
#' Each block lasted 3 weeks and corresponded to a sequence of all three treatments (placebo, sildenafil 40mg, or sildenafil 80mg), each treatment being allocated during 1 week.
#' The principal outcome were Raynaud Condition Score (RCS) and the frequency and cumulative duration of attacks over 24 hours.
#'
#' \itemize{
#' \item \code{id}: study participant.
#' \item \code{age}: age in years, between 18 and 74 years.
#' \item \code{male}: indicator of whether the participant was a male. 
#' \item \code{period}: number of the block for a given patient. From 1 to 5.
#' \item \code{time}: number of days since inclusion.
#' \item \code{treatment}: treatment allocated during the period (placebo, lowDose, highDose)
#' \item \code{rcs}: RCS score, ranging from 0 to 10.
#' \item \code{number}: daily number of attack, ranging from 0 to 9.
#' \item \code{duration}: cumulative daily duration of attacks in minutes, between 0 and 755 minutes.
#' }
#'
#' @author Brice Ozenne took a subset of the columns from the original dataset (available at \url{https://datadryad.org/dataset/doi:10.5061/dryad.c670tq2})
#' and translated their names into english.
#' @docType data
#' @usage data(profil, package = "BuyseTest")
#' @references Matthieu Roustit, Joris Giai, Olivier Gaget, et al. On-Demand Sildenafil as a Treatment for Raynaud Phenomenon: A Series of n-of-1 Trials. Ann Intern Med.2018;169:694-703. [Epub 30 October 2018]. doi:10.7326/M18-0517
#' @keywords datasets
"profil"

## * PRODIGE
#' @title RCT In Metastatic Pancreatic Cancer Comparing Two Chemoterapy. 
#' @name prodige
#' @rdname data-prodige
#'
#' @description Simulated data inspired from the PRODIGE trial comparing the survival of patients with metastatic pancreatic cancer
#' treated with FOLFIRINOX or gemcitabine . 
#'
#' \itemize{
#' \item \code{id}: study participant.
#' \item \code{treatment}: treatment arm: FOLFIRINOX (T) or gemcitabine (C).
#' \item \code{OS}: time from inclusion (say diagnosis) to end of follow-up.
#' \item \code{statusOS}: event triggering the end of follow-up: death (1), drop-out (0).
#' \item \code{PFS}: time from inclusion (say diagnosis) to progression of the disease or end of follow-up.
#' \item \code{statusPFS}: progression (1) or end of follow-up (0).
#' \item \code{toxicity}: most serious side effect observed during the follow-up (1 mild, 6 severe).
#' \item \code{sex}: male (M) or female (F)
#' }
#'
#' @author Brice Ozenne
#' @docType data
#' @usage data(prodige, package = "BuyseTest")
#' @references Conroy, Thierry, et al. "FOLFIRINOX versus gemcitabine for metastatic pancreatic cancer" New England Journal of Medicine (2011) 364(19):1817-25. doi: 10.1056/NEJMoa1011923.
#' @keywords datasets
"prodige"


##----------------------------------------------------------------------
### doc-data.R ends here

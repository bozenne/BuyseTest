### doc-data.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 13 2025 (11:30) 
## Version: 
## Last-Updated: Feb 16 2025 (18:05) 
##           By: Brice Ozenne
##     Update #: 34
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

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

##----------------------------------------------------------------------
### doc-data.R ends here

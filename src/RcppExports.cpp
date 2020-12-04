// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GPC_cpp
Rcpp::List GPC_cpp(arma::mat endpoint, arma::mat status, std::vector< arma::uvec > indexC, std::vector< arma::uvec > posC, std::vector< arma::uvec > indexT, std::vector< arma::uvec > posT, std::vector< double > threshold, arma::vec weight, arma::vec method, unsigned int D, unsigned int D_UTTE, unsigned int n_strata, arma::vec nUTTE_analyzedPeron_M1, std::vector<unsigned int> index_endpoint, std::vector<unsigned int> index_status, std::vector<int> index_UTTE, std::vector< std::vector< arma::mat > > list_survTimeC, std::vector< std::vector< arma::mat > > list_survTimeT, std::vector< std::vector< arma::mat > > list_survJumpC, std::vector< std::vector< arma::mat > > list_survJumpT, std::vector< arma::mat > list_lastSurv, arma::mat p_C, arma::mat p_T, std::vector< std::vector< arma::mat > > iid_survJumpC, std::vector< std::vector< arma::mat > > iid_survJumpT, double zeroPlus, int correctionUninf, bool hierarchical, int hprojection, bool neutralAsUninf, bool keepScore, bool precompute, int returnIID, int debug);
RcppExport SEXP _BuyseTest_GPC_cpp(SEXP endpointSEXP, SEXP statusSEXP, SEXP indexCSEXP, SEXP posCSEXP, SEXP indexTSEXP, SEXP posTSEXP, SEXP thresholdSEXP, SEXP weightSEXP, SEXP methodSEXP, SEXP DSEXP, SEXP D_UTTESEXP, SEXP n_strataSEXP, SEXP nUTTE_analyzedPeron_M1SEXP, SEXP index_endpointSEXP, SEXP index_statusSEXP, SEXP index_UTTESEXP, SEXP list_survTimeCSEXP, SEXP list_survTimeTSEXP, SEXP list_survJumpCSEXP, SEXP list_survJumpTSEXP, SEXP list_lastSurvSEXP, SEXP p_CSEXP, SEXP p_TSEXP, SEXP iid_survJumpCSEXP, SEXP iid_survJumpTSEXP, SEXP zeroPlusSEXP, SEXP correctionUninfSEXP, SEXP hierarchicalSEXP, SEXP hprojectionSEXP, SEXP neutralAsUninfSEXP, SEXP keepScoreSEXP, SEXP precomputeSEXP, SEXP returnIIDSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type endpoint(endpointSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type status(statusSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type indexC(indexCSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type posC(posCSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type indexT(indexTSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type posT(posTSEXP);
    Rcpp::traits::input_parameter< std::vector< double > >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type method(methodSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type D(DSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type D_UTTE(D_UTTESEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_strata(n_strataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nUTTE_analyzedPeron_M1(nUTTE_analyzedPeron_M1SEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type index_endpoint(index_endpointSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type index_status(index_statusSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type index_UTTE(index_UTTESEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type list_survTimeC(list_survTimeCSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type list_survTimeT(list_survTimeTSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type list_survJumpC(list_survJumpCSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type list_survJumpT(list_survJumpTSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::mat > >::type list_lastSurv(list_lastSurvSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p_C(p_CSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p_T(p_TSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type iid_survJumpC(iid_survJumpCSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type iid_survJumpT(iid_survJumpTSEXP);
    Rcpp::traits::input_parameter< double >::type zeroPlus(zeroPlusSEXP);
    Rcpp::traits::input_parameter< int >::type correctionUninf(correctionUninfSEXP);
    Rcpp::traits::input_parameter< bool >::type hierarchical(hierarchicalSEXP);
    Rcpp::traits::input_parameter< int >::type hprojection(hprojectionSEXP);
    Rcpp::traits::input_parameter< bool >::type neutralAsUninf(neutralAsUninfSEXP);
    Rcpp::traits::input_parameter< bool >::type keepScore(keepScoreSEXP);
    Rcpp::traits::input_parameter< bool >::type precompute(precomputeSEXP);
    Rcpp::traits::input_parameter< int >::type returnIID(returnIIDSEXP);
    Rcpp::traits::input_parameter< int >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(GPC_cpp(endpoint, status, indexC, posC, indexT, posT, threshold, weight, method, D, D_UTTE, n_strata, nUTTE_analyzedPeron_M1, index_endpoint, index_status, index_UTTE, list_survTimeC, list_survTimeT, list_survJumpC, list_survJumpT, list_lastSurv, p_C, p_T, iid_survJumpC, iid_survJumpT, zeroPlus, correctionUninf, hierarchical, hprojection, neutralAsUninf, keepScore, precompute, returnIID, debug));
    return rcpp_result_gen;
END_RCPP
}
// GPC2_cpp
Rcpp::List GPC2_cpp(arma::mat endpoint, arma::mat status, std::vector< arma::uvec > indexC, std::vector< arma::uvec > posC, std::vector< arma::uvec > indexT, std::vector< arma::uvec > posT, std::vector< double > threshold, arma::vec weight, arma::vec method, unsigned int D, unsigned int D_UTTE, unsigned int n_strata, arma::vec nUTTE_analyzedPeron_M1, std::vector<unsigned int> index_endpoint, std::vector<unsigned int> index_status, std::vector<int> index_UTTE, std::vector< std::vector< arma::mat > > list_survTimeC, std::vector< std::vector< arma::mat > > list_survTimeT, std::vector< std::vector< arma::mat > > list_survJumpC, std::vector< std::vector< arma::mat > > list_survJumpT, std::vector< arma::mat > list_lastSurv, arma::mat p_C, arma::mat p_T, std::vector< std::vector< arma::mat > > iid_survJumpC, std::vector< std::vector< arma::mat > > iid_survJumpT, double zeroPlus, int correctionUninf, bool hierarchical, int hprojection, bool neutralAsUninf, bool keepScore, bool precompute, int returnIID, int debug);
RcppExport SEXP _BuyseTest_GPC2_cpp(SEXP endpointSEXP, SEXP statusSEXP, SEXP indexCSEXP, SEXP posCSEXP, SEXP indexTSEXP, SEXP posTSEXP, SEXP thresholdSEXP, SEXP weightSEXP, SEXP methodSEXP, SEXP DSEXP, SEXP D_UTTESEXP, SEXP n_strataSEXP, SEXP nUTTE_analyzedPeron_M1SEXP, SEXP index_endpointSEXP, SEXP index_statusSEXP, SEXP index_UTTESEXP, SEXP list_survTimeCSEXP, SEXP list_survTimeTSEXP, SEXP list_survJumpCSEXP, SEXP list_survJumpTSEXP, SEXP list_lastSurvSEXP, SEXP p_CSEXP, SEXP p_TSEXP, SEXP iid_survJumpCSEXP, SEXP iid_survJumpTSEXP, SEXP zeroPlusSEXP, SEXP correctionUninfSEXP, SEXP hierarchicalSEXP, SEXP hprojectionSEXP, SEXP neutralAsUninfSEXP, SEXP keepScoreSEXP, SEXP precomputeSEXP, SEXP returnIIDSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type endpoint(endpointSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type status(statusSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type indexC(indexCSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type posC(posCSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type indexT(indexTSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::uvec > >::type posT(posTSEXP);
    Rcpp::traits::input_parameter< std::vector< double > >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type method(methodSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type D(DSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type D_UTTE(D_UTTESEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_strata(n_strataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nUTTE_analyzedPeron_M1(nUTTE_analyzedPeron_M1SEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type index_endpoint(index_endpointSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type index_status(index_statusSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type index_UTTE(index_UTTESEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type list_survTimeC(list_survTimeCSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type list_survTimeT(list_survTimeTSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type list_survJumpC(list_survJumpCSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type list_survJumpT(list_survJumpTSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::mat > >::type list_lastSurv(list_lastSurvSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p_C(p_CSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p_T(p_TSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type iid_survJumpC(iid_survJumpCSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::mat > > >::type iid_survJumpT(iid_survJumpTSEXP);
    Rcpp::traits::input_parameter< double >::type zeroPlus(zeroPlusSEXP);
    Rcpp::traits::input_parameter< int >::type correctionUninf(correctionUninfSEXP);
    Rcpp::traits::input_parameter< bool >::type hierarchical(hierarchicalSEXP);
    Rcpp::traits::input_parameter< int >::type hprojection(hprojectionSEXP);
    Rcpp::traits::input_parameter< bool >::type neutralAsUninf(neutralAsUninfSEXP);
    Rcpp::traits::input_parameter< bool >::type keepScore(keepScoreSEXP);
    Rcpp::traits::input_parameter< bool >::type precompute(precomputeSEXP);
    Rcpp::traits::input_parameter< int >::type returnIID(returnIIDSEXP);
    Rcpp::traits::input_parameter< int >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(GPC2_cpp(endpoint, status, indexC, posC, indexT, posT, threshold, weight, method, D, D_UTTE, n_strata, nUTTE_analyzedPeron_M1, index_endpoint, index_status, index_UTTE, list_survTimeC, list_survTimeT, list_survJumpC, list_survJumpT, list_lastSurv, p_C, p_T, iid_survJumpC, iid_survJumpT, zeroPlus, correctionUninf, hierarchical, hprojection, neutralAsUninf, keepScore, precompute, returnIID, debug));
    return rcpp_result_gen;
END_RCPP
}
// calcIntegralSurv_cpp
std::vector< double > calcIntegralSurv_cpp(const arma::mat& survival, double start, double lastSurv, double lastdSurv, bool returnDeriv, arma::colvec& derivSurv, arma::colvec& derivSurvD);
RcppExport SEXP _BuyseTest_calcIntegralSurv_cpp(SEXP survivalSEXP, SEXP startSEXP, SEXP lastSurvSEXP, SEXP lastdSurvSEXP, SEXP returnDerivSEXP, SEXP derivSurvSEXP, SEXP derivSurvDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type survival(survivalSEXP);
    Rcpp::traits::input_parameter< double >::type start(startSEXP);
    Rcpp::traits::input_parameter< double >::type lastSurv(lastSurvSEXP);
    Rcpp::traits::input_parameter< double >::type lastdSurv(lastdSurvSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDeriv(returnDerivSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type derivSurv(derivSurvSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type derivSurvD(derivSurvDSEXP);
    rcpp_result_gen = Rcpp::wrap(calcIntegralSurv_cpp(survival, start, lastSurv, lastdSurv, returnDeriv, derivSurv, derivSurvD));
    return rcpp_result_gen;
END_RCPP
}
// calcIntegralCif_cpp
double calcIntegralCif_cpp(const arma::mat& cifJump, double start_val, double stop_val, arma::rowvec cifTimeT, double lastCIF, int type, bool returnDeriv, arma::colvec& derivSurv, arma::colvec& derivSurvD);
RcppExport SEXP _BuyseTest_calcIntegralCif_cpp(SEXP cifJumpSEXP, SEXP start_valSEXP, SEXP stop_valSEXP, SEXP cifTimeTSEXP, SEXP lastCIFSEXP, SEXP typeSEXP, SEXP returnDerivSEXP, SEXP derivSurvSEXP, SEXP derivSurvDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type cifJump(cifJumpSEXP);
    Rcpp::traits::input_parameter< double >::type start_val(start_valSEXP);
    Rcpp::traits::input_parameter< double >::type stop_val(stop_valSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type cifTimeT(cifTimeTSEXP);
    Rcpp::traits::input_parameter< double >::type lastCIF(lastCIFSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type returnDeriv(returnDerivSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type derivSurv(derivSurvSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type derivSurvD(derivSurvDSEXP);
    rcpp_result_gen = Rcpp::wrap(calcIntegralCif_cpp(cifJump, start_val, stop_val, cifTimeT, lastCIF, type, returnDeriv, derivSurv, derivSurvD));
    return rcpp_result_gen;
END_RCPP
}
// calcIntegralSurv2_cpp
Rcpp::List calcIntegralSurv2_cpp(const std::vector<double>& time, const std::vector<double>& survival, const std::vector<double>& dSurvival, const std::vector<int>& index_survival, const std::vector<int>& index_dSurvival1, const std::vector<int>& index_dSurvival2, double lastSurv, double lastdSurv, bool iidNuisance, int nJump);
RcppExport SEXP _BuyseTest_calcIntegralSurv2_cpp(SEXP timeSEXP, SEXP survivalSEXP, SEXP dSurvivalSEXP, SEXP index_survivalSEXP, SEXP index_dSurvival1SEXP, SEXP index_dSurvival2SEXP, SEXP lastSurvSEXP, SEXP lastdSurvSEXP, SEXP iidNuisanceSEXP, SEXP nJumpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type survival(survivalSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dSurvival(dSurvivalSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type index_survival(index_survivalSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type index_dSurvival1(index_dSurvival1SEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type index_dSurvival2(index_dSurvival2SEXP);
    Rcpp::traits::input_parameter< double >::type lastSurv(lastSurvSEXP);
    Rcpp::traits::input_parameter< double >::type lastdSurv(lastdSurvSEXP);
    Rcpp::traits::input_parameter< bool >::type iidNuisance(iidNuisanceSEXP);
    Rcpp::traits::input_parameter< int >::type nJump(nJumpSEXP);
    rcpp_result_gen = Rcpp::wrap(calcIntegralSurv2_cpp(time, survival, dSurvival, index_survival, index_dSurvival1, index_dSurvival2, lastSurv, lastdSurv, iidNuisance, nJump));
    return rcpp_result_gen;
END_RCPP
}
// rowCumSum_cpp
arma::mat rowCumSum_cpp(const arma::mat X);
RcppExport SEXP _BuyseTest_rowCumSum_cpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(rowCumSum_cpp(X));
    return rcpp_result_gen;
END_RCPP
}
// rowCumProd_cpp
arma::mat rowCumProd_cpp(const arma::mat X);
RcppExport SEXP _BuyseTest_rowCumProd_cpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(rowCumProd_cpp(X));
    return rcpp_result_gen;
END_RCPP
}
// colCenter_cpp
arma::mat colCenter_cpp(const arma::mat X, const arma::colvec& center);
RcppExport SEXP _BuyseTest_colCenter_cpp(SEXP XSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(colCenter_cpp(X, center));
    return rcpp_result_gen;
END_RCPP
}
// rowCenter_cpp
arma::mat rowCenter_cpp(const arma::mat X, const arma::rowvec& center);
RcppExport SEXP _BuyseTest_rowCenter_cpp(SEXP XSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(rowCenter_cpp(X, center));
    return rcpp_result_gen;
END_RCPP
}
// colScale_cpp
arma::mat colScale_cpp(const arma::mat X, const arma::colvec& scale);
RcppExport SEXP _BuyseTest_colScale_cpp(SEXP XSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(colScale_cpp(X, scale));
    return rcpp_result_gen;
END_RCPP
}
// rowScale_cpp
arma::mat rowScale_cpp(const arma::mat X, const arma::rowvec& scale);
RcppExport SEXP _BuyseTest_rowScale_cpp(SEXP XSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rowScale_cpp(X, scale));
    return rcpp_result_gen;
END_RCPP
}
// colMultiply_cpp
arma::mat colMultiply_cpp(const arma::mat X, const arma::colvec& scale);
RcppExport SEXP _BuyseTest_colMultiply_cpp(SEXP XSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(colMultiply_cpp(X, scale));
    return rcpp_result_gen;
END_RCPP
}
// rowMultiply_cpp
arma::mat rowMultiply_cpp(const arma::mat X, const arma::rowvec& scale);
RcppExport SEXP _BuyseTest_rowMultiply_cpp(SEXP XSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rowMultiply_cpp(X, scale));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BuyseTest_GPC_cpp", (DL_FUNC) &_BuyseTest_GPC_cpp, 34},
    {"_BuyseTest_GPC2_cpp", (DL_FUNC) &_BuyseTest_GPC2_cpp, 34},
    {"_BuyseTest_calcIntegralSurv_cpp", (DL_FUNC) &_BuyseTest_calcIntegralSurv_cpp, 7},
    {"_BuyseTest_calcIntegralCif_cpp", (DL_FUNC) &_BuyseTest_calcIntegralCif_cpp, 9},
    {"_BuyseTest_calcIntegralSurv2_cpp", (DL_FUNC) &_BuyseTest_calcIntegralSurv2_cpp, 10},
    {"_BuyseTest_rowCumSum_cpp", (DL_FUNC) &_BuyseTest_rowCumSum_cpp, 1},
    {"_BuyseTest_rowCumProd_cpp", (DL_FUNC) &_BuyseTest_rowCumProd_cpp, 1},
    {"_BuyseTest_colCenter_cpp", (DL_FUNC) &_BuyseTest_colCenter_cpp, 2},
    {"_BuyseTest_rowCenter_cpp", (DL_FUNC) &_BuyseTest_rowCenter_cpp, 2},
    {"_BuyseTest_colScale_cpp", (DL_FUNC) &_BuyseTest_colScale_cpp, 2},
    {"_BuyseTest_rowScale_cpp", (DL_FUNC) &_BuyseTest_rowScale_cpp, 2},
    {"_BuyseTest_colMultiply_cpp", (DL_FUNC) &_BuyseTest_colMultiply_cpp, 2},
    {"_BuyseTest_rowMultiply_cpp", (DL_FUNC) &_BuyseTest_rowMultiply_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_BuyseTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

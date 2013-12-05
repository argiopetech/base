// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// initBase
void initBase(std::string modelDir, int msFilter, int msModel, int wdModel, int ifmr);
RcppExport SEXP BASE9_initBase(SEXP modelDirSEXP, SEXP msFilterSEXP, SEXP msModelSEXP, SEXP wdModelSEXP, SEXP ifmrSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< std::string >::type modelDir(modelDirSEXP );
        Rcpp::traits::input_parameter< int >::type msFilter(msFilterSEXP );
        Rcpp::traits::input_parameter< int >::type msModel(msModelSEXP );
        Rcpp::traits::input_parameter< int >::type wdModel(wdModelSEXP );
        Rcpp::traits::input_parameter< int >::type ifmr(ifmrSEXP );
        initBase(modelDir, msFilter, msModel, wdModel, ifmr);
    }
    return R_NilValue;
END_RCPP
}
// setClusterParameters
void setClusterParameters(double age, double feh, double distMod, double av, double y, double carbonicity);
RcppExport SEXP BASE9_setClusterParameters(SEXP ageSEXP, SEXP fehSEXP, SEXP distModSEXP, SEXP avSEXP, SEXP ySEXP, SEXP carbonicitySEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type age(ageSEXP );
        Rcpp::traits::input_parameter< double >::type feh(fehSEXP );
        Rcpp::traits::input_parameter< double >::type distMod(distModSEXP );
        Rcpp::traits::input_parameter< double >::type av(avSEXP );
        Rcpp::traits::input_parameter< double >::type y(ySEXP );
        Rcpp::traits::input_parameter< double >::type carbonicity(carbonicitySEXP );
        setClusterParameters(age, feh, distMod, av, y, carbonicity);
    }
    return R_NilValue;
END_RCPP
}
// changeModels
void changeModels(int msFilter, int msModel, int wdModel, int ifmr);
RcppExport SEXP BASE9_changeModels(SEXP msFilterSEXP, SEXP msModelSEXP, SEXP wdModelSEXP, SEXP ifmrSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type msFilter(msFilterSEXP );
        Rcpp::traits::input_parameter< int >::type msModel(msModelSEXP );
        Rcpp::traits::input_parameter< int >::type wdModel(wdModelSEXP );
        Rcpp::traits::input_parameter< int >::type ifmr(ifmrSEXP );
        changeModels(msFilter, msModel, wdModel, ifmr);
    }
    return R_NilValue;
END_RCPP
}
// setIFMRParameters
void setIFMRParameters(double intercept, double slope, double quadCoef);
RcppExport SEXP BASE9_setIFMRParameters(SEXP interceptSEXP, SEXP slopeSEXP, SEXP quadCoefSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type intercept(interceptSEXP );
        Rcpp::traits::input_parameter< double >::type slope(slopeSEXP );
        Rcpp::traits::input_parameter< double >::type quadCoef(quadCoefSEXP );
        setIFMRParameters(intercept, slope, quadCoef);
    }
    return R_NilValue;
END_RCPP
}
// evolve
std::array<double, 8> evolve(double mass);
RcppExport SEXP BASE9_evolve(SEXP massSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type mass(massSEXP );
        std::array<double, 8> __result = evolve(mass);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

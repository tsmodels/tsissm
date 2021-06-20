#ifndef _FILTERS_H
#define _FILTERS_H
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

//

Rcpp::List issestimation(Rcpp::NumericVector& , Rcpp::NumericVector& , Rcpp::NumericVector& ,
                         Rcpp::NumericVector& , Rcpp::NumericVector& , Rcpp::NumericVector& ,
                         Rcpp::IntegerVector& , Rcpp::NumericVector& , Rcpp::NumericVector& ,
                         Rcpp::NumericVector&, Rcpp::NumericVector&);
Rcpp::List issfilter(Rcpp::NumericVector& , Rcpp::NumericVector& , Rcpp::NumericVector& ,
                     Rcpp::NumericVector& , Rcpp::NumericVector& , Rcpp::NumericVector& ,
                     Rcpp::IntegerVector& , Rcpp::NumericVector& , Rcpp::NumericVector& ,
                     Rcpp::NumericVector& , Rcpp::NumericVector&, Rcpp::NumericVector&);
Rcpp::List isspredict(Rcpp::NumericVector& , Rcpp::NumericVector& , Rcpp::NumericVector& , Rcpp::NumericVector& ,
                      Rcpp::NumericVector& , Rcpp::NumericVector& , Rcpp::IntegerVector& ,
                      Rcpp::NumericVector& , Rcpp::NumericVector& , Rcpp::NumericVector&);

// inline/template functions

template<typename T>
inline
T
boxcox(const T& y, const double lambda)
{
    T x(y.n_elem);

    if (lambda < 1.0e-9) {
        x = arma::log(y);
    } else {
        x = (arma::sign(y) % arma::pow(arma::abs(y), lambda) - 1.0)/lambda;
    }

    return(x);
}

inline
double
loglikfun(const int time, const double lambda, const arma::vec& error, const arma::vec& y)
{
    return static_cast<double>(time) * std::log(arma::accu( error % error )) - 2.0 * (lambda - 1.0) * arma::accu(arma::log(y));
}

inline
void
initstate(const int time, const arma::vec& w, const arma::vec& ytrans, const arma::vec& g,
          const arma::vec& kappa, const arma::mat& F, const arma::mat& D, const arma::mat& xreg,
          arma::vec& yaux, arma::vec& eaux, arma::mat& waux, arma::mat& xaux, const arma::vec& good)
{
    for (int i = 1; i < time; i++) {
        yaux(i) = arma::as_scalar(xaux.row(i-1) * w + xreg.row(i) * kappa);
        if (good(i) > 0.5) {
            eaux(i) = ytrans(i) - yaux(i);
        } else {
            eaux(i) = 0.0;
        }
        xaux.row(i) = arma::trans(F * xaux.row(i-1).t() + g * eaux(i));
        waux.row(i) = waux.row(i-1) * D;
    }
}

inline
void
seedstate(const int stype, const int sstart, const int send,
          const arma::mat& waux, const arma::vec& eaux,
          arma::rowvec& xseed)
{
    // trigonometric
    switch(stype) {
        case 1: {
            xseed.head_cols(send) = arma::solve(waux.head_cols(send), eaux).t();
            break;
        }
        case 2: {
            double frequency = static_cast<double>(send - sstart + 1);
            xseed.head_cols(send - 1) = arma::solve(waux.head_cols(send - 1), eaux).t();

            // normalize the seasonal seed by the frequency
            double normseas = arma::accu(xseed.subvec(sstart - 1, send - 2))/frequency;

            xseed.subvec(sstart - 1, send - 2) = xseed.subvec(sstart - 1, send - 2) - normseas;
            xseed(send - 1) = -1.0 * normseas;

            break;
        }
        default: {
            xseed.head_cols(send) = arma::solve(waux.head_cols(send), eaux).t();
            break;
        }
    }
}

// macro functions

#define ISS_ESTIMATION_SETUP                                                                    \
const int states = mdim[0];                                                                     \
const int time  = mdim[1];                                                                      \
const int stype = mdim[2];                                                                      \
const int sstart = mdim[3];                                                                     \
const int send = mdim[4];                                                                       \
const int narma = mdim[5];                                                                      \
const int nxreg = mdim[6];                                                                      \
const double lambda = lambda_[0];                                                               \
arma::vec y = Rcpp::as<arma::vec>(y_);                                                          \
arma::vec good = Rcpp::as<arma::vec>(good_);                                                    \
arma::vec yhat = arma::zeros(time + 1);                                                         \
arma::vec ytrans = boxcox(y, lambda);                                                           \
arma::vec error = arma::zeros<arma::vec>(time + 1);                                             \
arma::mat F0 = arma::reshape(Rcpp::as<arma::vec>(f0_), states, states);                         \
arma::mat F1 = arma::reshape(Rcpp::as<arma::vec>(f1_), states, states);                         \
arma::mat F2 = arma::reshape(Rcpp::as<arma::vec>(f2_), states, states);                         \
arma::mat F = F0 % F1 % F2;                                                                     \
arma::vec w = Rcpp::as<arma::vec>(w_);                                                          \
arma::vec g = Rcpp::as<arma::vec>(g_);                                                          \
arma::mat x = arma::zeros<arma::mat>(time + 1, states);                                         \
arma::mat xreg = arma::reshape(Rcpp::as<arma::vec>(xreg_), time, nxreg);                        \
arma::vec kappa = Rcpp::as<arma::vec>(kappa_);                                                  \
arma::mat D = F -  g * w.t();                                                                   \
arma::cx_vec eigval;                                                                            \
arma::cx_mat eigvec;                                                                            \
arma::rowvec x0 = arma::zeros<arma::rowvec>(states);                                            \
arma::mat waux = arma::zeros<arma::mat>(time, states);                                          \
arma::mat xaux = arma::zeros<arma::mat>(time, states);                                          \
arma::vec eaux = arma::zeros<arma::vec>(time);                                                  \
arma::vec yaux = arma::zeros<arma::vec>(time);                                                  \
yaux(0) = arma::as_scalar(x0 * w + xreg.row(0) * kappa);                                        \
eaux(0) = ytrans(0) - yaux(0);                                                                  \
xaux.row(0) = arma::trans(F * x0.t() + g * eaux(0));                                            \
waux.row(0) = w.t();                                                                            \
arma::rowvec xseed(send + narma, arma::fill::zeros);                                            \
//                                                                                              \

// macro functions
#define ISS_FILTER_SETUP                                                                        \
const int states = mdim[0];                                                                     \
const int time  = mdim[1];                                                                      \
const int stype = mdim[2];                                                                      \
const int sstart = mdim[3];                                                                     \
const int send = mdim[4];                                                                       \
const int narma = mdim[5];                                                                      \
const int nxreg = mdim[6];                                                                      \
const double lambda = lambda_[0];                                                               \
arma::vec y = Rcpp::as<arma::vec>(y_);                                                          \
arma::vec good = Rcpp::as<arma::vec>(good_);                                                    \
arma::vec yhat = arma::zeros(time + 1);                                                         \
arma::vec ytrans = boxcox(y, lambda);                                                           \
arma::vec error = arma::zeros<arma::vec>(time + 1);                                             \
arma::mat F0 = arma::reshape(Rcpp::as<arma::vec>(f0_), states, states);                         \
arma::mat F1 = arma::reshape(Rcpp::as<arma::vec>(f1_), states, states);                         \
arma::mat F2 = arma::reshape(Rcpp::as<arma::vec>(f2_), states, states);                         \
arma::mat F = F0 % F1 % F2;                                                                     \
arma::vec w = Rcpp::as<arma::vec>(w_);                                                          \
arma::vec g = Rcpp::as<arma::vec>(g_);                                                          \
arma::mat x = arma::zeros<arma::mat>(time + 1, states);                                         \
arma::mat xreg = arma::reshape(Rcpp::as<arma::vec>(xreg_), time, nxreg);                        \
arma::vec kappa = Rcpp::as<arma::vec>(kappa_);                                                  \
arma::mat D = F -  g * w.t();                                                                   \
arma::rowvec xseed = Rcpp::as<arma::rowvec>(xseed_);                                            \
ytrans = join_cols(arma::zeros<arma::vec>(1), ytrans);                                          \
xreg = join_cols(arma::zeros<arma::mat>(1, nxreg), xreg);                                       \
x.row(0) = xseed;                                                                               \
//                                                                                              \

#define ISS_PREDICT_SETUP                                                                       \
const int states = mdim[0];                                                                     \
const int nsim  = mdim[1];                                                                      \
const int horizon  = mdim[2];                                                                   \
const int nxreg = mdim[3];                                                                      \
arma::rowvec xseed = Rcpp::as<arma::rowvec>(xseed_);                                            \
arma::mat ysim = arma::zeros<arma::mat>(nsim, horizon + 1);                                     \
arma::mat error = arma::reshape(Rcpp::as<arma::vec>(error_), nsim, horizon);                    \
error = join_rows(arma::zeros<arma::mat>(nsim, 1), error);                                      \
arma::mat F0 = arma::reshape(Rcpp::as<arma::vec>(f0_), states, states);                         \
arma::mat F1 = arma::reshape(Rcpp::as<arma::vec>(f1_), states, states);                         \
arma::mat F2 = arma::reshape(Rcpp::as<arma::vec>(f2_), states, states);                         \
arma::mat F = F0 % F1 % F2;                                                                     \
arma::vec w = Rcpp::as<arma::vec>(w_);                                                          \
arma::vec g = Rcpp::as<arma::vec>(g_);                                                          \
arma::cube xstates = arma::zeros<arma::cube>(horizon + 1,states, nsim);                         \
arma::mat xreg = arma::reshape(Rcpp::as<arma::vec>(xreg_), horizon, nxreg);                     \
xreg = join_cols(arma::zeros<arma::mat>(1, nxreg), xreg);                                       \
arma::vec kappa = Rcpp::as<arma::vec>(kappa_);                                                  \
arma::mat D = F -  g * w.t();                                                                   \
//                                                                                              \

#endif

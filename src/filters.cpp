#include "filters.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::export]]
Rcpp::List issestimation(NumericVector& f0_, NumericVector& f1_, NumericVector& f2_, NumericVector& w_,
                        NumericVector& g_, NumericVector& y_, IntegerVector& mdim, NumericVector& lambda_,
                        NumericVector& xreg_, NumericVector& kappa_, NumericVector& good_)
{

    try {
        ISS_ESTIMATION_SETUP
        
        initstate(time, w, ytrans, g, kappa, F, D, xreg, yaux, eaux, waux, xaux, good);
        
        // check if any waux are illegal values
        if(waux.has_nan()) {
            Rcpp::List output = Rcpp::List::create(Rcpp::Named("condition") = 1.0,
                                                   Rcpp::Named("loglik") = 1.0e6);
            return(output);
        }
        // solve for the initial states (beta coefficients error[t] ~ w[t-1] D)
        eig_gen(eigval, eigvec, D, "balance");
        arma::vec realeig = arma::abs(arma::real(eigval));
        arma::uvec test = find(realeig > 1.01);
        double lerr = 0.0;
        double logy = 0.0;
        double ngood = 0.0;
        //Rcout << "The value of realeig : " << realeig << "\n";
        
        if (test.n_elem > 0) {
            Rcpp::List output = Rcpp::List::create(Rcpp::Named("condition") = 1.0,
                                                   Rcpp::Named("loglik") = 1.0e6);
            return(output);
        } else {
            // calculate the initial seed states using the backwards recursion
            // if no xreg in equation then xreg is a one column matrix with zeros
            // and kappa is a vector with one element of zero.
            // filter model
            // append with extra time (0)
            seedstate(stype, sstart, send, waux, eaux, xseed);

            ytrans = join_cols(arma::zeros<arma::vec>(1), ytrans);
            good = join_cols(arma::ones<arma::vec>(1), good);
            xreg = join_cols(arma::zeros<arma::mat>(1, nxreg), xreg);
            x.row(0) = xseed;
            for (int i = 1; i <= time; i++) {
                yhat(i) = arma::as_scalar(x.row(i-1) * w + xreg.row(i) * kappa);
                if (good(i) > 0.5) {
                    error(i) = ytrans(i) - yhat(i);
                    lerr+= arma::as_scalar(error(i) * error(i));
                    logy+= std::log(arma::as_scalar(y(i-1)));
                    ngood+=1;
                } else {
                    error(i) = 0.0;
                }
                x.row(i) = arma::trans(F * x.row(i-1).t() + g * error(i));
            }
            double loglik = 0.0;
            if (flag == 1.0) {
                loglik = ngood * std::log(lerr);
            } else {
                loglik = ngood * std::log(lerr) - 2.0 * (lambda - 1.0) * logy;
            }
            Rcpp::List output = Rcpp::List::create(Rcpp::Named("xseed") = wrap(xseed),
                                                   Rcpp::Named("condition") = 0.0,
                                                   Rcpp::Named("loglik") = wrap(loglik));
            
            return(output);
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsissm--> iss estimation exception (unknown reason)" );
    }
    return R_NilValue;
}


// [[Rcpp::export]]
Rcpp::List issfilter(NumericVector& f0_, NumericVector& f1_, NumericVector& f2_, NumericVector& w_,
                     NumericVector& g_, NumericVector& y_, IntegerVector& mdim, NumericVector& lambda_,
                     NumericVector& xreg_, NumericVector& kappa_, NumericVector& xseed_, NumericVector& good_)
{

    try {
        ISS_FILTER_SETUP
        // filter model
        good = join_cols(arma::zeros<arma::vec>(1), good);
        double lerr = 0.0;
        double logy = 0.0;
        double ngood = 0.0;
        for (int i = 1; i <= time; i++) {
            yhat(i) = arma::as_scalar(x.row(i-1) * w + xreg.row(i) * kappa);
            if (good(i) > 0.5) {
                error(i) = ytrans(i) - yhat(i);
                lerr+= arma::as_scalar(error(i) * error(i));
                logy+= std::log(arma::as_scalar(y(i-1)));
                ngood+=1;
            } else {
                error(i) = 0.0;
            }
            x.row(i) = arma::trans(F * x.row(i-1).t() + g * error(i));
        }
        double loglik = 0.0;
        if (flag == 1.0) {
            loglik = ngood * std::log(lerr);
        } else {
            loglik = ngood * std::log(lerr) - 2.0 * (lambda - 1.0) * logy;
        }
        Rcpp::List output = Rcpp::List::create(Rcpp::Named("xseed") = wrap(xseed),
                                               Rcpp::Named("states") = wrap(x),
                                               Rcpp::Named("w") = wrap(w),
                                               Rcpp::Named("g") = wrap(g),
                                               Rcpp::Named("F") = wrap(F),
                                               Rcpp::Named("D") = wrap(D),
                                               Rcpp::Named("fitted") = wrap(yhat),
                                               Rcpp::Named("transformed") = wrap(ytrans),
                                               Rcpp::Named("error") = wrap(error),
                                               Rcpp::Named("loglik") = wrap(loglik));
        
        return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsissm--> iss filter exception (unknown reason)" );
    }
    return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List isspredict(NumericVector& f0_, NumericVector& f1_, NumericVector& f2_, NumericVector& w_,
                      NumericVector& g_, NumericVector& error_, IntegerVector& mdim,
                      NumericVector& xreg_, NumericVector& kappa_, NumericVector& xseed_)
{

    try {
        ISS_PREDICT_SETUP
        // filter model
        // append with extra time (0)

        // TODO: consider wrapping this (j) in an OpenMP loop
        for (int j = 0; j < nsim; j++) {
            arma::mat x = xstates.slice(j);
            x.row(0) = xseed;

            for (int i = 1; i <= horizon; i++) {
                ysim(j,i) = arma::as_scalar(x.row(i-1) * w + xreg.row(i) * kappa) + error(j,i);
                x.row(i) = arma::trans(F * x.row(i-1).t() + g * error(j,i));
            }

            xstates.slice(j) = x;
        }

        Rcpp::List output = Rcpp::List::create(Rcpp::Named("states") = wrap(xstates),
                                               Rcpp::Named("simulated") = wrap(ysim));
        
        return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsissm--> iss predict exception (unknown reason)" );
    }
    return R_NilValue;
}

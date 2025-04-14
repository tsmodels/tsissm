#include<RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp14)]]
inline Eigen::VectorXd boxcox(const Eigen::VectorXd& y, const double lambda) {
    Eigen::VectorXd x(y.size());
    if (lambda < 1e-12) {
        x = y.array().log();
    } else {
        x = (pow(y.array(), lambda) - 1.0) / lambda;
    }
    return x;
}

inline Eigen::MatrixXd vec2mat(const Eigen::VectorXd &vec, int nrow, int ncol) {
    // Eigen::Map is column-major by default.
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), nrow, ncol);
}

// [[Rcpp::export(.initialize_states)]]
Eigen::MatrixXd initialize_states(const Eigen::Map<Eigen::VectorXd>& y, 
                                  const Eigen::Map<Eigen::VectorXd>& good, 
                                  const Eigen::Map<Eigen::VectorXi>& valid_index,
                                  const Eigen::Map<Eigen::MatrixXd>& X, 
                                  const Eigen::Map<Eigen::VectorXd>& kappa, 
                                  Eigen::Map<Eigen::VectorXd> V,
                                  Eigen::Map<Eigen::VectorXd> allpars,
                                  const Eigen::Map<Eigen::VectorXi>& findex, 
                                  const Eigen::Map<Eigen::VectorXi>& fpindex, 
                                  const Eigen::Map<Eigen::VectorXi>& ppindex, 
                                  const Eigen::Map<Eigen::VectorXi>& fshape,
                                  const Eigen::Map<Eigen::VectorXi>& modeli,
                                  const Eigen::Map<Eigen::VectorXd>& issmpars,
                                  const double lambda) {
    
    int n_states = modeli(0);
    int n_seasonal = modeli(3);
    int n_arma = modeli(4);
    int timesteps = y.size();
    
    Eigen::VectorXd ytrans = boxcox(y, lambda);

    for (int i = 0; i < ppindex.size(); i++) {
        allpars(ppindex(i)) = issmpars(i);
    }
    for (int i = 0; i < findex.size(); i++) {
        V(findex(i)) = allpars(fpindex(i));
    }
    
    Eigen::VectorXd F0tmp = V.segment(fshape(0), fshape(1));
    Eigen::VectorXd F1tmp = V.segment(fshape(2), fshape(3));
    Eigen::VectorXd F2tmp = V.segment(fshape(4), fshape(5));
    
    Eigen::MatrixXd F0 = vec2mat(F0tmp, n_states, n_states);
    Eigen::MatrixXd F1 = vec2mat(F1tmp, n_states, n_states);
    Eigen::MatrixXd F2 = vec2mat(F2tmp, n_states, n_states);
    
    Eigen::MatrixXd F = (F0.array() * F1.array() * F2.array()).matrix();
    
    Eigen::VectorXd Wvec = V.segment(fshape(6), fshape(7));  // size = modeli(0)
    Eigen::VectorXd Gvec = V.segment(fshape(8), fshape(9));  // size = modeli(0)
    Eigen::MatrixXd Wmat = Eigen::Map<const Eigen::MatrixXd>(Wvec.data(), n_states, 1);
    Eigen::MatrixXd Gmat = Eigen::Map<const Eigen::MatrixXd>(Gvec.data(), n_states, 1);
    Eigen::MatrixXd B = Gmat * Wmat.transpose();
    Eigen::MatrixXd D = F - B;
    
    Eigen::VectorXd xreg = X * kappa;
    
    Eigen::MatrixXd xaux = Eigen::MatrixXd::Zero(timesteps, n_states);
    Eigen::MatrixXd waux = Eigen::MatrixXd::Zero(timesteps, n_states);
    Eigen::VectorXd eaux = Eigen::VectorXd::Zero(timesteps);
    Eigen::VectorXd yaux = Eigen::VectorXd::Zero(timesteps);

    yaux(0)      = xreg(0);
    waux.row(0)  = Wvec.transpose();
    eaux(0)      = ytrans(0) - yaux(0);
    xaux.row(0)  = (Gmat * eaux(0)).transpose();
    
    for (int i = 1; i < timesteps; i++) {
        Eigen::RowVectorXd x_prev = xaux.row(i - 1);
        double yhat_i = x_prev.dot(Wvec) + xreg(i); // same as (x_prev.array() * Wvec.array()).sum() + xreg(i)
        yaux(i) = yhat_i;
        // error
        if (good(i) > 0.5) {
            eaux(i) = ytrans(i) - yaux(i);
        } else {
            eaux(i) = eaux(i - 1);
        }
        Eigen::VectorXd Fx = F * x_prev.transpose();
        Eigen::VectorXd Gc = Gmat * eaux(i);
        xaux.row(i) = (Fx + Gc).transpose();
        waux.row(i) = waux.row(i - 1) * D;
    }
    Eigen::MatrixXd Bsub(valid_index.size(), n_states);
    for (int i = 0; i < valid_index.size(); i++) {
        Bsub.row(i) = waux.row(valid_index(i));
    }
    
    Eigen::MatrixXd A = Bsub.leftCols(n_seasonal);
    Eigen::VectorXd E(valid_index.size());
    for (int i = 0; i < valid_index.size(); i++) {
        E(i) = eaux(valid_index(i));
    }
    Eigen::VectorXd xseed = A.householderQr().solve(E);
    Eigen::MatrixXd init_states = Eigen::MatrixXd::Zero(1, n_seasonal + n_arma);
    init_states.block(0, 0, 1, n_seasonal) = xseed.transpose();
    return init_states;
}



Eigen::VectorXd garchrec(const Eigen::VectorXd& alpha, const Eigen::VectorXd& beta, const Eigen::VectorXd& error, 
                         const Eigen::VectorXd& initial_arch, const Eigen::VectorXd& initial_variance, const Eigen::VectorXi& vmodel, const double target_omega) {
    const int timesteps = error.size() + vmodel(0);
    Eigen::VectorXd residuals = Eigen::VectorXd::Zero(timesteps);
    Eigen::VectorXd sigma_squared = Eigen::VectorXd::Zero(timesteps);
    
    int start = vmodel(0);
    residuals.segment(start, error.size()) = error;
    Eigen::VectorXd residuals_squared = residuals.array().square();
    
    for (int j = 0; j < vmodel(0); ++j) {
        sigma_squared(j) += initial_variance(j);
        residuals_squared(j) += initial_arch(j);
    }
    
    for (int i = start; i < timesteps; ++i) {
        sigma_squared(i) += target_omega;
        for (int j = 0; j < vmodel(1); ++j) {
            int lag = i - j - 1;
            sigma_squared(i) += alpha(j) * residuals_squared(lag);
        }
        for (int j = 0; j < vmodel(2); ++j) {
            int lag = i - j - 1;
            sigma_squared(i) += beta(j) * sigma_squared(lag);
        }
    }
    Eigen::VectorXd sigma = sigma_squared.array().sqrt();
    return sigma;
}

inline void filter_loop(const Eigen::MatrixXd& F,
                        const Eigen::Ref<const Eigen::VectorXd>& g,
                        const Eigen::Ref<const Eigen::VectorXd>& w,
                        const Eigen::Ref<const Eigen::VectorXd>& kappa,
                        const Eigen::Ref<const Eigen::MatrixXd>& X,
                        const Eigen::Ref<const Eigen::VectorXd>& good,
                        const Eigen::Ref<const Eigen::VectorXd>& ytrans,
                        int time,
                        Eigen::MatrixXd& x,
                        Eigen::VectorXd& error,
                        Eigen::VectorXd& yhat,
                        double& ngood) {
    for (int i = 1; i <= time; ++i) {
        yhat(i) = x.row(i - 1).dot(w) + X.row(i).dot(kappa);
        if (good(i) > 0.5) {
            error(i) = ytrans(i) - yhat(i);
            ngood += 1.0;
        } else {
            error(i) = 0.0;
        }
        x.row(i) = F * x.row(i - 1).transpose() + g * error(i);
    }
}

// [[Rcpp::export(.issfilter_constant)]]
Rcpp::List issfilter_constant(Eigen::Map<Eigen::MatrixXd>& F0, Eigen::Map<Eigen::MatrixXd>& F1,
                              Eigen::Map<Eigen::MatrixXd>& F2, Eigen::Map<Eigen::VectorXd>& w,
                              Eigen::Map<Eigen::VectorXd>& g, Eigen::Map<Eigen::VectorXd>& y,
                              Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::VectorXd>& kappa,
                              Eigen::Map<Eigen::VectorXd>& xseed, Eigen::Map<Eigen::VectorXd>& good,
                              Rcpp::IntegerVector& mdim, const double lambda) {
    
    try {
        const int states = mdim[0];
        const int time  = mdim[1];
        Eigen::VectorXd yhat = Eigen::VectorXd::Zero(y.size());
        Eigen::VectorXd ytrans = Eigen::VectorXd::Zero(y.size());
        ytrans = boxcox(y, lambda);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(y.size());
        Eigen::MatrixXd F = F0.array() * F1.array() * F2.array();
        Eigen::MatrixXd x = Eigen::MatrixXd::Zero(time + 1, states);
        Eigen::MatrixXd D = F - g * w.transpose();
        x.row(0) = xseed;
        double ngood = 0.0;
        
        // Call the helper function to perform the loop
        filter_loop(F, g, w, kappa, X, good, ytrans, time, x, error, yhat, ngood);
        
        Rcpp::List output = Rcpp::List::create(
            Rcpp::Named("xseed") =  Rcpp::wrap(xseed),
            Rcpp::Named("states") =  Rcpp::wrap(x),
            Rcpp::Named("w") =  Rcpp::wrap(w),
            Rcpp::Named("g") =  Rcpp::wrap(g),
            Rcpp::Named("F") =  Rcpp::wrap(F),
            Rcpp::Named("D") =  Rcpp::wrap(D),
            Rcpp::Named("fitted") =  Rcpp::wrap(yhat),
            Rcpp::Named("transformed") =  Rcpp::wrap(ytrans),
            Rcpp::Named("error") =  Rcpp::wrap(error)
        );
        return output;
    } catch (std::exception& ex) {
        forward_exception_to_r(ex);
    } catch (...) {
        Rf_error("iss constant filter exception (unknown reason)");
    }
    return R_NilValue;
}


// [[Rcpp::export(.issfilter_dynamic)]]
Rcpp::List issfilter_dynamic(Eigen::Map<Eigen::MatrixXd>& F0, Eigen::Map<Eigen::MatrixXd>& F1,
                              Eigen::Map<Eigen::MatrixXd>& F2, Eigen::Map<Eigen::VectorXd>& w,
                              Eigen::Map<Eigen::VectorXd>& g, Eigen::Map<Eigen::VectorXd>& y,
                              Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::VectorXd>& kappa,
                              Eigen::Map<Eigen::VectorXd>& xseed, Eigen::Map<Eigen::VectorXd>& good,
                              Eigen::Map<Eigen::VectorXd>& eta, Eigen::Map<Eigen::VectorXd>& delta,
                              Rcpp::IntegerVector& mdim, Eigen::Map<Eigen::VectorXi>& vmodel, 
                              Eigen::Map<Eigen::VectorXd>& initial_arch, Eigen::Map<Eigen::VectorXd>& initial_variance,
                              const double omega, const double lambda) {
    
    try {
        const int states = mdim[0];
        const int time  = mdim[1];
        Eigen::VectorXd yhat = Eigen::VectorXd::Zero(y.size());
        Eigen::VectorXd ytrans = Eigen::VectorXd::Zero(y.size());
        ytrans = boxcox(y, lambda);
        Eigen::VectorXd error = Eigen::VectorXd::Zero(y.size());
        Eigen::MatrixXd F = F0.array() * F1.array() * F2.array();
        Eigen::MatrixXd x = Eigen::MatrixXd::Zero(time + 1, states);
        Eigen::MatrixXd D = F - g * w.transpose();
        x.row(0) = xseed;
        double ngood = 0.0;
        
        // Call the helper function to perform the ISSM loop
        filter_loop(F, g, w, kappa, X, good, ytrans, time, x, error, yhat, ngood);
        Eigen::VectorXd e = error.segment(1, error.size() - 1);
        // Call the helper function to perform the GARCH loop
        Eigen::VectorXd sigma = garchrec(eta, delta, e, initial_arch, initial_variance, vmodel, omega);
        
        Rcpp::List output = Rcpp::List::create(
            Rcpp::Named("xseed") =  Rcpp::wrap(xseed),
            Rcpp::Named("states") =  Rcpp::wrap(x),
            Rcpp::Named("w") =  Rcpp::wrap(w),
            Rcpp::Named("g") =  Rcpp::wrap(g),
            Rcpp::Named("F") =  Rcpp::wrap(F),
            Rcpp::Named("D") =  Rcpp::wrap(D),
            Rcpp::Named("fitted") =  Rcpp::wrap(yhat),
            Rcpp::Named("transformed") =  Rcpp::wrap(ytrans),
            Rcpp::Named("error") =  Rcpp::wrap(error),
            Rcpp::Named("sigma") =  Rcpp::wrap(sigma));
        return output;
    } catch (std::exception& ex) {
        forward_exception_to_r(ex);
    } catch (...) {
        Rf_error("iss dynamic filter exception (unknown reason)");
    }
    return R_NilValue;
}


// [[Rcpp::export(.isspredict_constant)]]
Rcpp::List isspredict_constant(Eigen::Map<Eigen::MatrixXd>& F0, Eigen::Map<Eigen::MatrixXd>& F1,
                               Eigen::Map<Eigen::MatrixXd>& F2, Eigen::Map<Eigen::VectorXd>& w,
                               Eigen::Map<Eigen::VectorXd>& g, Eigen::Map<Eigen::MatrixXd>& E,
                               Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::VectorXd>& kappa,
                               Eigen::Map<Eigen::VectorXd>& xseed, Rcpp::IntegerVector& mdim)
{
    
    try {
        const int states = mdim[0];
        const int nsim  = mdim[1];
        const int horizon  = mdim[2];
        Eigen::MatrixXd F = F0.cwiseProduct(F1).cwiseProduct(F2);
        // Initialize output structures
        Eigen::MatrixXd ysim = Eigen::MatrixXd::Zero(nsim, horizon + 1);
        Rcpp::List xstates(nsim);
        
        // Main simulation loop
        for(int j = 0; j < nsim; ++j) {
            Eigen::MatrixXd x = Eigen::MatrixXd::Zero(horizon + 1, states);
            x.row(0) = xseed;
            for(int i = 1; i <= horizon; ++i) {
                ysim(j, i) = x.row(i-1).dot(w) + X.row(i).dot(kappa) + E(j, i);
                // Update state
                x.row(i) = (F * x.row(i-1).transpose()).array() + (g * E(j, i)).array();
            }
            xstates[j] = Rcpp::wrap(x);
        }
        return Rcpp::List::create(
            Rcpp::Named("states") = xstates,
            Rcpp::Named("simulated") = Rcpp::wrap(ysim)
        );
        
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        Rf_error("tsissm--> iss predict (constant) exception (unknown reason)");
    }
    return R_NilValue;
}



Eigen::VectorXd garchsimvec(Eigen::VectorXd epsilon, Eigen::VectorXd z, 
                            Eigen::VectorXd init_arch, Eigen::VectorXd init_garch, Eigen::VectorXd alpha,
                            const Eigen::VectorXd beta, const Eigen::VectorXi order,
                            const double variance_intercept, const int h) {
    const int maxpq = (int) order.maxCoeff();
    Eigen::VectorXd sigma_sim = Eigen::VectorXd::Zero(h + maxpq);
    Eigen::VectorXd sigma_sqr_sim = Eigen::VectorXd::Zero(h + maxpq);
    int i,j = 0;
    for(i=maxpq; i<(maxpq+h); i++) {
        sigma_sqr_sim(i) = variance_intercept;
        if (order(0) > 0) {
            for(j=0; j<order(0); j++) {
                if((order(0) + j) >= i) {
                    sigma_sqr_sim(i) += alpha(j) * init_arch(j);
                } else {
                    sigma_sqr_sim(i) += alpha(j) * (epsilon(i - j - 1) * epsilon(i - j - 1));
                }
            }
        }
        if (order(1) > 0) {
            for(j=0; j<order(1); j++) {
                if((order(1) + j) >= i) {
                    sigma_sqr_sim(i) += beta(j) * init_garch(j);
                } else {
                    sigma_sqr_sim(i) += beta(j) * sigma_sqr_sim(i - j - 1);
                }
            }
        }
        sigma_sim(i) = sqrt(sigma_sqr_sim(i));
        epsilon(i) = z(i) * sigma_sim(i);
    }
    return sigma_sim;
}

// [[Rcpp::export(.isspredict_dynamic)]]
Rcpp::List isspredict_dynamic(Eigen::Map<Eigen::MatrixXd>& F0, Eigen::Map<Eigen::MatrixXd>& F1,
                               Eigen::Map<Eigen::MatrixXd>& F2, Eigen::Map<Eigen::VectorXd>& w,
                               Eigen::Map<Eigen::VectorXd>& g, Eigen::Map<Eigen::MatrixXd>& X, 
                               Eigen::Map<Eigen::VectorXd>& kappa, Eigen::Map<Eigen::VectorXd>& xseed, 
                               Eigen::Map<Eigen::VectorXd>& e, Eigen::Map<Eigen::MatrixXd>& Z, 
                               Eigen::Map<Eigen::VectorXd>& init_arch, Eigen::Map<Eigen::VectorXd>& init_garch, 
                               Eigen::Map<Eigen::VectorXd>& alpha, 
                               Eigen::Map<Eigen::VectorXd>& beta, 
                               const double variance_intercept, Eigen::Map<Eigen::VectorXi>& order,
                               Rcpp::IntegerVector& mdim)
{
    
    try {
        const int states = mdim[0];
        const int nsim  = mdim[1];
        const int horizon  = mdim[2];
        const int maxpq = order.maxCoeff();
        Eigen::MatrixXd F = F0.cwiseProduct(F1).cwiseProduct(F2);
        // Initialize output structures
        Eigen::MatrixXd ysim = Eigen::MatrixXd::Zero(nsim, horizon + 1);
        Rcpp::List xstates(nsim);
        Eigen::MatrixXd sigma_sim = Eigen::MatrixXd::Zero(nsim, horizon + 1);
        Eigen::MatrixXd sigma_tmp = Eigen::VectorXd::Zero(horizon + maxpq);
        Eigen::MatrixXd E = Eigen::MatrixXd::Zero(nsim, horizon + 1);
        Eigen::VectorXd epsilon = Eigen::VectorXd::Zero(horizon + maxpq);
        for (int j = 0; j< maxpq; j++) {
            epsilon(j) = e(j);
        }
        // Main simulation loop

        for(int j = 0; j < nsim; ++j) {
            // simulate sigma
            Eigen::VectorXd sig = garchsimvec(epsilon, Z.row(j), init_arch, init_garch, alpha, beta, order, variance_intercept, horizon);
            sigma_sim.row(j) = sig.tail(horizon + 1);
            Eigen::MatrixXd x = Eigen::MatrixXd::Zero(horizon + 1, states);
            x.row(0) = xseed;
            Eigen::VectorXd tmp_z = Z.row(j).tail(horizon + 1);
            for(int i = 1; i <= horizon; ++i) {
                // \varepsilon_t = z_t * \sigma_t|t-1
                E(j,i) = tmp_z(i) * sigma_sim(j,i);
                ysim(j, i) = x.row(i-1).dot(w) + X.row(i).dot(kappa) + E(j, i);
                // Update state
                x.row(i) = (F * x.row(i-1).transpose()).array() + (g * E(j, i)).array();
            }
            xstates[j] = Rcpp::wrap(x);
        }
        return Rcpp::List::create(
            Rcpp::Named("states") = xstates,
            Rcpp::Named("simulated") = Rcpp::wrap(ysim),
            Rcpp::Named("sigma") = Rcpp::wrap(sigma_sim),
            Rcpp::Named("Error") = Rcpp::wrap(E)
        );
        
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        Rf_error("tsissm--> iss predict (dynamic) exception (unknown reason)");
    }
    return R_NilValue;
}


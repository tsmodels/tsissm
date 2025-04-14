/// @file issmgarch.hpp
#ifndef issmgarch_hpp
#define issmgarch_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type issmgarch(objective_function<Type>* obj) {
    DATA_VECTOR(V);
    DATA_MATRIX(X);
    DATA_VECTOR(good);
    DATA_IVECTOR(valid_index);
    DATA_VECTOR(y);
    DATA_VECTOR(allpars);
    DATA_IVECTOR(findex);
    DATA_IVECTOR(fpindex);
    DATA_IVECTOR(ppindex);
    DATA_IVECTOR(fshape);
    DATA_IVECTOR(modeli);
    DATA_VECTOR(xseed);
    DATA_UPDATE(xseed);
    DATA_STRING(initmethod);
    DATA_IVECTOR(vmodel);
    PARAMETER_VECTOR(issmpars);
    PARAMETER(lambda);
    PARAMETER_VECTOR(kappa);
    PARAMETER_VECTOR(eta);
    PARAMETER_VECTOR(delta);
    PARAMETER_VECTOR(distribution);
    
    
    // Slot the parameters into the vectors
    for (int i = 0; i < ppindex.size(); i++) {
        allpars(ppindex(i)) = issmpars(i);
    }
    for (int i = 0; i < findex.size(); i++) {
        V(findex(i)) = allpars(fpindex(i));
    }
    vector<Type> F0tmp = V.segment(fshape(0), fshape(1));
    matrix<Type> F0 = asMatrix(F0tmp, modeli(0), modeli(0));
    vector<Type> F1tmp = V.segment(fshape(2), fshape(3));
    matrix<Type> F1 = asMatrix(F1tmp, modeli(0), modeli(0));
    vector<Type> F2tmp = V.segment(fshape(4), fshape(5));
    matrix<Type> F2 = asMatrix(F2tmp, modeli(0), modeli(0));
    matrix<Type> F = (F0.array() * F1.array() * F2.array()).matrix();
    vector<Type> xreg = X * kappa;
    vector<Type> W = V.segment(fshape(6), fshape(7));
    vector<Type> G = V.segment(fshape(8), fshape(9));
    matrix<Type> M = asMatrix(W,modeli(0),1).transpose();
    matrix<Type> B = asMatrix(G,modeli(0),1) * M;
    matrix<Type> D = F.array() -  B.array();
    int timesteps = y.size();
    vector<Type> ytrans = issmextra::box_cox_transform(y, good, lambda);
    matrix<Type> states(timesteps + 1, modeli(0));
    states.setZero();
    states.row(0) = xseed.transpose();
    
    // Remaining computation...
    vector<Type> tmp2(modeli(0));
    vector<Type> gtmp2 = F * states.row(0).transpose();
    vector<Type> yhat(timesteps + 1);
    vector<Type> error(timesteps + 1);
    vector<Type> ylog(timesteps + 1);
    yhat.setZero();
    error.setZero();
    ylog.setZero();
    for (int i = 1; i <= timesteps; i++) {
        tmp2 = states.row(i - 1);
        yhat(i) = (tmp2.array() * W.array()).sum() + xreg(i - 1);
        if (good(i - 1) > 0.5) {
            error(i) = ytrans(i - 1) - yhat(i);
            ylog(i) = log(y(i - 1));
        } else {
            error(i) = 0.0;
        }
        gtmp2 = F * states.row(i - 1).transpose();
        states.row(i) = (gtmp2 + G * error(i));
    }
    vector<Type> e = error.segment(1, error.size() - 1);
    
    Type initial_variance = garchfun::init_variance(e(valid_index), initmethod, vmodel(3));
    vector<Type> initial_arch(vmodel(0));
    
    for(int i = 0; i<vmodel(0);i++) {
        initial_arch(i) = initial_variance;
    }
    Type persistence = eta.sum() + delta.sum();
    Type sample_variance = e(valid_index).square().mean();
    Type ptmp = sample_variance * (Type(1.0) - persistence);
    Type target_omega = CppAD::CondExpLe(ptmp, Type(0.0), Type(1e-12), ptmp);
    
    vector<Type> sigma_tmp = garchfun::garchrec(eta, delta, e, initial_variance, target_omega, vmodel, initmethod);
    vector<Type> log_y = y(valid_index).array().log();
    Type good_timesteps = good.sum();
    vector<Type> sigma = sigma_tmp.segment(vmodel(0), sigma_tmp.size() - vmodel(0));
    vector<Type> std_residuals = e(valid_index).array() * (Type(1.0)/sigma(valid_index).array());
    vector<Type> nll = distfun::distlike(std_residuals, distribution(0), distribution(1), modeli(6)) - sigma(valid_index).array().log();
    nll = Type(-1.0) * nll;
    nll = nll.array() - (lambda - Type(1.0)) * log_y.array();
    REPORT(D);
    REPORT(G);
    REPORT(W);
    REPORT(F);
    REPORT(error);
    REPORT(yhat);
    REPORT(std_residuals);
    REPORT(ytrans);
    REPORT(states);
    REPORT(sigma);
    REPORT(initial_arch);
    REPORT(nll);
    REPORT(log_y);
    REPORT(target_omega);
    REPORT(persistence);
    ADREPORT(nll);
    Type loglik = (good_timesteps / timesteps) * nll.sum();
    return loglik;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
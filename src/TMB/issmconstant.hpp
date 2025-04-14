/// @file issmconstant.hpp
#ifndef issmconstant_hpp
#define issmconstant_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type issmconstant(objective_function<Type>* obj) {
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
    PARAMETER_VECTOR(issmpars);
    PARAMETER(lambda);
    PARAMETER_VECTOR(kappa);
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
    // Call the new function
    // matrix<Type> init_states = issmextra::initialize_states(ytrans, xreg, good, valid_index, F, W, G, D, timesteps, modeli);
    // Continue with the rest of the original function...
    matrix<Type> states(timesteps + 1, modeli(0));
    states.setZero();
    states.row(0) = xseed.transpose();
    
    // Remaining computation...
    vector<Type> tmp2(modeli(0));
    vector<Type> gtmp2 = F * states.row(0).transpose();
    vector<Type> yhat(timesteps + 1);
    vector<Type> error(timesteps + 1);
    yhat.setZero();
    error.setZero();
    for (int i = 1; i <= timesteps; i++) {
        tmp2 = states.row(i - 1);
        yhat(i) = (tmp2.array() * W.array()).sum() + xreg(i - 1);
        if (good(i - 1) > 0.5) {
            error(i) = ytrans(i - 1) - yhat(i);
        } else {
            error(i) = 0.0;
        }
        gtmp2 = F * states.row(i - 1).transpose();
        states.row(i) = (gtmp2 + G * error(i));
    }
    
    vector<Type> e = error.segment(1, error.size() - 1);
    vector<Type> log_y = y(valid_index).array().log();
    Type sigma = issmextra::sigma_calc(e(valid_index));
    Type good_timesteps = good.sum();
    vector<Type> std_residuals = e(valid_index).array() * (Type(1.0)/sigma);
    vector<Type> xnll = distfun::distlike(std_residuals, distribution(0), distribution(1), modeli(6)) - log(sigma);
    vector<Type> nll = Type(-1.0) * xnll - (lambda - Type(1.0)) * log_y;
    REPORT(D);
    REPORT(G);
    REPORT(W);
    REPORT(F);
    REPORT(error);
    REPORT(sigma);
    REPORT(yhat);
    REPORT(std_residuals);
    REPORT(states);
    REPORT(log_y);
    REPORT(good_timesteps);
    REPORT(ytrans);
    REPORT(nll);
    ADREPORT(nll);
    Type loglik = Type(good_timesteps / timesteps) * nll.sum();
    return loglik;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
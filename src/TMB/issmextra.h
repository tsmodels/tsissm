namespace issmextra {
     template <class Type>
     vector<Type> box_cox_transform(vector<Type> y, vector<Type> good, Type lambda) {
         const int n = y.size();
         vector<Type> ytrans(n);
         for (int i = 0; i < n; i++) {
             if (good(i) > 0.5) {
                ytrans(i) = CppAD::CondExpLe(lambda, Type(1e-12), log(y(i)), (pow(y(i), lambda) - Type(1.0)) / lambda);
             } else {
                 ytrans(i) = Type(0.0);
             }
         }
         return ytrans;
     }
     
     
     template <class Type>
     vector<Type> seed_transform(vector<Type> x, Type lambda) {
         const int n = x.size();
         vector<Type> xt(n);
         for (int i = 0; i < n; i++) {
             xt(i) = CppAD::CondExpLe(lambda, Type(1e-12), log(x(i)), (pow(x(i), lambda) - Type(1.0)) / lambda);
         }
         return xt;
     }
    
    template <class Type>
    Type sigma_calc(vector<Type> error) {
        Type sigma = sqrt(error.array().square().mean());
        return sigma;
    }
    
    
    template<class Type>
    matrix<Type> initialize_states(vector<Type> ytrans, vector<Type> xreg, vector<Type> good, 
                                   vector<int> valid_index, matrix<Type> F, vector<Type> W, vector<Type> G, 
                                   matrix<Type> D, int timesteps, vector<int> modeli) {
        // Initialize necessary matrices and vectors
        matrix<Type> xaux(timesteps, modeli(0));
        matrix<Type> waux(timesteps, modeli(0));
        vector<Type> eaux(timesteps);
        vector<Type> yaux(timesteps);
        xaux.setZero();
        waux.setZero();
        eaux.setZero();
        yaux.setZero();
        yaux(0) = xreg(0);
        waux.row(0) = W;
        eaux(0) = ytrans(0) - yaux(0);
        xaux.row(0) = G * eaux(0);
        vector<Type> tmp = xaux.row(0);
        vector<Type> gtmp = F * xaux.row(0).transpose();
        matrix<Type> init_states(1, modeli(3) + modeli(4));
        init_states.setZero();
        for (int i = 1; i < timesteps; i++) {
            tmp = xaux.row(i - 1);
            yaux(i) = (tmp.array() * W.array()).sum() + xreg(i);
            if (good(i) > 0.5) {
                eaux(i) = ytrans(i) - yaux(i);
            } else {
                eaux(i) = eaux(i - 1);
            }
            gtmp = F * xaux.row(i-1).transpose();
            xaux.row(i) = (gtmp + G * eaux(i));
            waux.row(i) = waux.row(i - 1) * D;
        }
        matrix<Type> B = waux(valid_index, Eigen::all);
        matrix<Type> A = B.leftCols(modeli(3));
        vector<Type> E = eaux(valid_index);
        matrix<Type> xseed = A.householderQr().solve(E.matrix());
        init_states.leftCols(modeli(3)) = xseed.transpose();
        return init_states;
    }
}
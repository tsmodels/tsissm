namespace garchfun{

    template <class Type>
    Type init_variance(vector<Type> x, std::string method, int samplen)
    {
        Type initv = 0.0;
        if (method == "unconditional") {
            initv = x.abs().pow(2.0).mean();
        } else  {
            initv = x.head(samplen).abs().pow(2.0).mean();
        }
        return initv;
    }
    
    template <class Type>
    vector<Type> garchrec(vector<Type> alpha, vector<Type> beta, vector<Type> error, Type initial_variance, Type target_omega, vector<int> vmodel, std::string initmethod) {
        const int timesteps = error.size() + vmodel(0);
        vector<Type> residuals(timesteps);
        vector<Type> sigma_squared(timesteps);
        sigma_squared.setZero();
        residuals.setZero();
        int start = vmodel(0);
        residuals.segment(start, residuals.size() - start) = error;
        vector<Type> residuals_squared = residuals.array().square();
        int j = 0;
        for(j = 0;j<vmodel(0);j++) {
            sigma_squared(j) += initial_variance;
            residuals_squared(j) = initial_variance;
        }
        
        for(int i = vmodel(0);i<timesteps;i++){
            sigma_squared(i) += target_omega;
            for(j = 0;j<vmodel(1);j++){
                sigma_squared(i) += alpha(j) * residuals_squared(i - j - 1);
            }
            for(j = 0;j<vmodel(2);j++){
                sigma_squared(i) += beta(j) * sigma_squared(i - j - 1);
            }
        }
        vector<Type> sigma = sigma_squared.sqrt();
        return sigma;
    }
}

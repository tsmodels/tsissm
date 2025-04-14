namespace distfun{
    template <class Type>
    Type mygammafn(Type x)
    {
        Type out = exp(lgamma(x));
        return out;
    }
    VECTORIZE1_t(mygammafn)
    
    template <class Type>
    Type norm_std(Type x, int give_log)
    {
        Type pdf;
        pdf = dnorm(x, Type(0.0), Type(1.0), give_log);
        return pdf;
    }
    VECTORIZE1_t(norm_std)
                    
    template <class Type>
    Type student_std(Type x, Type shape, int give_log)
    {
        Type pdf;
        Type s = sqrt(shape/(shape - Type(2.0)));
        pdf = dt(x * s, shape, 0) * s;
        if (give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE3_tti(student_std)
                        
                        
    template <class Type>
    Type jsu_std(Type x, Type skew, Type shape, int give_log)
    {
        Type rtau = Type(1.0)/shape;
        Type w = CppAD::CondExpLt(rtau, Type(0.0000001), Type(1.0), exp(rtau*rtau));
        Type omega = -skew * rtau;
        Type c = sqrt(Type(1.0)/(Type(0.5)*(w - Type(1.0))*(w*cosh(Type(2.0)*omega)+Type(1.0))));
        Type z = (x - (c * sqrt(w) * sinh(omega)))/c;
        Type r = -skew + log(z + sqrt(z*z+Type(1.0)))/rtau;
        Type pdf = -log(c) - log(rtau) - Type(0.5) * log( z * z + Type(1.0)) - Type(0.5) * log(Type(2.0) * M_PI) - Type(0.5) * r * r;
        if(give_log == 0) {
            pdf = exp(pdf);
        }
        return pdf;
    }
    VECTORIZE4_ttti(jsu_std)
                            
    template <class Type>
    Type sstd_std(Type x, Type skew, Type shape, int give_log)
    {
        Type a = Type(1.0)/Type(2.0);
        Type b = shape/Type(2.0);
        Type beta = exp((lgamma(a) - lgamma(a+b)) + lgamma(b));
        Type m1 = Type(2.0) * sqrt(shape - Type(2.0))/(shape - Type(1.0))/beta;
        Type mu = m1 * (skew - Type(1.0)/skew);
        Type sigma = sqrt((Type(1.0) - m1 * m1) * (skew * skew + Type(1.0)/(skew * skew)) + Type(2.0) * m1 * m1 - Type(1.0));
        Type z = x * sigma + mu;
        Type xxi_tmp = CppAD::CondExpLt(z, Type(0.0), Type(1.0)/skew, skew);
        Type xxi = CppAD::CondExpEq(z, Type(0.0), Type(1.0), xxi_tmp);
        Type g = Type(2.0)/(skew+Type(1.0)/skew);
        Type pdf = g * student_std(z/xxi,shape,0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE4_ttti(sstd_std)
                                
    template <class Type>
    Type ged_std(Type x, Type shape, int give_log)
    {
        Type lambda = sqrt(pow(Type(1.0)/Type(2.0),Type(2.0)/shape)*mygammafn(Type(1.0)/shape)/mygammafn(Type(3.0)/shape));
        Type g = shape/(lambda*(pow(Type(2.0),Type(1.0)+(Type(1.0)/shape)))*mygammafn(Type(1.0)/shape));
        Type pdf = g * exp(Type(-0.5)*pow(fabs(x/lambda),shape));
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE3_tti(ged_std)
                                    
    template <class Type>
    Type sged_std(Type x, Type skew, Type shape, int give_log)
    {
        Type lambda = sqrt(pow(Type(1.0)/Type(2.0),(Type(2)/shape))*mygammafn(Type(1.0)/shape)/mygammafn(Type(3.0)/shape));
        Type m1 = pow(Type(2.0), (Type(1.0)/shape))*lambda*mygammafn(Type(2.0)/shape)/mygammafn(Type(1.0)/shape);
        Type mu = m1*(skew-Type(1.0)/skew);
        Type skew2 = skew * skew;
        Type m12 = m1 * m1;
        Type sigma = sqrt((Type(1.0) - m12)*(skew2 + Type(1.0)/skew2) + Type(2.0) * m12 - Type(1.0));
        Type z = x*sigma + mu;
        Type xxi_tmp = CppAD::CondExpLt(z, Type(0.0), Type(1.0)/skew, skew);
        Type xxi = CppAD::CondExpEq(z, Type(0.0), Type(1.0), xxi_tmp);
        Type g = Type(2.0)/(skew + Type(1.0)/skew);
        Type pdf = g * ged_std(z/xxi, shape, 0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE4_ttti(sged_std)
                                        
    template <class Type>
    Type snorm_std(Type x, Type skew, int give_log)
    {
        Type m1 = Type(2.0)/sqrt(Type(2.0)*M_PI);
        Type m12 = m1 * m1;
        Type xi2 = skew*skew;
        Type mu = m1*(skew-Type(1.0)/skew);
        Type sigma = sqrt((Type(1.0)-m12) * (xi2 + Type(1.0)/xi2) + Type(2.0) * m12 - Type(1.0));
        Type z = x*sigma+mu;
        Type xxi = CppAD::CondExpLt(z, Type(0.0), Type(1.0)/skew, skew);
        Type g = Type(2.0)/(skew + Type(1.0)/skew);
        Type pdf = g * dnorm(z/xxi,Type(0), Type(1.0), 0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE3_tti(snorm_std)

    template <class Type>
    Type distlike(Type x, Type skew, Type shape, int dclass)
    {
        Type out = 0.0;
        const int give_log = 1;
        switch(dclass){
        case 1:
            out = norm_std(x,give_log);
            break;
        case 2:
            out = student_std(x, shape, give_log);
            break;
        case 3:
            out = snorm_std(x, skew, give_log);
            break;
        case 4:
            out = sstd_std(x, skew, shape, give_log);
            break;
        case 5:
            out = ged_std(x, shape, give_log);
            break;
        case 6:
            out = sged_std(x, skew, shape, give_log);
            break;
        case 7:
            out = jsu_std(x, skew, shape, give_log);
            break;
        default:
            out = norm_std(x,give_log);
        }
        return(out);
    }
    VECTORIZE4_ttti(distlike)
}

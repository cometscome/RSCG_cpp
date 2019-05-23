#include "rscg.h"
#include <algorithm>

RSCG::RSCG(int n,std::function<vector <complex <double> >(vector <complex <double> >&)> f){
        matvec = f;
        matvec_dble = nullptr;
        _n = n;
        _maximumsteps = 100;
        _eps = 1e-12;
    }

RSCG::RSCG(int n,std::function<vector <double>(vector <double>&)> f){
        matvec = nullptr;
        matvec_dble = f;
        _n = n;
        _maximumsteps = 100;
        _eps = 1e-12;
    }    

complex <double> dot_product(vector <complex <double> > &a,vector <complex <double> > &b){
    int n = a.size();
    complex <double> csum=0.0;
    for(int i = 0; i < n; i++)
    {
        csum += conj(a[i])*b[i];
    }
    return csum;
};

double dot_product(vector <double> &a,vector <double> &b){
    int n = a.size();
    double csum=0.0;
    for(int i = 0; i < n; i++)
    {
        csum += a[i]*b[i];
    }
    return csum;
};

void add_vec(complex <double> alpha,vector <complex <double> > &a,
    complex <double> beta,vector <complex <double> > &b){ //alpha*a + beta*b
    int n = a.size();
    for(int i = 0; i < n; i++)
    {
        a[i] = alpha*a[i];
        a[i] += beta*b[i];
    }
};

void add_vec(double alpha,vector <double> &a,
    double beta,vector <double> &b){ //alpha*a + beta*b
    int n = a.size();
    for(int i = 0; i < n; i++)
    {
        a[i] = alpha*a[i];
        a[i] += beta*b[i];
    }
};

double find_max(vector <complex <double> > &a){
    double amax;
    amax = abs(a[0]);
    if (a.size() > 1){
        for(int k = 1; k < a.size(); k++)
        {
            if(amax < abs(a[k])){
                amax = abs(a[k]);
            } 
        }
    }
    return amax;
};



vector <complex <double> > RSCG::calc_greenfunction_complexH(int i,int j,
    vector <complex <double> > vec_z //z_j
    ){
    int nsigma = vec_z.size();
    using Vector = vector<complex <double> >;
    using Vector_z = vector<complex <double> >;
    using Scalar= complex <double>;
    using Scalar_z = complex <double>;
    Vector_z Gij(nsigma,0.0);

    Vector b(_n,0.0);
    b[j] = 1.0;

    Vector x(_n,0.0);
    Vector r(b);
    Vector p(b);
    Vector Ap(_n,0.0);

    Scalar  alpham = 1.0;
    Scalar  betam = 1.0;


    Scalar Sigma = b[i];
    Vector_z theta(nsigma,0.0);
    Vector_z Pir(nsigma,Sigma);

    Vector_z rhok(nsigma,1.0);
    Vector_z rhokm(nsigma,1.0);
    Vector_z rhokp(nsigma,1.0);

    Scalar rnorm =0.0;
    Scalar alpha =0.0;
    Scalar beta =0.0;


    for(int k = 0; k < _maximumsteps+1; k++)
    {
        Ap = matvec(p);
        rnorm = dot_product(r,r);
        alpha = -rnorm/dot_product(p,Ap);
        add_vec(1.0,x,alpha,p); //x = x+ alpha*p
        add_vec(1.0,r,alpha,Ap); //r = r+alpha*Ap
        if (rnorm == 0.0){
            cout << k << " rnorm  " << rnorm  << std::endl;
            return theta;
        };


        beta = dot_product(r,r)/rnorm;
        add_vec(beta,p,1.0,r); //p = beta*p+r
        Sigma = r[i];

        for(int j = 0; j < nsigma; j++)
        {
            //cout << k << " " << j << " rhok " << rhok[j] << std::endl;
            
            if(abs(rhok[j]) > _eps){
                rhokp[j] = rhok[j]*rhokm[j]*alpham/(rhokm[j]*alpham*(1.0+alpha*vec_z[j])
                    +alpha*betam*(rhokm[j]-rhok[j])); //#Line 13
                Scalar_z alphakj = alpha*rhokp[j]/rhok[j]; //#Line 14
                theta[j] += alphakj*Pir[j]; //#Line 15
                Scalar_z betakj = ((rhokp[j]/rhok[j])*(rhokp[j]/rhok[j]))*beta; // #Line 16
                Pir[j] = rhokp[j]*Sigma+ betakj*Pir[j]; //#Line 17
                rhokm[j] = rhok[j];
                rhok[j] = rhokp[j];
            }
        }
        alpham = alpha;
        betam = beta;
        double hi = real(rnorm)*find_max(rhok);
        //cout << k << " " << hi << " " << rnorm << std::endl;
        
        if (hi < _eps && k > 0){
            cout << k << " " << hi << std::endl;
            return theta;
        };
    }

    return theta;

    };

vector <complex <double>> RSCG::calc_greenfunction_realH(int i,int j,
    vector <complex < double> > vec_z //z_j
    ){
    int nsigma = vec_z.size();
    using Vector = vector<double>;
    using Vector_z = vector<complex <double> >;
    using Scalar= double;
    using Scalar_z = complex <double>;
    Vector_z Gij(nsigma,0.0);

    Vector b(_n,0.0);
    b[j] = 1.0;

    Vector x(_n,0.0);
    Vector r(b);
    Vector p(b);
    Vector Ap(_n,0.0);

    Scalar  alpham = 1.0;
    Scalar  betam = 1.0;


    Scalar Sigma = b[i];
    Vector_z theta(nsigma,0.0);
    Vector_z Pir(nsigma,Sigma);

    Vector_z rhok(nsigma,1.0);
    Vector_z rhokm(nsigma,1.0);
    Vector_z rhokp(nsigma,1.0);

    Scalar rnorm =0.0;
    Scalar alpha =0.0;
    Scalar beta =0.0;


    for(int k = 0; k < _maximumsteps+1; k++)
    {
        Ap = matvec_dble(p);
        rnorm = dot_product(r,r);
        alpha = -rnorm/dot_product(p,Ap);
        add_vec(1.0,x,alpha,p); //x = x+ alpha*p
        add_vec(1.0,r,alpha,Ap); //r = r+alpha*Ap
        if (rnorm == 0.0){
            cout << k << " rnorm  " << rnorm  << std::endl;
            return theta;
        };


        beta = dot_product(r,r)/rnorm;
        add_vec(beta,p,1.0,r); //p = beta*p+r
        Sigma = r[i];

        for(int j = 0; j < nsigma; j++)
        {
            //cout << k << " " << j << " rhok " << rhok[j] << std::endl;
            
            if(abs(rhok[j]) > _eps){
                rhokp[j] = rhok[j]*rhokm[j]*alpham/(rhokm[j]*alpham*(1.0+alpha*vec_z[j])
                    +alpha*betam*(rhokm[j]-rhok[j])); //#Line 13
                Scalar_z alphakj = alpha*rhokp[j]/rhok[j]; //#Line 14
                theta[j] += alphakj*Pir[j]; //#Line 15
                Scalar_z betakj = ((rhokp[j]/rhok[j])*(rhokp[j]/rhok[j]))*beta; // #Line 16
                Pir[j] = rhokp[j]*Sigma+ betakj*Pir[j]; //#Line 17
                rhokm[j] = rhok[j];
                rhok[j] = rhokp[j];
            }
        }
        alpham = alpha;
        betam = beta;
        double hi = rnorm*find_max(rhok);
        //cout << k << " " << hi << " " << rnorm << std::endl;
        
        if (hi < _eps && k > 0){
            cout << k << " " << hi << std::endl;
            return theta;
        };
    }

    return theta;

    };

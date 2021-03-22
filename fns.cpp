#include "fns.hpp"
#include <math.h>
// #include "tensor.h"

// --------------------------------------------------------------
//                       MISC
// --------------------------------------------------------------

double qij(tensor2<double> xi,tensor2<double> xj,const double& h){
    return (xi-xj).norm()/h;
}

void scale_m(const int& lenX, double* Rho, const double& rho0, double& m){
    double sum = 0.0;
    for (int idx = 0; idx< lenX; idx++){
        sum+=Rho[idx];
    }
    m = sqrt(rho0*double(lenX)/sum);


}

// --------------------------------------------------------------
//                       PHI DERIVATIVES 
// --------------------------------------------------------------

double phi_dij(tensor2<double> xi,tensor2<double> xj,const double& h){
    double q = qij(xi,xj,h);
    if (q >=1.0){
        return 0.0;
    }
    return 4.0*pow(1.0-q*q,3)/(M_PI*h*h);
}

double grad2_phi_vij(tensor2<double> xi,tensor2<double> xj,const double& h){
    double q = qij(xi,xj,h);
    if (q >=1.0){
        return 0.0;
    }
    return 40.0*(1-q)/(M_PI*pow(h,4));
}

tensor2<double> grad_phi_pij(tensor2<double> xi, tensor2<double> xj,const double& h){
    tensor2<double> rij = xi-xj;
    double q = qij(xi,xj,h);
    if (q >=1.0){
        return rij*0.0;
    }
    return rij*(-30.0*pow(1.0-q,2)/(M_PI*h*h*h*q));
}

// --------------------------------------------------------------
//                       Pressure and Density
// --------------------------------------------------------------

double rho_i(const int& i, const int& lenX,  tensor2<double>* X, const double& m, const double& h){
    double sum = 0.0;
    for(int j = 0; j<lenX; j++){
        sum+=phi_dij(X[i],X[j],h);
    }
    return m*sum;

}



double P_i(const int& i, const int& lenRho,double* Rho, const double& k, const double& rho0){
    return (Rho[i]-rho0)*k;
}

tensor2<double> Fp_i(const int& i, const int& lenX,  tensor2<double>* X, double* P, double* Rho,const double& m, const double& h){
    tensor2<double> sum(0.0,0.0);
    for (int j =0;j<lenX;j++){
        if(j==i){
            continue;
        }
        sum-=grad_phi_pij(X[i],X[j],h)*(m/Rho[j]/2)*(P[i]+P[j]);
    }
    return sum;
}

// --------------------------------------------------------------
//                       FORCES
// --------------------------------------------------------------

tensor2<double> Fv_i(const int& i, const int& lenX,  tensor2<double>* X, tensor2<double>* V, double* Rho,const double& m, const double& h, const double& mu){
    tensor2<double> sum(0.0,0.0);
    tensor2<double> vij;
    for (int j =0;j<lenX;j++){
        if(j==i){
            continue;
        }
        vij = V[i]-V[j];
        sum-=vij*grad2_phi_vij(X[i],X[j],h) * (m/Rho[j]/2);
    }
    return sum*mu;

}


tensor2<double> Fg_i(const int& i, double* Rho,const double& g){
    tensor2<double> Force(0.0,-Rho[i]*g);
    return Force;
}






// --------------------------------------------------------------
//                       FILL FUNCS
// --------------------------------------------------------------

void Fill_P(const int& lenX,  tensor2<double>* X, double* Rho, double* P, const double& k, const double& rho0){
    for(int idx = 0; idx< lenX; idx++){
        P[idx] = P_i(idx,lenX,Rho,k,rho0);
    }

}

void Fill_Rho(const int& lenX,  tensor2<double>* X, double* Rho, const double& m, const double& h){
    for(int idx = 0; idx< lenX; idx++){
        Rho[idx] = rho_i(idx,lenX,X,m,h);
    }

}

void Fill_Fp(const int& lenX,  tensor2<double>* X, double* P, double* Rho,tensor2<double>* Fp,const double& m, const double& h){
    for(int idx = 0; idx< lenX; idx++){
        Fp[idx] = Fp_i(idx,lenX,X,P,Rho,m,h);
    }

}

void Fill_Fv(const int& lenX,  tensor2<double>* X, tensor2<double>* V, double* Rho,tensor2<double>* Fv,const double& m, const double& h, const double& mu){
    for(int idx = 0; idx< lenX; idx++){
        Fv[idx] = Fv_i(idx,lenX,X,V,Rho,m,h,mu);
    }

}

void Fill_Fg(const int& lenX,double* Rho,tensor2<double>* Fg,const double& g){
    for(int idx = 0; idx< lenX; idx++){
        Fg[idx] = Fg_i(idx,Rho,g);
    }

}

void Fill_A(const int& lenX, tensor2<double>* Fp,tensor2<double>* Fv,tensor2<double>* Fg,double* Rho, tensor2<double>* A){
    for(int idx = 0; idx< lenX; idx++){
        A[idx] = (Fp[idx]+Fv[idx]+Fg[idx])/Rho[idx];
    }
}


// --------------------------------------------------------------
//                       TIMESTEPS
// --------------------------------------------------------------
 
void Tstep_init(const int& lenX, tensor2<double>* A, tensor2<double>* X, tensor2<double>* V, const double& tstep){
    for (int idx = 0; idx < lenX; idx++){
        V[idx] +=A[idx]*tstep/2.0;
        X[idx] += V[idx]*tstep;
    }
}


void Tstep(const int& lenX, tensor2<double>* A, tensor2<double>* X, tensor2<double>* V, const double& tstep){
    for (int idx = 0; idx < lenX; idx++){
        V[idx] +=A[idx]*tstep;
        X[idx] += V[idx]*tstep;
    }
}

void EnforceBC(const int& lenX, tensor2<double>* X, tensor2<double>* V, const double& e, const double& h){
     for (int idx = 0; idx < lenX; idx++){

        //BOUNDARY CONDITIONS X-axis
        if(X[idx].x1() < h){
            X[idx].setx1(h);
            V[idx].setx1(-e*V[idx].x1());
        }else if(X[idx].x1() >(1-h)){
            X[idx].setx1(1-h);
            V[idx].setx1(-e*V[idx].x1());

        }

        //BOUNDARY CONDITIONS Y-axis
        if(X[idx].x2() < h){
            X[idx].setx2(h);
            V[idx].setx2(-e*V[idx].x2());
        }else if(X[idx].x2() >(1-h)){
            X[idx].setx2(1-h);
            V[idx].setx2(-e*V[idx].x2());

        }
    }
    
}
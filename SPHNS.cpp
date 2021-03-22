
// -------------------------
//         HEADERS
// -------------------------

#include "fns.hpp"
#include "tensor.cpp"



#include <iostream>
#include <math.h>
#include <iomanip>  

// double tt(tensor2<double>& xj){
//     return xj.norm();
// }


int main(){

    std::cout<<"\n\n\n\n\n";


    const double h = sqrt(2)*0.5;
    double m = 1.0;
    const int len = 6;
    const double k = 1.0;
    const double rho0 = 1000.0;
    const double g = 9.81;
    const double mu = 1.0;


    
    tensor2<double> p0(0.0,0.0);
    tensor2<double> p1(0.5,0.0);
    tensor2<double> p2(1.0,0.0);
    tensor2<double> p3(1.0,0.5);
    tensor2<double> p4(0.5,0.5);
    tensor2<double> p5(0.0,0.5);

    tensor2<double> v0(1.0,1.0);
    tensor2<double> v1(0.0,1.0);
    tensor2<double> v2(-1.0,1.0);
    tensor2<double> v3(-1.0,0.0);
    tensor2<double> v4(0.0,0.0);
    tensor2<double> v5(1.0,0.0);







    tensor2<double>* X = new tensor2<double>[len];
    tensor2<double>* V = new tensor2<double>[len];
    double* P = new double[len];
    double* Rho = new double[len];
    tensor2<double>* Fp = new tensor2<double>[len];
    tensor2<double>* Fv = new tensor2<double>[len];
    tensor2<double>* Fg = new tensor2<double>[len];
    tensor2<double>* A = new tensor2<double>[len];

    X[0] = p0;
    X[1] = p1;
    X[2] = p2;
    X[3] = p3;
    X[4] = p4;
    X[5] = p5;

    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    V[3] = v3;
    V[4] = v4;
    V[5] = v5;










    
    Fill_Rho(len,X,Rho,m,h);
    scale_m(len,Rho,rho0,m);
    std::cout<<"scaled mass: "<<m<<std::endl;

    Fill_P(len,X,Rho,P,k,rho0);
    Fill_Fp(len,X,P,Rho,Fp,m,h);
    Fill_Fv(len,X,V,Rho,Fv,m,h,mu);
    Fill_Fg(len,Rho,Fg,g);
    Fill_A(len,Fp,Fv,Fg,Rho,A);


    std::cout<<"VELOCITIES:\n";
    print_arr(len,V);
    std::cout<<"POSITIONS:\n";
    print_arr(len,X);
    std::cout<<"DENSITIES:\n";
    print_arr(len,Rho);
    std::cout<<"PRESSURES:\n";
    print_arr(len,P);
    std::cout<<"PRESSURE FORCE:\n";
    print_arr(len,Fp);
    std::cout<<"VISCOUS FORCE:\n";
    print_arr(len,Fv);
    std::cout<<"GRAVITATIONAL FORCE:\n";
    print_arr(len,Fg);
   
    



    delete[] X;
    delete[] P;
    delete[] Rho;
    delete[] Fp;
    delete[] Fv;
    delete[] Fg;
    

}
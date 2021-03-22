
// -------------------------
//         HEADERS
// -------------------------

#include "fns.hpp"
#include "tensor.cpp"



#include <iostream>
#include <math.h>
#include <iomanip>  

//-------------------------
//       CONSTANTS
// ------------------------
const double h = 0.01;
double m = 1.0;
const int len = 6;
const double k = 2000.0;
const double rho0 = 1000.0;
const double g = 9.81;
const double mu = 1.0;
const double step = pow(10,-4);
const double e = 0.5;
const double sim_time = 100.00; //simulation time in seconds
const int ntimesteps = int((sim_time+1)/step);



/**
 * @brief Runs SPH simulation
 * 
 * @param Rho 
 * @param P 
 * @param X 
 * @param V 
 * @param Fp 
 * @param Fv 
 * @param Fg 
 * @param A 
 */
void Simulate(double* Rho, double* P, tensor2<double>* X,
                 tensor2<double>* V,tensor2<double>* Fp,
                 tensor2<double>* Fv,tensor2<double>* Fg, tensor2<double>* A){

    //-----------------------
    //     INITIAL CONDITIONS
    //-----------------------


    Fill_Rho(len,X,Rho,m,h);
    scale_m(len,Rho,rho0,m);
    std::cout<<"scaled mass: "<<m<<std::endl;
    Fill_P(len,X,Rho,P,k,rho0);
    Fill_Fp(len,X,P,Rho,Fp,m,h);
    Fill_Fv(len,X,V,Rho,Fv,m,h,mu);
    Fill_Fg(len,Rho,Fg,g);
    Fill_A(len,Fp,Fv,Fg,Rho,A);

    std::cout<<"\nINITITIAL CONDITIONS\n";
    std::cout<<"VELOCITIES:\n";
    print_arr(len,V);
    std::cout<<"POSITIONS:\n";
    print_arr(len,X);



    //-----------------------
    //       STEP 1
    //-----------------------
    std::cout<<"\nINIT TSTEP\n";
    Tintegrate_init(len,A,X,V,step);
    EnforceBC(len,X,V,e,h);
    std::cout<<"VELOCITIES:\n";
    print_arr(len,V);
    std::cout<<"POSITIONS:\n";
    print_arr(len,X);



    //-----------------------
    //       STEP 2
    //-----------------------

    for(int t = 1; t<ntimesteps;t++){
        Fill_Rho(len,X,Rho,m,h);
        Fill_P(len,X,Rho,P,k,rho0);
        Fill_Fp(len,X,P,Rho,Fp,m,h);
        Fill_Fv(len,X,V,Rho,Fv,m,h,mu);
        Fill_Fg(len,Rho,Fg,g);
        Fill_A(len,Fp,Fv,Fg,Rho,A);


        
        Tintegrate(len,A,X,V,step);
        EnforceBC(len,X,V,e,h);
        if(t%10000 == 0){
            std::cout<<"\n\nt =  "<<t/10000<<"s"<<std::endl;
            std::cout<<"VELOCITIES:\n";
            print_arr(len,V);
            std::cout<<"POSITIONS:\n";
            print_arr(len,X);
            std::cout<<"NET SCALED FORCE:\n";
            print_arr(len,A);
        }


        //------------------
        //  PRINT INTERACTIONS
        //-------------------

        // for(int p = 0; p<len; p++){
        //     if(A[p].x1() != 0){
        //     std::cout<<"\n\n***INTERACTION***:\n";
        //     std::cout<<"t =  "<<t/10000<<"s"<<std::endl;
        //     std::cout<<"VELOCITIES:\n";
        //     print_arr(len,V);
        //     std::cout<<"POSITIONS:\n";
        //     print_arr(len,X);
        //     std::cout<<"NET SCALED FORCE:\n";
        //     print_arr(len,A);
        //     break;

        //     }
        // }
        

    }

}

int main(){

    std::cout<<"\n\n\n\n\n";

    //-----------------------
    //     MEMORY ALLOC
    //-----------------------
    tensor2<double> p0(0.0,0.0);
    tensor2<double> p1(0.5,0.0);
    tensor2<double> p2(1.0,0.0);
    tensor2<double> p3(0.99,0.5);
    tensor2<double> p4(0.4999,0.5);
    tensor2<double> p5(0.000001,0.5);

    tensor2<double> v0(0.0,0.0);
    tensor2<double> v1(0.0,0.0);
    tensor2<double> v2(0.0,0.0);
    tensor2<double> v3(0.0,0.0);
    tensor2<double> v4(0.0,0.0);
    tensor2<double> v5(0.0,0.0);



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


    //-----------------------
    //     SIM
    //-----------------------

    Simulate(Rho,P,X,V,Fp,Fv,Fg,A);

    // ----------------
    //    TIDY UP
    // ----------------
    delete[] X;
    delete[] P;
    delete[] Rho;
    delete[] Fp;
    delete[] Fv;
    delete[] Fg;
    delete[] A;
    

}
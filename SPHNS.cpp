
// -------------------------
//         HEADERS
// -------------------------
#include <mpi.h>
#include "fns.hpp"
#include "tensor.cpp"
#include<tuple>



#include <iostream>
#include <math.h>
#include <iomanip>  

// --------------------
//    GLOBAL VARS
// --------------------
auto communicator = MPI_COMM_WORLD; // All processes
//---------------------

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
                 tensor2<double>* Fv,tensor2<double>* Fg, tensor2<double>* A, const bool& verbose = false){

    //-----------------------
    //     INITIAL CONDITIONS
    //-----------------------

    std::string delim = ",";//output file deliminator


    Fill_Rho(len,X,Rho,m,h);
    scale_m(len,Rho,rho0,m);
    std::cout<<"scaled mass: "<<m<<std::endl;
    Fill_P(len,X,Rho,P,k,rho0);
    Fill_Fp(len,X,P,Rho,Fp,m,h);
    Fill_Fv(len,X,V,Rho,Fv,m,h,mu);
    Fill_Fg(len,Rho,Fg,g);
    Fill_A(len,Fp,Fv,Fg,Rho,A);

    if(verbose){
        std::cout<<"\nINITITIAL CONDITIONS\n";
        std::cout<<"VELOCITIES:\n";
        print_arr(len,V);
        std::cout<<"POSITIONS:\n";
        print_arr(len,X);
    }

    //SAVE DATA AND ADD HEAD
    save_position(len,X,0.0,delim,true,false);
    


    //-----------------------
    //       STEP 1
    //-----------------------
    
    Tintegrate_init(len,A,X,V,step);
    EnforceBC(len,X,V,e,h);
    if(verbose){
        std::cout<<"\nINIT TSTEP\n";
        std::cout<<"VELOCITIES:\n";
        print_arr(len,V);
        std::cout<<"POSITIONS:\n";
        print_arr(len,X);
        std::cout<<"POTENTIAL ENERGY:\n";
        std::cout<<Ep(len,X,m,g)<<std::endl;
    }

    // SAVE ENERGY AND ADD HEAD
    save_energy(Ek(len,V,m),Ep(len,X,m,g),step,delim,true,false);

    

    



    //-----------------------
    //       STEP 2
    //-----------------------

    for(int t = 2; t<ntimesteps;t++){
        Fill_Rho(len,X,Rho,m,h);
        Fill_P(len,X,Rho,P,k,rho0);
        Fill_Fp(len,X,P,Rho,Fp,m,h);
        Fill_Fv(len,X,V,Rho,Fv,m,h,mu);
        Fill_Fg(len,Rho,Fg,g);
        Fill_A(len,Fp,Fv,Fg,Rho,A);


        
        Tintegrate(len,A,X,V,step);
        EnforceBC(len,X,V,e,h);


        if(t%1000 == 0){
            save_energy(Ek(len,V,m),Ep(len,X,m,g),t*step,delim,false,false);
        }
        if(t%10000 == 0){
            if (verbose){
                std::cout<<"\nt =  "<<t*step<<"s"<<std::endl;
                std::cout<<"VELOCITIES:\n";
                print_arr(len,V);
                std::cout<<"POSITIONS:\n";
                print_arr(len,X);
                // std::cout<<"NET SCALED FORCE:\n";
                // print_arr(len,A);
                // std::cout<<"KINETIC ENERGY:\n";
                // std::cout<<Ek(len,V,m)<<std::endl;
                std::cout<<"POTENTIAL ENERGY:\n";
                std::cout<<Ep(len,X,m,g)<<std::endl;
            }

            //SAVE DATA
            save_position(len,X,t*step,delim,false,false);



        }
        



        

    }
    save_position(len,X,sim_time,",",false,true);

}



int main(int argc, char* argv[]){

    // ----------------------
    // INITIALISE MPI
    // ----------------------
    int rank = 0;
    int size = 0;
    int confirmation = 0;

    int err = MPI_Init (&argc, &argv);
    if(err != MPI_SUCCESS){
        std::cout<<"Error initialising MPI"<<std::endl;
        return -1;
    }

    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    if(rank == 0){
        std::cout<<"\n\n\n\n\n";
        std::cout<<"MPI INIT SUCCESSFUL\nRUNNING ON "<<size<<" PROCESSES\n";
        
    }
    MPI_Bcast(&confirmation, 1, MPI_INT, 0, MPI_COMM_WORLD);

    

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


    //segment work
    auto [start,end] = segment_work(len,size,rank);
    std::cout<<"RANK "<<rank<<"\n SLICE: "<<start<<"-"<<end<<std::endl;

    
    if(rank == 0){
        std::cout<<"RANK "<<rank<<" BUSY\n";
        Simulate(Rho,P,X,V,Fp,Fv,Fg,A);
        std::cout<<"RANK "<<rank<<" DONE\n";
    }else{
        std::cout<<"RANK "<<rank<<" IDLE\n";
    }

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

    //-----------------------
    // CLOSE THREADS/PROCESSES
    // -----------------------

    MPI_Finalize();
    

}
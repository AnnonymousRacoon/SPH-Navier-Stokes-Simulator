
#include<iostream>
#include <mpi.h>
#include <math.h>
#include "SPH.hpp"






// ----------------------------
//            MAIN
// ----------------------------
int main(int argc, char** argv) {

    
    double* X = new double[5];


    MPI_Comm world = MPI_COMM_WORLD;

    SPH sph;
    sph.InitMPI(&argc,argv,world);
    sph.set_delim(" ");

  
    // sph.Model_Block_Droplet(10);

    sph.Model_Droplet(100);


    

   

    // X[0] = 0.01;
    // X[10] = 0.011;
    // X[2] = 0.99;
    // X[8] = 0.5111;
    // X[4] = 0.99;
    // X[6] = 0.989;

    // X[1] = X[3] = X[5] = 0.0;
    // X[7] = X[9] = X[11] = 0.5;

    // X[12] = 0.5;
    // X[13] = 0.25;

    // X[14] = 0.35;
    // X[15] = 0.25;

    // X[16] = 0.35;
    // X[17] = 0.25;

    // sph.Val_Single_Particle();

    // sph.Simulate_SC(1,X);
    // sph.Simulate_SC(7,X);

 

    




}
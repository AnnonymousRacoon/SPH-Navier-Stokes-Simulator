#pragma once
#include <math.h>
#include <mpi.h>
#include<tuple>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

class SPH
{
private:

    // ---------------------------------------------------------
    // *********************************************************
    //                    PRIVATE VARS
    // *********************************************************
    // ---------------------------------------------------------


    // PHYSICAL CONSTANTS
    const double __g = 9.81; //gravitational constant
    double __mu = 1.0; //viscosity
    double __e = 0.5; //coefficient of restitution
    double __k = 2000.0;// Gas constant
    double __rho0 = 1000.0;//resting density
    double __h = 0.01;// Radius of influence
    double __m = 1.0;//scaled mass
    
    
    

    // SIMULATION PARAMS
    double __step = pow(10,-4); //timestep
    double __sim_time = 100.00; //simulation time in seconds
    int __ntimesteps = int((__sim_time+1)/__step); //maximum number of simulation timesteps
    std::string __delim = ",";//save file delimeter

    //MPI PARAMS
    MPI_Comm __world; //MPI world
    int __rank = 0; //MPI COMM SIZE
    int __size = 0; //MPI PROCESS RANK

    // ---------------------------------------------------------
    // *********************************************************
    //                    PRIVATE FUNCS
    // *********************************************************
    // ---------------------------------------------------------
  

    // -------------------------------------
    //             MISC FUNCTIONS
    // -------------------------------------

    double __nrm(const double& rx, const double& ry);
    double __qij(const double& rx, const double& ry);
    void __scale_m(const int& lenX, double* Rho);
    std::tuple<int,int> __segment_work(const int& lenX);

    // -------------------------------------
    //         PHI AND DERRIVATIVES
    // -------------------------------------

    double __phi_dij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj);
    std::tuple<double,double> __grad_phi_pij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj);
    double __grad2_phi_vij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj);

    double __rho_i(const int& i, const int& lenX,  double* X);
    

    // --------------------------------------
    //              FORCES
    // --------------------------------------
    double __P_i(const int& i, const int& lenX,double* Rho);
    std::tuple<double,double> __Fp_i(const int& i, const int& lenX,  double* X, double* Rho);
    std::tuple<double,double> __Fv_i(const int& i, const int& lenX,  double* X, double* V, double* Rho);
    double __Fg_i(const int& i, double* Rho);


    // --------------------------------------
    //          SINGLE CORE JOBS
    // --------------------------------------
    void __Fill_Rho_SC(const int& lenX,  double* X, double* Rho,const bool& initial_step = false);
    void __Fill_Fp_SC(const int& lenX,  double* X, double* Rho,double* Fp);
    void __Fill_Fv_SC(const int& lenX,  double* X, double* V, double* Rho,double* Fv);
    void __Fill_Fg_SC(const int& lenX,double* Rho, double* Fg);
    void __Fill_A_SC(const int& lenX, double* Fp,double* Fv,double* Fg,double* Rho, double* A);

    // -------------------------------------
    //               MPI JOBS
    // -------------------------------------

    void __Fill_Rho_MC(const int& start, const int& end, const int& lenX,  double* X, double* Rho, double* buff,const int& buffsize,const bool& initial_step = false);
    void __Fill_Fp_MC(const int& start, const int& end,const int& lenX,  double* X, double* Rho,double* Fp, double* buff,const int& buffsize);
    void __Fill_Fv_MC(const int& start, const int& end,const int& lenX,  double* X, double *V, double* Rho,double* Fv, double* buff,const int& buffsize);
    void __Fill_Fg_MC(const int& start, const int& end,const int& lenX,double* Rho, double* Fg,double* buff, const int& buffsize);
    void __Fill_A_MC(const int& start, const int& end,const int& lenX, double* Fp,double* Fv,double* Fg,double* Rho, double* A,double* buff, const int& buffsize);

    // -------------------------------------
    //             TIMESTEPS
    // -------------------------------------
    void __Tintegrate(const int& lenX, double* A, double* X, double* V,const bool& initial);
    void __EnforceBC(const int& lenX, double* X, double* V);

    // -------------------------------------
    //            ENERGIES
    // -------------------------------------

    void __Ek_SC(const int& lenX, double* V, double& Energy);
    void __Ep_SC(const int& lenX, double* X, double& Energy);

    void __Ek_MC(const int& lenX, double* V,double& Energy);
    void __Ep_MC(const int& lenX, double* X, double& Energy);

    // -------------------------------------
    //         SAVING DATA
    // -------------------------------------


    void __save_position(const int& lenX, double* X,const double& timestamp, const std::string& delim, const bool& Open, const bool& Close);
    void __save_energy(const double& KE, const double& PE, const double& timestamp, const std::string& delim, const bool& Open, const bool& Close);
    void __write_buff(std::string& buff, const std::string& fname);

        
public:

    // ---------------------------------------------------------
    // *********************************************************
    //                    PUBLIC FUNCS
    // *********************************************************
    // ---------------------------------------------------------

    // ----------------------------------------
    //          CONSTRUCTOR/DESTRUCTOR
    // ----------------------------------------
    
    SPH(){srand((unsigned) time(0));}//set random seed on init
    ~SPH(){ MPI_Finalize();}// close threads on destruction


    // -----------------------------------------
    //            DRIVER FUNCTIONS
    // -----------------------------------------

    int InitMPI(int* argc,char** argv,MPI_Comm world);
    void Simulate_MC(const int& nParticles, double* X, const bool& keepXp = true);
    void Simulate_SC(const int& nParticles, double* X, const bool& keepXp = true);



    void Val_Single_Particle();
    void Val_Two_Particles();
    void Val_Four_Particles();
    void Model_Damn_Break(const int& resolution, const double& noise_amplitude = 0.001);
    void Model_Block_Droplet(const int& resolution, const double& noise_amplitude = 0.001);
    void Model_Droplet(const int& nParticles);



    // ------------------------------------------
    //          GETTERS AND SETTERS
    // ------------------------------------------

    int rank();
    int size();
    

    void set_tstep(const double & step);
    void set_max_sim_time(const double & time);
    void set_delim(const std::string& delim);
    void set_mu(const double& mu);
    void set_e(const double& e);
    void set_k(const double& k);
    void set_rho0(const double& rho);
    void set_h(const double& h);

    // ------------------------------------------
    //             HELPER FUNCS
    // ------------------------------------------

    double gen_noise(const double& amplitude);
    void gen_grid(const double& px1,const double& py1,const double& px2,const double& py2, const int& resolution, const double& noise_amplitude, const int& buffsize, double* buff);
    int gen_sphere(const double& radius, const double& rx,const double& ry, const int& buffsize, double* buff);
};
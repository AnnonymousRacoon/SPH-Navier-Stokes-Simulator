#pragma once
#include <math.h>
#include <mpi.h>
#include<tuple>
#include <string>

class SPH
{
private:
    // ----------------------------
    // ----------------------------
    // PRIVATE VARS
    // ----------------------------
    // ----------------------------

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

    



    // ----------------------------
    // ----------------------------
    // PRIVATE FUNCS
    // ----------------------------
    // ----------------------------
    
public:
    SPH(){}
    ~SPH(){ MPI_Finalize();}


    // --------------------------------------------------------
    //                    DRIVER FUNCTIONS
    // ---------------------------------------------------------

    int InitMPI(int* argc,char** argv,MPI_Comm world);
    void Simulate_MC(const int& nParticles, double* X);
    void Simulate_SC(const int& nParticles, double* X);

    // --------------------------------------------------------
    //                    MISC FUNCTIONS
    // ---------------------------------------------------------

    double nrm(const double& rx, const double& ry);
    double qij(const double& rx, const double& ry);
    void scale_m(const int& lenX, double* Rho);
    std::tuple<int,int> segment_work(const int& lenX);

    // --------------------------------------------------------
    //                   PHI AND DERRIVATIVES
    // ---------------------------------------------------------

    double phi_dij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj);
    std::tuple<double,double> grad_phi_pij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj);
    double grad2_phi_vij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj);

    double rho_i(const int& i, const int& lenX,  double* X);
    

    // --------------------------------------------------------
    //                      FORCES
    // ---------------------------------------------------------
    double P_i(const int& i, const int& lenX,double* Rho);
    std::tuple<double,double> Fp_i(const int& i, const int& lenX,  double* X, double* Rho);
    std::tuple<double,double> Fv_i(const int& i, const int& lenX,  double* X, double* V, double* Rho);
    double Fg_i(const int& i, double* Rho);


    // --------------------------------------------------------
    //                   SINGLE CORE JOBS
    // ---------------------------------------------------------
    void Fill_Rho_SC(const int& lenX,  double* X, double* Rho,const bool& initial_step = false);
    void Fill_Fp_SC(const int& lenX,  double* X, double* Rho,double* Fp);
    void Fill_Fv_SC(const int& lenX,  double* X, double* V, double* Rho,double* Fv);
    void Fill_Fg_SC(const int& lenX,double* Rho, double* Fg);
    void Fill_A_SC(const int& lenX, double* Fp,double* Fv,double* Fg,double* Rho, double* A);

    // --------------------------------------------------------
    //                      MPI JOBS
    // ---------------------------------------------------------

    void Fill_Rho_MC(const int& start, const int& end, const int& lenX,  double* X, double* Rho, double* buff,const int& buffsize,const bool& initial_step = false);
    void Fill_Fp_MC(const int& start, const int& end,const int& lenX,  double* X, double* Rho,double* Fp, double* buff,const int& buffsize);
    void Fill_Fv_MC(const int& start, const int& end,const int& lenX,  double* X, double *V, double* Rho,double* Fv, double* buff,const int& buffsize);
    void Fill_Fg_MC(const int& start, const int& end,const int& lenX,double* Rho, double* Fg,double* buff, const int& buffsize);
    void Fill_A_MC(const int& start, const int& end,const int& lenX, double* Fp,double* Fv,double* Fg,double* Rho, double* A,double* buff, const int& buffsize);

    // --------------------------------------------------------------
    //                       TIMESTEPS
    // --------------------------------------------------------------
    void Tintegrate(const int& lenX, double* A, double* X, double* V,const bool& initial);
    void EnforceBC(const int& lenX, double* X, double* V);

    // --------------------------------------------------------------
    //                       ENERGIES
    // --------------------------------------------------------------

    void Ek_SC(const int& lenX, double* V, double& Energy);
    void Ep_SC(const int& lenX, double* X, double& Energy);

    void Ek_MC(const int& lenX, double* V,double& Energy);
    void Ep_MC(const int& lenX, double* X, double& Energy);

    // --------------------------------------------------------------
    //                       SAVING DATA
    // --------------------------------------------------------------


    void save_position(const int& lenX, double* X,const double& timestamp, const std::string& delim, const bool& Open, const bool& Close);
    void save_energy(const double& KE, const double& PE, const double& timestamp, const std::string& delim, const bool& Open, const bool& Close);

    // --------------------------------------------------------------
    //                       GETTERS AND SETTERS
    // --------------------------------------------------------------

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
};
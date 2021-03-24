#include "SPH.hpp"
#include <iostream>
#include <iomanip>  
#include <fstream>
#include <time.h>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

// --------------------------------------------------------
//                        DRIVER FUNCTIONS
// ---------------------------------------------------------
/**
 * @brief Initialise the SPH object for use with MPI
 * 
 * @param argc 
 * @param argv 
 */
int SPH::InitMPI(int* argc,char** argv,MPI_Comm world){
    __world = world;
    int err = MPI_Init (argc, &argv);
    if(err != MPI_SUCCESS){
        std::cout<<"Error initialising MPI"<<std::endl;
        return -1;
    }

    MPI_Comm_rank(world, &__rank);
    MPI_Comm_size(world, &__size);
   

    if (__size< 2){
        std::cout<<"WARNING MPIC++ CONFIGURED TO RUN ON A SINGLE CORE"<<std::endl;
        return -1;
    }

    return 0;




}

/**
 * @brief Runs an SPH Simulation using MPI on multiple cores
 * 
 * @param nParticles the number of particles represented by the system
 * @param Xp array of initial positions of the particles. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 */
void SPH::Simulate_MC(const int& nParticles, double* Xp, const bool& keepXp){

    if(__rank ==0){
        std::cout<<"SIMULATING ON "<<__size<<" CORES\n";
        
    }

    //START TIMER
    auto t1 = Clock::now();

    int vector_len = 2*nParticles;
    int scalarbuffsize = (nParticles %__size == 0)? nParticles/__size : 1+(nParticles/__size);
  

    
    int vectorbuffsize = 2*scalarbuffsize;


    // ------------------------------------------
    //      ALLOCATE MEMORY
    // ------------------------------------------
    double* X = new double[vectorbuffsize*__size]; //X, Y Coordinate
    double* V = new double[vectorbuffsize*__size]; //U,V Velocity
    double* RHO= new double[scalarbuffsize*__size];//Density
    double* buff = new double[vectorbuffsize*__size];//Memory Buffer
    double* FP = new double[vectorbuffsize*__size]; //Pressure Force
    double* FV = new double[vectorbuffsize*__size]; //Viscous Force
    double* FG = new double[scalarbuffsize*__size]; //Gravitational Force
    double* A = new double[vectorbuffsize*__size]; //A Matrix

    int start = __rank*scalarbuffsize;// rank start 
    int end = (__rank+1)*scalarbuffsize;// rank end



   
    for (int i = 0; i< vector_len; i++){
        V[i] = 0.0; //Particles initially at rest
        X[i] = Xp[i]; //Copy Xp into P (perhaps the user would like to keep X)
    }

    //clear XP if not required
    if (not keepXp){
        delete[] Xp;
    }


 

    __Fill_Rho_MC(start,end,nParticles,X,RHO,buff,scalarbuffsize,true);
    __Fill_Fp_MC(start,end,nParticles,X,RHO,FP,buff,vectorbuffsize);
    __Fill_Fv_MC(start,end,nParticles,X,V,RHO,FV,buff,vectorbuffsize);
    __Fill_Fg_MC(start,end,nParticles,RHO,FG,buff,scalarbuffsize);
    __Fill_A_MC(start,end,nParticles,FP,FV,FG,RHO,A,buff,vectorbuffsize);


   
    __Tintegrate(nParticles,A,X,V,true);
    __EnforceBC(nParticles,X,V);

    // for (int i = 0; i< nParticles; i++){
    //     std::cout<<"("<<X[2*i]<<","<<X[1+2*i]<<")\n";
    // }

    double Ek;
    double Ep;
    __Ek_MC(nParticles,V,Ek);
    __Ep_MC(nParticles,X,Ep);

    if (__rank==0){
        __save_energy(Ek,Ep,__step,__delim,true,false);
        __save_position(nParticles,X,__step,__delim,true,false);
    }


    for(int t = 2; t<__ntimesteps;t++){
        __Fill_Rho_MC(start,end,nParticles,X,RHO,buff,scalarbuffsize);
        __Fill_Fp_MC(start,end,nParticles,X,RHO,FP,buff,vectorbuffsize);
        __Fill_Fv_MC(start,end,nParticles,X,V,RHO,FV,buff,vectorbuffsize);
        __Fill_Fg_MC(start,end,nParticles,RHO,FG,buff,scalarbuffsize);
        __Fill_A_MC(start,end,nParticles,FP,FV,FG,RHO,A,buff,vectorbuffsize);

        __Tintegrate(nParticles,A,X,V,false);
        __EnforceBC(nParticles,X,V);

        __Ek_MC(nParticles,V,Ek);
        __Ep_MC(nParticles,X,Ep);

        
        //SAVE STATES AT MICROSECOND INTERVALS
        if(t%100 == 0 and __rank == 0){
            __save_energy(Ek,Ep,t*__step,__delim,false,false);
            __save_position(nParticles,X,t*__step,__delim,false,false);
        }

       

        // BREAK EARLY IF AT REST AND SIM OVER 5s
        if( Ek <= nParticles*1.3485118e-07 and Ep <=nParticles*0.028 and t*__step > 5){
            break;
        }


        
 

    }


    // ------------------------------------------
    //      CLEAN UP
    // ------------------------------------------
    delete[] X;
    delete[] V;
    delete[] RHO;
    delete[] buff;
    delete[] FP;
    delete[] FV;
    delete[] FG;
    delete[] A;


    if(__rank == 0){
        auto t2 = Clock::now();
        std::cout<<"JOB COMPLETED\nCOMPUTE TIME: "<< std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()*pow(10,-9)
              << " seconds" << std::endl;
    }
    

}

/**
 * @brief Runs an SPH Simulation using a single core
 * 
 * @param nParticles the number of particles represented by the system
 * @param Xp array of initial positions of the particles. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 */
void SPH::Simulate_SC(const int& nParticles, double* Xp, const bool& keepXp){

    if(__rank ==0){
        std::cout<<"SIMULATING ON ONE CORE\n";
        
    }else{
        std::cout<<"PROCESS "<<__rank<<" IDLE\n";
        //Exit on non root processes
        return;
    }
    auto t1 = Clock::now();

    int vector_len = 2*nParticles;
    int scalarbuffsize = vector_len;

    int vectorbuffsize = 2*vector_len;

    
    double* X;
    // Program will overwrite Xp if the user requires
    if (not keepXp){
        X = Xp;
    //Else copy into a new array
    }else{
        X = new double[vectorbuffsize*__size];

    }

 

    
    double* V = new double[vectorbuffsize*__size];
    double* RHO= new double[scalarbuffsize*__size];
  
    double* FP = new double[vectorbuffsize*__size];
    double* FV = new double[vectorbuffsize*__size];
    double* FG = new double[scalarbuffsize*__size];
    double* A = new double[vectorbuffsize*__size];



    for (int i = 0; i< vector_len; i++){
        V[i] = 0.0;
        X[i] = Xp[i];
    }

    


    __Fill_Rho_SC(nParticles,X,RHO,true);
    __Fill_Fp_SC(nParticles,X,RHO,FP);
    __Fill_Fv_SC(nParticles,X,V,RHO,FV);
    __Fill_Fg_SC(nParticles,RHO,FG);
    __Fill_A_SC(nParticles,FP,FV,FG,RHO,A);

    __Tintegrate(nParticles,A,X,V,true);
    __EnforceBC(nParticles,X,V);

    double Ek;
    double Ep;
    __Ek_SC(nParticles,V,Ek);
    __Ep_SC(nParticles,X,Ep);

   
    if (__rank==0){
        __save_energy(Ek,Ep,__step,__delim,true,false);
        __save_position(nParticles,X,__step,__delim,true,false);
    }

    

    for(int t = 2; t<__ntimesteps;t++){
        __Fill_Rho_SC(nParticles,X,RHO);
        __Fill_Fp_SC(nParticles,X,RHO,FP);
        __Fill_Fv_SC(nParticles,X,V,RHO,FV);
        __Fill_Fg_SC(nParticles,RHO,FG);
        __Fill_A_SC(nParticles,FP,FV,FG,RHO,A);

        __Tintegrate(nParticles,A,X,V,false);
        __EnforceBC(nParticles,X,V);

        __Ek_SC(nParticles,V,Ek);
        __Ep_SC(nParticles,X,Ep);

        

        //SAVE STATES AT MICROSECOND INTERVALS
        if(t%100 == 0 and __rank == 0){
            __save_energy(Ek,Ep,t*__step,__delim,false,false);
            __save_position(nParticles,X,t*__step,__delim,false,false);
        }


        // BREAK EARLY IF AT REST AND SIM OVER 5s
        if( Ek <= nParticles*1.3485118e-07 and Ep <=nParticles*0.028 and t*__step > 5){
            break;
        }
 

    }


    delete[] X;
    delete[] V;
    delete[] RHO;
    
    delete[] FP;
    delete[] FV;
    delete[] FG;
    delete[] A;


    if(__rank == 0){
        auto t2 = Clock::now();
        std::cout<<"JOB COMPLETED\nCOMPUTE TIME: "<< std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()*pow(10,-9)
              << " seconds" << std::endl;
    }

}

// --------------------------------------------------------
//                        MISC FUNCTIONS
// ---------------------------------------------------------

/**
 * @brief computes the eucledian norm of a 2D tensor
 * 
 * @param rx component 1
 * @param ry component 2
 * @return double 
 */
double SPH::__nrm(const double& rx, const double& ry){
    return sqrt(pow(rx,2)+pow(ry,2));
}

/**
 * @brief returns qij of an interaction radius
 * 
 * @param rx radius x component
 * @param ry radius y component
 * @return double 
 */
double SPH::__qij(const double& rx, const double& ry){
    return __nrm(rx,ry)/__h;
}



/**
 * @brief scales particle mass and the initial unit density array
 * 
 * @param lenX number of particles in the system
 * @param Rho array of densities
 */
void SPH::__scale_m(const int& lenX, double* Rho){
    double sum = 0.0;
    for (int idx = 0; idx< lenX; idx++){
        sum+=Rho[idx];
    }
  
    __m = sqrt(__rho0*lenX/sum);


    for (int idx = 0; idx< lenX; idx++){
        Rho[idx]*=__m;
    }

    std::cout<<"SCALED MASS: "<<__m<<"\n";


}


/**
 * @brief fetches the required interval for an efficiently segmented process
 * 
 * @param lenX number of particles in the system
 * @return std::tuple<int,int> 
 */
std::tuple<int,int> SPH::__segment_work(const int& lenX){
    int start;
    int end;
    
    int switchpoint = lenX%__size;
    int small_chunk= (lenX - switchpoint)/__size;
    int big_chunk = small_chunk+1;

    if (__rank < switchpoint){
        start = __rank*big_chunk;
        end = (__rank+1)*big_chunk;
    }else{
        start = (__rank-switchpoint)*small_chunk + switchpoint*big_chunk;
        end = (__rank+1-switchpoint)*small_chunk + switchpoint*big_chunk;
    }

    return {start,end};
}


// --------------------------------------------------------
//                        PHI
// ---------------------------------------------------------

/**
 * @brief Computes Phi for density for particle i
 * 
 * @param Pxi x component of P_i
 * @param Pyi y component of P_i
 * @param Pxj x component of P_j
 * @param Pyj y component of P_j
 * @return double 
 */


double SPH::__phi_dij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj){
    double q = __qij(Pxi-Pxj,Pyi-Pyj);
    if (q >=1.0){
        return 0.0;
    }

    
    double ret = 4.0*pow((1.0-q*q),3)/(M_PI*__h*__h);

    return (isnan(ret))? 0.0 : ret;
}

/**
 * @brief computes grad phi of pressure for particle i
 * 
 * @param Pxi x component of P_i
 * @param Pyi y component of P_i
 * @param Pxj x component of P_j
 * @param Pyj y component of P_j
 * @return std::tuple<double,double> a tuple <phi_x, phi_y>
 */
std::tuple<double,double> SPH::__grad_phi_pij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj){
    double rx = Pxi - Pxj;
    double ry = Pyi - Pyj;
    double q = __qij(rx,ry);
  
    if (q >=1.0 or q== 0.0){
        return {0.0,0.0};
    }
    double multiplier = (-30.0*pow(1.0-q,2)/(M_PI*__h*__h*__h*q));
    
    // std::cout<<multiplier<<"\n";

    return{multiplier*rx,multiplier*ry};
}


/**
 * @brief computes grad^2 of phi of viscosity for particle i
 * 
 * @param Pxi x component of P_i
 * @param Pyi y component of P_i
 * @param Pxj x component of P_j
 * @param Pyj y component of P_j
 * @return double 
 */
double SPH::__grad2_phi_vij(const double& Pxi, const double& Pyi, const double& Pxj, const double& Pyj){
    double q = __qij(Pxi-Pxj,Pyi-Pyj);
    if (q >=1.0){
        return 0.0;
    }
    return 40.0*(1-q)/(M_PI*pow(__h,4));
}


/**
 * @brief Computes the scaled density for particle i
 * 
 * @param i index of particle in array storage
 * @param lenX total number of particles described by the system
 * @param X array of particle x positions
 * @param Y array of particle y positions
 * @return double 
 */
double SPH::__rho_i(const int& i, const int& lenX,  double* X){
   
    double sum = 0.0;
    for(int j = 0; j<lenX; j++){
        sum+=__phi_dij(X[2*i],X[(2*i)+1],X[2*j],X[(2*j)+1]);
    }
    return __m*sum;

}



// --------------------------------------------------------
//                      FORCES
// ---------------------------------------------------------


/**
 * @brief Computes the pressure of particle i
 * 
 * @param i index of particle in array storage 
 * @param lenX total number of particle represented by the system
 * @param Rho array of particle densities
 * @return double 
 */
double SPH::__P_i(const int& i, const int& lenX,double* Rho){
    return (Rho[i]-__rho0)*__k;
}



/**
 * @brief Calulates the Pressure Force of particle i
 * 
 * @param i index of particle in array storage 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Rho array of particle densities
 * @return std::tuple<double,double> 
 */
std::tuple<double,double> SPH::__Fp_i(const int& i, const int& lenX,  double* X, double* Rho){
    double sumx = 0.0;
    double sumy = 0.0;
    for (int j =0;j<lenX;j++){
        if(j==i){
           
            continue;
        }
    
        auto [cx, cy] = __grad_phi_pij(X[2*i],X[(2*i)+1],X[2*j],X[(2*j)+1]);
  


        double multiplier = -0.5*(__m/Rho[j])*(__P_i(i,lenX,Rho)+__P_i(j,lenX,Rho));

        sumx+=cx*multiplier;
        sumy+=cy*multiplier;



    }

    return {sumx, sumy};
}



/**
 * @brief calculates the viscous force of particle i
 * 
 * @param i index of particle in array storage 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param V 2D tensor array of particle velocities. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Rho array of particle densities
 * @return std::tuple<double,double> 
 */
std::tuple<double,double> SPH::__Fv_i(const int& i, const int& lenX,  double* X, double* V, double* Rho){
    double sumx = 0.0;
    double sumy = 0.0;
    double uij;
    double vij;
    for (int j =0;j<lenX;j++){
        if(j==i){
            continue;
        }
        uij = V[(2*i)]-V[(2*j)];
        vij = V[(2*i)+1]-V[(2*j)+1];
        double multiplier = __grad2_phi_vij(X[2*i],X[(2*i)+1],X[2*j],X[(2*j)+1]) * (__m/Rho[j]/2);
        sumx-=uij*multiplier;
        sumy-=vij*multiplier;
    }
    return {sumx*__mu,sumy*__mu};

}


/**
 * @brief Returns the vertical component of the gravitational force 
 * 
 * @param i index of particle in array storage 
 * @param Rho array of particle densities
 * @return double 
 */
double SPH::__Fg_i(const int& i, double* Rho){
    return -Rho[i]*__g;
}



// --------------------------------------------------------
//                   SINGLE CORE JOBS
// ---------------------------------------------------------


/**
 * @briefSingle-core. Fills Array Rho
 * 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Rho total number of particle represented by the system
 * @param initial_step boolean value indicating whether the function is being called on the inital half timestep
 */
void SPH::__Fill_Rho_SC(const int& lenX,  double* X, double* Rho,const bool& initial_step){
    for(int idx = 0; idx< lenX; idx++){
        Rho[idx] = __rho_i(idx,lenX,X);
    }

    if (initial_step){
        __scale_m(lenX,Rho);
    }

}
/**
 * @brief Single-core. Fills Array Fp
 * 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Rho total number of particle represented by the system
 * @param Fp 2D tensor array of particlePressure forces. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 */
void SPH::__Fill_Fp_SC(const int& lenX,  double* X, double* Rho,double* Fp){
    for(int idx = 0; idx< lenX; idx++){
        for(int idx = 0; idx< lenX; idx++){
        
        auto [Fx, Fy] = __Fp_i(idx,lenX,X,Rho);

        int idxx = 2*idx;
        int idxy = idxx+1;
        
        Fp[idxx] = Fx;
        Fp[idxy] = Fy;
    }
    }

}

/**
 * @brief Single-core. Fills Array Fv
 * 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param V 2D tensor array of particle velocities. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Rho total number of particle represented by the system
 * @param Fv 2D tensor array of particle viscous forces. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 */
void SPH::__Fill_Fv_SC(const int& lenX,  double* X, double* V, double* Rho,double* Fv){
    for(int idx = 0; idx< lenX; idx++){
        
        auto [Fx, Fy] = __Fv_i(idx,lenX,X,V,Rho);

        int idxx = 2*idx;
        int idxy = idxx+1;
        
        Fv[idxx] = Fx;
        Fv[idxy] = Fy;
    }

}

/**
 * @brief Single-core. Fills array Fg
 * 
 * @param lenX total number of particle represented by the system
 * @param Rho total number of particle represented by the system
 * @param Fg 1D tensor array of particle viscous gravitaional Forces acting in the Y direction. 
 */
void SPH::__Fill_Fg_SC(const int& lenX,double* Rho, double* Fg){
    for(int idx = 0; idx< lenX; idx++){
        Fg[idx] = __Fg_i(idx,Rho);
       
    }

}

/**
 * @brief Single-core. fills array A, the tensor Sum of particle forces Normalised by density. 
 * 
 * @param lenX total number of particle represented by the system
 * @param Fp 2D tensor array of particlePressure forces. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Fv 2D tensor array of particle viscous forces. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Fg 1D tensor array of particle viscous gravitaional Forces acting in the Y direction. 
 * @param Rho total number of particle represented by the system
 * @param A 2D tensor Sum of particle forces Normalised by density. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 */
void SPH::__Fill_A_SC(const int& lenX, double* Fp,double* Fv,double* Fg,double* Rho, double* A){
    
    for(int idx = 0; idx< lenX; idx++){

        int idxx = 2*idx;
        int idxy = idxx+1;

        A[idxx] = (Fp[idxx]+Fv[idxx])/Rho[idx];
        A[idxy] = (Fp[idxy]+Fv[idxy] + Fg[idx])/Rho[idx];
    }
}




// --------------------------------------------------------
//                      MPI JOBS
// ---------------------------------------------------------

/**
 * @brief Multi-core. Fills Array Rho
 * 
 * @param start rank start index
 * @param end rank end index
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Rho total number of particle represented by the system
 * @param buff memory buffer for MPI
 * @param buffsize size of buff
 * @param initial_step boolean value indicating whether the function is being called on the inital half timestep
 */
void SPH::__Fill_Rho_MC(const int& start, const int& end, const int& lenX,  double* X, double* Rho, double* buff,const int& buffsize, const bool& initial_step){
    for(int i = start; i<end; i++){
        if (i >= lenX){
            buff[i-start] =-1.0;
            continue;
        }
        buff[i-start] = __rho_i(i,lenX,X);
     
    }


    MPI_Gather(buff, buffsize, MPI_DOUBLE, Rho,buffsize, MPI_DOUBLE, 0, __world);
    MPI_Barrier(__world);
    if (__rank==0 and initial_step){
        __scale_m(lenX,Rho);
    }
    if(initial_step){
        MPI_Bcast(&__m, 1, MPI_DOUBLE, 0, __world);
    }

    MPI_Bcast(Rho, buffsize*__size, MPI_DOUBLE, 0, __world);

}


/**
 * @brief Multi-core. Fills Array Fp
 * 
 * @param start rank start index
 * @param end rank end index
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Rho total number of particle represented by the system
 * @param Fp 2D tensor array of particle Pressure forces. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param buff memory buffer for MPI. 
 * @param buffsize size of buff
 */
void SPH::__Fill_Fp_MC(const int& start, const int& end,const int& lenX,  double* X, double* Rho,double* Fp, double* buff,const int& buffsize){
    
    
    for(int i = start; i<end; i++){
        int idx = i-start;
        if (i >= lenX){
            buff[2*idx] =-1.0;
            continue;
        }
        auto [Fx, Fy] = __Fp_i(i,lenX,X,Rho);
        buff[2*idx] = Fx;
        buff[(2*idx)+1] = Fy;
        
    }



   
    MPI_Gather(buff, buffsize,MPI_DOUBLE, Fp,buffsize, MPI_DOUBLE, 0, __world);
   
    MPI_Barrier(__world);
    MPI_Bcast(Fp, buffsize*__size, MPI_DOUBLE, 0, __world);
}


/**
 * @brief Multi-core. Fills Array Fv
 * 
 * @param start rank start index
 * @param end rank end index
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param V 2D tensor array of particle velocities. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Rho total number of particle represented by the system
 * @param Fv 2D tensor array of particle viscous forces. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param buff memory buffer for MPI
 * @param buffsize size of buff
 */
void SPH::__Fill_Fv_MC(const int& start, const int& end,const int& lenX,  double* X, double *V, double* Rho,double* Fv, double* buff,const int& buffsize){
    for(int i = start; i<end; i++){
        int idx = i-start;
        if (i >= lenX){
            buff[2*idx] =-1.0;
            continue;
        }
        auto [Fx, Fy] = __Fv_i(i,lenX,X,V,Rho);
        buff[2*idx] = Fx;
        buff[(2*idx)+1] = Fy;
        
    }



  
    MPI_Gather(buff, buffsize, MPI_DOUBLE, Fv,buffsize, MPI_DOUBLE, 0, __world);
  
    MPI_Barrier(__world);
    MPI_Bcast(Fv, buffsize*__size, MPI_DOUBLE, 0, __world);
    MPI_Barrier(__world);


}

/**
 * @brief Multi-core. Fills array Fg
 * 
 * @param start rank start index
 * @param end rank end index
 * @param lenX total number of particle represented by the system
 * @param Rho total number of particle represented by the system
 * @param Fg Target array. 1D tensor array of particle viscous gravitaional Forces acting in the Y direction. 
 * @param buff memory buffer for MPI
 * @param buffsize size of buff
 */
void SPH::__Fill_Fg_MC(const int& start, const int& end,const int& lenX,double* Rho, double* Fg,double* buff, const int& buffsize){
    for(int i = start; i<end; i++){
        int idx = i-start;
        if (i >= lenX){
            buff[idx] =-1.0;
            continue;
        }
        buff[idx] = __Fg_i(i,Rho); 
    }

 
    MPI_Gather(buff,buffsize, MPI_DOUBLE, Fg,buffsize, MPI_DOUBLE, 0, __world);
  
    MPI_Barrier(__world);
    MPI_Bcast(Fg, buffsize*__size, MPI_DOUBLE, 0, __world);
    MPI_Barrier(__world);

}

/**
 * @brief Multi-core. fills array A, the tensor Sum of particle forces Normalised by density. 
 * 
 * @param start rank start index
 * @param end rank end index
 * @param lenX total number of particle represented by the system
 * @param Fp 2D tensor array of particlePressure forces. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Fv 2D tensor array of particle viscous forces. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Fg 1D tensor array of particle viscous gravitaional Forces acting in the Y direction. 
 * @param Rho total number of particle represented by the system
 * @param A Target Array. 2D tensor Sum of particle forces Normalised by density. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param buff memory buffer for MPI
 * @param buffsize size of buff
 */
void SPH::__Fill_A_MC(const int& start, const int& end,const int& lenX, double* Fp,double* Fv,double* Fg,double* Rho, double* A,double* buff, const int& buffsize){
    
    
    for(int i = start; i<end; i++){
        int idx = i-start;
        if (i >= lenX){
            buff[2*idx] =-1.0;
            continue;
        }

        int idxx = 2*idx;
        int idxy = idxx+1;

       
        buff[idxx] = (Fv[idxx] + Fp[idxx])/Rho[idx];
        buff[idxy] = (Fv[idxy] + Fp[idxy] + Fg[idx])/Rho[idx];


        
    }

    MPI_Gather(buff, buffsize, MPI_DOUBLE, A,buffsize, MPI_DOUBLE, 0, __world);
   
    MPI_Barrier(__world);
    MPI_Bcast(A, buffsize*__size, MPI_DOUBLE, 0, __world);
    MPI_Barrier(__world);



}


// --------------------------------------------------------------
//                       TIMESTEPS
// --------------------------------------------------------------
 

 /**
  * @brief Increments the simulation 1 timestep in time
  * 
  * @param lenX total number of particle represented by the system
  * @param A 2D tensor Sum of particle forces Normalised by density. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
  * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
  * @param V 2D tensor array of particle velocities. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
  * @param initial boolean value indicating whether the function is being called on the inital half timestep
  */
void SPH::__Tintegrate(const int& lenX, double* A, double* X, double* V,const bool& initial){

    double adjustor = (initial)? 2.0 : 1.0;//Adjustment to be applied in the initial time step

    for (int idx = 0; idx < lenX; idx++){

        int idxx = 2*idx;
        int idxy = idxx+1;

        
        V[idxx] +=A[idxx]*__step/adjustor;
        V[idxy] +=A[idxy]*__step/adjustor;

        X[idxx] += V[idxx]*__step;
        X[idxy] += V[idxy]*__step;
    }
}




/**
 * @brief Enforce boundary Conditions
 * 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param V 2D tensor array of particle velocities. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 */
void SPH::__EnforceBC(const int& lenX, double* X, double* V){
     for (int idx = 0; idx < lenX; idx++){

        int idxx = 2*idx;
        int idxy = idxx+1;

        //BOUNDARY CONDITIONS X-axis
        if(X[idxx] < __h){
            X[idxx] = __h;
            V[idxx]*=__e;
        }else if(X[idxx] >(1-__h)){
            X[idxx] = (1-__h);
            V[idxx]*=__e;

        }

        //BOUNDARY CONDITIONS Y-axis
        if(X[idxy] < __h){
            X[idxy] = __h;
            V[idxy]*=__e;
        }else if(X[idxy] >(1-__h)){
            X[idxy] = (1-__h);
            V[idxy]*=__e;

        }
    }
    
}


// --------------------------------------------------------------
//                       ENERGIES
// --------------------------------------------------------------


/**
 * @brief Perform Single-Core Kinetic energy calculation
 * 
 * @param lenX total number of particle represented by the system
 * @param V 2D tensor array of particle velocities. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Energy Variable to write into
 */
void SPH::__Ek_SC(const int& lenX, double* V, double& Energy){

    double sum = 0.0;
    for (int idx = 0; idx < lenX; idx++){
        int idxx = 2*idx;
        int idxy = idxx+1;

        sum+=pow(__nrm(V[idxx],V[idxy]),2);
    }

    Energy = sum*0.5*__m;
}


/**
 * @brief Perform Multi-Core Kinetic energy calculation
 * 
 * @param lenX total number of particle represented by the system
 * @param V 2D tensor array of particle velocities. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Energy Variable to write into
 */
void SPH::__Ek_MC(const int& lenX, double* V, double& Energy){

    auto [start, end] = __segment_work(lenX);
    double sum = 0.0;
    for (int idx = start; idx < end; idx++){

        int idxx = 2*idx;
        int idxy = idxx+1;


        sum+=pow(__nrm(V[idxx],V[idxy]),2);
    }

    sum*=0.5*__m;


    MPI_Allreduce(&sum,&Energy,1,MPI_DOUBLE,MPI_SUM,__world);


}


/**
 * @brief Perform Single-Core Potential energy calculation
 * 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Energy Variable to write into
 */
void SPH::__Ep_SC(const int& lenX, double* X, double& Energy){

    double sum = 0.0;
    for (int idx = 0; idx < lenX; idx++){
        int idxy = (2*idx)+1;
        sum+=X[idxy];
    }
    Energy = sum*__m*__g;
}

/**
 * @brief Perform Multi-Core Potential energy calculation
 * 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param Energy Variable to write into
 */
void SPH::__Ep_MC(const int& lenX, double* X, double& Energy){

    auto [start, end] = __segment_work(lenX);
    double sum = 0.0;
    for (int idx = start; idx < end; idx++){
        int idxy = (2*idx)+1;
        sum+=X[idxy];
    }
    sum*=__m*__g;

    MPI_Allreduce(&sum,&Energy,1,MPI_DOUBLE,MPI_SUM,__world);


}



// --------------------------------------------------------------
//                       SAVING DATA
// --------------------------------------------------------------


/**
 * @brief save particle positions to a txt or csv file
 * 
 * @param lenX total number of particle represented by the system
 * @param X 2D tensor array of particle postions. Formatted using row-major monvention (e.g [X1,Y1,X2,Y2......Xn,Yn])
 * @param timestamp timestamp
 * @param delim delimiter. use  " " for plain text, use "," for .csv
 * @param Open True for Open and start new. False for append
 * @param Close Close File (obselete)
 */
void SPH::__save_position(const int& lenX, double* X,const double& timestamp, const std::string& delim, const bool& Open, const bool& Close){
    std::string ext = (delim == ",")? ".csv" : ".txt";//file extension
    std::string pfname = "output";//position filename
    // create new file
    std::ofstream PositionFile;
    
    //-----------------------
    //   OPEN FILE AND FORMAT
    //----------------------

    if(Open){//OVERWRITE AND ADD HEADER
        PositionFile.open(pfname+ext);
        PositionFile<<"Timestamp(s)"<<delim;

        for(int particle_idx = 0; particle_idx < lenX; particle_idx++){
            PositionFile <<"px"<<particle_idx << delim;
        }
        for(int particle_idx = 0; particle_idx < lenX; particle_idx++){
            PositionFile <<"py"<<particle_idx << ((particle_idx<lenX-1)? delim: "\n");

        }
    }else{//PARSE TO END BEFORE BEGINING WRITE
        PositionFile.open(pfname+ext,std::ios_base::app);
    }

    //-----------------
    //    WRITE DATA
    //-----------------
    PositionFile<<timestamp<<delim;

    for(int particle_idx = 0; particle_idx < lenX; particle_idx++){
        int idxx = 2*particle_idx;
        PositionFile << X[idxx] << delim;

    }
    for(int particle_idx = 0; particle_idx < lenX; particle_idx++){
        int idxy = (2*particle_idx)+1;
        PositionFile << X[idxy] << ((particle_idx<lenX-1)? delim: "\n");

    }


    //-----------------
    //    CLOSE
    //-----------------
    if (Close){
        PositionFile.close();
    }
    
}


/**
 * @brief saves particle energies to a txt or csv file
 * 
 * @param KE Kinetic Enery
 * @param PE Potential Energy
 * @param timestamp timestamp
 * @param delim delimiter. use  " " for plain text, use "," for .csv
 * @param Open True for Open and start new. False for append
 * @param Close Close File (obselete)
 */
void SPH::__save_energy(const double& KE, const double& PE, const double& timestamp, const std::string& delim, const bool& Open, const bool& Close){
    std::string ext = (delim == ",")? ".csv" : ".txt";//file extension
    std::string pfname = "energy";//position filename
    // create new file
    std::ofstream NRGFile;
    

    if(Open){//OVERWRITE AND ADD HEADER
        NRGFile.open(pfname+ext);
        NRGFile<<"Timestamp(s)"<<delim<<"KE"<<delim<<"PE\n";

    }else{//PARSE TO END BEFORE BEGINING WRITE
        NRGFile.open(pfname+ext,std::ios_base::app);
    }
    //WRITE DATA
    NRGFile<<timestamp<<delim<<std::setprecision(8)<<KE<<delim<<std::setprecision(8)<<PE<<"\n";

    //close
    if (Close){
        NRGFile.close();
    }
}

// --------------------------------------------------------------
//                       PARTICLE MODELLING
// --------------------------------------------------------------



/**
 * @brief Validates case for a single particle released at (0.5,0.5). Runs on a single core for efficiency
 * 
 */
void SPH::Val_Single_Particle(){
    double* Xbuff = new double[2]; //Position tensor memory buffer
    //PLACE PARTICLES
    Xbuff[0] = 0.5;
    Xbuff[1] = 0.5;
    Simulate_SC(1,Xbuff,false);
}

/**
 * @brief Validates case for two particles released at (0.5,0.5) and (0.5,h). Runs on a single core for efficiency
 * 
 */
void SPH::Val_Two_Particles(){
    double* Xbuff = new double[4]; //Position tensor memory buffer
    //PLACE PARTICLES
    Xbuff[0] = 0.5;
    Xbuff[1] = 0.5;
    Xbuff[2] = 0.5;
    Xbuff[3] = __h;
    Simulate_SC(2,Xbuff,false);

}

/**
 * @brief Validates case for four particles released from (0.505, 0.5), (0.515, 0.5), (0.51, 0.45) and (0.5, 0.45). Runs on a single core for efficiency
 * 
 */
void SPH::Val_Four_Particles(){
    double* Xbuff = new double[8]; //Position tensor memory buffer
    //PLACE PARTICLES
    Xbuff[0] = 0.505;
    Xbuff[1] = 0.5;
    Xbuff[2] = 0.515;
    Xbuff[3] = 0.5;
    Xbuff[4] = 0.51;
    Xbuff[5] = 0.45;
    Xbuff[6] = 0.5;
    Xbuff[7] = 0.45;
    Simulate_SC(4,Xbuff,false);

}


/**
 * @brief Models a a grid of particles initially occupying the region [0, 0.2]^2. Runs on multiple cores.
 * 
 * @param resolution number of spanwise blocks to discretise the block into (e.g a resolution of 10 would create 10x10 pixels or an 11x11 grid of particles)
 * @param noise_amplitude noise amplitude to give when placing particles
 */
void SPH::Model_Damn_Break(const int& resolution, const double & noise_amplitude){

    int buffsize =2*(resolution+1)*(resolution+1);

    //GENERATE BLOCK SIMULATION
    double* Xbuff = new double[buffsize];
    if (__rank == 0){
        gen_grid(0.0,0.0,0.2,0.2,resolution,noise_amplitude,buffsize,Xbuff);
    }

    //BROADCAST TO OTHER RANKS
    MPI_Bcast(Xbuff, buffsize, MPI_DOUBLE, 0, __world);
   
    //SIMULATE USING MULTIPLE CORES
    Simulate_MC(buffsize/2,Xbuff,false);



}

/**
 * @brief Models a block droplet defined in the region [0.1, 0.3]x[0.3, 0.6]. Runs on multiple cores.
 * 
 * @param resolution number of spanwise blocks to discretise the block into (e.g a resolution of 10 would create 10x10 pixels or an 11x11 grid of particles)
 * @param noise_amplitude noise amplitude to give when placing particles
 */
void SPH::Model_Block_Droplet(const int& resolution, const double & noise_amplitude){

    int buffsize =2*(resolution+1)*(resolution+1);

    //GENERATE BLOCK SIMULATION
    double* Xbuff = new double[buffsize];
    if (__rank == 0){
        gen_grid(0.1,0.3,0.3,0.6,resolution,noise_amplitude,buffsize,Xbuff);
    }

    //BROADCAST TO OTHER RANKS
    MPI_Bcast(Xbuff, buffsize, MPI_DOUBLE, 0, __world);
   
    //SIMULATE USING MULTIPLE CORES
    Simulate_MC(buffsize/2,Xbuff,false);



}

/**
 * @brief Models particles occupying a circle of radius 0.1, centred at the point [0.5,0.7].
 * 
 * @param nParticles number of particles in the system
 */
void SPH::Model_Droplet(const int& nParticles){

    int buffsize =2*nParticles;

    //GENERATE SPHERE SIMULATION
    double* Xbuff = new double[buffsize];
    if (__rank == 0){
        gen_sphere(0.1,0.5,0.7,buffsize,Xbuff);
    }

    //BROADCAST TO OTHER RANKS
    MPI_Bcast(Xbuff, buffsize, MPI_DOUBLE, 0, __world);
   
    //SIMULATE USING MULTIPLE CORES
    Simulate_MC(buffsize/2,Xbuff,false);

}

// --------------------------------------------------------------
//                       GETTERS AND SETTERS
// --------------------------------------------------------------

/**
 * @brief returns the process rank
 * 
 * @return int 
 */
int SPH::rank(){
    return __rank;
}
/**
 * @brief returns the MPI COMM size
 * 
 * @return int 
 */
int SPH::size(){
    return __size;
}
    
/**
 * @brief sets the simulation time step
 * 
 * @param step 
 */
void SPH::set_tstep(const double & step){
    __step = step;
    __ntimesteps = int((__sim_time+1)/__step);
}

/**
 * @brief sets the maximum simulation time
 * 
 * @param time 
 */
void SPH::set_max_sim_time(const double & time){
    __sim_time = time;
    __ntimesteps = int((__sim_time+1)/__step);

}

/**
 * @brief sets the save file delimter. "," will output a .csv whilst " " will ouput a space separated text file
 * 
 * @param delim 
 */
void SPH::set_delim(const std::string& delim){
    __delim = delim;

}

/**
 * @brief sets the system viscosity
 * 
 * @param mu 
 */
void SPH::set_mu(const double& mu){
    __mu = mu;

}

/**
 * @brief sets the coefficient of restitution
 * 
 * @param e 
 */
void SPH::set_e(const double& e){
    __e = e;

}

/**
 * @brief sets the gas constant
 * 
 * @param k 
 */
void SPH::set_k(const double& k){
    __k = k;

}

/**
 * @brief sets the resting density
 * 
 * @param rho 
 */
void SPH::set_rho0(const double& rho){
    __rho0 = rho;

}

/**
 * @brief sets the radius of influence
 * 
 * @param h 
 */
void SPH::set_h(const double& h){
    __h = h;

}




// --------------------------------------------------------------
//                       HELPER FUNCS
// --------------------------------------------------------------

/**
 * @brief Generates oscillatory noise about 0
 * 
 * @param amplitude noise amplitude
 * @return double 
 */
double SPH::gen_noise(const double& amplitude){
    return(1.0-(2.0*(double)rand()/RAND_MAX))*amplitude;
}


/**
 * @brief Generates an evenly spaced grid of particles in a box defined by [px1,px2]x[py1,py2]
 * 
 * @param px1 
 * @param py1 
 * @param px2 
 * @param py2 
 * @param resolution 
 * @param noise_amplitude 
 * @param buffsize 
 * @param buff 
 */
void SPH::gen_grid(const double& px1,const double& py1,const double& px2,const double& py2, const int& resolution, const double& noise_amplitude, const int& buffsize, double* buff){
    int minbuffsize =2*(resolution+1)*(resolution+1);

    //CHECK BUFFER IS CORRECTLY SIZED
    try{
        if (buffsize < minbuffsize){
            throw(buffsize);
        }
    }catch(int memsize){
        std::cout << "\nMEMORY ERROR IN SPH::gen_grid\nMinimum buffer size is "
            <<minbuffsize<<" but "<<buffsize<< " was provided\n";
        return;

    }
    
 
    //CALCULATE POINT SPACING AND GRID ENDPOINTS
    double xstart = px1;
    double xend = px2;
    double xspacing =(xend-xstart)/double(resolution);//spacing between x points
    double ystart = py1;
    double yend = py2;
    double yspacing =(yend-ystart)/double(resolution);//spacing between y points


    //GRID CURSOR
    double px = xstart;
    double py = ystart;
    

    //LOOP OVER EACH GRIDPOINT
    int idx = 0;//buffer index
    for(int i = 0; i < resolution+1;i++){
        for(int j = 0; j < resolution+1;j++){
            //append gridpoints + noise
            buff[idx] = py+gen_noise(noise_amplitude);
            buff[idx+1] = px+gen_noise(noise_amplitude);
            //increment buffer index
            idx+=2;
            //increment x pointer
            px+=xspacing;
        }
        //return x pointer to start
        px = xstart;
        //increment y pointer
        py+=yspacing;
    }
}

 
/**
 * @brief generates points in a circle using random numbers
 * 
 * @param radius circle radius
 * @param rx circle x centre
 * @param ry circle y centre
 * @param buffsize buffer size
 * @param buff buffer to write into. Writes coordinates in row major tensor form (e.g [X1,Y1,X2,Y2....Xn,Yn])
 * @return int number of particles generated
 */
int SPH::gen_sphere(const double& radius, const double& rx,const double& ry, const int& buffsize, double* buff){

    double px;
    double py;

    //GENERATE RANDOM UNIT CIRCLE POINTS ABOUT 0 UNTIL MEMORY BUFFER IS FULL
    int idx = 0;//buffer index
    while(idx < floor(buffsize/2)){
        px = gen_noise(1.0);
        py = gen_noise(1.0);

        // IF POINTS ARE WITHIN THE CIRCLE SCALE AND PLACE
        if ((px*px)+(py*py)<=1){
            buff[2*idx] = rx+radius*px;
            buff[1+(2*idx)] = ry+radius*py;
            idx++;
        }
    }

    return buffsize/2;



}
#pragma once
#include "tensor.h"
#include <iostream>
#include <iomanip>  



/**
 * @brief returns the normalised interaction radius ||rij ||/h, of two points in 2 dimensional state space
 * 
 * @param xi coordinate of point 1
 * @param xj coordinate of point 2
 * @param h Radius of influence
 * @return double 
 */
double qij(tensor2<double> xi,tensor2<double> xj, const double& h);


// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
//                       PHI DERIVATIVES
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------

/**
 * @brief returns the value of phi_density for two points in 2D state space.
 * 
 * @param xi coordinate of point 1
 * @param xj coordinate of point 2
 * @param h Radius of influence
 * @return double 
 */
double phi_dij(tensor2<double> xi,tensor2<double> xj,const double& h);


/**
 * @brief returns grad^2 phi_viscous for two points in 2D state space
 * 
 * @param xi coordinate of point 1
 * @param xj coordinate of point 2
 * @param h Radius of influence
 * @return double 
 */
double grad2_phi_vij(tensor2<double> xi,tensor2<double> xj,const double& h);


/**
 * @brief returns the value of grad phi_pressure for two points in 2D state space.
 * 
 * @param xi coordinate of point 1
 * @param xj coordinate of point 2
 * @param h Radius of influence
 * @return double 
 */
tensor2<double> grad_phi_pij(tensor2<double> xi, tensor2<double> xj,const double& h);



// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
//                       PRESSURES AND DENSITY
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------


/**
 * @brief returns the fluid density associated with particle i
 * 
 * @param i index of particle i in array X
 * @param lenX length of array x
 * @param X pointer to array X
 * @param m scaled mass
 * @param h Radius of influence
 * @return double 
 */
double rho_i(const int& i, const int& lenX,  tensor2<double>* X, const double& m, const double& h);



/**
 * @brief returns the pressure associated with particle i using the ideal gas law
 * 
 * @param i index of the particle in the density array Rho
 * @param lenRho the number of particles represented in array Rho
 * @param Rho the density array
 * @param k Gas constant
 * @param rho0 resting density
 * @return double 
 */
double P_i(const int& i, const int& lenRho,double* Rho, const double& k, const double& rho0);
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
//                       FORCES
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------

/**
 * @brief Calculates the pressure force associated with particle i as a 2D tensor
 * 
 * @param i index of the particle in the arrays X,P, and Rho
 * @param lenX the number of particles represented by arrays X,P and Rho
 * @param X position array (2D tensor)
 * @param P Pressure Array
 * @param Rho Density Array
 * @param m scaled mass
 * @param h Radius of influence
 * @return tensor2<double> 
 */
tensor2<double> Fp_i(const int& i, const int& lenX,  tensor2<double>* X, double* P, double* Rho,const double& m, const double& h);

/**
 * @brief Calculates the viscous force associated with particle i as a 2D tensor
 * 
 * @param i index of the particle in the arrays X,P, and V
 * @param lenX the number of particles represented by arrays X,P and V
 * @param X position array (2D tensor)
 * @param V velocity array(2D tensor)
 * @param Rho Densiity Array
 * @param m scaled mass
 * @param h Radius of influence
 * @param mu viscosity
 * @return tensor2<double> 
 */
tensor2<double> Fv_i(const int& i, const int& lenX,  tensor2<double>* X, tensor2<double>* V, double* Rho,const double& m, const double& h, const double& mu);

/**
 * @brief Calculates the gravitational force associated with particle i as a 2D tensor
 * 
 * @param i index of the particle in the density array Rho
 * @param Rho the density array
 * @param g gravitational constant
 * @return tensor2<double> 
 */
tensor2<double> Fg_i(const int& i, double* Rho,const double& g);





// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
//                       FILL FUNCS
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
/**
 * @brief Populates the array of Pressures (P) using the ideal gas law
 * 
 * @param lenX length of array X
 * @param X position array
 * @param Rho density array
 * @param P target Pressure array
 * @param k Gas constant
 * @param rho0 resting density
 * @return double 
 */
void Fill_P(const int& lenX,  tensor2<double>* X, double* Rho, double* P, const double& k, const double& rho0);


/**
 * @brief Populates the array of Densities (Rho) using the ideal gas law
 * 
 * @param lenX length of array X
 * @param X position array
 * @param Rho target density array
 * @param m scaled mass
 * @param h Radius of influence
 */
void Fill_Rho(const int& lenX,  tensor2<double>* X, double* Rho, const double& m, const double& h);





/**
 * @brief Populates the Array of pressure forces (Fp) as a 2D tensor
 * 
 * @param lenX the number of particles represented by Fp
 * @param X position array (2D tensor)
 * @param P Pressure array 
 * @param Rho Density array 
 * @param Fp target array to write into
 * @param m scaled mass
 * @param h Radius of influence
 */
void Fill_Fp(const int& lenX,  tensor2<double>* X, double* P, double* Rho,tensor2<double>* Fp,const double& m, const double& h);


/**
 * @brief Populates the Array of viscous forces (Fv) as a 2D tensor
 * 
 * @param lenX the number of particles represented by Fp
 * @param X position array (2D tensor)
 * @param V velocity array (2D tensor)
 * @param Rho Density array 
 * @param Fp target array to write into
 * @param m scaled mass
 * @param h Radius of influence
 * @param mu Viscosity
 */
void Fill_Fv(const int& lenX,  tensor2<double>* X, tensor2<double>* V, double* Rho,tensor2<double>* Fv,const double& m, const double& h, const double& mu);


/**
 * @brief Populates the Array of Gravitational forces (Fg) as a 2D tensor
 * 
 * @param lenX the number of particles represented by Fg
 * @param Rho Density array 
 * @param Fg target array to write into
 * @param g gravitational constant
 */
void Fill_Fg(const int& lenX, double* Rho,tensor2<double>* Fg,const double& g);


/**
 * @brief Populates the A array for use in explicit time integration
 * 
 * @param lenX number of particle represented by the system
 * @param Fp An array of Pressure forces pertaining to each particle
 * @param Fv An array of viscous forces pertaining to each particle
 * @param Fg An array of gravotational forces pertaining to each particle
 * @param Rho An array of denisities pertaining to each particle
 * @param A The target array to write A
 */
void Fill_A(const int& lenX, tensor2<double>* Fp,tensor2<double>* Fv,tensor2<double>* Fg,double* Rho, tensor2<double>* A);
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
//                      TIMESTEP FUNCS
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------







// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
//                      HELPER FUNCS
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------

/**
 * @brief prints an array in the form "[x0, x1,...xn]"
 * 
 * @tparam T 
 * @param LenA length of array A
 * @param A pointer to array A
 */
template <class T>
void print_arr(const int& LenA, T* A){
    std::cout<<"[";
    for(int idx = 0; idx< LenA;idx++){
        std::cout<<A[idx];
        std::string endl = (idx<LenA-1)? ", ": "";
        std::cout<<endl;
    }
    std::cout<<"]\n";

}
#pragma once
#include <string>



// ----------------------------------//
//              CONFIG               //
// ----------------------------------//
// double g = 9.81
// double mu = 1.0
// double e = 0.5
// double k = 2000.0
// double rho0 = 1000.0
// double h = 0.01
// double step = 0.0001
// double max_sim_time = 100.00
// string delim = ,
// ----------------------------------//
/**
 * @brief SPH configuration
 * @param g double: gravitational constant
 * @param mu double: viscosity
 * @param e double: coefficient of restitution
 * @param k double: gas constant
 * @param rho0 double: resting density
 * @param h double: radius of influence
 * @param step double: time step
 * @param n_particles int: approx number of particles to use in simulations
 * @param max_sim_time double: maximum simulation time
 * @param delim double: delimiter, "," outputs .csv, " " outputs a space separated .txt file
 */
struct Config {
    double g;//double: gravitational constant
    double mu;//double: viscosity
    double e;//double: coefficient of restitution
    double k;//double: gas constant
    double rho0;//double: resting density
    double h;//double: radius of influence
    double step;//double: time step
    double max_sim_time;//int: approx number of particles to use in simulations
    int n_particles;//double: maximum simulation time
    std::string delim;// double: delimiter, "," outputs .csv, " " outputs a space separated .txt file
};
void load_Config(Config& config,const std::string& path = "config.txt");
void save_as_Config(Config& config,const std::string& fname);
void save_Config(Config& config,const bool& overwrite = true);
void display_Config(Config& config);


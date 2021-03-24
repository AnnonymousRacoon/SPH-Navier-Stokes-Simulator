#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "config.hpp"
using namespace std;

// // ----------------------------------//
// //              CONFIG               //
// // ----------------------------------//

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
//            LOAD CONFIG            //
// ----------------------------------//
// Loads saved configuration into memory from config.txt
// inputs: struct Config&
// outputs: None
// ----------------------------------//
void load_Config(Config& config,const string& path) {
    ifstream fin(path);
    string line;
    while (getline(fin, line)) {
        istringstream sin(line.substr(line.find("=") +1));
        if (line.find("g") != -1){
            sin >> config.g;
        }
        else if (line.find("mu") != -1){
            sin >> std::setprecision(6)>>config.mu;
        }
        
        else if (line.find("k") != -1)
            sin >> config.k;
        else if (line.find("rho0") != -1)
            sin >> config.rho0;
        else if (line.find("h") != -1)
            sin >> config.h;
        else if (line.find("dt") != -1)
            sin >> config.step;
        else if (line.find("n_particles") != -1)
            sin >> config.n_particles;
        else if (line.find("T") != -1)
            sin >>std::setprecision(6)>> config.max_sim_time;
        else if (line.find("delim") != -1){
            sin >> config.delim;
            //ENSURE ',' or ""
            if (config.delim != ","){
                config.delim = " ";
            }
        }
        else if (line.find("e") != -1)
            sin >> config.e;
    }
}


// ----------------------------------//
//          SAVE CONFIG AS           //
// ----------------------------------//
// Saves new configuration. Default into config.txt unless overwrite specified to false
// inputs: struct Config&, string& filename
// outputs: None
// ----------------------------------//
void save_as_Config(Config& config,const string& fname){

    // create new file
    ofstream outFile;
    outFile.open(fname+".txt");

    // write
    outFile << "g= "<<config.g<< endl;
    outFile << "mu = "<<config.mu<< endl;
    outFile << "e = "<<config.e<< endl;
    outFile << "k = "<<config.k<< endl;
    outFile << "rho0 = "<<config.rho0<< endl;
    outFile << "h = "<<config.h<< endl;
    outFile << "dt = "<<config.step<< endl;
    outFile << "n_particles = "<<config.n_particles<< endl;
    outFile << "T = "<<config.max_sim_time<< endl;
    outFile << "delim = "<<config.delim<< endl;
    outFile.close();


}

// ----------------------------------//
//            SAVE CONFIG            //
// ----------------------------------//
// Saves new configuration. Default into config.txt unless overwrite specified to false
// inputs: struct Config&, optional bool overwrire
// outputs: None
// ----------------------------------//
void save_Config(Config& config,const bool& overwrite){

     // default filename
    string fname = "temp";

    // prompt user for filename if new file
    if (not overwrite){
        std::cout<<"SAVE CONFIG FILE AS: ";
        std::cin>>fname;
        std::cout<<endl;
    }


    // save_as
    save_as_Config(config,fname);


    // overwrite file
    if (overwrite){
        remove("config.txt");
        rename("temp.txt", "config.txt");
        
    }


}


// ----------------------------------//
//          DISPLAY CONFIG           //
// ----------------------------------//
// Displays the current configuration settings
// inputs: struct Config&
// outputs: None
// ----------------------------------//
void display_Config(Config& config){

    std::cout << "g= "<<config.g<< endl;
    std::cout << "mu = "<<config.mu<< endl;
    std::cout << "e = "<<config.e<< endl;
    std::cout << "k = "<<config.k<< endl;
    std::cout << "rho0 = "<<config.rho0<< endl;
    std::cout << "h = "<<config.h<< endl;
    std::cout << "dt = "<<config.step<< endl;
    std::cout << "n_particles = "<<config.n_particles<< endl;
    std::cout << "T = "<<config.max_sim_time<< endl;
    std::cout << "delim = "<<config.delim<< endl;

}
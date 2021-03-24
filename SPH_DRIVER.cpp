
#include<iostream>
#include <mpi.h>
#include <math.h>
#include "SPH.hpp"
#include "config.hpp"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

// using namespace std;
namespace po = boost::program_options;



// ----------------------------
//            MAIN
// ----------------------------
int main(int argc, char** argv) {


    // --------------------------//
    //     INIT SPH OBJECT       //
    // --------------------------//
    MPI_Comm world = MPI_COMM_WORLD;
    SPH sph;
    sph.InitMPI(&argc,argv,world);
    int rank = sph.rank();//process rank


    // --------------------------//
    //    LOAD CONFIG FILE       //
    // --------------------------//
    Config config;
    try{
        load_Config(config);
    }catch(...){
        //HANDLE MISSING FILE
        std::cout<<"ERROR LOADING CONFIG FILE\nLOADING DEFAULTS";
        config.g = 9.81;
        config.mu = 1.0;
        config.e = 0.5;
        config.k = 2000.0;
        config.rho0 = 1000.0;
        config.h = 0.01;
        config.step = 0.0001;
        config.n_particles = 100;
        config.max_sim_time = 100.00;
        config.delim = " ";

        save_Config(config,true);
    }

        // SET SPH PARAMS
        sph.set_tstep(config.step);
        sph.set_max_sim_time(config.max_sim_time);
        sph.set_e(config.e);
        sph.set_mu(config.mu);
        sph.set_k(config.k);
        sph.set_rho0(config.rho0);
        sph.set_h(config.h);
        sph.set_delim(config.delim);


    // --------------------------//
    //           CLI             //
    // --------------------------//

    try {

        // allowed options
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("config","display current configuration")
            ("load_config",po::value<std::string>(),"load a configuration configuration file from $FILENAME$.txt")
            ("save_config_as", po::value<std::string>(),"save current configuration to $FILENAME$.txt")
            ("ic-dam-break", "Use dam-break initial condition.")
            ("ic-block-drop", "Use block-drop initial condition.")
            ("ic-droplet", "Use droplet initial condition.")
            ("ic-one-particle", "Use one particle validation case initial condition.")
            ("ic-two-particles", "Use two particle validation case initial condition.")
            ("ic-four-particles", "Use four particle validation case initial condition.")
            ("dt",po::value<double>(), "Time-step to use.")
            ("T",po::value<double>(), "Total integration time.")
            ("h",po::value<double>(), "Radius of influence of each particle.")
            ("n_particles",po::value<double>(), "Approximate number of particles to simulate.")
            ("delim",po::value<std::string>(), "sets the delimiter for output files. can be \",\" or \" \"")
            ("edit_multi", po::value<std::vector<std::string>>()->multitoken(),"change settings. takes two arguments $SETTING$ and $VALUE$. accepts multiple arguments\ne.g --edit dt 0.001 h 0.5")
            
        ;

        // parse command line
        po::variables_map vm;        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        // ------------------//
        //       help        //
        // ------------------//
        if (vm.count("help")) {
            if(rank == 0){
                std::cout << desc << "\n";
                return 0;
            }
        }


        // ------------------//
        //   display config  //
        // ------------------//
        if (vm.count("config")) {
            std::cout << "Current configuration settings:\n";
            display_Config(config);
            return 0;
        }

        // ------------------//
        //    load config    //
        // ------------------//
        if (vm.count("load_config")) {
            std::cout << "Loading "
                << vm["load_config"].as<std::string>()<<".txt"<<std::endl;
            try{
                load_Config(config,vm["load_config"].as<std::string>()+".txt");
                save_Config(config);
            }catch(std::exception& e) {
                std::cerr << "error: " << e.what() << "\n";
                return 1;
            }
            return 0;
        } 

        // ------------------//
        //    save config    //
        // ------------------//
        if (vm.count("save_config_as")) {
            std::cout << "Saving compression to... "
                << vm["save_config"].as<std::string>()<<".txt"<<std::endl;
            try{
                save_as_Config(config,vm["save_config"].as<std::string>());
            }catch(std::exception& e) {
                std::cerr << "error: " << e.what() << "\n";
                return 1;
            }
            return 0;
        } 

        // ------------------//
        //     dam break     //
        // ------------------//
        if (vm.count("ic-dam-break")) {
            
            //COMPUTE RESOLUTION
            int res = int(floor(sqrt(config.n_particles)))-1;//resolution
            //RUN SPH ROUTINE
            sph.Model_Block_Droplet(res,0.001);
            return 0;
        }

        // ------------------//
        //    block drop     //
        // ------------------//
        if (vm.count("ic-block-drop")) {
            //COMPUTE RESOLUTION
            int res = int(floor(sqrt(config.n_particles)))-1;//resolution
            //RUN SPH ROUTINE
            sph.Model_Block_Droplet(res,0.001);
            return 0;
        }

        // ------------------//
        //       droplet     //
        // ------------------//
        if (vm.count("ic-droplet")) {
            //RUN SPH ROUTINE
            sph.Model_Droplet(100);
            return 0;
        }

        // ------------------//
        //     1 particle    //
        // ------------------//
        if (vm.count("ic-one-particle")) {
            //RUN SPH ROUTINE
            sph.Val_Single_Particle();
            return 0;
        }

        // ------------------//
        //     2 particle    //
        // ------------------//
        if (vm.count("ic-two-particles")) {
            //RUN SPH ROUTINE
            sph.Val_Two_Particles();
            return 0;
        }

        // ------------------//
        //     4 particle    //
        // ------------------//
        if (vm.count("ic-four-particles")) {
            //RUN SPH ROUTINE
            sph.Val_Four_Particles();
            return 0;
        }

        // ------------------//
        //         dt        //
        // ------------------//
        if (vm.count("dt")) {
            try{
                config.step = vm["dt"].as<double>();
                save_Config(config);
                display_Config(config);
                std::cout<<"configuration saved successfully!";
            }catch(std::exception& e) {
                std::cerr << "error: " << e.what() << "\n";
                return 1;
            }
            return 0;
        }

        // ------------------//
        //         T         //
        // ------------------//
        if (vm.count("T")) {
            try{
                config.max_sim_time = vm["T"].as<double>();
                save_Config(config);
                display_Config(config);
                std::cout<<"configuration saved successfully!";
            }catch(std::exception& e) {
                std::cerr << "error: " << e.what() << "\n";
                return 1;
            }
            return 0;
        }

        // ------------------//
        //         h         //
        // ------------------//
        if (vm.count("h")) {
            try{
                config.h = vm["h"].as<double>();
                save_Config(config);
                display_Config(config);
                std::cout<<"configuration saved successfully!";
            }catch(std::exception& e) {
                std::cerr << "error: " << e.what() << "\n";
                return 1;
            }
            return 0;
        }

        // ------------------//
        //        delim       //
        // ------------------//
        if (vm.count("delim")) {
            try{
                config.delim = vm["delim"].as<std::string>();
                save_Config(config);
                display_Config(config);
                std::cout<<"configuration saved successfully!";
            }catch(std::exception& e) {
                std::cerr << "error: " << e.what() << "\n";
                return 1;
            }
            return 0;
        }

        // ------------------//
        //   n_particles     //
        // ------------------//
        if (vm.count("n_particles")) {
            try{
                config.n_particles = vm["n_particles"].as<double>();
                save_Config(config);
                display_Config(config);
                std::cout<<"configuration saved successfully!";
            }catch(std::exception& e) {
                std::cerr << "error: " << e.what() << "\n";
                return 1;
            }
            return 0;
        }
        // ------------------//
        //   edit config     //
        // ------------------//
        if (vm.count("edit_multi")) {
            // get args
            std::vector<std::string> options = vm["edit_multi"].as< std::vector<std::string>>();
            
            // PARSE ARGS
            // --------
            auto it = options.begin();
            while (std::next(it,1) != options.end()){
                // g
                if (*it == "g"){
                    try{
                        config.g = stod(*std::next(it,1));

                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }

                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"num\n";
                    }
                }
                // h
                if (*it == "h"){
                    try{
                        config.h = stod(*std::next(it,1));

                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }

                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"num\n";
                    }
                }
                // mu
                if (*it == "mu"){
                    try{
                        config.mu = stod(*std::next(it,1));

                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }

                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"num\n";
                    }
                }
                // e
                if (*it == "e"){
                    try{
                        config.e = stod(*std::next(it,1));

                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }

                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"num\n";
                    }
                }
                //rho0
                if (*it == "rho0"){
                    try{
                        config.rho0 = stod(*std::next(it,1));

                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }

                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"num\n";
                    }
                }
                //dt
                if (*it == "dt"){
                    try{
                        config.step = stod(*std::next(it,1));

                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }

                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"num\n";
                    }
                }
                //T
                if (*it == "T"){
                    try{
                        config.max_sim_time = stod(*std::next(it,1));

                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }

                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"num\n";
                    }
                }
                //n_particles
                if (*it == "n_particles"){
                    try{
                        config.n_particles = stoi(*std::next(it,1));

                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }

                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"num\n";
                    }
                }
                //delim
                if (*it == "delim"){
                    try{
                        config.delim = *std::next(it,1);

                        //guard against strange delims
                        if (config.delim != ","){
                            config.delim = " ";
                        }
                        if (std::next(it,2) != options.end()){
                            std::advance(it,2);
                            continue;
                        }
                        break;
                    // handle error
                    }catch(std::exception& e) {
                        std::cerr << "error: " << e.what() << "\n"<<"skipping "<<"str\n";
                    }
                }

                std::cout<<*it<<std::endl;
                std::advance(it, 1);





            }

            // CONFIRM SETTINGS
            // ------------
            std::cout << "New configuration settings:\n";
            display_Config(config);
            std::cout<<"Confirm\n[Y/N]\n";

            std::string confirmation;
            std::cin>>confirmation;

            // save config if yes else abort
            if (confirmation == "Y" or confirmation =="y"){
                save_Config(config);
                std::cout<<"parameters changed succesfully!"<<std::endl;

            }else{
                std::cout<<"change aborted\n";

            }
            
        } 


    }

    //HANDLE ERRORS
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }

    return 0;


}
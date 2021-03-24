CXXFLAGS = -Wall -O3 -pedantic -std=c++17
CC=mpicxx
CXX=mpicxx
HDRS = SPH.hpp config.hpp
OBJS = SPH_DRIVER.o SPH.o config.o
LDLIBS = -lboost_program_options

#default compile and execute
default: build_help

# build object files for all .cpp files
.cpp.o: $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

# compile object files
SPH_SIM: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDLIBS) -o SPH_SIM

# execute and delete objs
build_help: SPH_SIM
	mpiexec -np 2 ./SPH_SIM --help
	rm *.o

# execute on one process
.PHONY: droplet
droplet: SPH_SIM
	mpiexec -np 6 ./SPH_SIM --ic-droplet
	rm *.o

# execute on two processes
.PHONY: validate
validate: SPH_SIM
	mpiexec -np 1 ./SPH_SIM --ic-one-particle
	rm *.o

.PHONY: dam
dam: SPH_SIM
	mpiexec -np 6 ./SPH_SIM --ic-dam-break
	rm *.o

# purge object files
.PHONY:clean
clean:
	rm SPH_SIM *.o
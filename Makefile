CXXFLAGS = -Wall -O3 -pedantic -std=c++17
CC=mpicxx
CXX=mpicxx
HDRS = SPH.hpp
OBJS = SPH_DRIVER.o SPH.o
LDLIBS = -lboost_program_options

#default compile and execute
default: run3

# build object files for all .cpp files
.cpp.o: $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

# compile object files
SPH_SIM: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDLIBS) -o SPH_SIM

# execute and delete objs
exec_delete: SPH_SIM
	./SPH_SIM 
	rm *.o

# execute on one process
.PHONY: run1
run1: SPH_SIM
	mpiexec -np 1 ./SPH_SIM

# execute on two processes
.PHONY: run2
run2: SPH_SIM
	mpiexec -np 2 ./SPH_SIM

.PHONY: run3
run3: SPH_SIM
	mpiexec -np 3 ./SPH_SIM

.PHONY: run4
run4: SPH_SIM
	mpiexec -np 4 ./SPH_SIM


# purge object files
.PHONY:clean
clean:
	rm SPH_SIM *.o
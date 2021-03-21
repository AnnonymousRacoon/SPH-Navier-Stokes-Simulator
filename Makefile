CXXFLAGS = -Wall -O3 -pedantic -std=c++17
CC=mpicxx
CXX=mpicxx
HDRS = tensor.h, fns.hpp
OBJS = SPHNS.o tensor.o fns.o
LDLIBS = -lboost_program_options

#default compile and execute
default: run1

# build object files for all .cpp files
.cpp.o: $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

# compile object files
sph: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDLIBS) -o sph

# execute and delete objs
exec_delete: sph
	./sph 
	rm *.o

# execute on one process
.PHONY: run1
run1: sph
	mpiexec -np 1 ./sph

# execute on two processes
.PHONY: run2
run2: sph
	mpiexec -np 2 ./sph


# purge object files
.PHONY:clean
clean:
	rm sph *.o


# CXXFLAGS = -Wall -O3 -pedantic -std=c++17 -Wextra
# CC=mpicxx
# CXX=mpicxx
# EXENAME = prog1    
# OBJS = SPHNS.o tensor.o fns.o


# all : $(EXENAME)

# $(EXENAME) : $(OBJS)
#   $(CXX) $(CXXFLAGS) $(OBJS) -o $(EXENAME)

# # SPHNS.o : tensor.h fns.hpp
# # link.o : link.h

# clean :
#   -rm -f *.o $(EXENAME)


# Makefile for Shallow Water Equation Solver

objs = SWE_solver.o flowMath.o fvm.o integration.o mesh.o postprec.o preproc.o solve.o variable.o
    
SWE : $(objs)
	g++-9 -g -O3 -o SWE $(objs) -DARMA_DONT_USE_WRAPPER -std=c++11 -larmadillo -lhdf5 -headerpad_max_install_names \
	-L/usr/local/Cellar/openblas/0.3.9/lib -lopenblas -L/usr/local/Cellar/lapack/3.9.0/lib -llapack \
	-I/usr/local/Cellar/armadillo/9.850.1_1/include -L/usr/local/Cellar/armadillo/9.850.1_1/lib \
	-I/usr/local/hdf5/HDF5/1.10.5/include/ -L/usr/local/hdf5/HDF5/1.10.5/lib -lhdf5_cpp 

SWE_solver.o       : global.hpp preproc.hpp solve.hpp postprec.hpp
	g++-9 -g -c SWE_solver.cpp  -std=c++11  -I/usr/local/Cellar/armadillo/9.850.1_1/include -L/usr/local/Cellar/armadillo/9.850.1_1/lib -I/usr/local/hdf5/HDF5/1.10.5/include/ -L/usr/local/hdf5/HDF5/1.10.5/lib -lhdf5_cpp 

flowMath.o  : flowMath.hpp
	g++-9 -g -c flowMath.cpp

fvm.o   : fvm.hpp global.hpp
	g++-9 -g -c -std=c++11 fvm.cpp -I/usr/local/Cellar/armadillo/9.850.1_1/include -L/usr/local/Cellar/armadillo/9.850.1_1/lib

integration.o  : integration.hpp
	g++-9 -g -c -std=c++11 integration.cpp

mesh.o  : global.hpp mesh.hpp
	g++-9 -g -c -std=c++11 mesh.cpp -I/usr/local/Cellar/armadillo/9.850.1_1/include -L/usr/local/Cellar/armadillo/9.850.1_1/lib

postprec.o    : postprec.hpp
	g++-9 -g -c -std=c++11 postprec.cpp -I/usr/local/Cellar/armadillo/9.850.1_1/include -L/usr/local/Cellar/armadillo/9.850.1_1/lib

preproc.o    : preproc.hpp
	g++-9 -g -c preproc.cpp  -std=c++11 -I/usr/local/hdf5/HDF5/1.10.5/include/ -L/usr/local/hdf5/HDF5/1.10.5/lib -lhdf5_cpp

solve.o	: global.hpp solve.hpp
	g++-9 -g -c -std=c++11 solve.cpp -I/usr/local/Cellar/armadillo/9.850.1_1/include -L/usr/local/Cellar/armadillo/9.850.1_1/lib

variable.o : variable.hpp
	g++-9 -g -c -std=c++11 variable.cpp

clean:
	rm $(objs)

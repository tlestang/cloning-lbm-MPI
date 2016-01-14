CXXFLAGS=-O3 

all: main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o
	mpic++ -o TLGK_LBM main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o 

# --- MAIN FUNCTIONS

main.o: mpi_LBM_GK.cpp
	mpic++ -o main.o -c mpi_LBM_GK.cpp $(CXXFLAGS)

# --- CORE FUNCTIONS ---

initialize_lattice_arrays.o: src/initialize_lattice_arrays.cpp
	g++ -o initialize_lattice_arrays.o -c src/initialize_lattice_arrays.cpp $(CXXFLAGS)
streamCollCompute.o: src/streamCollCompute.cpp
	g++ -o streamCollCompute.o -c src/streamCollCompute.cpp $(CXXFLAGS)
domain_noSlipWalls.o: src/domain_noSlipWalls.cpp
	g++ -o domain_noSlipWalls.o -c src/domain_noSlipWalls.cpp $(CXXFLAGS)
domain_InletOutlet.o: src/domain_InletOutlet.cpp
	g++ -o domain_InletOutlet.o -c src/domain_InletOutlet.cpp $(CXXFLAGS)
square.o: src/square.cpp
	g++ -o square.o -c src/square.cpp $(CXXFLAGS)
force.o: src/force.cpp
	g++ -o force.o -c src/force.cpp $(CXXFLAGS)
write_vtk.o: src/write_vtk.cpp
	g++ -o write_vtk.o -c src/write_vtk.cpp $(CXXFLAGS)
write_vtk_mean.o: src/write_vtk_mean.cpp
	g++ -o write_vtk_mean.o -c src/write_vtk_mean.cpp $(CXXFLAGS)
clean:
	rm -rf *.o
mrproper: clean
	rm -rf all;



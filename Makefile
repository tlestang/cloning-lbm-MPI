CXXFLAGS=-O3 
LPFLAGS=-lfftw3 -lm
all: main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o perturbation.o generate_random_field.o take_curl.o covariance_fctn.o randNormal.o
	mpic++ -o TLGK_LBM main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o perturbation.o generate_random_field.o take_curl.o covariance_fctn.o randNormal.o $(LPFLAGS)

# --- MAIN FUNCTIONS --- #

main.o: mpi_LBM_GK.cpp
	mpic++ -o main.o -c mpi_LBM_GK.cpp $(CXXFLAGS)

# --- LBM FUNCTIONS --- #

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

# --- PERTURBATION OF INITIAL CONDITION --- #
perturbation.o: perturbation_TLGK/perturbation.cpp
	g++ -o perturbation.o -c perturbation_TLGK/perturbation.cpp
generate_random_field.o: perturbation_TLGK/circulant_embedding/generate_random_field.cpp
	g++ -o generate_random_field.o -c perturbation_TLGK/circulant_embedding/generate_random_field.cpp 
covariance_fctn.o: perturbation_TLGK/circulant_embedding/covariance_fctn.cpp
	g++ -o covariance_fctn.o -c perturbation_TLGK/circulant_embedding/covariance_fctn.cpp
take_curl.o: perturbation_TLGK/circulant_embedding/take_curl.cpp
	g++ -o take_curl.o -c perturbation_TLGK/circulant_embedding/take_curl.cpp
randNormal.o: perturbation_TLGK/circulant_embedding/randNormal.cpp
	g++ -o randNormal.o -c perturbation_TLGK/circulant_embedding/randNormal.cpp

# --- CLEANING --- #
clean:
	rm -rf *.o
mrproper: clean
	rm -rf all;



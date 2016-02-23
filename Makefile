CXXFLAGS=-O3 
LPFLAGS=-lfftw3 -lm

all: main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o generate_random_field.o take_curl.o covariance_fctn.o randNormal.o generate_mask.o
	mpic++ -o TLGK_LBM main.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o generate_random_field.o take_curl.o covariance_fctn.o randNormal.o generate_mask.o $(LPFLAGS)

forcing: main_forcing.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o generate_random_field.o take_curl.o covariance_fctn.o randNormal.o generate_mask.o
	mpic++ -o TLGK_LBM_FORCING main_forcing.o initialize_lattice_arrays.o streamCollCompute.o domain_noSlipWalls.o square.o force.o write_vtk.o generate_random_field.o take_curl.o covariance_fctn.o randNormal.o generate_mask.o $(LPFLAGS)

# --- MAIN FUNCTIONS --- #

main.o: src/mpi_LBM_GK.cpp
	mpic++ -o main.o -c src/mpi_LBM_GK.cpp $(CXXFLAGS)
main_forcing.o: src/mpi_LBM_GK_forcing.cpp
	mpic++ -o main_forcing.o -c src/mpi_LBM_GK_forcing.cpp $(CXXFLAGS)

# --- LBM FUNCTIONS --- #

initialize_lattice_arrays.o: src/lbm_src/initialize_lattice_arrays.cpp
	g++ -o initialize_lattice_arrays.o -c src/lbm_src/initialize_lattice_arrays.cpp $(CXXFLAGS)
streamCollCompute.o: src/lbm_src/streamCollCompute.cpp
	g++ -o streamCollCompute.o -c src/lbm_src/streamCollCompute.cpp $(CXXFLAGS)
domain_noSlipWalls.o: src/lbm_src/domain_noSlipWalls.cpp
	g++ -o domain_noSlipWalls.o -c src/lbm_src/domain_noSlipWalls.cpp $(CXXFLAGS)
domain_InletOutlet.o: src/lbm_src/domain_InletOutlet.cpp
	g++ -o domain_InletOutlet.o -c src/lbm_src/domain_InletOutlet.cpp $(CXXFLAGS)
square.o: src/lbm_src/square.cpp
	g++ -o square.o -c src/lbm_src/square.cpp $(CXXFLAGS)
force.o: src/lbm_src/force.cpp
	g++ -o force.o -c src/lbm_src/force.cpp $(CXXFLAGS)
write_vtk.o: src/lbm_src/write_vtk.cpp
	g++ -o write_vtk.o -c src/lbm_src/write_vtk.cpp $(CXXFLAGS)
write_vtk_mean.o: src/lbm_src/write_vtk_mean.cpp
	g++ -o write_vtk_mean.o -c src/lbm_src/write_vtk_mean.cpp $(CXXFLAGS)

# --- PERTURBATION OF INITIAL CONDITION --- #
perturbation.o: src/perturbation_routines/perturbation_routines/perturbation.cpp
	g++ -o perturbation.o -c src/perturbation_routines/perturbation_routines/perturbation.cpp
generate_random_field.o: src/perturbation_routines/generate_random_field.cpp
	g++ -o generate_random_field.o -c src/perturbation_routines/generate_random_field.cpp 
covariance_fctn.o: src/perturbation_routines/covariance_fctn.cpp
	g++ -o covariance_fctn.o -c src/perturbation_routines/covariance_fctn.cpp
take_curl.o: src/perturbation_routines/take_curl.cpp
	g++ -o take_curl.o -c src/perturbation_routines/take_curl.cpp
randNormal.o: src/perturbation_routines/randNormal.cpp
	g++ -o randNormal.o -c src/perturbation_routines/randNormal.cpp
generate_mask.o: src/perturbation_routines/generate_mask.cpp
	g++ -o generate_mask.o -c src/perturbation_routines/generate_mask.cpp

# --- CLEANING --- #
clean:
	rm -rf *.o
mrproper: clean
	rm -rf all; rm -rf forcing;



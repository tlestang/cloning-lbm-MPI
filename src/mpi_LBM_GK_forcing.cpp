#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <malloc.h>
#include <unistd.h>
#include "stdio.h"
#include <sys/time.h>
#include<cmath>
#include<mpi.h>

#ifndef __global__
#define __global__
#include "lbm_src/global.h"
#endif

#include "lbm_src/initialize_lattice_arrays.h"
#include "lbm_src/streamCollCompute.h"
#include "lbm_src/boundaryConditions.h"
#include "lbm_src/force.h"
#include "lbm_src/write_vtk.h"

#include "perturbation.h"

#define MASTER 0

using namespace std;

//GLOBAL PARAMETERS FOR LBM
int c[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
int Dx, Dy, xmin, xmax, ymin, ymax;

int main()
{
  
  // --- PARAMETERS FOR TLGK ALGO. ---
  int Nc = 4; // Number of clones
  double T = 8; // Total simulation time
  double dT = 2; // Cloning timestep
  double dT0 = 2.0/10.0;
  //------------------------

  // --- PARAMETERS FOR LBM ---
  double tau = 1.0, beta = 1.0, Ma = 1.0, t0 = 1.0;
  int Lx = 0, Ly = 0;
  //READ INPUT FILE
  ifstream input_file("input_LBM.datin");
  input_file >> Lx; Ly = Lx;
  input_file >> tau;
  input_file >> Ma;
  input_file >> t0; //t0 IS THE TURN AROUND TIME (GIVEN IN LBM TIMESTEP)
  input_file.close();

  //------------------

  //VARIABLES FOR TLGK
  double phi_alpha, phi_theor;
  double alpha=0.4, alphaMin = -0.5, alphaIncr = 0.004, alphaMax = 0.5;
  int NcPrime, deltaN, copyIdx, k;
  int nbrTimeSteps = floor(T/dT);
  int l= 0; int idx;

  int my_rank, p, local_Nc, cloneIdx, cloneIdxMin, tag = 0;
  int ctm, cte, nbComm;
  int sender, dest;
  bool flag;
  //-------------------

  //VARIABLES FOR LBM
  double *fin, *fout, *pivot, *rho, *ux, *uy, *map;
  Dy = 4*Ly + 1, Dx = Dy;
  xmin = (Dx-1)/2; xmax = xmin + Lx;
  ymin = (Dy-1)/2 - Ly/2; ymax = ymin + Ly;
  double cs = 1./sqrt(3); double rho0 = 1.0;
  double nu = 1./3.*(tau-0.5);
  double u0 = cs*cs*Ma; double uxSum, uxMean;
  double beta0 = 8*nu*u0/((Dy-1)/2)/((Dy-1)/2); double a=1.0;
  double omega = 1.0/tau;
  double F;
  double delta_t = 1.0/t0;
  int error;
  int lbmTimeSteps1 = floor(dT0/delta_t);
  int lbmTimeSteps2 = floor((dT-dT0)/delta_t);
   
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);
  local_Nc = Nc/p;
  cloneIdxMin = my_rank*local_Nc;
  
  double **state;
  state = new double*[local_Nc];
    for(int j=0;j<local_Nc;j++)
    {
      // ONE POP ARRAY PER CLONE
      state[j] = (double *) memalign(getpagesize(), Dx*Dy*9*sizeof(double)); 
    }
  fout = (double *) memalign(getpagesize(), Dx*Dy*9*sizeof(double));
  rho = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  ux = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  uy = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  map = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));


  int temp[2*Nc]; 
  int comm_instru_Send[local_Nc]; int comm_instru_Recv[local_Nc];
  double s[Nc]; double s_;
  int nbCreatedCopies[Nc];
  double R_record[nbrTimeSteps]; double R, total_R;

  string path_to_control_run;
  path_to_control_run = "/home/thibault/tailleur_lecomte/lbm/L32_square_domain_pops.dat";
  ifstream crFileID;

  MPI_Status status;

  //Set up variables and containers for output
  string folderName[local_Nc], fileName, instru;
  string masterFolderName = "output/";
  ofstream output_file, weightsFile;
  instru = "mkdir " + masterFolderName;
  for(int i=0;i<local_Nc;i++)
    {
      stringstream buf;
      buf << "clone_" << i + my_rank*local_Nc;
      folderName[i] = buf.str();
      instru = "mkdir " + masterFolderName + folderName[i];
      system(instru.c_str());
    }
   
  if(my_rank==MASTER)
    {
      crFileID.open(path_to_control_run.c_str(), ios::binary);
      for(int j=0;j<local_Nc;j++)
	{
	  crFileID.seekg(0, std::ios::beg); //Set cursor to beginning of file
	  crFileID.read((char*)&state[j][0], Dx*Dy*9*sizeof(double));
	}
      crFileID.close();
    }

  //BROADCAST OF INITIAL POPULATIONS FROM MASTER TO OTHER PROCESSES
  for(int j=0;j<local_Nc;j++)
    {
      MPI_Bcast(&state[j][0], Dx*Dy*9, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    }

   //TIME EVOLUTION OVER TOTAL TIME T (T/dT CLONING STEPS)
  for(int t=0;t<nbrTimeSteps;t++)
    {

      R = 0.0;

      if(my_rank==MASTER)
	{
	  stringstream weightsFileName;
	  weightsFileName << "weights_evolution_" << t;
	  instru = masterFolderName + weightsFileName.str();
	  weightsFile.open(instru.c_str(), ios::binary);
	}
      
      
      //SIMULATE THE SYSTEM DURING dT AND COMPUTE WEIGHT
      for(int j=0;j<local_Nc;j++) // Loop on clones
	{
	  
	  s_ = 0;

	  stringstream buf;
	  buf << "/evolution_" << t << "_" << "clone_" << j;
	  fileName = folderName[j] + buf.str();
	  output_file.open(fileName.c_str(), ios::binary);
	  
	  generate_random_field(Dx, map, error);	  
	  for(int t=0;t<lbmTimeSteps1;t++)
	    {
	      streamingAndCollisionComputeMacroBodyForceSpatial(state[j], fout, rho, ux, uy, beta0, map, tau);
	      computeDomainNoSlipWalls_BB(fout, state[j]);
	      computeSquareBounceBack_TEST(fout, state[j]);
	      // RESET NODES INSIDE THE SQUARE TO EQUILIBRIUM DISTRIBUTION
	      for(int x=xmin+1;x<xmax;x++)
		{
		  for(int y=ymin+1;y<ymax;y++)
		    {
		      for(int k=0;k<9;k++)
			{
			  fout[IDX(x,y,k)] = w[k];
			}
		    }
		}
	      //SWAP POINTERS ON POPULATIONS FOR NEXT ITERATION OF LBM
	      pivot = state[j];
	      state[j] = fout;
	      fout = pivot;
	      // COMPUTE FORCE ON SQUARE
	      F = computeForceOnSquare(state[j], omega);
	      output_file.write((char*)&F, sizeof(double));
	      // COMPUTE WEIGHT
	      s_ += F;
	    } //END LOOP ON TIMESTEPS

	  for(int t=0;t<lbmTimeSteps2;t++)
	    {
	      streamingAndCollisionComputeMacroBodyForce(state[j], fout, rho, ux, uy, beta0, tau);
	      computeDomainNoSlipWalls_BB(fout, state[j]);
	      computeSquareBounceBack_TEST(fout, state[j]);
	      // RESET NODES INSIDE THE SQUARE TO EQUILIBRIUM DISTRIBUTION
	      for(int x=xmin+1;x<xmax;x++)
		{
		  for(int y=ymin+1;y<ymax;y++)
		    {
		      for(int k=0;k<9;k++)
			{
			  fout[IDX(x,y,k)] = w[k];
			}
		    }
		}
	      //SWAP POINTERS ON POPULATIONS FOR NEXT ITERATION OF LBM
	      pivot = state[j];
	      state[j] = fout;
	      fout = pivot;
	      // COMPUTE FORCE ON SQUARE
	      F = computeForceOnSquare(state[j], omega);
	      output_file.write((char*)&F, sizeof(double));
	      // COMPUTE WEIGHT
	      s_ += F;
	    } //END LOOP ON TIMESTEPS

	  output_file.close();
	  
	  s_ *= delta_t;
	  //STORE WEIGHT IN WEIGHTS ARRAY s[local_Nc]
	  s[j] = exp(alpha*s_);
	  //UPDATE LOCAL AVERAGE WEIGHT
	  R += s[j];
	}


      // ---------------------------------------------------------------------------------------------
      // ---------------------------------------------------------------------------------------------

      //NOW TIME EVOLUTION OF COPIES IS DONE AND WE MUST :
      // -- DETERMINE HOW MANY CLONES ARE GENERATED BY EACH COPY
      // -- CLONE/PRUNE TO STICK TO A CONSTANT NUMBER OF COPIES AT T + dT
      // -- CREATE INITIAL CONDITION FOR NEXT ITERATION (TIME EVOLUTION)
       

      //COMMUNICATIONS BETWEEN MASTER AND OTHER PROCESSES PRIOR TO CLONING STEP
      if(my_rank==MASTER)
	{
	  total_R = R; //COMPUTE MASTER's CONTRIB. TO TOTAL AVERAGE WEIGHT TOTAL_R
	  //GATHER WEIGHTS AND LOCAL AVERAGE WEIGHTS FROM OTHER PROCESSES
	  for(int source=1;source<p;source++)
	    {
	      tag = 1;
	      MPI_Recv(&s[source*local_Nc], local_Nc, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	      tag = 3;
	      MPI_Recv(&R, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
	      total_R += R;
	    }
	}
      else //SEND s[] VALUES FOR THE PROCESS'S CLONES AND LOCAL AVERAGE WEIGHT TO MASTER
	{
	  tag = 1;
	  MPI_Send(&s[0], local_Nc, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD);
	  tag = 3;
	  MPI_Send(&R, 1, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD);
	}

      //MASTER POST-PROCESSES EVOLUTION OF COPIES AND DO THE CLONING
      if(my_rank==MASTER)
	{
	  //WRITE WEIGHTS ON DISK
	  weightsFile.write((char*)&s[0], Nc*sizeof(double));
	  weightsFile.close();

	  total_R /= Nc; //NORMALIZATION OF THE AVERAGE WEIGHT

	  NcPrime = 0; //NcPRIME IS THE NUMBER OF COPIES AFTER CLONING
	  for(int j=0;j<Nc;j++)
	    {
	      cout << s[j] << endl;
	      //EACH COPY j LEADS TO nbCreatedCopies[j] CLONES (INCLUDING ITSELF)
	      nbCreatedCopies[j] = floor(s[j]/total_R + drand48());
	      
	      NcPrime += nbCreatedCopies[j]; //COMPUTE NcPrime BY COUNTING
	    }
	  //deltaN IS THE DIFFERENCE BETWEEN NcPRIME AND THE IMPOSED NB OF CLONES Nc
	  deltaN = NcPrime - Nc;
	  //FILLS THE TEMP[] ARRAY FOR UNIFORM SAMPLING OVER THE NEW CLONES (SEE MANUAL)
	  k=0;
	  for(int j=0;j<Nc;j++)
	    {
	      for(int i=0;i<nbCreatedCopies[j];i++)
		{
		  temp[k] = j;
		  k++;
		}
	    }

	  //  -------------------------------------------------------------------------------------------
	  //  -------------------------------------------------------------------------------------------

	  
	  // CLONING/PRUNING PHASE
	  if(deltaN > 0) // IF NcPrime > Nc KILL deltaN CLONES
	    {
	      for(int i=0;i<deltaN;i++)
		{
		  //CHOOSE A CLONE AT RAND. UNIFORMLY THROUGH NEWLY CREATED CLONES
		  idx = rand()%(NcPrime-i); //random number between [0:NcPrime-i-1]
		  //EXTRACT CLONE ABSOLUTE INDEX FROM TEMP[] ARRAY
		  cloneIdx = temp[idx];
		  //RECORD THE KILL
		  nbCreatedCopies[cloneIdx]--;
		  //UPDATE TEMP[] ARRAY FOR NEXT ITERATION OF THE PRUNING PROCESS
		  //(SEE MANUAL)
		  temp[idx] = temp[NcPrime-i-1];
		}
	    }
		
	  else if(deltaN < 0) // IF NcPrime < Nc ADD deltaN CLONE OF RANDOM COPIES
	    {
	      for(int i=0;i<-deltaN;i++)
		{
		  //CHOOSE A CLONE AT RAND. UNIFORMLY THROUGH NEWLY CREATED CLONES
		  idx = rand()%NcPrime;
		  //EXTRACT CLONE ABSOLUTE INDEX FROM TEMP[] ARRAY
		  cloneIdx = temp[idx];
		  //RECORD THE KILL
		  nbCreatedCopies[cloneIdx]++;
		}
	    }

	  // NOW CREATE COMMUNICATION TABLE (TEMP[] IS RECYCLED)
	  nbComm = 0; //nbComm IS THE NUMBER OF POINT TO POINT COMM. (SENDER,DEST)
	  //LOOP ON ALL Nc CLONES
	  for(int i=0;i<Nc;i++) 
	    {
	      //IF COPY GAVE BIRTH TO CLONES, LOOP ON THEM
	      while(nbCreatedCopies[i] > 1)
		{
		  flag = true; k=0;
		  //SEARCH FOR A KILLED COPY TO REPLACE IT BY THE CLONE OF COPY i
		  while(flag)
		    {
		      if(nbCreatedCopies[k]==0) //IF COPY k WERE KILLED
			{
			  //RECORD THE COMMUNICATION IN TEMP[]
			  temp[2*nbComm] = i; // AN INSTANCE OF COPY i MUST REPLACE COPY k
			  temp[2*nbComm+1] = k; // COPY k must be replaced by A CLONE OF COPY i
			  //NEXT LINE SO THAT COPY k IS NOT CHOSEN AGAIN IN FURTHER ITERATIONS 
			  nbCreatedCopies[k] = -1; 
			  nbComm++; //COMPUTE THE NB OF COMM. BY COUNTING
			  flag = false; //EXIT WHILE LOOP WHEN FOUND A PLACE FOR CLONE OF COPY i
			}
		      k++;
		    }
		  nbCreatedCopies[i]--; //RECORD THAT 1 CLONE OF COPY i HAS BEEN TAKEN CARE OF
		}
	    }
	  R_record[t] = total_R; // STORE AVERAGE VALUE FOR SCGF CALCULATION AT A LATER STAGE
	} // IF MASTER

	  //BROADCAST OF NB OF COMM. FROM MASTER
      MPI_Bcast(&nbComm, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
      //MAKE SURE EACH AND EVERY PROCESS KNOWs THE NB OF COMM. BEFORE RECEVING COMM TABLE
      MPI_Barrier(MPI_COMM_WORLD);
      //BROADCAST OF COMMUNICATION TABLE TEMP[] FROM MASTER
      MPI_Bcast(&temp[0], 2*nbComm, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

      //EACH PROCESS RUNS THE FOLLOWING LOOP ON COMMUNICATIONS AND CHECK IF IT MUST DO
      // SOMETHING
      for(int i=0;i<nbComm;i++)
	{
	  //COMPUTE THE ABSOLUTE INDEXES OF CLONES CONCERNED BY COMMUNICATION
	  ctm = temp[2*i]; cte = temp[2*i+1]; //ctm == clone_to_move | cte == clone_to_erase
	  //COMPUTE RANK OF PROC. SENDING AND PROC. RECEIVING
	  //ALL VARIABLES ARE INTEGERS SO DIVISIONS ARE EUCLIDIAN.
	  sender = ctm/local_Nc; dest = cte/local_Nc; 
	  tag = i;
	  if(dest == sender) //IF COMM. IS INTERNAL
	    {
	      if(sender == my_rank) //IF CLONES RESIDE IN PROC. my_rank
		{
		  //DO THE COPY
		  state[cte%local_Nc] = state[ctm%local_Nc];
		}
	    }else{
	    if(sender == my_rank) //IF PROC. MUST SEND A CLONE
	      {
		//COMPUTE THE LOCAL INDEX OF THE CLONE TO BE MOVED
		cloneIdx = ctm%local_Nc;
		MPI_Send(&state[cloneIdx], 3, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
	      }
	    else if(dest == my_rank) //IF PROC. MUST RECEIVE A CLONE
	      {
		//COMPUTE THE LOCAL INDEX OF THE CLONE TO BE ERASED
		cloneIdx = cte%local_Nc;
		MPI_Recv(&state[cloneIdx], 3, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, &status);
	      }
	  }
	  //SYNCHRONIZE PROC. NOT SURE ITS NEEDED. MIGHT SEVERELY HARM PERFORMANCE.
	  MPI_Barrier(MPI_COMM_WORLD);
	}

    } //TIMESTEPS

  MPI_Finalize();
} // MAIN()
	  

     

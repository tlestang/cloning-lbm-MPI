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
#include "src/global.h"
#endif

#include "src/initialize_lattice_arrays.h"
#include "src/streamCollCompute.h"
#include "src/boundaryConditions.h"
#include "src/force.h"
#include "src/write_vtk.h"

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
  double T = 10; // Total simulation time
  double dT = 2.5; // Cloning timestep
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
  //SPECIFIC PARAMETERS FOR PROGRESSIVE FORCING
  double tau0;
  
  //------------------

  //VARIABLES FOR TLGK
  double phi_alpha, phi_theor;
  double alpha=0.04, alphaMin = -0.5, alphaIncr = 0.004, alphaMax = 0.5;
  int NcPrime, deltaN, copyIdx, k;
  int nbrTimeSteps = floor(T/dT);
  int l= 0; int idx;

  int my_rank, p, local_Nc, cloneIdx, cloneIdxMin, tag = 0;
  int ctm, cte, nbComm;
  int sender, dest;
  bool flag;
  //-------------------

  //VARIABLES FOR LBM
  double *fin, *fout, *pivot, *rho, *ux, *uy;
  Dy = 4*Ly + 1, Dx = 2*(Dy-1) + 1;
  xmin = (Dx-1)/2; xmax = xmin + Lx;
  ymin = (Dy-1)/2 - Ly/2; ymax = ymin + Ly;
  double cs = 1./sqrt(3); double rho0 = 1.0;
  double nu = 1./3.*(tau-0.5);
  double u0 = cs*cs*Ma; double uxSum, uxMean;
  double beta0 = 8*nu*u0/((Dy-1)/2)/((Dy-1)/2); double a=1.0;
  double omega = 1.0/tau;
  double F;
  double delta_t = 1.0/t0; 
  int lbmTimesteps = floor(dT/delta_t);
  int transientLength = 10*tau0; //tau0 IS THE CHARACTERISTIC TIME OF THE FORCING
  
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);
  local_Nc = Nc/p;
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

  int temp[2*Nc]; 
  int comm_instru_Send[local_Nc]; int comm_instru_Recv[local_Nc];
  double s[Nc]; double s_;
  int nbCreatedCopies[Nc];
  double R_record[nbrTimeSteps]; double R, total_R;

  string path_to_control_run;
  path_to_control_run = "/home/thibault/lbm_code/source_code/L32_FOR_TRANSIENT/pops.datout";
  ifstream crHandle;

  MPI_Status status;

  ofstream result;

  cloneIdxMin = my_rank*local_Nc;
  //OPEN FILE FOR WRITING SCGF ON DISK
  if(my_rank==MASTER){result.open("phi_alpha_mpi_V2.datout");}
  
  //START LOOP ON VALUES OF ALPHA
  //while(alpha<alphaMax)
  //{
  //alpha += alphaIncr;
      //INITIALIZE CLONES TO STATE POPS.DATOUT FROM path_to_control_run
      crHandle.open(path_to_control_run.c_str());
      for(int j=0;j<local_Nc;j++)
	{
	  crHandle.seekg(0, std::ios::beg); //Set cursor to beginning of file
	  for(int x=0;x<Dx;x++)
	    {
	      for(int y=0;y<Dy;y++)
	  	{
	  	  for(int k=0;k<9;k++)
	  	    {
	  	      crHandle >> state[j][IDX(x,y,k)];
		    }
	  	}
	    }
	}
      crHandle.close();
      
      //TIME EVOLUTION OVER TOTAL TIME T (T/dT CLONING STEPS)
      for(int t=0;t<nbrTimeSteps;t++)
	{
	  // if(my_rank==MASTER)
	  //   {
	  //     cout << "Going for " << lbmTimesteps << "lbm timesteps" << endl;
	  // cout << " " << endl;
	  // cout << "Timestep nb " << t <<"/"<<nbrTimeSteps<<endl;
	  // cout << " " << endl;
	  //   }
	  R = 0.0;  
          
	  //LATTICE BOLTZMANN DYNAMICS
	  //W8 FOR DYNAMICS TO RELAX DUE TO PROGRESSIVE FORCING
	  for(int j=0;j<local_Nc;j++) // Loop on clones
	    {
	      tau0 = t0 + t0*(2.0*drand48()-1)/2.0; /*Random time between 0.5*t0 and 2.5*t0*/
	      transientLength = floor(10*tau0);
	      for(int t=0;t<transientLength;t++)
		{
		  a = -t/tau0;
		  beta = beta0*(1.0-exp(a));
		  streamingAndCollisionComputeMacroBodyForce(state[j], fout, rho, ux, uy, beta, tau);
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
		}
	    }
	  //SIMULATE THE SYSTEM THROUGH STATIONARY DYNAMICS DURING dT AND COMPUTE WEIGHT
	  for(int j=0;j<local_Nc;j++) // Loop on clones
	    {
	      s_ = 0;
	      for(int t=0;t<lbmTimesteps;t++)
		{
		  streamingAndCollisionComputeMacroBodyForce(state[j], fout, rho, ux, uy, beta, tau);
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
		  //if(my_rank==MASTER){cout << "F = " << F << endl;}
		  // COMPUTE WEIGHT
		  s_ += F;
		}
	      
	      s_ *= delta_t;
	      //STORE WEIGHT IN WEIGHTS ARRAY s[local_Nc]
	      s[j] = exp(alpha*s_);
	      //if(my_rank==1){cout << "clone " << j << " : s = " << s_ << endl;}
	      //UPDATE LOCAL AVERAGE WEIGHT
	      R += s[j];
	    }
       

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
	      cout << "MASTER STARTS PREPROCESSING EVOLUTION" << endl;
	      total_R /= Nc; //NORMALIZATION OF THE AVERAGE WEIGHT

	      NcPrime = 0; //NcPRIME IS THE NUMBER OF COPIES AFTER CLONING
	      for(int j=0;j<Nc;j++)
		{
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
	      cout << "CLONING/PRUNING WITH Ncprime = " << NcPrime << endl;
	      // CLONNG/PRUNING PHASE
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

	      //----------------------------------
	      // for(int j=0;j<Nc;j++)
	      // 	{
	      // 	  cout << "Clone " << j << " created " << nbCreatedCopies[j] << endl;
	      // 	}
	      
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
      //} // ALPHA
      MPI_Finalize();
} // MAIN()
	  

     

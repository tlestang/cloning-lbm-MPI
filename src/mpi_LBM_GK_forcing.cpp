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
#define MAX_CHARS 50

using namespace std;

//GLOBAL PARAMETERS FOR LBM
int c[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
int Dx, Dy, xmin, xmax, ymin, ymax;

int main(int argc, char *argv[])
//int main()
{
  int error;
  ofstream ascii_output;
  
  // --- PARAMETERS FOR TLGK ALGO. ---
  int Nc = 0;
  double T = 1.0;
  double dT = 1.0;
  double dT0 = 1.0;
  double alpha = 1.0;
  string path_to_folder, masterFolderName;
  //------------------------

  // --- PARAMETERS FOR LBM ---
  double tau = 1.0, beta = 1.0, Ma = 1.0, t0 = 1.0, U0=1.0, eps=1.0;
  int Lx = 0, Ly = 0;
  
  //MPI_INIT
  int my_rank, p, tag=0;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);
  
  //READ INPUT FILE
  ifstream input_file(argv[1]);
  if(input_file.is_open())
    {
      input_file >> Nc;
      input_file >> T;
      input_file >> dT;
      input_file >> dT0;
      input_file >> eps;
      input_file >> Lx; Ly = Lx;
      input_file >> tau;
      input_file >> U0;
      input_file >> alpha;
      input_file >> path_to_folder;
      input_file >> masterFolderName;
      input_file.close();
    }
  else{error=1;}

  input_file.close();

    if(my_rank==MASTER)
    {
      ascii_output.open("output_TLGK.out");
      if(error){
	cout << "COULD NOT OPEN INPUT FILE IN MASTER" << endl;
	ascii_output << "COULD NOT OPEN INPUT FILE IN MASTER" << endl;
      }
      for (int source=1;source<p;source++)
	{
	  tag = source;
	  MPI_Recv(&error, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
	  if(error){
	    cout << "COULD NOT OPEN INPUT FILE IN PROC " << source << endl;
	    ascii_output << "COULD NOT OPEN INPUT FILE IN PROC " << source << endl;
	  }
	}
    }
  else
    {
      MPI_Send(&error, 1, MPI_INT, MASTER, my_rank, MPI_COMM_WORLD);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank==MASTER)
    {
      cout << "All processes sucessfully read input file " << argv[1] << endl;
    }
  //------------------

  //VARIABLES FOR TLGK
  int NcPrime, deltaN, copyIdx, k;
  int nbrTimeSteps = floor(T/dT);
  int l= 0; int idx;

  int local_Nc, cloneIdx, cloneIdxMin;
  local_Nc = Nc/p;
  cloneIdxMin = my_rank*local_Nc;
  int ctm, cte, nbComm;
  int sender, dest;
  bool flag;
  //-------------------

  //VARIABLES FOR LBM
  double *fin, *fout, *pivot, *rho, *ux, *uy, *map;
  Dy = 4*Ly + 1, Dx = Dy;
  xmin = (Dx-1)/2; xmax = xmin + Lx;
  ymin = (Dy-1)/2 - Ly/2; ymax = ymin + Ly;
  int N = Dx*Dy*9;
  double cs = 1./sqrt(3); double rho0 = 1.0;
  double nu = 1./3.*(tau-0.5);
  double u0 = cs*cs*Ma; double uxSum, uxMean;
  double a=1.0;
  double omega = 1.0/tau;
  double F, beta0, T0, F0, oneOvF0;

    //COMPUTE CHARESTICTC VELOCITY AND TIME
  beta0 = (1./(Dx-1))*((double)Lx/(Dy-1))*U0*U0;
  T0 = Lx/U0;
  F0 = (U0*U0)*Lx*0.5;
  oneOvF0 = 1./F0;   
   
  double delta_t = 1.0/T0; //LBM time steps in units of physical time T0
  int lbmTimeSteps1 = floor(dT0*T0);
  int lbmTimeSteps2 = floor((dT-dT0)*T0);

  if(my_rank==MASTER)
    {

      cout << "READ input file" << endl;
      cout << "---------------" << endl;
      cout << "  Nc = " << Nc << endl;
      cout << "  T = " << T << " dT = " << dT << endl;
      cout << "  L = " << Lx << " Dx = " << Dx << " Dy = " << Dy << endl;
      cout << "  beta0 = " << beta0 << endl;
      cout << "---------------------" << endl;
      cout << "---------------------" << endl;

      ascii_output << "READ input file" << endl;
      ascii_output << "---------------" << endl;
      ascii_output << "  Nc = " << Nc << endl;
      ascii_output << "  T = " << T << " dT = " << dT << endl;
      ascii_output << "  L = " << Lx << " Dx = " << Dx << " Dy = " << Dy << endl;
      ascii_output << "  beta0 = " << beta0 << endl;
      ascii_output << "---------------------" << endl;
      ascii_output << "---------------------" << endl;
    }

  //INIT SEED TO RANK SO THAT EACH PORC HAS ITS OWN RANDOM SEQUENCE
  srand(my_rank+1);
  //EQUILIBRATE RANDOM GENRATOR (TO BE PROVEN RELEVANT)
  for(int i=0;i<1000;i++){rand();}
  
  double **state;
  state = new double*[local_Nc];
  for(int j=0;j<local_Nc;j++)
    {
      // ONE POP ARRAY PER CLONE
      state[j] = (double *) memalign(getpagesize(), N*sizeof(double)); 
    }
  //double *state;
  
  //state = (double *) memalign(getpagesize(), local_Nc*N*sizeof(double));
  //state = new double[local_Nc*N*sizeof(double)];
  fout = (double *) memalign(getpagesize(), N*sizeof(double));
  rho = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  ux = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  uy = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  map = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));


  int temp[2*Nc]; 
  int comm_instru_Send[local_Nc]; int comm_instru_Recv[local_Nc];
  double s[Nc]; double s_;
  int nbCreatedCopies[Nc];
  double R_record[nbrTimeSteps]; double R, total_R;

    //Set up variables and containers for output
  string folderName[local_Nc], instru;
  stringstream weightsFileName, copiesFileName, buf;
  ofstream weightsFile, copiesFile;
#ifdef FORCE_IO
  ofstream output_file;
#endif
  instru = "mkdir " + masterFolderName;
  if(my_rank==MASTER)
    {
      cout << "Created parent folder " << masterFolderName.c_str() << endl;
      ascii_output << "Created parent folder " << masterFolderName.c_str() << endl;
    }
  if(my_rank==MASTER){system(instru.c_str());}
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef FORCE_IO
  for(int i=0;i<local_Nc;i++)
    {
      buf << "clone_" << i + my_rank*local_Nc;
      folderName[i] = masterFolderName + buf.str();
      buf.str(string());
      buf.clear();
      instru = "mkdir " + folderName[i];
      system(instru.c_str());
    }
#endif

  // INITIALIZATION PROCEDURE ------------------------------------------------------------------------
  ifstream popFile;
  string path_to_file, fileName;

  path_to_file = path_to_folder+"popfiles_list.dat";
  ifstream fileList;
  char* buffer;
  if(my_rank==MASTER)
    {
      fileList.open(path_to_file.c_str());
      if(fileList.is_open())
	{
	  cout << "Index file " << path_to_file.c_str() << " successfully opened" << endl;
	  cout << "-------------" << endl;
	  ascii_output << "Index file " << path_to_file.c_str() << " successfully opened" << endl;
	  ascii_output << "-------------" << endl;
	}
      else
	{
	  cout << "ERROR : COULD NOT OPEN INDEX FILE " << path_to_file << endl;
	  ascii_output << "ERROR : COULD NOT OPEN INDEX FILE " << path_to_file << endl;
	}
      for(int j=local_Nc;j<Nc;j++)
  	{
  	  fileList >> fileName;
	  cout << "    Copy " << j << " init. on " << fileName.c_str() << endl;
  	  tag = j;
  	  dest = j/local_Nc;
  	  //BEWARE : MUST SEND STRING.LENGTH()+1 CHARACTERS BECAUSE OF THE FINAL \0 !
  	  //OTHERWISE STRING IS CORRUPTED
  	  MPI_Send((char*)fileName.c_str(), fileName.length()+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  	}
      for(int j=0;j<local_Nc;j++)
  	{
  	  fileList >> fileName;
	  cout << "    Copy " << j << " init. on " << fileName.c_str() << endl;
  	  instru = path_to_folder + fileName;
  	  popFile.open(instru.c_str(), ios::binary);
  	  popFile.read((char*)&state[j][0], N*sizeof(double));
  	  popFile.close();
  	}
      fileList.close();
    }
  else
    {
      buffer = new char[MAX_CHARS];
      for(int j=0;j<local_Nc;j++)
  	{
  	  tag = j + my_rank*local_Nc;
  	  MPI_Recv(buffer, MAX_CHARS, MPI_CHAR, MASTER, tag, MPI_COMM_WORLD, &status);
  	  instru = path_to_folder + string(buffer);
    	  popFile.open(instru.c_str(), ios::binary);
  	  popFile.read((char*)&state[j][0], N*sizeof(double));
  	  popFile.close();
  	}
      delete[] buffer;
    }
  MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank==MASTER)
    {
      cout << "--------------------- " << endl;
      cout << "--------------------- " << endl;
      ascii_output << "--------------------- " << endl;
      ascii_output << "--------------------- " << endl;
    }

    if(my_rank==MASTER){
      cout << "Looking good, about to enter timestep loop : " << endl;
      ascii_output << "Looking good, about to enter timestep loop : " << endl;
    }
    

   //TIME EVOLUTION OVER TOTAL TIME T (T/dT CLONING STEPS)
  for(int t=0;t<nbrTimeSteps;t++)
    {

      R = 0.0;

      if(my_rank==MASTER)
	{
	  weightsFileName << "weights_evolution_" << t;
	  copiesFileName << "copies_evolution_" << t;
	  instru = masterFolderName + weightsFileName.str();
	  weightsFile.open(instru.c_str(), ios::binary);
	  instru = masterFolderName + copiesFileName.str();
	  copiesFile.open(instru.c_str(), ios::binary);
	  cout << "    This is cloning step " << t << endl;
	  cout << "    Created files " << weightsFileName.str().c_str() << " and " << copiesFileName.str().c_str() << " in " << masterFolderName.c_str() << endl;
	  weightsFileName.str(string());
	  weightsFileName.clear();
	  copiesFileName.str(string());
	  copiesFileName.clear();
	}
      
      //SIMULATE THE SYSTEM DURING dT AND COMPUTE WEIGHT
      for(int j=0;j<local_Nc;j++) // Loop on clones
	{
	  
	  s_ = 0;
	  
#ifdef FORCE_IO
	  buf << "/evolution_" << t << "_" << "clone_" << j+my_rank*local_Nc;
	  fileName = folderName[j] + buf.str();
	  buf.str(string());
	  buf.clear();
	  output_file.open(fileName.c_str(), ios::binary);
#endif

	 
	  generate_random_field(Dx, map, error);
	 
	  for(int tt=0;tt<lbmTimeSteps1;tt++)
	    {
	      streamingAndCollisionComputeMacroBodyForceSpatial(state[j], fout, rho, ux, uy, beta0, map, tau, eps);
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
	      //	      memcpy(state[j], fout, N*sizeof(double));
	      // COMPUTE FORCE ON SQUARE
	      F = computeForceOnSquare(state[j], omega);
	      F = F*oneOvF0;
#ifdef FORCE_IO
	      output_file.write((char*)&F, sizeof(double));
#endif
	      // COMPUTE WEIGHT
	      s_ += F;
	    } //END LOOP ON TIMESTEPS
	  
	  for(int tt=0;tt<lbmTimeSteps2;tt++)
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
	      //memcpy(state[j], fout, N*sizeof(double));
	      pivot = state[j];
	      state[j] = fout;
	      fout = pivot;
	      // COMPUTE FORCE ON SQUARE
	      F = computeForceOnSquare(state[j], omega);
	      F = F*oneOvF0;
#ifdef FORCE_IO
	      output_file.write((char*)&F, sizeof(double));
#endif
	      // COMPUTE WEIGHT
	      s_ += F;
	    } //END LOOP ON TIMESTEPS
#ifdef FORCE_IO
	  output_file.close();
#endif
	  
	  s_ *= delta_t;
	  //STORE WEIGHT IN WEIGHTS ARRAY s[local_Nc]
	  s[j] = exp(alpha*s_);
	  //UPDATE LOCAL AVERAGE WEIGHT
	  R += s[j];
	}
      MPI_Barrier(MPI_COMM_WORLD);



      

      // ---------------------------------------------------------------------------------------------
      // ---------------------------------------------------------------------------------------------

      // //NOW TIME EVOLUTION OF COPIES IS DONE AND WE MUST :
      // // -- DETERMINE HOW MANY CLONES ARE GENERATED BY EACH COPY
      // // -- CLONE/PRUNE TO STICK TO A CONSTANT NUMBER OF COPIES AT T + dT
      // // -- CREATE INITIAL CONDITION FOR NEXT ITERATION (TIME EVOLUTION)
       

      //COMMUNICATIONS BETWEEN MASTER AND OTHER PROCESSES PRIOR TO CLONING STEP
      if(my_rank==MASTER)
      	{
      	  total_R = R; //COMPUTE MASTER's CONTRIB. TO TOTAL AVERAGE WEIGHT TOTAL_R
      	  //GATHER WEIGHTS AND LOCAL AVERAGE WEIGHTS FROM OTHER PROCESSES
      	  for(int source=1;source<p;source++)
      	    {
      	      tag = source;
      	      MPI_Recv(&s[source*local_Nc], local_Nc, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      	      tag = 2*source;
      	      MPI_Recv(&R, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      	      total_R += R;
      	    }
      	}
      else //SEND s[] VALUES FOR THE PROCESS'S CLONES AND LOCAL AVERAGE WEIGHT TO MASTER
      	{
      	  tag = my_rank;
      	  MPI_Send(&s[0], local_Nc, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD);
      	  tag = 2*my_rank;
      	  MPI_Send(&R, 1, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD);
      	}
      MPI_Barrier(MPI_COMM_WORLD);
      
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
      	      //EACH COPY j LEADS TO nbCreatedCopies[j] CLONES (INCLUDING ITSELF)
      	      nbCreatedCopies[j] = floor(s[j]/total_R + drand48());
      	      NcPrime += nbCreatedCopies[j]; //COMPUTE NcPrime BY COUNTING
      	    }
      	  //WRITE NB OF CREATED COPIES
      	  copiesFile.write((char*)&nbCreatedCopies[0], Nc*sizeof(int));
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

      	  copiesFile.write((char*)&nbCreatedCopies[0], Nc*sizeof(int));
      	  copiesFile.close();
	  
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
      	  // for(int i=0;i<nbComm;i++)
      	  //   {
          //     cout << temp[2*i] << " ---> " << temp[2*i+1] << endl;
          //   }
	} // IF MASTER

      //MPI_Barrier(MPI_COMM_WORLD);
      //BROADCAST OF NB OF COMM. FROM MASTER
      MPI_Bcast(&nbComm, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
      //MAKE SURE EACH AND EVERY PROCESS KNOWs THE NB OF COMM. BEFORE RECEVING COMM TABLE
      //MPI_Barrier(MPI_COMM_WORLD);
      //BROADCAST OF COMMUNICATION TABLE TEMP[] FROM MASTER
      MPI_Bcast(&temp[0], 2*nbComm, MPI_INT, MASTER, MPI_COMM_WORLD);

      //if(my_rank==MASTER){cout << "Now doing communications" << endl;}

      // //----------------------------
      // // for(int j=0;j<local_Nc;j++)
      // // 	{
      // // 	  cout << "Proc " << my_rank << " clone " << j << " : " << computeForceOnSquare(state[j], omega) << endl;
      // // 	}
      
      

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
      	      if(my_rank == sender) //IF CLONES RESIDE IN PROC. my_rank
      		{
      		  //DO THE COPY
      		  //state[cte%local_Nc] = state[ctm%local_Nc];
      		  memcpy(state[cte%local_Nc], state[ctm%local_Nc], N*sizeof(double));
      		}
      	    }else{
      	    if(my_rank == sender) //IF PROC. MUST SEND A CLONE
      	      {
      		//COMPUTE THE LOCAL INDEX OF THE CLONE TO BE MOVED
      		cloneIdx = ctm%local_Nc;
      		MPI_Send(state[cloneIdx], N, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      	      }
      	    else if(my_rank == dest) //IF PROC. MUST RECEIVE A CLONE
      	      {
      		//COMPUTE THE LOCAL INDEX OF THE CLONE TO BE ERASED
      		cloneIdx = cte%local_Nc;
      		MPI_Recv(state[cloneIdx], N, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, &status);
      	      }
      	  }
      	  //SYNCHRONIZE PROC. NOT SURE ITS NEEDED. MIGHT SEVERELY HARM PERFORMANCE.
      	}
      MPI_Barrier(MPI_COMM_WORLD);
        
    } //TIMESTEPS
  if(my_rank==MASTER)
    {
      instru = masterFolderName + "rvalues.dat";
      //WRITE SEQUENCE OF R's ON DISK FOR SCGF CALCULATION LATER
      ofstream rvalues(instru.c_str(), ios::binary);
      rvalues.write((char*)&R_record[0], nbrTimeSteps*sizeof(double));
      rvalues.close();
    }

    // FREE MEMORY -------------------------
  for(int j=0;j<local_Nc;j++)
    {
      delete[] state[j];
    }
  delete[] state;
  
  // -------------------------------------
  
  MPI_Finalize();
} // MAIN()
	  

     

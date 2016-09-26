#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;

int main()
{
  stringstream fileName, ssbuf;
  string dirName;
  int Nc, facquVTK;
  double T, dT;
  string buf, masterFolderName = "outputFolder/", folderName, vtkName;
  
  cout << "Input output directory (/ AT THE END)" << endl;
  //cin >> dirName;
  dirName = "/home/thibault/tailleur_lecomte/lbm/runs_from_psmn/L32_T_1500_dT_2_k_20_eps_0.01_128_clones_IO_ZERO_MEAN/";
  cout << "Input number of replica" << endl;
  //cin >> Nc; 
  Nc = 128;
  cout << "Input total time of run and cloning time" << endl;
  //cin >> T; cin >> dT;
  cout << "Input VTK file sampling period (in timesteps)" << endl;
  //cin >> facquVTK;
  T = 1500; dT = 2;

  int nt = floor(T/dT);
  int L; double U, T0, dt;
  L = 32; U = 0.03; T0 = L/U;
  dt = 1/T0;
  int ndt = dT/dt;
  
  buf = dirName + "labels.dat";
  ifstream labels(buf.c_str(), ios::binary);
  int *labelArray;
  labelArray = new int[Nc];
  int cloneIdx;
  ofstream reconstructed; ifstream forceFile;
  double *f = new double[nt*ndt];

  /*Create master folder*/
  buf = "rm -rf outputFolder/ && mkdir outputFolder";
  system(buf.c_str());
  
  /* LOOP ON CLONES */

  for (int j=0;j<1;j++)
    {
      //cout << "Reconstructing trajetory and drag signal for clone " << j+1 << "/" << Nc << endl;
      
      labelArray[j] = j;

      /*CREATE REPO. FOR CLONE J IN OUTPUT FOLDER FOR RECONSTRUCTED TRAJ.*/
      ssbuf << "clone_" << j << "/";
      folderName = masterFolderName + ssbuf.str();
      buf = "mkdir " + masterFolderName + ssbuf.str();
      ssbuf.str(string()); ssbuf.clear();
      system(buf.c_str());
      /*OPEN RECONSTRCUTED TRAJ. FILE*/
      buf = folderName + "data_force.datout";
      reconstructed.open(buf.c_str(), ios::binary);



      /* LOOP ON TIME STEPS BACKWARDS IN TIME */
      cloneIdx = j;
      int lbmTimeSteps2 = floor(dT*T0);
      int nbrVTK = floor(lbmTimeSteps2/facquVTK);
      int vtkIdxStart = nbrVTK*nt -1;
      for(int t=nt-1;t>-1;t--)
	{
	  cout<<"Evolution : "<<t<<" | reading clone " << cloneIdx << " ";
	  /* OPEN TRAJ. FILE FOR CLONE J*/
	  ssbuf << "clone_" << cloneIdx << "/data_force.datout";
	  buf = dirName + ssbuf.str();
	  ssbuf.str(string()); ssbuf.clear();
	  forceFile.open(buf.c_str(), ios::binary);
	  
	  if(t>0)
	    {
	      cout << "Setting cursor to " << Nc*(t-1) << endl;
	      labels.seekg(Nc*(t-1)*sizeof(int), ios::beg);
	      forceFile.seekg(t*ndt*sizeof(double), ios::beg);
	    }
	  else{labels.seekg(0, ios::beg); forceFile.seekg(0, ios::beg);}
	  
	  forceFile.read((char*)&f[t*ndt], ndt*sizeof(double));
	  forceFile.close();
	  
	  // for(int vtkIdx=vtkIdxStart;vtkIdx>vtkIdxStart-nbrVTK;vtkIdx--)
	  //   {
	  //     ssbuf << "clone_" << cloneIdx << "/" << "fluid_t" << vtkIdx << ".vtk";
	  //     vtkName = dirName + ssbuf.str(); ssbuf.str(string()); ssbuf.clear(); 

	  //     ssbuf << "clone_" << j << "/" << "fluid_t" << vtkIdx << ".vtk";
	  //     buf = "cp " + vtkName + " " + masterFolderName + ssbuf.str();
	  //     ssbuf.str(string()); ssbuf.clear();
	  //     system(buf.c_str());

	  //   }
	  // vtkIdxStart -= nbrVTK;


	  labels.read((char*)&labelArray[0], Nc*sizeof(int));
	  cloneIdx = labelArray[cloneIdx];
	}
      reconstructed.write((char*)&f[0], nt*ndt*sizeof(double));
      reconstructed.close();
	
    }
  delete[] f;
}
      

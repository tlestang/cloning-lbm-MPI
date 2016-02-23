#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

void generate_random_field(int n, double* map_, int &error);
void take_curl(int n, double* psi, double* curlx, double* curly, int &error);

int main()
{
  int n = 257;
  double *grid, *vx, *vy;
  grid = new double[n*n];
  vx = new double[n*n];
  vy = new double[n*n];
  
  int error;

  generate_random_field(n, grid, error);

  //WRITE NOISE MAP ON FILE MATRIX.DAT TO BE READ WITH GNUPLOT (IMAGE MODE)
  // ofstream result("result_matrix.dat");
  // for(int i=0;i<n;i++)
  //   {
  //     for(int j=0;j<n;j++)
  // 	{
  // 	  result << grid[j + n*i] << " ";
  // 	}
  //     result << endl;
  //   }

  //  result.close();
  //system("gnuplot plot_RF.gnu");
  cout << "DONE with generation of RF." << endl;
  cout << "Procesding to curl calculation NOW" << endl;

  take_curl(n,grid,vx,vy,error);
  cout << "DONE with curl" << endl;
  cout << "write result on disk NOW" << endl;
  //WRITE velocity field ON FILE ux(uy)_matrix.DAT TO BE READ WITH GNUPLOT (IMAGE MODE)
  ofstream resultx("ux_matrix.dat");
  ofstream resulty("uy_matrix.dat");
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
  	{
  	  resultx << vx[j + n*i] << " ";
	  resulty << vy[j + n*i] << " ";
  	}
      resultx << endl;
      resulty << endl;
    }

  cout << "DONE" << endl;
  cout << vx[0] << endl;
  delete[] grid; delete[] vx; delete[] vy;
   
}
    

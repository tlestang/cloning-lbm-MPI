#include <cstdlib>
#include <cmath>

int generatePerturbedState(double *fin, double **popsForPerturb, double eps, int NN, int N)
  
{
  double *fprim; double norm_0, norm_prim, randomPrefactor;
  double oneOvA, a, weight;
  int error = 0;
  fprim = new double[N];
  double sum_prefactors = 0.0;
  for (int nn=0;nn<NN;nn++)
    {
      randomPrefactor = drand48();
      sum_prefactors += randomPrefactor;
      for(int i=0;i<N;i++)
	{
  	  fprim[i] += randomPrefactor * popsForPerturb[nn][i];
  	}
    }
  norm_0 = 0.0; norm_prim = 0.0;
  for(int i=0;i<N;i++)
    {
      norm_0 += fin[i]*fin[i];
      norm_prim += fprim[i]*fprim[i];
    }
  norm_0 = sqrt(norm_0); norm_prim = sqrt(norm_prim);
  weight = norm_0 / norm_prim;
  a = 1.0 + sum_prefactors*eps*weight;
  oneOvA = 1./a;
  for(int i=0;i<N;i++)
    {
      fin[i] += eps*fprim[i]*weight;
      fin[i] = fin[i]*oneOvA;
    }
  delete[] fprim;
}

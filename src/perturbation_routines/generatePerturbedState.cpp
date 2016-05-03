#include <cstdlib>
#include <cmath>

int generatePerturbedState(double *fin, double **popsForPerturb, double eps, int NN, int N)
  
{
  double *fprim; double norm, randomPrefactor;
  double oneOvNorm, oneOvA, a;
  int error = 1;
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
  norm = 0.0;
  for(int i=0;i<N;i++)
    {
      norm += fprim[i]*fprim[i];
    }
  norm = sqrt(norm);
  oneOvNorm = 1./norm;
  a = 1.0 + sum_prefactors*eps*oneOvNorm;
  oneOvA = 1./a;
  for(int i=0;i<N;i++)
    {
      fin[i] += eps*fprim[i]*oneOvNorm*oneOvA;
    }
}

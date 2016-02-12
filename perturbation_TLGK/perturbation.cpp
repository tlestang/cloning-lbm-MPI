#ifndef __global__
#define __global__
#include "../src/global.h"
#endif

#include "perturbation.h"

int perturbPopulationField(int n, double *f)
{
  double *psi, *du, *dv;
  double u,v,rho;
  double eueu, eu, u2, ftemp;
  double eps = 1.0e-2;
  int error;
  psi = new double[n*n];
  du = new double[n*n];
  dv = new double[n*n];
  
  generate_random_field(n, psi, error);
  apply_mask_TLGK(n, psi, 4.0);
  take_curl(n,psi,du,dv,error);

  for(int x=0;x<n;x++)
    {
      for(int y=0;y<n;y++)
	{
	  rho = 0.0;
	  u = v = 0.0;
	  for (int k=0;k<9;k++)
	    {
	      ftemp = f[IDX(x,y,k)];
	      rho += ftemp;
	      u += ftemp*c[k][0];
	      v += ftemp*c[k][1];
	    }
	  u = u/rho; v = v/rho;
	  u += eps*du[idx(x,y)]; v += eps*dv[idx(x,y)];
	  u2 = -1.5*(u*u + v*v);
	  for(int k=0;k<9;k++)
	    {
	      eu = c[k][0]*u + c[k][1]*v;
	      eueu = 4.5*eu*eu;
	      f[IDX(x,y,k)] = w[k]*rho*(1.0+3.0*eu+eueu+u2);
	    }
	}
    }

  delete[] psi; delete[] du; delete[] dv;
}
	  
				    
  

  

#include <cmath>

#ifndef __global__
#define __global__
#include "../lbm_src/global.h"
#endif

using namespace std;

double ramp(double x);

void generate_mask(int n, double *mask, double d0)
{
  double dist;
  int ii;
  //NORTH AND SOUTH WALLS
  for (int x=0;x<n;x++)
    {
      for(int y=0;y<n;y++)
	{
	  if(n-1-y < y){dist = n-1-y;}
	  else{dist = y;}
	  if(dist < d0)
	    {
	      mask[idx(x,y)]= ramp(dist/d0);
	    }
	  else{mask[idx(x,y)] = 1.0;}
	}
    }
  double distSQ;
  //SQUARE
  for(int x=0;x<n;x++)
    {
      for(int y=0;y<n;y++)
	{
	  ii=0;
	  distSQ = 2*n;
	  for(int xsq=xmin;xsq<xmax+1;xsq++)
	    {
	      for(int ysq=ymin;ysq<ymax+1;ysq++)
		{
		  dist = sqrt((x-xsq)*(x-xsq) + (y-ysq)*(y-ysq));
		  if(dist < distSQ){distSQ = dist;}
		}
	    }
	  if(distSQ < d0)
	    {
	      mask[idx(x,y)] = ramp(distSQ/d0);
	    }
	  else{mask[idx(x,y)] = 1.0;
	}
	}
    }
}
    

double ramp(double x)
{
  double a,b,c;
  a = 15.0/8.0;
  b = 10./8.;
  c = 3./8.;

  return a*x - b*x*x*x + c*x*x*x*x*x;
}

#define idx(i,j) ((j) + n*(i))

void take_curl(int n, double* psi, double* curlx, double* curly, int &error)
{
  double diffx, diffy;
  error = 0;
  // Span the bulk to do 2nd order central differences
  for (int x=1;x<n-1;x++)
    {
      for (int y=1;y<n-1;y++)
	{
	  diffy = psi[idx(x,y+1)] - psi[idx(x,y-1)];
	  diffx = psi[idx(x+1,y)] - psi[idx(x-1,y)];
	  curlx[idx(x,y)] = - diffy/2.0;
	  curly[idx(x,y)] = diffx/2.0;
	}
    }
  // Compute gradients along x axis at west and east boundaries
  // Backward 2nd order for x = 0
  // Forward 2nd order for x =  n-1
  for (int y=1;y<n-1;y++)
    {
      //WEST BOUNDARY
      diffx = - psi[idx(2,y)]+4.0*psi[idx(1,y)]-3.0*psi[idx(0,y)];
      diffy = psi[idx(0,y+1)] - psi[idx(0,y-1)];
      curly[idx(0,y)] = -(diffx/2.0);
      curlx[idx(0,y)] = diffy/2.0;
      //EAST BOUNDARY
      diffx = 3.0*psi[idx(n-1,y)]-4.0*psi[idx(n-2,y)]+psi[idx(n-3,y)];
      diffy = psi[idx(n-1,y+1)] - psi[idx(n-1,y-1)];
      curly[idx(n-1,y)] = -(diffx/2.0);
      curlx[idx(n-1,y)] = diffy/2.0;
    }
  
  // Compute gradients along y axis at south and east boundaries
  // Backward 2nd order for y = 0
  // Forward 2nd order for y =  n-1
  for (int x=1;x<n-1;x++)
    {
      //SOUTH BOUNDARY
      diffy =  - psi[idx(x,2)]+4.0*psi[idx(x,1)]-3.0*psi[idx(x,0)];
      diffx = psi[idx(x+1,0)] - psi[idx(x-1,0)];
      curlx[idx(x,0)] = diffy/2.0;
      curly[idx(x,0)] = -(diffx/2.0);
      //NORTH BOUNDARY
      diffy = 3.0*psi[idx(x,n-1)]-4.0*psi[idx(x,n-2)]+psi[idx(x,n-3)];
      diffx = psi[idx(x+1,n-1)] - psi[idx(x-1,n-1)];
      curlx[idx(x,n-1)] = diffy/2.0;
      curly[idx(x,n-1)] = -(diffx/2.0);
    }
  //CORNERS
  //SOUTH WEST (FORWARDY FORWARDX)
  diffy = -psi[idx(0,2)]+4.0*psi[idx(0,1)]-3.0*psi[idx(0,0)];
  diffx = -psi[idx(2,0)]+4.0*psi[idx(1,0)]-3.0*psi[idx(0,0)];
  curlx[idx(0,0)] = diffy/2.0;
  curly[idx(0,0)] = -(diffx/2.0);
  //SOUTH EAST (FORWARDY BACKWARDX)
  diffy = -psi[idx(n-1,2)]+4.0*psi[idx(n-1,1)]-3.0*psi[idx(n-1,0)];
  diffx = 3.0*psi[idx(n-1,0)]-4.0*psi[idx(n-2,0)]+psi[idx(n-3,0)];
  curlx[idx(n-1,0)] = diffy/2.0;
  curly[idx(n-1,0)] = -(diffx/2.0);
  //NORTH WEST (BACKWARDY FORWARDX)
  diffy = 3.0*psi[idx(0,n-1)]-4.0*psi[idx(0,n-2)]+psi[idx(0,n-3)];
  diffx = -psi[idx(2,n-1)]+4.0*psi[idx(1,n-1)]-3.0*psi[idx(0,n-1)];
  curlx[idx(0,n-1)] = diffy/2.0;
  curly[idx(0,n-1)] = -(diffx/2.0);
  //NORTH EAST (BACKWARDY BACKWARDX)
  diffy = 3.0*psi[idx(n-1,n-1)]-4.0*psi[idx(n-1,n-2)]+psi[idx(n-1,n-3)];
  diffx = 3.0*psi[idx(n-1,n-1)]-4.0*psi[idx(n-2,n-1)]+psi[idx(n-3,n-1)];
  curlx[idx(n-1,n-1)] = diffy/2.0;
  curly[idx(n-1,n-1)] = -(diffx/2.0);
}

#include <cmath>

double f(int x, int y)
{
  double d, a;
  double sigma = 32.0;
  d = x*x + y*y;
  a = exp(-d/(2*sigma*sigma));
  return a;
}

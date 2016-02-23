#include <cmath>
#include <cstdlib>

double randNormal(const double mean_, const double sigma_)
{
  /* Return a random number sampled in N(mean_, sigma_).
     Box-Muller method.
  */

  double x1, x2, w;
  do {
    x1 = 2.0 * (rand () / double (RAND_MAX)) - 1.0;
    x2 = 2.0 * (rand () / double (RAND_MAX)) - 1.0;
    w = x1 * x1 + x2 * x2;
  } while (w >= 1.0);

  w = sqrt (-2.0 * log (w)/w);
  const double y1 = x1 * w;
  const double y2 = x2 * w;

  return mean_ + y1 * sigma_;
}

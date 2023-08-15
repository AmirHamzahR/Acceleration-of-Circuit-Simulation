#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
double t_k = 0;
double t_stop = 10;
double h_k = 0.1;
double t_br = 1;
int maxord = 3;
double h_o = 0.1;
int n = 1;
double h;
int k = 0;
double t_o = 0;

while (t_k < t_stop)
{
if (t_k + h_k > t_br)
{
h_k = t_br - t_k;
}

// Solve circuit at t_k + h_k
bool convergedSlowlyOrFailed = false;
if (convergedSlowlyOrFailed)
{
  n = 1;
  double h = h_k / 2;
  h_k = h;
}
else
{
  n = min(n + 1, maxord);
  double h = min(2 * h_k, h);
  if (h < 1e-6)
  {
    cout << "Timestep too small, ending simulation." << endl;
    return 0;
  }
  if (t_k + h_k == t_br)
  {
    n = 1;
  }
  t_k = t_k + h_k;
  h_k = h;
  k++;
}
}

return 0;
}




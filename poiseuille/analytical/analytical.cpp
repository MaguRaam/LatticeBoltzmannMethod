// analytical solution for laminar flow in a square duct
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>

using namespace std;

inline vector<double> linspace(double start, double stop, int n)
{
    vector<double> result(n);
    double step = (stop - start) / double(n - 1);
    for (int i = 0; i < n; ++i)
        result[i] = start + i * step;
    return result;
}


struct channel_flow
{

  double operator()(double y, double z)
  {
    // compute z/c
    double zbyc = z/c;

    // lambda to compute kth term in the sum
    auto term = [=](int k)
    {
        double minus_one_power_k = (k%2 == 0) ? 1 : -1;
    
        double alpha_k = (2.0*k - 1.0)*0.5*M_PI;
        double alpha_k_cube = alpha_k*alpha_k*alpha_k;

        return (minus_one_power_k/alpha_k_cube)*( cosh(alpha_k*y/c)/cosh(alpha_k*b/c) ) * cos(alpha_k*z/c) ;
	
    };
    
    // compute sum of n terms
    double sum = 0.0;
      
    for (int k = 1; k < n; ++k) sum += term(k);
    
    return   u0*(1.0 - zbyc*zbyc + 4.0*sum); 
  
  }
  

  int n = 200;                     // no of terms in the series
  double dpdx = 1;                 // pressure gradient
  double mu = 0.2286;             // dynamic viscosity
  double c = 49;                  // half length along z 
  double b = 49;                  // half length along y
  double u0 = (4.0*c*c*dpdx)/(8.0*mu);
};


int main()
{
  channel_flow u;

  double z = 0.0;
  auto y_vec = linspace(-49, 49, 1000);
  std::vector<double> u_vec(y_vec.size());

  // compute velocity along y
  std::transform(begin(y_vec), end(y_vec), begin(u_vec), [&u, &z](double y){return u(y, z);});

  // compute maximum of velocity
  double max = *std::max_element(begin(u_vec), end(u_vec));

  ofstream file("plot.dat", ios::out) ;
  file.flags( ios::dec | ios::scientific );
  file.precision(16) ;

  for (int i = 0; i < y_vec.size(); ++i)
    file << y_vec[i] + 48.5 << " " << u_vec[i]/max << "\n";
}








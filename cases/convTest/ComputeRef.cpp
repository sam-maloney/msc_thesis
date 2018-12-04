#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <cctype>
#include <vector>

double computeRho(double, double, double);
void writeHeader(std::fstream&, std::string, std::string, int, int, int, int);
void writeFooter(std::fstream&, std::string, std::string);

//int NX, NY;

int main(int argc, char** argv)
{
  // read # of grid cells from cmd line if present, otherwise use default value
//  NX = ( argc > 1 ? std::stoi(argv[1]) : 400 );
//  NY = ( argc > 2 ? std::stoi(argv[2]) : 400 );
//  double eps = ( argc > 3 ? std::stod(argv[3]) : 0.0 );

  std::fstream rhoFile;

  rhoFile.open("./ref_rho0_calc", std::fstream::out | std::fstream::trunc);

  if ( !rhoFile.is_open() ) {
    std::cerr << "Unable to open rho_ref0_calc file!\n" << std::endl;
    return -1;
  }

  rhoFile << std::setprecision(15);

  double sigma = 0.2;
  double mu = 0.0;

  double rhoPrev = 1.0;
  double Phi = 0.0;

  std::vector<int> NYs = {25, 50, 100, 200, 400, 800, 1600, 3200, 6400};

  for (int NY : NYs) {
    for (int i = 0; i < NY; ++i) {
      double y = -2.0 + (i+0.5)*4.0/NY;

      /// quintic
      if        ( y < -1.5 ) {
        Phi = 0.0; // left constant region
      } else if ( y <  1.5 ) {
        Phi = 0.024691358024691358*y*y*y*y*y - 0.185185185185185185*y*y*y + 0.625*y + 0.5; // left cubic
      } else {
        Phi = 1.0;// right constant region
      }

      double rho = computeRho(y, Phi, rhoPrev);
//      rho += eps/(sigma*2*M_PI)*std::exp(-0.5*std::pow((y-mu)/sigma,2));
  //    double p   = 0.75*std::pow(rho, 4.0/3.0);
  //    double T   = (7.0/5.0)*(p/rho);
  //    double U   = -0.394688351072542/rho;

      rhoFile << y << " " << rho << "\n";

      rhoPrev = rho;
    }
  }

  rhoFile.close();

  return 0;
}

double computeRho(const double y, double Phi, double rhoPrev)
{
  double B = 4.077889447236181; // Bernoulli constant at inlet
//  double B = 4.013862648201948;
  double rho = rhoPrev;

  unsigned int iter = 0;

  do {
    rhoPrev = rho;
    rho -= (0.077889447236181/(rho*rho) + 3.0*std::cbrt(rho) + Phi - B)/
           (-0.155778894472362/(rho*rho*rho) + std::pow(rho,-2.0/3.0));
    ++iter;
  } while ( ( std::abs(rho-rhoPrev) > 2e-15 ) && ( iter < 1000 ) );

  return rho;
}

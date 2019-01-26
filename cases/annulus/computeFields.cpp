#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <cctype>

double computeRho(double, double, double);
void writeHeader(std::fstream&, std::string, std::string, int, int, int, int);
void writeFooter(std::fstream&, std::string, std::string);

int Ntheta, Nr;
constexpr double pi = 3.141592653589793238462643383279;

int main(int argc, char** argv)
{
  // read # of grid cells from cmd line if present, otherwise use default value
  Ntheta = ( argc > 1 ? std::stoi(argv[1]) : 50 );
  Nr = ( argc > 2 ? std::stoi(argv[2]) : 20 );

  std::fstream pFile, TFile, UFile, PhiFile, rhoFile;

//  PhiFile.open("./0/potential", std::fstream::out | std::fstream::trunc);
  rhoFile.open("./0/rho"      , std::fstream::out | std::fstream::trunc);
    pFile.open("./0/p"        , std::fstream::out | std::fstream::trunc);
    TFile.open("./0/T"        , std::fstream::out | std::fstream::trunc);
    UFile.open("./0/U"        , std::fstream::out | std::fstream::trunc);

  if ( !( rhoFile.is_open() && pFile.is_open()
         && TFile.is_open() &&   UFile.is_open() ) ) {
    std::cerr << "Unable to open field files!\n" << std::endl;
    return -1;
  }

//  PhiFile << std::setprecision(15);
  rhoFile << std::setprecision(15);
    pFile << std::setprecision(15);
    TFile << std::setprecision(15);
    UFile << std::setprecision(15);

//  writeHeader(PhiFile, "potential", "Scalar", 0,  2, -2, 0);
  writeHeader(rhoFile, "rho"      , "Scalar", 1, -3,  0, 0);
  writeHeader(  pFile, "p"        , "Scalar", 1, -1, -2, 0);
  writeHeader(  TFile, "T"        , "Scalar", 0,  0,  0, 1);
  writeHeader(  UFile, "U"        , "Vector", 0,  1, -1, 0);

for (int k = 0; k < 4; ++k) {
  double rhoPrev = 1.0;
  for (int i = 0; i < Nr; ++i) {
    double r = 2.0 - (i+0.5)/Nr;
    double Phi = 0.935973200965767;
//    double Phi = 0.467986600482883*(std::tanh(20.0*y)+1.0);
    rhoPrev = 1.0;
    double rho = computeRho(r, Phi, rhoPrev);
    double p   = 0.75*std::pow(rho, 4.0/3.0);
    double T   = (7.0/5.0)*(p/rho);
    double U   = 0.592032526608812/(rho*r);

    int theta_i = Ntheta*k;

    for (int j = 0; j < Ntheta; ++j) {
      double theta = 2*pi*(theta_i+0.5)/(4*Ntheta);
      ++theta_i;
//      PhiFile << Phi << "\n";
      rhoFile << rho << "\n";
        pFile <<  p  << "\n";
        TFile <<  T  << "\n";
        UFile <<  "("<< U*cos(theta) << " " << U*sin(theta) << " 0)\n";
    }

    rhoPrev = rho;
  }
}

//  writeFooter(PhiFile, "fixedValue;\n        value           uniform 0.935973200965766",
//                       "fixedValue;\n        value           uniform 0");
  writeFooter(rhoFile, "calculated;\n        value           uniform 1.037662503736345",
                       "calculated;\n        value           uniform 0.837929055827603");
  writeFooter(  TFile, "fixedValue;\n        value           uniform 1.063019766808279",
                       "zeroGradient");
//  writeFooter(  UFile, "fixedValue;\n        value           uniform (0 -0.394688351072542 0)",
//                       "zeroGradient");
  writeFooter(  UFile, "groovyBC;\n         value           uniform (0 0 0);\n         variables (\n             \"Ur=-0.285272198078403;\"\n             \"r=sqrt(pow(pos().x,2)+pow(pos().y,2));\"\n             \"val=vector(Ur*pos().x/r, Ur*pos().y/r, 0);\"\n         );\n         valueExpression \"val\"",
                       "zeroGradient");
  writeFooter(  pFile, "fixedValue;\n        value           uniform 0.787896966248218",
                       "zeroGradient");

//  PhiFile.close();
  rhoFile.close();
    pFile.close();
    TFile.close();
    UFile.close();

  return 0;
}

double computeRho(const double r, double Phi, double rhoPrev)
{
  double B = 4.013862648201948; // Bernoulli constant at inlet
  double rho = rhoPrev;

  do {
    rhoPrev = rho;
//    rho -= (0.077889447236181/(rho*rho) + 3.0*std::cbrt(rho) + Phi - B)/
//           (-0.155778894472362/(rho*rho*rho) + std::pow(rho,-2.0/3.0));
    rho -= (0.175251256281407/(rho*rho*r*r) + 3.0*std::cbrt(rho) + Phi - B)/
           (-0.350502512562814/(rho*rho*rho*r*r) + std::pow(rho,-2.0/3.0));
  } while ( std::abs(rho-rhoPrev) > 2e-15 );

  return rho;
}

void writeHeader(std::fstream& file, std::string name, std::string type,
                 int mass, int length, int time, int temp)
{
  file << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
       << "| =========                 |                                                 |\n"
       << "| \\\\      /  F ield         | foam-extend: Open Source CFD                    |\n"
       << "|  \\\\    /   O peration     | Version:     4.0                                |\n"
       << "|   \\\\  /    A nd           | Web:         http://www.foam-extend.org         |\n"
       << "|    \\\\/     M anipulation  |                                                 |\n"
       << "\\*---------------------------------------------------------------------------*/\n"
       << "FoamFile\n"
       << "{\n"
       << "    version     2.0;\n"
       << "    format      ascii;\n"
       << "    class       vol" << type << "Field;\n"
       << "    object      " << name << ";\n"
       << "}\n"
       << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
       << "\n"
       << "dimensions      [" << mass << " " << length << " " << time << " " << temp << " 0 0 0];\n"
       << "\n";
  type[0] = std::tolower(type[0]);
  file << "internalField   nonuniform List<" << type << "> \n"
       << 4*Ntheta*Nr << "\n"
       << "(\n";
}

void writeFooter(std::fstream& file, std::string inBC, std::string outBC)
{
  file << ")\n"
       << ";\n"
       << "\n"
       << "boundaryField\n"
       << "{\n"
       << "    inlet\n"
       << "    {\n"
       << "        type            " << inBC << ";\n"
       << "    }\n"
       << "    outlet\n"
       << "    {\n"
       << "        type            " << outBC << ";\n"
       << "    }\n"
       << "    defaultFaces\n"
       << "    {\n"
       << "        type            empty;\n"
       << "    }\n"
       << "}\n"
       << "\n"
       << "\n"
       << "// ************************************************************************* //\n";
}

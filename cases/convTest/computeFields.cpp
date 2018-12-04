#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <cctype>

double computeRho(double, double, double);
void writeHeader(std::fstream&, std::string, std::string, int, int, int, int);
void writeFooter(std::fstream&, std::string, std::string);

int NX, NY;

int main(int argc, char** argv)
{
  // read # of grid cells from cmd line if present, otherwise use default value
  NX = ( argc > 1 ? std::stoi(argv[1]) : 400 );
  NY = ( argc > 2 ? std::stoi(argv[2]) : 400 );
  double eps = ( argc > 3 ? std::stod(argv[3]) : 0.0 );

  std::fstream pFile, TFile, UFile, PhiFile, rhoFile;

  PhiFile.open("./0/potential", std::fstream::out | std::fstream::trunc);
  rhoFile.open("./0/rho"      , std::fstream::out | std::fstream::trunc);
    pFile.open("./0/p"        , std::fstream::out | std::fstream::trunc);
    TFile.open("./0/T"        , std::fstream::out | std::fstream::trunc);
    UFile.open("./0/U"        , std::fstream::out | std::fstream::trunc);

  if ( !( PhiFile.is_open() && rhoFile.is_open() && pFile.is_open()
         && TFile.is_open() &&   UFile.is_open() ) ) {
    std::cerr << "Unable to open field files!\n" << std::endl;
    return -1;
  }

  PhiFile << std::setprecision(15);
  rhoFile << std::setprecision(15);
    pFile << std::setprecision(15);
    TFile << std::setprecision(15);
    UFile << std::setprecision(15);

  writeHeader(PhiFile, "potential", "Scalar", 0,  2, -2, 0);
  writeHeader(rhoFile, "rho"      , "Scalar", 1, -3,  0, 0);
  writeHeader(  pFile, "p"        , "Scalar", 1, -1, -2, 0);
  writeHeader(  TFile, "T"        , "Scalar", 0,  0,  0, 1);
  writeHeader(  UFile, "U"        , "Vector", 0,  1, -1, 0);

  double sigma = 0.2;
  double mu = 0.0;

  double rhoPrev = 1.0;
  double Phi = 0.0;

  for (int i = 0; i < NY; ++i) {
    double y = -2.0 + (i+0.5)*4.0/NY;

//    /// quadratic connected linear
//    if        ( y < -1.5 ) {
//      Phi = 0.0; // left constant region
//    } else if ( y < -1.0 ) {
//      Phi = 0.4*y*y + 1.2*y + 0.9; // left quadratic transition
//    } else if ( y <  1.0 ) {
//      Phi = 0.4*y + 0.5; // linear centre region
//    } else if ( y <  1.5 ) {
//      Phi = -0.4*y*y + 1.2*y + 0.1; // right quadratic transition
//    } else {
//      Phi = 1.0;// right constant region
//    }

//    /// opposing cubics
//    if        ( y < -1.5 ) {
//      Phi = 0.0; // left constant region
//    } else if ( y < 0 ) {
//      Phi = (148.0/999.0)*y*y*y + (2.0/3.0)*y*y + y + 0.5; // left cubic
//    } else if ( y < 1.5 ) {
//      Phi = (148.0/999.0)*y*y*y - (2.0/3.0)*y*y + y + 0.5; // right cubic
//    } else {
//      Phi = 1.0;// right constant region
//    }

    /// quintic
    if        ( y < -1.5 ) {
      Phi = 0.0; // left constant region
    } else if ( y <  1.5 ) {
      Phi = (2.0/81.0)*y*y*y*y*y - (5.0/27.0)*y*y*y + 0.625*y + 0.5; // left cubic
    } else {
      Phi = 1.0;// right constant region
    }

//    double Phi = 0.467986600482883*(std::tanh(2.0*y)+1.0) + offset;
    double rho = computeRho(y, Phi, rhoPrev);
    rho += eps/(2*sigma*std::sqrt(2*M_PI))*std::exp(-0.5*std::pow((y-mu)/sigma,2));
    double p   = 0.75*std::pow(rho, 4.0/3.0);
    double T   = (7.0/5.0)*(p/rho);
    double U   = -0.394688351072542/rho;

    for (int j = 0; j < NX; ++j) {
      PhiFile << Phi << "\n";
      rhoFile << rho << "\n";
        pFile <<  p  << "\n";
        TFile <<  T  << "\n";
        UFile <<  "(0 " << U << " 0)\n";
    }

    rhoPrev = rho;
  }

//  std::string waveOutCond = "waveTransmissive;\n";
//  waveOutCond += "        correctSupercritical off;\n";
//  waveOutCond += "        inletOutlet     off;\n";
//  waveOutCond += "        gamma           4.0/3.0;\n";
//  waveOutCond += "        phi             phi;\n";
//  waveOutCond += "        rho             rho;\n";
//  waveOutCond += "        psi             psi;\n";
//  waveOutCond += "        lInf            0.05;\n";
//
//  std::string pOutCond = waveOutCond;
//  pOutCond += "        field           p;\n";
//  pOutCond += "        fieldInf        2.529013440114932;\n";
//  pOutCond += "        value           uniform 2.529013440114932";
//
//  std::string UOutCond = waveOutCond;
//  UOutCond += "        field           U;\n";
//  UOutCond += "        fieldInf        (0 -0.158612314779744 0);\n";
//  UOutCond += "        value           uniform (0 -0.158612314779744 0)";
//
//  std::string TOutCond = waveOutCond;
//  TOutCond += "        field           T;\n";
//  TOutCond += "        fieldInf        1.422858679912701;\n";
//  TOutCond += "        value           uniform 1.422858679912701";

  writeFooter(PhiFile, "fixedValue;\n        value           uniform 1",
                       "fixedValue;\n        value           uniform 0");
  writeFooter(rhoFile, "calculated;\n        value           uniform 1",
                       "calculated;\n        value           uniform 2.488384030083816");
  writeFooter(  TFile, "fixedValue;\n        value           uniform 1.05",
                       "zeroGradient");
  writeFooter(  UFile, "fixedValue;\n        value           uniform (0 -0.394688351072542 0)",
                       "zeroGradient");
  writeFooter(  pFile, "zeroGradient",
                       "zeroGradient");

  PhiFile.close();
  rhoFile.close();
    pFile.close();
    TFile.close();
    UFile.close();

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
       << NX*NY << "\n"
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
//       << "    sides\n"
//       << "    {\n"
//       << "        type            cyclic;\n"
//       << "    }\n"
       << "    defaultFaces\n"
       << "    {\n"
       << "        type            empty;\n"
       << "    }\n"
       << "}\n"
       << "\n"
       << "\n"
       << "// ************************************************************************* //\n";
}
